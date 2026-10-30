"""
Interactive HTML report of a ResistNet run (runResistnet.py -o DIR/NAME):
DIR/NAME/NAME.report.html, with Plotly figures of each model's results
and, with --model both, of their comparison.

The report is built only from the run's output files (each figure from a
table the run writes), so it can be rebuilt with reportResistnet.py, and
the same tables can feed other reports (e.g. MultiQC).
"""
import html
import os
import pickle
import re
from datetime import datetime

import numpy as np
import pandas as pd

from resistnet.plotting import (AXIS, BLUES_FOR_LINES, COLOURS,
                                FAMILY_NAMES, GRID, MUTED, TEXT, TEXT2,
                                plotly_colorscale, plotly_layout)

TRANSFORMATIONS = {
    0: "linear", 1: "Ricker", 2: "reverse Ricker", 3: "inverse Ricker",
    4: "reverse-inverse Ricker", 5: "monomolecular",
    6: "reverse monomolecular", 7: "inverse monomolecular",
    8: "reverse-inverse monomolecular"
}
NUMBER = re.compile(r"-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?")
NP_SCALAR = re.compile(r"np\.\w+\(")
REFERENCE_NAMES = {
    "null": "Population effects only (null)",
    "distance_only": "Distance only (river km)",
    "samc_distance": "SAMC distance alone",
    "samc_flow": "SAMC flow alone",
}


# -------------------------------------------------------------- run tables
def river_distance_matrix(network):
    """
    River distance between every pair of sites: the sum of the reaches'
    lengths (length_col) along the path, as for the distance_only
    reference model (isolation by distance).

    Args:
        network: A network with its sites and incidence matrix.

    Returns:
        numpy.ndarray: Distances, ordered as network._points_labels.
    """
    import resistnet.resist_dist as rd
    import resistnet.utils as utils
    df = utils.nx_to_df(network._K)
    id_col = "EDGE_ID" if network.minimize else network.reachid_col
    lengths = df.groupby(id_col)[network.length_col].sum()
    lengths = lengths.loc[list(network._edge_order)]
    return rd.effectiveResistanceMatrix(network._points_snapped,
                                        network._inc, lengths)


def write_pairwise_table(network, prefix):
    """
    Writes {prefix}.pairwise.tsv: for every pair of sites, the genetic
    distance, the river distance (km; isolation by distance), and the best
    and model-averaged models' distances (for SAMC, the passage times of
    the two directions combined as fitted).

    Args:
        network: The fitted network (with genetic distances).
        prefix (str): Output prefix of the run.
    """
    samc = os.path.exists(f"{prefix}.Model-Best.CFPTMatrix.tsv")
    kind = "CFPTMatrix" if samc else "ResistanceMatrix"
    labels = [str(x) for x in network._points_labels.values()]
    iu = np.triu_indices(len(labels), 1)
    names = np.array(labels)
    table = pd.DataFrame({"site1": names[iu[0]], "site2": names[iu[1]],
                          "genetic_distance": network._gendist[iu],
                          "river_km": river_distance_matrix(network)[iu]})
    for which in ("Best", "Average"):
        path = f"{prefix}.Model-{which}.{kind}.tsv"
        if not os.path.exists(path):
            continue
        M = pd.read_csv(path, sep="\t", index_col=0)
        M.index = M.index.astype(str)
        M.columns = M.columns.astype(str)
        M = M.loc[labels, labels].to_numpy(dtype=float)
        if samc:
            M = (network.pair_distance(M) if hasattr(network,
                                                     "pair_distance")
                 else (M + M.T) / 2.0)
        table[f"{which.lower()}_model_distance"] = M[iu]
    table.to_csv(f"{prefix}.pairwise.tsv", sep="\t", index=False)


def _read(path, **kw):
    return (pd.read_csv(path, sep="\t", **kw) if os.path.exists(path)
            else None)


def describe(row, variables, samc):
    """A Hall of Fame model's description, from its row alone."""
    parts = []
    several = sum(row.get(v) == 1 for v in variables) > 1
    for v in variables:
        if row.get(v) != 1:
            continue
        t = int(row[f"{v}_trans"])
        label = TRANSFORMATIONS[t]
        if t:
            label += f", shape {float(row[f'{v}_shape']):g}"
        if samc and int(row.get(f"{v}_asym", 0) or 0) == 1:
            label += ", directional"
        if several:
            label += f", weight {float(row[f'{v}_weight']):.2f}"
        parts.append(f"{v} ({label})")
    flow = float(row.get("flow", 0) or 0) if samc else 0.0
    if flow:
        parts.append(f"flow {flow:+.3f}/km")
    text = " + ".join(parts) or "(no variables)"
    # The range of its resistance per km (SAMC: with the flow term, if any)
    max_r = row.get("max_resistance")
    has_var = any(row.get(v) == 1 for v in variables)
    if has_var and max_r is not None and pd.notna(max_r):
        text += f"; resistance 1-{float(max_r):.3g} per km"
    return text


class RunOutputs:
    """The output tables of one model's run (DIR/NAME/{model}/NAME.*)."""

    def __init__(self, prefix, model):
        self.prefix, self.model = prefix, model
        self.samc = model == "samc"
        self.hof = _read(f"{prefix}.HallOfFame.tsv", na_values="-")
        # Standard model: resistances per km (runs of earlier versions,
        # per reach, have no max_resistance column)
        self.per_km = not self.samc
        if (self.per_km and self.hof is not None
                and "max_resistance" not in self.hof.columns):
            raise ValueError(
                f"{prefix}.HallOfFame.tsv is a standard run of an earlier "
                "version (resistances per reach), which is no longer "
                "supported; rerun it")
        self.variables = [c[:-6] for c in self.hof.columns
                          if c.endswith("_trans")]
        self.hof["model"] = [describe(r, self.variables, self.samc)
                             for _, r in self.hof.iterrows()]
        self.models = self._unique_models()
        self.null = _read(f"{prefix}.Null-Model.tsv", index_col=0,
                          keep_default_na=False, na_values=[""])
        self.best = _read(f"{prefix}.bestModel.tsv", na_values="-")
        self.rvi = _read(f"{prefix}.varImportance.tsv")
        self.maw = _read(f"{prefix}.modavgWeights.tsv", na_values="-")
        self.pairwise = _read(f"{prefix}.pairwise.tsv")
        self.flow_profile = _read(f"{prefix}.flowProfile.tsv")
        self.points = _read(f"{prefix}.pointsTable.txt")
        kind = "TransitionEdges" if self.samc else "ResistanceEdges"
        self.edges = {w: _read(f"{prefix}.Model-{w}.{kind}.tsv",
                               na_values="NA")
                      for w in ("Average", "Best")}
        self.reach_values = pd.DataFrame()
        self.geometry = self._geometry()

    def _unique_models(self):
        """
        The Hall of Fame with each model once: the best model of several
        chains can be the same (per km, with maximum resistances differing
        only by search noise), splitting its Akaike weight.
        """
        from resistnet.compare import duplicate_key
        hof = self.hof.copy()
        hof["key"] = [duplicate_key(m, a)
                      for m, a in zip(hof["model"], hof["aic"])]
        hof["kept"] = hof["keep"].astype(str).str.lower() == "true"
        models = hof.groupby("key", sort=False).agg(
            model=("model", "first"), aic=("aic", "first"),
            r2m=("r2m", "first"), akaike_weight=("akaike_weight", "sum"),
            chains=("model", "size"), kept=("kept", "any"),
            **{c: (c, "first") for c in ("delta_aic_distance",
                                         "delta_aic_flow") if c in hof})
        return models.sort_values("aic").reset_index(drop=True)

    def _geometry(self):
        """Reach lines (lon, lat) by reach ID, from the run's network."""
        path = f"{self.prefix}.subgraph.net"
        spec = f"{self.prefix}.modelSpec.json"
        if not (os.path.exists(path) and os.path.exists(spec)):
            return {}
        reachid = pd.read_json(spec, typ="series")["reachid_col"]
        with open(path, "rb") as fh:
            G = pickle.load(fh)
        lines, values = {}, []
        for u, v, d in G.edges(data=True):
            geom = d.get("geometry")
            if geom is not None and hasattr(geom, "coords"):
                xy = np.asarray(geom.coords)[:, :2]
            else:
                xy = np.array([u[:2], v[:2]], dtype=float)
            lines[str(d.get(reachid))] = xy
            values.append({k: x for k, x in d.items()
                           if k in self.variables or k == reachid})
        # Each reach's values of the variables (mean over its edges)
        if values:
            self.reach_values = (pd.DataFrame(values)
                                 .groupby(reachid).mean(numeric_only=True))
        return lines

    def sites(self):
        """Site labels and coordinates (lon, lat)."""
        if self.points is None:
            return pd.DataFrame(columns=["label", "lon", "lat"])
        # Nodes are written as "(x, y)", or with numpy's scalar reprs
        # ("(np.float64(x), np.float64(y))")
        xy = [[float(v) for v in NUMBER.findall(NP_SCALAR.sub("", str(n)))]
              for n in self.points["Node"]]
        return pd.DataFrame({"label": self.points["Label"].astype(str),
                             "lon": [p[0] for p in xy],
                             "lat": [p[1] for p in xy]})


# ------------------------------------------------------------------ figures
def _go():
    import plotly.graph_objects as go
    return go


def fig_aic_profile(run):
    """AIC of each Hall of Fame model, best first."""
    go = _go()
    hof = run.models
    colour = COLOURS[run.model]
    kept = hof["kept"].to_numpy()
    weight = hof["akaike_weight"].fillna(0).to_numpy()
    size = 9 + 22 * np.sqrt(weight)
    hover = [f"<b>{html.escape(m)}</b><br>AIC {a:.2f} (Δ {d:.2f})<br>"
             f"Akaike weight {w:.3f}<br>marginal R² {r:.3f}<br>"
             f"best of {c} chain{'s' if c > 1 else ''}"
             for m, a, d, w, r, c in zip(hof["model"], hof["aic"],
                                         hof["aic"] - hof["aic"].min(),
                                         weight, hof["r2m"],
                                         hof["chains"])]
    fig = go.Figure()
    for flag, name, col in ((True, "kept for model averaging", colour),
                            (False, "not kept", MUTED)):
        sel = kept == flag
        if not sel.any():
            continue
        fig.add_trace(go.Scatter(
            x=hof.index[sel] + 1, y=hof["aic"][sel], mode="markers",
            name=name, hovertext=np.array(hover)[sel], hoverinfo="text",
            marker=dict(size=size[sel], color=col,
                        line=dict(color="white", width=1.5))))
    fig.add_hline(y=hof["aic"].min() + 2, line=dict(color=TEXT2, width=1,
                                                    dash="dot"),
                  annotation_text="best + 2", annotation_font_color=TEXT2,
                  annotation_position="top left")
    fig.update_layout(**plotly_layout(
        title="AIC of the Hall of Fame models (marker size: Akaike weight)",
        height=360))
    fig.update_xaxes(title="Model (best first)", dtick=1)
    fig.update_yaxes(title="AIC")
    return fig


def fig_importance(run):
    """Relative importance and model-averaged weight of each variable."""
    go = _go()
    from plotly.subplots import make_subplots
    colour = COLOURS[run.model]
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.18,
                        subplot_titles=("Relative importance (RVI)",
                                        "Model-averaged weight"))
    for col, table, column in ((1, run.rvi, "RVI"), (2, run.maw, "MAW")):
        if table is None:
            continue
        t = table.dropna(subset=[column]).sort_values(column)
        fig.add_trace(go.Bar(
            x=t[column], y=t["variable"], orientation="h",
            marker=dict(color=colour), text=[f"{v:.2f}" for v in t[column]],
            textposition="outside", textfont=dict(color=TEXT2),
            hovertemplate="%{y}: %{x:.3f}<extra></extra>", showlegend=False),
            row=1, col=col)
    fig.update_layout(**plotly_layout(height=max(260, 46 * len(run.variables)
                                                 + 120)))
    fig.update_xaxes(range=[0, 1.15], gridcolor=GRID, linecolor=AXIS)
    fig.update_yaxes(gridcolor="rgba(0,0,0,0)", linecolor=AXIS,
                     tickfont=dict(color=TEXT))
    for a in fig.layout.annotations:
        a.update(font=dict(color=TEXT, size=12), x=a.x - 0.2,
                 xanchor="left")
    return fig


def fig_contributions(run):
    """
    Standard model per km: the share of the best model's resistance per km
    due to distance (the baseline of 1) and to each of its variables,
    averaged over the network's km (bestModel.tsv, Contribution), so that
    distance can be compared with the variables.
    """
    b = run.best
    if not run.per_km or b is None or "Contribution" not in b:
        return None
    t = b.dropna(subset=["Contribution"])
    if t.empty:
        return None
    t = t.sort_values("Contribution")
    go = _go()
    colours = [COLOURS["reference"] if v == "distance" else COLOURS[run.model]
               for v in t["Variable"]]
    fig = go.Figure(go.Bar(
        x=t["Contribution"], y=t["Variable"], orientation="h",
        marker=dict(color=colours), text=[f"{v:.0%}" for v in
                                          t["Contribution"]],
        textposition="outside", textfont=dict(color=TEXT2),
        hovertemplate="%{y}: %{x:.1%} of the resistance per km"
                      "<extra></extra>", showlegend=False))
    fig.update_layout(**plotly_layout(
        height=max(200, 46 * len(t) + 110),
        title=dict(text="Share of the best model's resistance per km",
                   font=dict(size=14, color=TEXT), x=0, xanchor="left")))
    fig.update_xaxes(range=[0, 1.15], tickformat=".0%", gridcolor=GRID,
                     linecolor=AXIS)
    fig.update_yaxes(gridcolor="rgba(0,0,0,0)", linecolor=AXIS,
                     tickfont=dict(color=TEXT))
    return fig


def fig_pairwise(run):
    """Genetic distance against the model distances of each pair."""
    go = _go()
    pw = run.pairwise
    colour = COLOURS[run.model]
    fig = go.Figure()
    unit = "passage time (steps)" if run.samc else "resistance distance"
    for column, name, style in (
            ("average_model_distance", "Model average",
             dict(color="white", size=8, line=dict(color=colour, width=1.5))),
            ("best_model_distance", "Best model",
             dict(color=colour, size=8, line=dict(color="white",
                                                  width=1)))):
        if column not in pw:
            continue
        x, y = pw[column].to_numpy(), pw["genetic_distance"].to_numpy()
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="markers", name=name, marker=style,
            customdata=np.stack([pw["site1"], pw["site2"]], axis=1),
            hovertemplate=("%{customdata[0]} – %{customdata[1]}<br>"
                           "model distance %{x:.3g}<br>genetic distance "
                           "%{y:.4f}<extra>" + name + "</extra>")))
        if column == "best_model_distance":
            slope, intercept = np.polyfit(x, y, 1)
            xs = np.array([x.min(), x.max()])
            fig.add_trace(go.Scatter(
                x=xs, y=intercept + slope * xs, mode="lines",
                line=dict(color=TEXT2, width=1.5), hoverinfo="skip",
                name="Least-squares line (best model)"))
    r2m = float(run.hof.sort_values("aic")["r2m"].iloc[0])
    fig.update_layout(**plotly_layout(
        title=f"Genetic distance against model distance (best model: "
              f"marginal R² {r2m:.2f})", height=420))
    fig.update_xaxes(title=f"Model {unit}")
    fig.update_yaxes(title="Genetic distance")
    return fig


def fig_map(run, layers, colorbar_title, fixed_range=None, log=False):
    """
    Map of the network's reaches coloured by a value, with the sites; one
    layer per entry of layers ({name: Series by reach ID}), chosen from a
    menu. With log, the colours follow log10 of the values (for skewed
    positive values, such as resistances per km).
    """
    go = _go()
    lines = run.geometry
    fig = go.Figure()
    colours = BLUES_FOR_LINES
    n = len(colours)
    menu, vis_all = [], []
    for li, (name, raw) in enumerate(layers.items()):
        raw = raw.dropna()
        raw.index = raw.index.astype(str)
        if log:
            raw = raw[raw > 0]
        values = np.log10(raw) if log else raw
        lo, hi = (fixed_range if fixed_range is not None
                  else (float(values.min()), float(values.max())))
        span = hi - lo if hi > lo else 1.0
        colorbar = dict(title=dict(text=colorbar_title, side="right",
                                   font=dict(color=TEXT2)),
                        tickfont=dict(color=TEXT2), thickness=14,
                        outlinewidth=0)
        if log:
            ticks = _log_ticks(10 ** lo, 10 ** hi)
            colorbar.update(tickvals=np.log10(ticks).tolist(),
                            ticktext=[f"{t:g}" for t in ticks])
        vis = []
        # Reaches without a value, in grey
        missing = [r for r in lines if r not in values.index]
        xs, ys = _line_coords([lines[r] for r in missing])
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", hoverinfo="skip",
                                 line=dict(color=AXIS, width=2),
                                 showlegend=False, visible=li == 0))
        vis.append(len(fig.data) - 1)
        # Lines in bins of the colour ramp
        bins = np.clip(((values - lo) / span * n).astype(int), 0, n - 1)
        for b in range(n):
            ids = [r for r in values.index[bins == b] if r in lines]
            if not ids:
                continue
            xs, ys = _line_coords([lines[r] for r in ids])
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode="lines", hoverinfo="skip",
                line=dict(color=colours[b], width=3.2), showlegend=False,
                visible=li == 0))
            vis.append(len(fig.data) - 1)
        # Invisible markers at reach midpoints: hover and colour bar
        ids = [r for r in values.index if r in lines]
        mid = np.array([lines[r][len(lines[r]) // 2] for r in ids])
        fig.add_trace(go.Scatter(
            x=mid[:, 0] if len(mid) else [], y=mid[:, 1] if len(mid) else [],
            mode="markers", showlegend=False, visible=li == 0,
            marker=dict(size=10, opacity=0, color=values[ids],
                        colorscale=plotly_colorscale(colours), cmin=lo,
                        cmax=hi, colorbar=colorbar),
            customdata=np.column_stack([ids, raw[ids].to_numpy()])
            if ids else [],
            hovertemplate="reach %{customdata[0]}<br>" + colorbar_title
                          + " %{customdata[1]:.4g}<extra></extra>"))
        vis.append(len(fig.data) - 1)
        vis_all.append(vis)
        menu.append(name)

    # Sites, on every layer
    sites = run.sites()
    fig.add_trace(go.Scatter(
        x=sites["lon"], y=sites["lat"], mode="markers+text",
        text=sites["label"], textposition="top center",
        textfont=dict(color=TEXT, size=10), name="Sampling sites",
        marker=dict(size=9, color=COLOURS["samc"] if not run.samc
                    else COLOURS["resistance"],
                    line=dict(color="white", width=1.5)),
        hovertemplate="site %{text}<extra></extra>"))
    site_trace = len(fig.data) - 1

    buttons = []
    for li, name in enumerate(menu):
        visible = [False] * len(fig.data)
        for i in vis_all[li]:
            visible[i] = True
        visible[site_trace] = True
        buttons.append(dict(label=name, method="update",
                            args=[{"visible": visible}]))
    lat = float(sites["lat"].mean()) if len(sites) else 0.0
    layout = plotly_layout(height=560, margin=dict(l=60, r=24, t=72, b=48))
    if len(buttons) > 1:
        layout["updatemenus"] = [dict(
            buttons=buttons, direction="down", x=0, xanchor="left", y=1.12,
            yanchor="top", bgcolor="white", bordercolor=AXIS,
            font=dict(color=TEXT, size=11), showactive=True)]
    else:
        layout["title"] = dict(text=menu[0], font=dict(size=14, color=TEXT),
                               x=0, xanchor="left")
    fig.update_layout(**layout)
    fig.update_xaxes(title="Longitude", showgrid=False)
    fig.update_yaxes(title="Latitude", showgrid=False, scaleanchor="x",
                     scaleratio=1 / max(np.cos(np.radians(lat)), 0.1))
    return fig


def _log_ticks(lo, hi):
    """Colour bar ticks (1, 2, 5 x powers of 10) from lo to hi."""
    ticks = [m * 10 ** e for e in range(int(np.floor(np.log10(lo))),
                                        int(np.ceil(np.log10(hi))) + 1)
             for m in (1, 2, 5) if lo <= m * 10 ** e <= hi]
    if len(ticks) < 2:
        ticks = [float(f"{t:.2g}") for t in np.geomspace(lo, hi, 3)]
    return ticks


def fig_response(run):
    """
    How the best model's resistance responds to each of its variables (as
    ResistanceGA's Plot.trans): the resistance per km due to that variable
    (1 + (K - 1) x its weight share x its 0-1 value; the same for the
    standard model and SAMC), against its raw values over its range in the
    network, from the model's specification (modelSpec.json); the
    network's reaches are marked on each curve.
    """
    path = f"{run.prefix}.modelSpec.json"
    if not os.path.exists(path):
        return None
    from resistnet.predict import load_spec, variable_values
    model = load_spec(path)["models"][0]
    variables = model.get("variables", {})
    if not variables:
        return None
    go = _go()
    from plotly.subplots import make_subplots
    names = list(variables)
    ncol = min(3, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=names,
                        horizontal_spacing=0.1, vertical_spacing=0.22)
    total = sum(var["weight"] for var in variables.values())
    per_km = "max_resistance" in model and total > 0
    label = "Resistance per km" if per_km else "Transformed value (0-1)"
    colour = COLOURS[run.model]

    def curve(x, var):
        f = variable_values(x, var)
        if per_km:
            share = var["weight"] / total
            return 1 + (model["max_resistance"] - 1) * share * f
        return f

    for i, (v, var) in enumerate(variables.items()):
        row, col = i // ncol + 1, i % ncol + 1
        x = np.linspace(var["min"], var["max"], 201)
        line = go.Scatter(
            x=x, y=curve(x, var), mode="lines", showlegend=False,
            line=dict(color=colour, width=2.5),
            hovertemplate=(v + " %{x:.4g}<br>" + label
                           + " %{y:.4g}<extra></extra>"))
        fig.add_trace(line, row=row, col=col)
        if v in run.reach_values:
            r = run.reach_values[v].dropna()
            dots = go.Scatter(
                x=r.to_numpy(), y=curve(r.to_numpy(), var), mode="markers",
                showlegend=False, customdata=r.index.astype(str),
                marker=dict(size=7, color="white",
                            line=dict(color=colour, width=1.2)),
                hovertemplate=("reach %{customdata}<br>" + v
                               + " %{x:.4g}<br>" + label
                               + " %{y:.4g}<extra></extra>"))
            fig.add_trace(dots, row=row, col=col)
        fig.update_xaxes(title_text=f"{v} (raw value)", row=row, col=col)
    fig.update_layout(**plotly_layout(height=300 * nrow + 60))
    fig.update_xaxes(gridcolor=GRID, linecolor=AXIS)
    fig.update_yaxes(gridcolor=GRID, linecolor=AXIS)
    fig.update_yaxes(title_text=label, col=1)
    for a in fig.layout.annotations:
        a.update(font=dict(color=TEXT, size=12))
    return fig


def _line_coords(arrays):
    """Concatenates lines for one Plotly trace (None between lines)."""
    xs, ys = [], []
    for xy in arrays:
        xs.extend(xy[:, 0].tolist() + [None])
        ys.extend(xy[:, 1].tolist() + [None])
    return xs, ys


def fig_flow_profile(run):
    """
    AIC of flow alone at each rate, both directions (signed log scale):
    the evidence for gene flow biased downstream or upstream.
    """
    go = _go()
    t = run.flow_profile.sort_values("rate")
    rate = t["rate"].to_numpy()
    low = max(np.abs(rate[rate != 0]).min(), 1e-6) if (rate != 0).any() \
        else 1e-3

    def position(r):
        r = np.asarray(r, dtype=float)
        return np.where(r == 0, 0.0, np.sign(r) * (
            np.log10(np.maximum(np.abs(r), low) / low) + 1))

    x = position(rate)
    best = t["delta_aic"].to_numpy().argmin()
    zero = t[t["rate"] == 0]
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x, y=t["delta_aic"], mode="lines+markers", name="Flow alone",
        line=dict(color=COLOURS["samc"], width=2),
        marker=dict(size=7, color=COLOURS["samc"],
                    line=dict(color="white", width=1)),
        customdata=np.stack([rate, t["aic"]], axis=1),
        hovertemplate=("rate %{customdata[0]:+.4f}/km<br>AIC "
                       "%{customdata[1]:.2f} (Δ %{y:.2f})<extra></extra>")))
    fig.add_trace(go.Scatter(
        x=[x[best]], y=[0], mode="markers", name="Best rate",
        marker=dict(size=13, color=COLOURS["samc"],
                    line=dict(color=TEXT, width=1.5)),
        hovertemplate=f"best rate {rate[best]:+.4f}/km<extra></extra>"))
    if len(zero):
        fig.add_hline(y=float(zero["delta_aic"].iloc[0]),
                      line=dict(color=MUTED, width=1, dash="dot"),
                      annotation_text="distance alone (no flow)",
                      annotation_font_color=TEXT2,
                      annotation_position="top left")
    ticks = [r for r in (-1, -0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1)
             if r == 0 or (low <= abs(r) <= np.abs(rate).max() * 1.0001)]
    top = float(t["delta_aic"].max())
    for xpos, text in ((x.min(), "← gene flow downstream"),
                       (x.max(), "gene flow upstream →")):
        fig.add_annotation(x=xpos, y=top, text=text, showarrow=False,
                           xanchor="left" if xpos < 0 else "right",
                           yanchor="bottom", font=dict(color=TEXT2,
                                                       size=11))
    fig.update_layout(**plotly_layout(
        title="Flow alone: AIC difference from the best rate", height=380,
        margin=dict(l=60, r=24, t=72, b=48)))
    fig.update_xaxes(title="Flow rate (per km; signed log scale)",
                     tickvals=position(ticks).tolist(),
                     ticktext=[f"{r:g}" for r in ticks])
    fig.update_yaxes(title="ΔAIC (lower is better)", rangemode="tozero")
    return fig


def fig_ranking(table, top=20):
    """The joint ranking of both families and the reference models."""
    go = _go()
    t = table.head(top)
    fig = go.Figure()
    for family in ("resistance", "samc", "reference"):
        sel = t[t["family"] == family]
        if sel.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sel["delta_aic"], y=sel["model"], mode="markers",
            name=FAMILY_NAMES[family],
            marker=dict(size=12, color=COLOURS[family],
                        line=dict(color="white", width=1.5)),
            customdata=np.stack([sel["aic"], sel["akaike_weight"],
                                 sel["k"], sel["rank"]], axis=1),
            hovertemplate=("<b>%{y}</b><br>rank %{customdata[3]} · ΔAIC "
                           "%{x:.2f} (AIC %{customdata[0]:.2f}, k "
                           "%{customdata[2]})<br>Akaike weight "
                           "%{customdata[1]:.3f}<extra></extra>")))
    fig.add_vline(x=2, line=dict(color=TEXT2, width=1, dash="dot"),
                  annotation_text="ΔAIC 2", annotation_font_color=TEXT2,
                  annotation_position="bottom right")
    fig.update_layout(**plotly_layout(
        title="AIC difference from the best model (best at the top)",
        height=max(340, 32 * len(t) + 140)))
    fig.update_xaxes(title="ΔAIC", rangemode="tozero")
    fig.update_yaxes(categoryorder="array",
                     categoryarray=list(t["model"])[::-1], automargin=True,
                     tickfont=dict(size=11, color=TEXT))
    return fig


def fig_joint_importance(importance, families):
    go = _go()
    imp = importance.iloc[::-1]
    fig = go.Figure()
    for family in families:
        column = f"RVI_{family}"
        if column not in imp:
            continue
        fig.add_trace(go.Bar(
            x=imp[column], y=imp["variable"], orientation="h",
            name=FAMILY_NAMES[family],
            marker=dict(color=COLOURS[family],
                        line=dict(color="white", width=1.5)),
            hovertemplate="%{y}: %{x:.3f} from " + FAMILY_NAMES[family]
                          + " models<extra></extra>"))
    fig.update_layout(**plotly_layout(
        barmode="stack", height=max(240, 44 * len(imp) + 120),
        title="Variable importance over both families"))
    fig.update_xaxes(title="Relative importance (summed Akaike weight of "
                     "the models with the variable)", range=[0, 1.05])
    fig.update_yaxes(gridcolor="rgba(0,0,0,0)")
    return fig


def fig_best_fits(runs):
    """
    Genetic distance against river distance (isolation by distance) and
    against each family's best model, side by side.
    """
    go = _go()
    from plotly.subplots import make_subplots
    panels = []
    first = next(iter(runs.values()))
    if first.pairwise is not None and "river_km" in first.pairwise:
        ref = first.null.loc["distance_only"] if first.null is not None \
            else None
        title = "Isolation by distance (river km)"
        if ref is not None:
            title += (f"<br><sup>AIC {float(ref['aic']):.1f} · marginal "
                      f"R² {float(ref['r2m']):.2f}</sup>")
        panels.append(("reference", first.pairwise, "river_km", title))
    for family, run in runs.items():
        if run.pairwise is None or "best_model_distance" not in run.pairwise:
            continue
        best = run.models.iloc[0]
        panels.append((family, run.pairwise, "best_model_distance",
                       f"{FAMILY_NAMES[family]} best model<br><sup>AIC "
                       f"{float(best['aic']):.1f} · marginal R² "
                       f"{float(best['r2m']):.2f}</sup>"))
    fig = make_subplots(rows=1, cols=len(panels), shared_yaxes=True,
                        horizontal_spacing=0.05,
                        subplot_titles=[p[3] for p in panels])
    for i, (family, pw, column, _) in enumerate(panels, start=1):
        raw = pw[column].to_numpy()
        x = (raw - raw.mean()) / raw.std()
        y = pw["genetic_distance"].to_numpy()
        unit = "km" if column == "river_km" else "model distance"
        points = go.Scatter(
            x=x, y=y, mode="markers", name=FAMILY_NAMES[family],
            marker=dict(size=8, color=COLOURS[family],
                        line=dict(color="white", width=1)),
            customdata=np.stack([pw["site1"], pw["site2"], raw], axis=1),
            hovertemplate=("%{customdata[0]} – %{customdata[1]}<br>"
                           + unit + " %{customdata[2]:.3g}<br>genetic "
                           "distance %{y:.4f}<extra></extra>"))
        fig.add_trace(points, row=1, col=i)
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.array([x.min(), x.max()])
        fig.add_trace(go.Scatter(x=xs, y=intercept + slope * xs,
                                 mode="lines", hoverinfo="skip",
                                 line=dict(color=TEXT2, width=1.5),
                                 showlegend=False), row=1, col=i)
        fig.update_xaxes(title=("River distance (standardised)"
                                if column == "river_km"
                                else "Best model's distance (standardised)"),
                         row=1, col=i)
    fig.update_yaxes(title="Genetic distance", row=1, col=1)
    fig.update_layout(**plotly_layout(height=440, showlegend=False,
                                      margin=dict(l=60, r=24, t=80, b=48)))
    for a in fig.layout.annotations:
        a.update(font=dict(color=TEXT, size=12))
    return fig


# ------------------------------------------------------------------- report
CSS = """
:root { --ink: #0b0b0b; --ink2: #52514e; --muted: #898781; --rule: #e1e0d9;
  --axis: #c3c2b7; --bg: #f9f9f7; --surface: #ffffff; --accent: #2a78d6;
  --samc: #eb6834; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
  font-family: "IBM Plex Sans", "Helvetica Neue", Arial, sans-serif;
  font-size: 15px; line-height: 1.55; }
header { background: var(--surface); border-bottom: 1px solid var(--rule);
  padding: 22px 24px 14px; }
header h1 { margin: 0 0 4px; font-size: 1.55rem; font-weight: 600; }
header .meta { color: var(--ink2); font-size: 0.9rem; }
nav { position: sticky; top: 0; z-index: 10; background: var(--surface);
  border-bottom: 1px solid var(--rule); padding: 8px 24px; display: flex;
  gap: 18px; flex-wrap: wrap; font-size: 0.9rem; }
nav a { color: var(--ink2); text-decoration: none; }
nav a:hover { color: var(--accent); }
main { max-width: 1180px; margin: 0 auto; padding: 8px 16px 64px; }
section { margin-top: 30px; }
section > h2 { font-size: 1.3rem; font-weight: 600; margin: 0 0 6px;
  padding-bottom: 6px; border-bottom: 2px solid var(--rule); }
section > h2 .tag { display: inline-block; width: 10px; height: 10px;
  border-radius: 50%; margin-right: 8px; }
h3 { font-size: 1.02rem; font-weight: 600; margin: 22px 0 6px; }
p { margin: 6px 0; color: var(--ink2); max-width: 78ch; }
.cards { display: grid; gap: 10px; margin: 12px 0;
  grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); }
.card { background: var(--surface); border: 1px solid var(--rule);
  border-radius: 8px; padding: 10px 14px; }
.card .label { color: var(--muted); font-size: 0.75rem; letter-spacing: .04em;
  text-transform: uppercase; }
.card .value { font-size: 1.15rem; font-weight: 600; margin-top: 2px; }
.card .sub { color: var(--ink2); font-size: 0.82rem; }
.panel { background: var(--surface); border: 1px solid var(--rule);
  border-radius: 8px; padding: 8px; margin: 10px 0; }
.tablewrap { overflow-x: auto; background: var(--surface);
  border: 1px solid var(--rule); border-radius: 8px; margin: 10px 0; }
table { border-collapse: collapse; width: 100%; font-size: 0.85rem; }
th, td { padding: 6px 10px; text-align: left; border-bottom: 1px solid
  var(--rule); vertical-align: top; }
th { color: var(--ink2); font-weight: 600; background: #f3f2ee; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums;
  white-space: nowrap; }
tr:last-child td { border-bottom: 0; }
footer { color: var(--muted); font-size: 0.8rem; margin-top: 40px; }
"""


def _fig_html(fig):
    return fig.to_html(full_html=False, include_plotlyjs=False, config={
        "displaylogo": False, "responsive": True,
        "toImageButtonOptions": {"format": "svg"}})


def _panel(fig):
    return f'<div class="panel">{_fig_html(fig)}</div>'


def _fmt(v, d=2):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "–"
    if isinstance(v, (int, np.integer)):
        return str(v)
    if isinstance(v, (float, np.floating)):
        if v != 0 and abs(v) < 0.001:
            return f"{v:.1e}"
        return f"{v:,.{d}f}"
    return html.escape(str(v))


def _table(df, columns, numeric=()):
    """An HTML table of the given columns ({column: header})."""
    head = "".join(
        f'<th class="{"num" if c in numeric else ""}">{html.escape(h)}</th>'
        for c, h in columns.items())
    rows = []
    for _, r in df.iterrows():
        cells = "".join(
            f'<td class="{"num" if c in numeric else ""}">{_fmt(r.get(c))}'
            "</td>" for c in columns)
        rows.append(f"<tr>{cells}</tr>")
    return (f'<div class="tablewrap"><table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def _cards(items):
    return '<div class="cards">' + "".join(
        f'<div class="card"><div class="label">{html.escape(a)}</div>'
        f'<div class="value">{b}</div><div class="sub">{c}</div></div>'
        for a, b, c in items) + "</div>"


HOW_IT_WORKS = {
    "km": (
        "How this model works: each selected variable is transformed and "
        "rescaled to 0-1, and their weighted mean gives each reach a "
        "resistance per km of 1 + (K - 1) × mean, from 1 (distance alone) "
        "to the maximum resistance K, which is fitted. A reach costs its "
        "resistance per km times its length, and the distance between two "
        "sites is the sum along the river between them. This is "
        "ResistanceGA's formulation (resistance 1 + max × f(x), times the "
        "distance between cells): distance is always part of the model, "
        "K = 1 is distance alone, and a large K means the variables "
        "dominate."),
    "samc": (
        "How this model works: samc's random walker on a lattice of the "
        "river (points about 1 km apart along each reach, and the nodes "
        "where reaches meet and where the sites are). It moves between "
        "neighbouring points with samc's conductance, 1 / (resistance per "
        "km × distance), the resistance per km being the standard "
        "model's (1 to a fitted maximum K with the weighted variables), "
        "times exp(-rate × km moved upstream) with a flow term (by "
        "default nothing absorbs it). The distance between two sites "
        "combines the passage times of the two directions (by default, "
        "the larger). In simulations, the walker moved like a lineage "
        "traced back in time: a negative flow rate meant gene flow biased "
        "downstream."),
}


def _distance_share(run):
    """Distance's share of the best model's resistance per km (from
    bestModel.tsv), or None."""
    b = run.best
    if (not run.per_km or b is None or "Contribution" not in b
            or "distance" not in set(b["Variable"])):
        return None
    return float(b.set_index("Variable").loc["distance", "Contribution"])


def _kind(run):
    """The model's formulation: "samc" or "km"."""
    return "samc" if run.samc else "km"


def _model_section(run):
    """The section of one model's run."""
    hof = run.models
    best = hof.iloc[0]
    kept = int(hof["kept"].sum())
    cards = [("Best model", html.escape(best["model"]), ""),
             ("AIC", _fmt(float(best["aic"])),
              f"marginal R² {_fmt(float(best['r2m']), 3)}")]
    if "delta_aic_distance" in hof:
        cards.append(("vs distance alone",
                      _fmt(float(best["delta_aic_distance"]), 1),
                      "AIC better (positive) or worse"))
    share = _distance_share(run)
    if share is not None:
        cards.append(("Distance's share", f"{share:.0%}",
                      "of the best model's resistance per km (the rest "
                      "from its variables)"))
    if run.samc and "delta_aic_flow" in hof:
        cards.append(("vs flow alone", _fmt(float(best["delta_aic_flow"]),
                                            1),
                      "AIC better (positive) or worse"))
    cards.append(("Hall of Fame", f"{len(hof)} models",
                  f"{kept} kept for model averaging; from "
                  f"{int(hof['chains'].sum())} chains"))
    parts = [f"<p>{HOW_IT_WORKS[_kind(run)]}</p>", _cards(cards)]

    if run.null is not None:
        ref = run.null.copy()
        ref["name"] = [REFERENCE_NAMES.get(i, i) for i in ref.index]
        ref["delta"] = ref["aic"] - float(best["aic"])
        parts.append("<h3>Reference models</h3><p>Fitted to the same genetic "
                     "distances. ΔAIC is the reference's AIC minus the best "
                     "model's: positive means the best model fits better."
                     "</p>")
        cols = {"name": "Model", "aic": "AIC", "delta": "ΔAIC",
                "r2m": "Marginal R²", "loglik": "Log-likelihood"}
        if "flow" in ref:
            cols["flow"] = "Flow rate (/km)"
        parts.append(_table(ref, cols, numeric=set(cols) - {"name"}))

    if run.samc and run.flow_profile is not None:
        parts.append("<h3>Flow: gene flow biased downstream or "
                     "upstream?</h3><p>"
                     "The flow-alone reference model at each rate tried "
                     "when fitting it: negative rates make the walker's "
                     "upstream moves easier, positive rates harder, and "
                     "rate 0 is distance alone. The walker moves like a "
                     "lineage traced back in time, against gene flow, so a "
                     "negative rate means gene flow biased downstream (in "
                     "simulations, with no absorption and the larger of the "
                     "two passage times). Gradients along the river (e.g. "
                     "temperature) and barriers the model does not "
                     "represent can mimic or hide a direction, and gains of "
                     "a few AIC units are weak evidence.</p>")
        parts.append(_panel(fig_flow_profile(run)))
    parts.append("<h3>Hall of Fame</h3>")
    parts.append(_panel(fig_aic_profile(run)))
    cols = {"model": "Model", "aic": "AIC", "delta": "ΔAIC",
            "akaike_weight": "Akaike weight", "r2m": "Marginal R²",
            "chains": "Chains"}
    if "delta_aic_distance" in hof:
        cols["delta_aic_distance"] = "vs distance alone"
    if run.samc and "delta_aic_flow" in hof:
        cols["delta_aic_flow"] = "vs flow alone"
    hof = hof.copy()
    hof["delta"] = hof["aic"] - hof["aic"].min()
    parts.append(_table(hof.head(15), cols, numeric=set(cols) - {"model"}))

    if run.rvi is not None or run.maw is not None:
        parts.append("<h3>Variables</h3><p>Relative importance is the summed "
                     "Akaike weight of the models with the variable; "
                     "model-averaged weights are averaged over those models."
                     "</p>")
        parts.append(_panel(fig_importance(run)))
    shares = fig_contributions(run)
    if shares is not None:
        parts.append("<p>Distance is part of every model per km, so it "
                     "has no relative importance of its own. Instead, this "
                     "is how much of the best model's resistance per km "
                     "comes from distance (the baseline of 1) and from "
                     "each variable, averaged over the river's km, as "
                     "ResistanceGA's p.contribution.</p>")
        parts.append(_panel(shares))
    response = fig_response(run)
    if response is not None:
        parts.append(
            "<h3>The best model's response to each variable</h3><p>"
            + ("The resistance per km that each variable adds over its range "
               "in this network (1 = distance alone; the variables' "
               "contributions add up), as ResistanceGA's Plot.trans. "
               if run.per_km else
               "Each variable's transformed value (0-1) over its range in "
               "this network, as ResistanceGA's Plot.trans. ")
            + "Dots are this network's reaches.</p>")
        parts.append(_panel(response))
    if run.pairwise is not None:
        parts.append("<h3>Fit to the genetic distances</h3><p>Each point is "
                     "a pair of sites. Hover for the pair.</p>")
        parts.append(_panel(fig_pairwise(run)))
    if run.geometry and any(e is not None for e in run.edges.values()):
        layers = {}
        for which in ("Average", "Best"):
            e = run.edges[which]
            if e is not None and "Resistance" in e:
                label = "Model average" if which == "Average" else \
                    "Best model"
                layers[f"{label}: resistance"] = \
                    e.set_index("EDGE_ID")["Resistance"]
        fig = fig_map(run, layers, "Resistance per km", log=True)
        note = ("Resistance per km of each reach (1 = distance alone, up "
                "to the model's maximum K; colours on a log scale; choose "
                "the model average or the best model in the menu).")
        if run.samc:
            note += (" The flow term, if any, acts on top of it, the same "
                     "per km everywhere.")
        parts.append(f"<h3>Map</h3><p>{note} Hover over a reach for its "
                     "value; scroll to zoom.</p>")
        parts.append(_panel(fig))
    return "\n".join(parts)


def _comparison_section(runs, table, summary, importance):
    best = summary.set_index("family")
    cards = []
    for family in summary["family"]:
        r = best.loc[family]
        cards.append((FAMILY_NAMES[family], html.escape(r["best_model"]),
                      f"AIC {_fmt(r['best_aic'])} · ΔAIC "
                      f"{_fmt(r['delta_aic'], 1)} · the family's summed "
                      f"Akaike weight {_fmt(r['akaike_weight'], 3)}"))
    parts = [
        "<p>Both models were fitted with the same MLPE to the same genetic "
        "distances, so their AICs are comparable. Akaike weights are "
        "computed over the models of both and the reference models. AIC "
        "ranks fits, not mechanisms: path resistance and a walker's "
        "passage time are different ideas.</p>",
        _cards(cards), _panel(fig_ranking(table))]
    if any(r.pairwise is not None for r in runs.values()):
        parts.append("<h3>Each family's best model, and isolation by "
                     "distance</h3><p>Genetic distance against river "
                     "distance (the distance-only reference model) and "
                     "against each family's best model. Distances are "
                     "standardised, as fitted; lines are least-squares fits "
                     "for display.</p>")
        parts.append(_panel(fig_best_fits(runs)))
    if importance is not None and len(importance):
        parts.append("<h3>Variables</h3>")
        parts.append(_panel(fig_joint_importance(importance, list(runs))))
    parts.append("<h3>Ranking</h3>")
    parts.append(_table(table.head(15), {
        "rank": "Rank", "family": "Family", "model": "Model", "k": "k",
        "aic": "AIC", "delta_aic": "ΔAIC", "akaike_weight": "Akaike weight"},
        numeric={"rank", "k", "aic", "delta_aic", "akaike_weight"}))
    return "\n".join(parts)


def build_report(out, verbose=True):
    """
    Builds DIR/NAME/NAME.report.html from the outputs of
    runResistnet.py -o DIR/NAME (one or both models).

    Args:
        out (str): The run's -o value (DIR/NAME).
        verbose (bool): Print the report's path.

    Returns:
        str: Path of the report.
    """
    from plotly.offline import get_plotlyjs
    name = os.path.basename(os.path.normpath(out))
    runs = {}
    for model in ("resistance", "samc"):
        prefix = os.path.join(out, model, name)
        if os.path.exists(f"{prefix}.HallOfFame.tsv"):
            runs[model] = RunOutputs(prefix, model)
    if not runs:
        raise FileNotFoundError(f"No ResistNet outputs found in {out}")

    sections = []
    table = _read(os.path.join(out, f"{name}.modelComparison.tsv"))
    summary = _read(os.path.join(out, f"{name}.familySummary.tsv"))
    if len(runs) > 1 and table is not None and summary is not None:
        importance = _read(os.path.join(out, f"{name}.jointImportance.tsv"))
        sections.append(("comparison", "Model comparison", MUTED,
                         _comparison_section(runs, table, summary,
                                             importance)))
    for model, run in runs.items():
        sections.append((model, FAMILY_NAMES[model], COLOURS[model],
                         _model_section(run)))

    first = next(iter(runs.values()))
    n_sites = len(first.points) if first.points is not None else "?"
    nav = "".join(f'<a href="#{sid}">{html.escape(title)}</a>'
                  for sid, title, _, _ in sections)
    body = "".join(
        f'<section id="{sid}"><h2><span class="tag" style="background:'
        f'{colour}"></span>{html.escape(title)}</h2>{content}</section>'
        for sid, title, colour, content in sections)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ResistNet report: {html.escape(name)}</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script>
</head><body>
<header><h1>ResistNet report: {html.escape(name)}</h1>
<div class="meta">{len(runs)} model{"s" if len(runs) > 1 else ""}
({", ".join(FAMILY_NAMES[m] for m in runs)}) · {n_sites} sites ·
{len(first.variables)} variables ({", ".join(first.variables)}) ·
{datetime.now().strftime("%Y-%m-%d %H:%M")}</div></header>
<nav>{nav}</nav>
<main>{body}
<footer>Built from the run's output files in {html.escape(out)}; each
figure's data is in those tables (e.g. *.pairwise.tsv, *.HallOfFame.tsv).
Hover for values; drag to zoom; double-click to reset; the camera icon saves
a figure as SVG.</footer></main></body></html>
"""
    path = os.path.join(out, f"{name}.report.html")
    with open(path, "w") as fh:
        fh.write(page)
    if verbose:
        print(f"\nReport: {path}")
    return path
