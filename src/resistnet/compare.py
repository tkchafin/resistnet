"""
Compares the standard (resistance) and SAMC models fitted to the same
genetic distances (runResistnet.py --model both): one ranking of the models
of both, and of the reference models, by AIC with joint Akaike weights.

Both models are fitted with the same MLPE to the same genetic distances, so
their AICs are comparable; the reference models computed by both (null and
distance_only) must therefore be the same, which is checked.
"""
import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from resistnet.plotting import COLOURS, FAMILY_NAMES, GRID, TEXT, TEXT2
from resistnet.samc_network import ResistanceNetworkSAMC

TRANSFORMATIONS = {
    0: "linear", 1: "Ricker", 2: "reverse Ricker", 3: "inverse Ricker",
    4: "reverse-inverse Ricker", 5: "monomolecular",
    6: "reverse monomolecular", 7: "inverse monomolecular",
    8: "reverse-inverse monomolecular"
}

# Reference models: (family, description, number of parameters)
REFERENCES = {
    "null": ("reference", "population effects only (null)", 3),
    "distance_only": ("reference", "distance only (km; standard model)", 4),
    "samc_distance": ("reference", "SAMC distance alone", 4),
    "samc_flow": ("reference", "SAMC flow alone", 5),
}


def read_null_table(prefix):
    """Reads {prefix}.Null-Model.tsv (keeping its "null" row)."""
    return pd.read_csv(f"{prefix}.Null-Model.tsv", sep="\t", index_col=0,
                       keep_default_na=False, na_values=[""])


def describe(row, variables, network, runner):
    """
    Describes a Hall of Fame model: its variables with their transformation
    and, for SAMC, their direction in raw values and the flow rate; for the
    standard model per km, the range of its resistance per km (1 to K).

    Args:
        row (pd.Series): Hall of Fame row (as read from HallOfFame.tsv).
        variables (list): The variables, in the network's order.
        network: The network the model was fitted to.
        runner: The ModelRunner that fitted it.

    Returns:
        str: The description.
    """
    samc = isinstance(network, ResistanceNetworkSAMC)
    directions = (network.effective_directions(runner._row_individual(row))
                  if samc else {})
    parts = []
    for v in variables:
        if row[v] != 1:
            continue
        t = int(row[f"{v}_trans"])
        label = TRANSFORMATIONS[t]
        if t:
            label += f", shape {float(row[f'{v}_shape']):g}"
        if samc:
            if int(row[f"{v}_asym"]) == 1:
                share = directions.get(v)
                label += (
                    ", directional" if share is None else
                    ", directional: towards higher values harder"
                    if share >= 0.95 else
                    ", directional: towards lower values harder"
                    if share <= 0.05 else
                    ", directional: mixed"
                )
            else:
                label += ", symmetric"
        parts.append(f"{v} ({label})")
    flow = float(row["flow"]) if samc and "flow" in row else 0.0
    if flow:
        # The walker moves like a lineage, against gene flow
        parts.append(f"flow {flow:+.3f}/km (gene flow "
                     f"{'upstream' if flow > 0 else 'downstream'})")
    text = " + ".join(parts)
    # Standard model per km: the range of its resistance per km
    if not samc and parts and "max_resistance" in row:
        text += f"; resistance 1-{float(row['max_resistance']):.3g} per km"
    return text


def duplicate_key(model, aic):
    """
    Models that are the same: the same description apart from the maximum
    resistance per km (which chains converging on the same optimum find
    with small differences), and the same AIC to 0.01.
    """
    return model.split("; resistance 1-")[0], round(float(aic), 2)


def comparison_tables(runs, tol=1e-6):
    """
    Ranks the models of several fitted model families together.

    Args:
        runs (dict): For each family ("resistance", "samc"), a tuple of the
                     network, the ModelRunner, and the output prefix of its
                     run (with HallOfFame.tsv and Null-Model.tsv written).
        tol (float): Largest allowed difference between the families'
                     shared reference models (AIC).

    Returns:
        tuple: The ranking, the summary per family, and the joint relative
               importance of each variable (DataFrames).

    Raises:
        ValueError: If the families were not fitted to the same genetic
                    distances (their shared reference models differ).
    """
    rows, nulls, selected, owners = [], {}, [], []
    for family, (network, runner, prefix) in runs.items():
        hof = pd.read_csv(f"{prefix}.HallOfFame.tsv", sep="\t",
                          na_values="-")
        variables = list(network._predictors.columns)
        seen = set()
        for _, row in hof.iterrows():
            model = describe(row, variables, network, runner)
            # Each chain's best model is kept, so a model can appear more
            # than once (for SAMC without variables, with different
            # maximum resistances that have no effect; per km, with
            # maximum resistances differing only by search noise): list it
            # once
            key = duplicate_key(model, row["aic"])
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "family": family, "model": model,
                "k": int(round((float(row["aic"]) + 2 * float(row["loglik"]))
                               / 2)),
                "loglik": float(row["loglik"]), "aic": float(row["aic"])})
            selected.append({v for v in variables if row[v] == 1})
            owners.append(family)
        nulls[family] = read_null_table(prefix)

    # The families share the null and distance-only models
    first = next(iter(nulls.values()))
    for family, table in nulls.items():
        for ref in ("null", "distance_only"):
            if abs(float(table.loc[ref, "aic"])
                   - float(first.loc[ref, "aic"])) > tol:
                raise ValueError(
                    f"The {ref} models of the fitted families differ, so "
                    "they were not fitted to the same genetic distances")

    # Reference models, each once
    refs = {}
    for table in nulls.values():
        for ref in table.index:
            if ref in REFERENCES and ref not in refs:
                refs[ref] = table.loc[ref]
    for ref, r in refs.items():
        family, model, k = REFERENCES[ref]
        if ref == "samc_flow":
            model += f" ({float(r['flow']):+.3f}/km)"
        rows.append({"family": family, "model": model, "k": k,
                     "loglik": float(r["loglik"]), "aic": float(r["aic"])})
        selected.append(set())
        owners.append(family)

    table = pd.DataFrame(rows)
    table["delta_aic"] = table["aic"] - table["aic"].min()
    w = np.exp(-table["delta_aic"] / 2)
    table["akaike_weight"] = w / w.sum()
    order = np.argsort(table["aic"].to_numpy(), kind="stable")
    table = table.iloc[order].reset_index(drop=True)
    selected = [selected[i] for i in order]
    owners = [owners[i] for i in order]
    table.insert(0, "rank", np.arange(1, len(table) + 1))

    summary = []
    for family in list(runs) + ["reference"]:
        part = table[table["family"] == family]
        if part.empty:
            continue
        best = part.iloc[0]
        summary.append({
            "family": family, "best_model": best["model"],
            "best_aic": best["aic"], "delta_aic": best["delta_aic"],
            "akaike_weight": part["akaike_weight"].sum()})
    summary = pd.DataFrame(summary)

    # Relative importance over both families, and the part of it from
    # each family's models
    variables = sorted(set().union(*selected))
    importance = pd.DataFrame({"variable": variables})
    weights = table["akaike_weight"].to_numpy()
    importance["RVI"] = [sum(w for w, s in zip(weights, selected) if v in s)
                         for v in variables]
    for family in runs:
        importance[f"RVI_{family}"] = [
            sum(w for w, s, f in zip(weights, selected, owners)
                if v in s and f == family) for v in variables]
    importance = importance.sort_values("RVI", ascending=False,
                                        kind="stable")
    return table, summary, importance


def write_comparison(runs, prefix, verbose=True, plot=True):
    """
    Writes the joint ranking of the fitted model families (see
    comparison_tables()) to {prefix}.modelComparison.tsv, the summary per
    family to {prefix}.familySummary.tsv, and the joint relative variable
    importance to {prefix}.jointImportance.tsv; and plots them
    ({prefix}.modelComparison.pdf, {prefix}.jointImportance.pdf), with each
    family's best model against the genetic distances
    ({prefix}.bestFits.pdf).

    Args:
        runs (dict): As for comparison_tables().
        prefix (str): Output prefix.
        verbose (bool): Print the summary.
        plot (bool): Write the plots.

    Returns:
        tuple: The three tables.
    """
    table, summary, importance = comparison_tables(runs)
    table.to_csv(f"{prefix}.modelComparison.tsv", sep="\t", index=False)
    summary.to_csv(f"{prefix}.familySummary.tsv", sep="\t", index=False)
    importance.to_csv(f"{prefix}.jointImportance.tsv", sep="\t",
                      index=False)
    if plot:
        plot_ranking(table, summary, f"{prefix}.modelComparison.pdf")
        plot_importance(importance, list(runs),
                        f"{prefix}.jointImportance.pdf")
        plot_best_fits(runs, summary, f"{prefix}.bestFits.pdf")
    if verbose:
        print("\nModels of both families and the reference models, ranked "
              "together by AIC:")
        with pd.option_context("display.max_rows", 20,
                               "display.max_colwidth", 80,
                               "display.width", 200):
            print(table.head(10).to_string(index=False))
            print("\nBy family:")
            print(summary.to_string(index=False))
    return table, summary, importance


def _style(ax, grid="x"):
    """Quiet axes: no box, light grid, secondary-ink ticks."""
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=TEXT2, labelsize=8, length=0)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)


def _legend(ax, families, **kw):
    handles = [Line2D([], [], marker="o", linestyle="", markersize=7,
                      markerfacecolor=COLOURS[f], markeredgecolor="white",
                      label=FAMILY_NAMES[f]) for f in families]
    leg = ax.legend(handles=handles, frameon=False, fontsize=8, **kw)
    for text in leg.get_texts():
        text.set_color(TEXT2)
    return leg


def plot_ranking(table, summary, path, top=15):
    """
    Plots the best models of the joint ranking as their AIC difference
    from the best (symmetric log scale), coloured by family, with each
    family's summed Akaike weight alongside.
    """
    rows = table.head(top).iloc[::-1].reset_index(drop=True)
    labels = [textwrap.fill(m, 62) for m in rows["model"]]
    lines = np.array([lab.count("\n") + 1 for lab in labels])
    gaps = 0.55 + 0.32 * lines
    y = np.cumsum(gaps) - gaps / 2
    height = max(2.6, 0.28 * gaps.sum() + 1.4)
    fig, (ax, axw) = plt.subplots(
        1, 2, figsize=(11.5, height),
        gridspec_kw={"width_ratios": [4.2, 1.1], "wspace": 0.55})

    delta = rows["delta_aic"].to_numpy()
    colours = [COLOURS[f] for f in rows["family"]]
    ax.hlines(y, 0, delta, color=GRID, linewidth=1.2, zorder=1)
    ax.scatter(delta, y, s=48, c=colours, edgecolors="white",
               linewidths=1.5, zorder=3)
    ax.set_xscale("symlog", linthresh=2, linscale=0.6)
    top_x = max(10.0, delta.max() * 1.6)
    ticks = [t for t in (0, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
             if t <= top_x]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.set_xlim(-0.3, top_x)
    ax.axvline(2, color=TEXT2, linewidth=0.8, linestyle=(0, (3, 3)),
               zorder=2)
    ax.text(2, 0.05, " ΔAIC = 2", color=TEXT2, fontsize=7, va="bottom")
    for yi, d, w in zip(y, delta, rows["akaike_weight"]):
        if w >= 0.005:
            ax.text(d, yi, f"   w = {w:.2f}", color=TEXT2, fontsize=7,
                    va="center")
    ax.set_yticks(y, labels, fontsize=7.5, color=TEXT)
    ax.set_ylim(0, gaps.sum())
    ax.set_xlabel("AIC difference from the best model (symmetric log "
                  "scale)", color=TEXT2, fontsize=8)
    ax.set_title("Models of both families and the reference models, "
                 "by AIC", loc="left", fontsize=10, color=TEXT, pad=20)
    _style(ax)
    _legend(ax, [f for f in COLOURS if f in set(table["family"])],
            loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3,
            borderaxespad=0.2)

    fam = summary.iloc[::-1].reset_index(drop=True)
    yb = np.arange(len(fam))
    axw.barh(yb, fam["akaike_weight"], height=0.55,
             color=[COLOURS[f] for f in fam["family"]], edgecolor="white",
             linewidth=1.5)
    for yi, w in zip(yb, fam["akaike_weight"]):
        axw.text(min(w, 1.0) + 0.03, yi, f"{w:.2f}" if w >= 0.005
                 else "<0.01", va="center", color=TEXT2, fontsize=8)
    axw.set_yticks(yb, [FAMILY_NAMES[f] for f in fam["family"]],
                   fontsize=8, color=TEXT)
    axw.set_xlim(0, 1.25)
    axw.set_xticks([0, 0.5, 1])
    axw.set_ylim(-0.8, len(fam) - 0.2)
    axw.set_title("Summed Akaike weight", loc="left", fontsize=9,
                  color=TEXT, pad=20)
    _style(axw)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_importance(importance, families, path):
    """
    Plots each variable's relative importance over both families, split
    into the parts from each family's models.
    """
    imp = importance.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(imp))
    fig, ax = plt.subplots(figsize=(7, max(2.0, 0.42 * len(imp) + 1.2)))
    left = np.zeros(len(imp))
    for family in families:
        part = imp[f"RVI_{family}"].to_numpy()
        ax.barh(y, part, left=left, height=0.55, color=COLOURS[family],
                edgecolor="white", linewidth=1.5)
        left += part
    for yi, total in zip(y, imp["RVI"]):
        ax.text(total + 0.02, yi, f"{total:.2f}" if total >= 0.005
                else "<0.01", va="center", color=TEXT2, fontsize=8)
    ax.set_yticks(y, imp["variable"], fontsize=8, color=TEXT)
    ax.set_xlim(0, 1.12)
    ax.set_xlabel("Relative importance (summed Akaike weight of the models "
                  "with the variable)", color=TEXT2, fontsize=8)
    ax.set_title("Variable importance over both families", loc="left",
                 fontsize=10, color=TEXT, pad=20)
    _style(ax)
    _legend(ax, families, loc="lower right", bbox_to_anchor=(1.0, 1.0),
            ncol=2, borderaxespad=0.2)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_best_fits(runs, summary, path):
    """
    Plots the genetic distances against river distance (isolation by
    distance; the distance-only reference model) and against each family's
    best model distance (standardised, as fitted by MLPE), side by side.
    """
    tables = {f: pd.read_csv(f"{p}.pairwise.tsv", sep="\t")
              for f, (_, _, p) in runs.items()
              if os.path.exists(f"{p}.pairwise.tsv")}
    if not tables:
        return
    best = summary.set_index("family")
    panels = []
    first = next(iter(tables))
    if "river_km" in tables[first]:
        ref = read_null_table(runs[first][2]).loc["distance_only"]
        panels.append(("reference", tables[first], "river_km",
                       f"Isolation by distance: AIC {float(ref['aic']):.1f}, "
                       f"marginal R² {float(ref['r2m']):.2f}\nDistance only "
                       "(river km)", "River distance (standardised)"))
    for family, table in tables.items():
        r2m = pd.read_csv(f"{runs[family][2]}.HallOfFame.tsv", sep="\t",
                          na_values="-")["r2m"].iloc[0]
        panels.append((
            family, table, "best_model_distance",
            f"{FAMILY_NAMES[family]}: AIC "
            f"{best.loc[family, 'best_aic']:.1f}, marginal R² {r2m:.2f}\n"
            + textwrap.fill(best.loc[family, "best_model"], 52),
            "Best model's distance (standardised)"))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.3),
                             sharey=True, squeeze=False)
    for ax, (family, table, column, title, xlabel) in zip(axes[0], panels):
        x = table[column].to_numpy(dtype=float)
        x = (x - x.mean()) / x.std()
        y = table["genetic_distance"].to_numpy(dtype=float)
        ax.scatter(x, y, s=16, color=COLOURS[family], alpha=0.8,
                   edgecolors="white", linewidths=0.5)
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.array([x.min(), x.max()])
        ax.plot(xs, intercept + slope * xs, color=TEXT2, linewidth=1.2)
        ax.set_title(title, loc="left", fontsize=8, color=TEXT)
        ax.set_xlabel(xlabel, color=TEXT2, fontsize=8)
        _style(ax, grid="both")
    axes[0][0].set_ylabel("Genetic distance", color=TEXT2, fontsize=8)
    fig.suptitle("Genetic distance against river distance and each family's "
                 "best model", x=0.01, ha="left", fontsize=10, color=TEXT,
                 y=1.1)
    fig.text(0.01, -0.04, "Each point is a pair of sites. Lines are "
             "least-squares fits, for display: the models are fitted with "
             "MLPE (population effects).", color=TEXT2, fontsize=7)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
