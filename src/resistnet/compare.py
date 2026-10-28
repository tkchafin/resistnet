"""
Compares the standard (resistance) and SAMC models fitted to the same
genetic distances (runResistnet.py --model both): one ranking of the models
of both, and of the reference models, by AIC with joint Akaike weights.

Both models are fitted with the same MLPE to the same genetic distances, so
their AICs are comparable; the reference models computed by both (null and
distance_only) must therefore be the same, which is checked.
"""
import numpy as np
import pandas as pd

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
    and, for SAMC, their direction in raw values and the flow rate.

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
        parts.append(f"flow {flow:+.3f}/km (upstream "
                     f"{'harder' if flow > 0 else 'easier'})")
    return " + ".join(parts)


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
    rows, nulls, selected = [], {}, []
    for family, (network, runner, prefix) in runs.items():
        hof = pd.read_csv(f"{prefix}.HallOfFame.tsv", sep="\t",
                          na_values="-")
        variables = list(network._predictors.columns)
        seen = set()
        for _, row in hof.iterrows():
            model = describe(row, variables, network, runner)
            # Each chain's best model is kept, so a model can appear more
            # than once (for SAMC without variables, with different
            # maximum resistances that have no effect): list it once
            key = (model, round(float(row["aic"]), 6))
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "family": family, "model": model,
                "k": int(round((float(row["aic"]) + 2 * float(row["loglik"]))
                               / 2)),
                "loglik": float(row["loglik"]), "aic": float(row["aic"])})
            selected.append({v for v in variables if row[v] == 1})
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

    table = pd.DataFrame(rows)
    table["delta_aic"] = table["aic"] - table["aic"].min()
    w = np.exp(-table["delta_aic"] / 2)
    table["akaike_weight"] = w / w.sum()
    order = np.argsort(table["aic"].to_numpy(), kind="stable")
    table = table.iloc[order].reset_index(drop=True)
    selected = [selected[i] for i in order]
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

    variables = sorted(set().union(*selected))
    importance = pd.DataFrame({
        "variable": variables,
        "RVI": [sum(w for w, s in zip(table["akaike_weight"], selected)
                    if v in s) for v in variables]
    }).sort_values("RVI", ascending=False, kind="stable")
    return table, summary, importance


def write_comparison(runs, prefix, verbose=True):
    """
    Writes the joint ranking of the fitted model families (see
    comparison_tables()) to {prefix}.modelComparison.tsv, the summary per
    family to {prefix}.familySummary.tsv, and the joint relative variable
    importance to {prefix}.jointImportance.tsv.

    Args:
        runs (dict): As for comparison_tables().
        prefix (str): Output prefix.
        verbose (bool): Print the summary.

    Returns:
        tuple: The three tables.
    """
    table, summary, importance = comparison_tables(runs)
    table.to_csv(f"{prefix}.modelComparison.tsv", sep="\t", index=False)
    summary.to_csv(f"{prefix}.familySummary.tsv", sep="\t", index=False)
    importance.to_csv(f"{prefix}.jointImportance.tsv", sep="\t",
                      index=False)
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
