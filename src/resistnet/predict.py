"""
Applies a fitted model to the reaches of another network (e.g. another
watershed), using the constants saved with the model ({out}.modelSpec.json;
see ResistanceNetwork.model_spec()). Every rescaling uses the training
network's ranges, so a reach with the same values as a training reach gets
the same result; values outside the training range are clamped to it. The
variables' combined value is a weighted mean of their 0-1 values, so it is
always within range. Standard and SAMC models share the resistance per km
of each reach; a SAMC model adds its flow term.
"""
import json

import numpy as np
import pandas as pd

import resistnet.transform as trans


def load_spec(path):
    """
    Reads a model specification file ({out}.modelSpec.json).

    Args:
        path (str): Path of the file.

    Returns:
        dict: The specification; its "models" are best first.
    """
    with open(path) as fh:
        spec = json.load(fh)
    # Standard models of earlier versions had resistances per reach, and
    # SAMC models of earlier versions moved between reaches
    for m in spec.get("models", []):
        if m.get("type") == "standard" and m.get("units") != "km":
            raise ValueError(
                f"{path} has standard models of an earlier version "
                "(resistances per reach), which are no longer supported; "
                "refit them")
        if m.get("type") == "samc" and "cell_km" not in m:
            raise ValueError(
                f"{path} has SAMC models of an earlier version (moves "
                "between reaches), which are no longer supported; refit "
                "them")
    return spec


def _rescale(x, lo, hi):
    """Rescales x so that lo -> 0 and hi -> 1 (all 0 if hi == lo)."""
    x = np.asarray(x, dtype=float)
    if hi > lo:
        return (x - lo) / (hi - lo)
    return np.zeros_like(x)


def out_of_range(raw, var):
    """Share of values outside the training range of a variable."""
    raw = np.asarray(raw, dtype=float)
    return float(np.mean((raw < var["min"]) | (raw > var["max"])))


def variable_values(raw, var):
    """
    A variable's values on the training network's 0-1 scale: rescaled to
    0-10 with the training range (and clamped to it), transformed, and
    rescaled with the training range of the transformed values.

    Args:
        raw (array-like): Raw values of the variable.
        var (dict): The variable's entry in a model specification.

    Returns:
        numpy.ndarray: Values from 0 to 1.
    """
    z = np.clip(10.0 * _rescale(raw, var["min"], var["max"]), 0.0, 10.0)
    # Some transformations use the range of the data (e.g. to reverse it);
    # the training values spanned exactly 0-10, so adding those two values
    # gives the same constants
    padded = pd.Series(np.concatenate([z, [0.0, 10.0]]))
    y = trans.apply_transformation(padded, var["transform"], var["shape"])
    y = np.asarray(y, dtype=float)[:len(z)]
    return np.clip(_rescale(y, var["tmin"], var["tmax"]), 0.0, 1.0)


def standard_composite(raw, model):
    """
    Combined value of each reach under a standard (path-sum) model: the
    weighted mean of the variables (each clamped to 0-1, so it is within
    0-1).

    Args:
        raw (DataFrame): Raw values of the model's variables, one row per
                         reach.
        model (dict): A model from a specification (type "standard").

    Returns:
        numpy.ndarray: Combined values (None if the model has no variables).
    """
    multi = None
    total = 0.0
    for v, var in model["variables"].items():
        x = variable_values(raw[v], var) * var["weight"]
        total += var["weight"]
        multi = x if multi is None else multi + x
    if multi is None:
        return None
    return multi / total if total > 0 else multi * 0.0


def standard_resistance(raw, model):
    """
    Resistance per km of each reach under a standard (path-sum) model,
    1 + (K - 1) * composite (from 1, distance alone, to the maximum
    resistance K).

    Args:
        raw (DataFrame): Raw values of the model's variables, one row per
                         reach.
        model (dict): A model from a specification (type "standard").

    Returns:
        numpy.ndarray: Resistance of each reach.
    """
    x = standard_composite(raw, model)
    if x is None:
        return np.ones(len(raw))
    x = np.clip(x, 0.0, 1.0)
    return 1.0 + (model["max_resistance"] - 1.0) * x


def samc_downstream(lengths, model):
    """
    Probability of moving downstream from a state inside each reach under a
    SAMC model: its two neighbours are on the same reach, the same distance
    away, so only the flow term makes them differ,
    1 / (1 + exp(-2 x flow x link length)); NaN for reaches too short to
    have inner states (one link).

    Args:
        lengths (array-like): Length of each reach (km).
        model (dict): A model from a specification (type "samc").

    Returns:
        numpy.ndarray: Probability for each reach.
    """
    lengths = np.asarray(lengths, dtype=float)
    # As the lattice (resistnet.lattice): round(L / cell_km) links, halves
    # to even as Python's round()
    pieces = np.maximum(1, np.rint(lengths / model["cell_km"]))
    link = lengths / pieces
    p = 1.0 / (1.0 + np.exp(-2.0 * model.get("flow", 0.0) * link))
    return np.where(pieces > 1, p, np.nan)


def predict_reaches(reaches, spec, model=0, reachid_col=None,
                    length_col=None):
    """
    Applies a model to each reach of a reach layer.

    Each reach's resistance: per km (from 1), and times the reach's length
    (reach_resistance). For a SAMC model, also the probability of moving
    downstream from a state inside the reach (p_downstream; see
    samc_downstream()).

    Args:
        reaches (DataFrame): One row per reach, with the model's variables.
        spec (dict): Model specification (see load_spec()).
        model (int): Index of the model in the specification (0 = best).
        reachid_col, length_col (str): Column names (default: as when the
                                       model was fitted).

    Returns:
        tuple: Per-reach DataFrame, and the share of reaches with values
               outside the training range (clamped to it) for each
               variable.
    """
    m = spec["models"][model]
    reachid_col = reachid_col or spec.get("reachid_col", "HYRIV_ID")
    length_col = length_col or spec.get("length_col", "LENGTH_KM")
    missing = [v for v in m["variables"] if v not in reaches.columns]
    if missing:
        raise ValueError(f"Variables missing from the reaches: {missing}")
    raw = reaches[list(m["variables"])].astype(float)
    clamped = {v: out_of_range(raw[v], var)
               for v, var in m["variables"].items()}
    out = pd.DataFrame({reachid_col: reaches[reachid_col].to_numpy()})
    out["resistance"] = standard_resistance(raw, m)
    # Also each reach's resistance over its length (its contribution to the
    # distances between sites)
    if length_col in reaches:
        lengths = reaches[length_col].to_numpy(dtype=float)
        out["reach_resistance"] = out["resistance"].to_numpy() * lengths
        if m["type"] == "samc":
            out["p_downstream"] = samc_downstream(lengths, m)
    elif m["type"] == "samc":
        raise ValueError(f"SAMC models need the reach lengths "
                         f"({length_col})")
    return out, clamped
