"""
Applies a fitted model to the reaches of another network (e.g. another
watershed), using the constants saved with the model ({out}.modelSpec.json;
see ResistanceNetwork.model_spec()). Every rescaling uses the training
network's ranges, so a reach with the same values as a training reach gets
the same result; values outside the training range are clamped to it. The
variables' combined value is clamped to its training range as well, so
reaches beyond it all get the same value.
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
        return json.load(fh)


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


def _outside(x):
    """Share of rescaled values outside 0-1 (beyond rounding)."""
    x = np.asarray(x, dtype=float)
    if not x.size:
        return 0.0
    return float(np.mean((x < -1e-9) | (x > 1.0 + 1e-9)))


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
    Combined value of each reach under a standard (path-sum) model, on the
    training network's 0-1 scale before clamping.

    Args:
        raw (DataFrame): Raw values of the model's variables, one row per
                         reach.
        model (dict): A model from a specification (type "standard").

    Returns:
        numpy.ndarray: Combined values (None if the model has no variables).
    """
    multi = None
    for v, var in model["variables"].items():
        x = variable_values(raw[v], var) * var["weight"]
        multi = x if multi is None else multi + x
    if multi is None:
        return None
    return _rescale(multi, model["composite_min"], model["composite_max"])


def standard_resistance(raw, model):
    """
    Resistance (0-1) of each reach under a standard (path-sum) model.

    Args:
        raw (DataFrame): Raw values of the model's variables, one row per
                         reach.
        model (dict): A model from a specification (type "standard").

    Returns:
        numpy.ndarray: Resistance of each reach.
    """
    x = standard_composite(raw, model)
    if x is None:
        return np.zeros(len(raw))
    return np.clip(x, 0.0, 1.0)


def samc_composite(raw, rows, cols, model):
    """
    Combined value of each move between adjacent reaches under a SAMC
    model, on the training network's 0-1 scale before clamping.

    Args:
        raw (DataFrame): Raw values of the model's variables, one row per
                         reach.
        rows, cols (array-like): Indices of the reaches moved from and to.
        model (dict): A model from a specification (type "samc").

    Returns:
        tuple: Combined values (None if the model has no variables), and
               whether any variable's link value was outside its training
               range (and clamped to it), for each move.
    """
    rows, cols = np.asarray(rows), np.asarray(cols)
    multi = None
    link_out = np.zeros(len(rows), dtype=bool)
    for v, var in model["variables"].items():
        x = variable_values(raw[v], var)
        if var["asym"]:
            link = x[cols] - x[rows]
        else:
            link = (x[rows] + x[cols]) / 2.0
        link = _rescale(link, var["link_min"], var["link_max"])
        link_out |= (link < -1e-9) | (link > 1.0 + 1e-9)
        link = np.clip(link, 0.0, 1.0) * var["weight"]
        multi = link if multi is None else multi + link
    if multi is None:
        return None, link_out
    return (_rescale(multi, model["composite_min"], model["composite_max"]),
            link_out)


def _move_weights(composite, upstream_km, model, n):
    """Movement weights from combined values (see samc_weights())."""
    if composite is None:
        weights = np.ones(n)
    else:
        x = np.clip(composite, 0.0, 1.0)
        weights = 1.0 / (1.0 + (model["max_resistance"] - 1.0) * x)
    if model.get("flow", 0.0):
        weights = weights * np.exp(-model["flow"]
                                   * np.asarray(upstream_km, dtype=float))
    return weights


def samc_weights(raw, rows, cols, upstream_km, model):
    """
    Movement weight of each move between adjacent reaches under a SAMC
    model (before normalising over each reach's moves).

    Args:
        raw (DataFrame): Raw values of the model's variables, one row per
                         reach.
        rows, cols (array-like): Indices of the reaches moved from and to.
        upstream_km (array-like): Distance moved upstream (negative when
                                  moving downstream), between the reaches'
                                  midpoints.
        model (dict): A model from a specification (type "samc").

    Returns:
        numpy.ndarray: Weight of each move.
    """
    composite, _ = samc_composite(raw, rows, cols, model)
    return _move_weights(composite, upstream_km, model, len(np.asarray(rows)))


def reach_links(reaches, reachid_col, next_down_col, length_col):
    """
    Moves between adjacent reaches of a reach layer (each reach and its next
    downstream reach, as in the SAMC model), and the distance moved upstream
    between the reaches' midpoints.

    Args:
        reaches (DataFrame): One row per reach.
        reachid_col, next_down_col, length_col (str): Column names.

    Returns:
        tuple: rows, cols (reach positions) and km moved upstream.
    """
    ids = reaches[reachid_col].to_numpy()
    pos = {r: i for i, r in enumerate(ids)}
    length = reaches[length_col].to_numpy(dtype=float)
    up, down = [], []
    for i, nd in enumerate(reaches[next_down_col].to_numpy()):
        if nd in pos:
            up.append(i)
            down.append(pos[nd])
    up, down = np.array(up, dtype=int), np.array(down, dtype=int)
    km = (length[up] + length[down]) / 2.0
    rows = np.concatenate([up, down])
    cols = np.concatenate([down, up])
    return rows, cols, np.concatenate([-km, km])


def predict_reaches(reaches, spec, model=0, reachid_col=None,
                    next_down_col="NEXT_DOWN", length_col=None):
    """
    Applies a model to each reach of a reach layer.

    For a standard model, each reach's resistance (0-1). For a SAMC model,
    the weights of moving to the downstream reach and to the upstream
    reaches, and the share of moves that go downstream (without
    absorption); and a table of the moves.

    Args:
        reaches (DataFrame): One row per reach, with the model's variables.
        spec (dict): Model specification (see load_spec()).
        model (int): Index of the model in the specification (0 = best).
        reachid_col, length_col (str): Column names (default: as when the
                                       model was fitted).
        next_down_col (str): Column giving the next downstream reach (SAMC).

    Returns:
        tuple: Per-reach DataFrame, moves DataFrame (SAMC; else None), and
               the share of out-of-range values for each variable; under
               "composite", the share of reaches (standard) or moves (SAMC)
               whose combined value was outside the training range, and was
               clamped to it.
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
    if m["type"] == "standard":
        x = standard_composite(raw, m)
        clamped["composite"] = 0.0 if x is None else _outside(x)
        out["resistance"] = (np.zeros(len(raw)) if x is None
                             else np.clip(x, 0.0, 1.0))
        return out, None, clamped

    rows, cols, upstream = reach_links(reaches, reachid_col, next_down_col,
                                       length_col)
    composite, link_out = samc_composite(raw, rows, cols, m)
    if composite is not None:
        out_any = (link_out | (composite < -1e-9)
                   | (composite > 1.0 + 1e-9))
        clamped["composite"] = float(np.mean(out_any)) if len(rows) else 0.0
    else:
        clamped["composite"] = 0.0
    w = _move_weights(composite, upstream, m, len(rows))
    n = len(reaches)
    is_down = upstream < 0
    w_down = np.bincount(rows[is_down], weights=w[is_down], minlength=n)
    w_up = np.bincount(rows[~is_down], weights=w[~is_down], minlength=n)
    has_down = np.bincount(rows[is_down], minlength=n) > 0
    has_up = np.bincount(rows[~is_down], minlength=n) > 0
    out["weight_downstream"] = np.where(has_down, w_down, np.nan)
    out["weight_upstream"] = np.where(has_up, w_up, np.nan)
    total = w_down + w_up
    with np.errstate(invalid="ignore", divide="ignore"):
        out["p_downstream"] = np.where(total > 0, w_down / total, np.nan)
    ids = reaches[reachid_col].to_numpy()
    moves = pd.DataFrame({"from": ids[rows], "to": ids[cols],
                          "direction": np.where(is_down, "downstream",
                                                "upstream"),
                          "km_upstream": upstream, "weight": w})
    return out, moves, clamped
