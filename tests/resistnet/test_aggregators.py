import pytest
import numpy as np
import scipy.stats
import resistnet.aggregators as aggr


# Test cases for the aggregateDist function
@pytest.mark.parametrize("stuff, method, expected_result,\
                         expected_exception", [
    # Valid cases for different aggregation methods
    ([1, 2, 3], "ARITH", np.mean([1, 2, 3]), None),
    ([1, 2, 3], "GEOM", scipy.stats.mstats.gmean([1, 2, 3]), None),
    ([1, 2, 3], "MEDIAN", np.median([1, 2, 3]), None),
    ([1, 2, 3], "MAX", np.max([1, 2, 3]), None),
    ([1, 2, 3], "MIN", np.min([1, 2, 3]), None),
    ([1, 2, 3], "SD", np.std([1, 2, 3]), None),
    ([1, 2, 3], "VAR", np.var([1, 2, 3]), None),
    ([1, 2, 3], "SUM", np.sum([1, 2, 3]), None),
    ([1, 2, 3, 0], "HARM", None, ValueError),
    ([1, 2, 3], "ADJHARM", aggr.adjustedHarmonicMean([1, 2, 3]), None),
    ([], "FIRST", None, ValueError),
    ([1, 2, 3], "CV", np.std([1, 2, 3]) / np.mean([1, 2, 3]), None),

    # Invalid cases
    ("not a list", "ARITH", None, ValueError),
    ([1, 2, 3], "unknown method", None, ValueError),
    ([1, 2, 3, 0], "HARM", None, ValueError),
    ([1, 2, 3, -1], "GEOM", None, ValueError)
])
def test_aggregateDist(stuff, method, expected_result, expected_exception):
    if expected_exception:
        with pytest.raises(expected_exception):
            aggr.aggregateDist(stuff, method)
    else:
        assert aggr.aggregateDist(stuff, method) == \
            pytest.approx(expected_result)


# Test cases for the adjustedHarmonicMean function
@pytest.mark.parametrize("stuff, expected_result", [
    # Valid cases
    ([1, 2, 3], aggr.adjustedHarmonicMean([1, 2, 3])),
    ([0, 1, 2, 3], aggr.adjustedHarmonicMean([0, 1, 2, 3])),
    ([], 0.0),

    # Case with all non-positive values
    ([-1, -2, -3], 0.0)
])
def test_adjustedHarmonicMean(stuff, expected_result):
    assert aggr.adjustedHarmonicMean(stuff) == pytest.approx(expected_result)
