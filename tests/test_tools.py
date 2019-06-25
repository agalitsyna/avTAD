from __future__ import division, print_function
import numpy as np
import pandas as pd
from avTAD.tools import *

def test_loading():
    """
    Example test:
      pytest tests/test_tools.py::test_loading
    """
    d1 = read_cooler('data/OSC_dm3.cool', balance=False)
    d2 = read_hiclib_heatmap('data/OSC_dm3_rawmap.hdf5', balance=False)
    d3 = read_hiclib_heatmap('data/OSC_dm3_rawmap.hdf5', balance=False)

    assert d1[2]==d2[2]
    assert d2[2]==d3[2]

    diagonals_to_remove = 2
    for ch in d2[1]:
        mtx2 = d2[0][ch].astype(float)
        mtx3 = d2[0][ch].astype(float)
        for i in range(1, diagonals_to_remove):
            np.fill_diagonal(mtx2[i:, :-i], np.nan)
            np.fill_diagonal(mtx2[:-i, i:], np.nan)
        if diagonals_to_remove:
            np.fill_diagonal(mtx2, np.nan)
        for i in range(1, diagonals_to_remove):
            np.fill_diagonal(mtx3[i:, :-i], np.nan)
            np.fill_diagonal(mtx3[:-i, i:], np.nan)
        if diagonals_to_remove:
            np.fill_diagonal(mtx3, np.nan)

        assert np.nansum(mtx2)==np.nansum(mtx3)