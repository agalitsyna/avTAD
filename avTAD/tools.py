import hiclib
import mirnylib
from mirnylib import numutils
import cooler

import glob

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

import pandas as pd
import numpy as np

import cooler

def read_cooler(in_cooler, balance=True):
    datasets = {}
    c = cooler.Cooler(in_cooler)
    chrms = c.chromnames
    for ch in chrms:
        datasets[ch] = c.matrix(as_pixels=False, balance=balance).fetch(f'{ch}')

    return datasets, chrms, c.binsize

def read_hiclib_heatmap(infile):
    pass

def read_hiclib_bychr(infile):
    pass


def shuffle_segmentation(segmentation, seed=None):

    if not seed is None:
        np.random.seed(seed)

    tads_lens      = segmentation[:,1] - segmentation[:,0]
    tad_idxs       = np.arange(len(tads_lens))
    tads_idxs_shuf = np.random.permutation(tad_idxs)
    tads_lens_shuf = tads_lens[tads_idxs_shuf]

    intertads_lens      = np.append(segmentation[0,0], segmentation[1:,0] - segmentation[:-1,1])
    intertads_lens_shuf = np.random.permutation(intertads_lens)

    ends = np.cumsum(intertads_lens_shuf+tads_lens_shuf)
    starts = ends-tads_lens_shuf

    segmentation_reconstructed = np.array([starts, ends]).T

    return segmentation_reconstructed, tads_idxs_shuf


def snipper(segmentations, dataset, window=1, key_bgn='bgn_bin', key_end='end_bin'):
    ret = []
    for i, r in segmentations.iterrows():
        size = int(window * (r[key_end] - r[key_bgn]))
        input_mtx = dataset[r.ch]
        bgn = max(0, r[key_bgn] - size)
        end = min(r[key_end] + size, len(input_mtx))
        mtx = input_mtx[bgn:end, bgn:end]
        mtx[np.isinf(mtx)] = np.nan
        ret.append(np.log2(mtx))
    return np.array(ret)

def zoom(snippets, finalShape=(30, 30), saveSum=True, order=1):
    ret = []
    for i, snip in enumerate(snippets):
        mtx = snip
        mtx = numutils.zoomArray(mtx, finalShape, saveSum, order=order)
        ret.append(mtx)
    return ret

def compute_enrichment(mtx, window=1, normalize=False):
    # TODO design a proper test
    size = len(mtx)
    bgn = size*window/(2*window+1)
    end = size-bgn
    center_mtx = mtx[bgn:end, bgn:end]
    if normalize:
        # not is used currently, might be used for calculation of local enrichment
        enrichment = (np.sum(center_mtx)/np.sum(mtx),
                      np.mean(center_mtx)/np.mean(mtx),
                      np.median(center_mtx)/np.median(mtx),
                      np.sum(np.isfinite(center_mtx)))
    else:
        enrichment = np.nansum(center_mtx), np.nanmean(center_mtx), np.nanmedian(center_mtx), np.sum(np.isfinite(center_mtx))
    return enrichment


def plot_heatmap(s, pngname, cmap='Reds', center=None, vmax=None, title='', mask_diags=2):
    plt.figure(figsize=[7, 7])
    mtx = s.copy()
    for i in range(1, mask_diags):
        np.fill_diagonal(mtx[i:, :-i], np.nan)
        np.fill_diagonal(mtx[:-i, i:], np.nan)
    if mask_diags:
        np.fill_diagonal(mtx, np.nan)

    if not center is None and not vmax is None:
        sns.heatmap(mtx, square=True, cmap=cmap,
                    center=center, vmin=-vmax, vmax=vmax)
    elif not center is None:
        sns.heatmap(mtx, square=True, cmap=cmap,
                    center=center)
    else:
        sns.heatmap(mtx, square=True, cmap=cmap)
    # plt.plot([0, 120], [120 / 3, 120 / 3], '--', color='grey', alpha=0.5)
    # plt.plot([0, 120], [2 * 120 / 3, 2 * 120 / 3], '--', color='grey', alpha=0.5)
    # plt.plot([120 / 3, 120 / 3], [0, 120], '--', color='grey', alpha=0.5)
    # plt.plot([2 * 120 / 3, 2 * 120 / 3], [0, 120], '--', color='grey', alpha=0.5)

    plt.xticks([])
    plt.yticks([])

    # enrichment = np.nanmean(mtx[40:80, 40:80])
    # plt.title(f'{title}\nRelative enrichment: {enrichment:.2f}')

    plt.savefig(pngname)
