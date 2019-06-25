import mirnylib
from mirnylib import numutils
import cooler
import pickle
import h5py

import glob

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np

def read_cooler(in_cooler, balance=True):
    datasets = {}
    c = cooler.Cooler(in_cooler)
    chrms = c.chromnames
    for ch in chrms:
        datasets[ch] = c.matrix(as_pixels=False, balance=balance).fetch(f'{ch}')

    return datasets, chrms, c.binsize

def read_hiclib_heatmap(infile, balance=False):
    datasets = {}
    f = h5py.File(infile)
    map = f['heatmap'][()]
    resolution = pickle.loads(f['resolution'][()])
    chrms_sizes = f['chromosomeStarts'][()]
    idx2chrms = pickle.loads(f['genomeIdxToLabel'][()])
    for idx in idx2chrms.keys():
        ch = idx2chrms[idx]
        ch = ch if 'chr' in ch else f'chr{ch}'
        bgn = chrms_sizes[idx]
        end = chrms_sizes[idx+1] if idx+1<len(chrms_sizes) else len(map)
        mtx = map[bgn:end, bgn:end]
        if balance:
            mtx = numutils.iterativeCorrection(mtx)
        datasets[ch] = mtx

    return datasets, sorted(datasets.keys()), resolution

def read_hiclib_bychr(infile, balance=False):
    datasets = {}
    f = h5py.File(infile)
    resolution = pickle.loads(f['resolution'][()])
    idx2chrms = pickle.loads(f['genomeIdxToLabel'][()])
    for idx in idx2chrms.keys():
        ch = idx2chrms[idx]
        mtx = f[f'{idx} {idx}'][()]
        ch = ch if 'chr' in ch else f'chr{ch}'
        if balance:
            mtx = numutils.iterativeCorrection(mtx)
        datasets[ch] = mtx

    return datasets, sorted(datasets.keys()), resolution

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
        mtx = np.log2(input_mtx[bgn:end, bgn:end])
        mtx[np.isinf(mtx)] = np.nan
        ret.append(mtx)
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


def plot_heatmap(s, pngname, cmap='Reds', center=None, vmax=None, title='', cbar=True, figsize=[7,7], window=1):

    plt.figure(figsize=figsize)
    mtx = s.copy()

    if not center is None and not vmax is None:
        sns.heatmap(mtx, square=True, cmap=cmap,
                    center=center, vmin=-vmax, vmax=vmax, cbar=cbar)
    elif not center is None:
        sns.heatmap(mtx, square=True, cmap=cmap,
                    center=center, cbar=cbar)
    else:
        sns.heatmap(mtx, square=True, cmap=cmap, cbar=cbar)

    size = len(mtx)
    bgn = size * window / (2 * window + 1)
    end = size - bgn
    plt.plot([0, size], [bgn, bgn], '--', color='gray', alpha=0.5)
    plt.plot([0, size], [end, end], '--', color='gray', alpha=0.5)
    plt.plot([bgn, bgn], [0, size], '--', color='gray', alpha=0.5)
    plt.plot([end, end], [0, size], '--', color='gray', alpha=0.5)

    plt.xticks([])
    plt.yticks([])

    plt.title(title)

    plt.savefig(pngname)
