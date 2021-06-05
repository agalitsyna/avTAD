# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli, get_logger
import click

from ..tools import *
from mirnylib import numutils
import pickle

@cli.command()
@click.argument(
    "segmentation",
    metavar="TAD_SEGMENTATION_BED")
@click.argument(
    "map",
    metavar="INPUT_MAP")
@click.argument(
    "output_prefix",
    metavar="OUTPUT_PREFIX")
@click.option(
    "--format", "-f",
    help="Input file format (cool, hiclib_heatmap, hiclib_bychr).",
    is_flag=False,
    default="cool",
    show_default=True)
@click.option(
    "--balance/--no-balance",
    help="Balance the map with iterative correction before snipping. "
         "For cool file it will read the file with option --balance, "
         "for hiclib it will perform deafault balancing with mirnylib.numutils.",
    is_flag=True,
    default=True,
    show_default=True)
@click.option(
    "--niter", "-n",
    help="Number of iterations for segmentation shuffling control.",
    is_flag=False,
    default=0,
    type=int,
    show_default=True)
@click.option(
    "--window", "-w",
    help="Size of window in TAD units. Default is +-1 TAD.",
    is_flag=False,
    default=1,
    type=float,
    show_default=True)
@click.option(
    "--diagonals-to-remove", "-d",
    help="Number of diagonals to remove from map.",
    is_flag=False,
    default=1,
    type=int,
    show_default=True)
@click.option(
    "--enrichment-only",
    help="Compute only dataframe with enrichment of TAD interactions (reduces the time if you don't need the average TAD plot).",
    is_flag=True,
    default=False,
    show_default=True)
def snip(segmentation, map, output_prefix, format, balance, niter, window, enrichment_only, diagonals_to_remove):
    """
    Create snips for TADs and calculate enrichment with shuffled control.
    OUTPUT_PREFIX: The prefix for writing output files (pickle with snips and tsv file with TAD info).

    Output files to be created:
      {OUTPUT_PREFIX}.TADmetadata.tsv
    if not --enrichment-only:
      {OUTPUT_PREFIX}.TADsnips.pickle
      {OUTPUT_PREFIX}.TADsnips_shuf0.pickle etc.

    Example run:
      avTAD snip data/OSC_TADS.bed data/OSC_dm3.cool tmp_results --format cool --diagonals-to-remove 2 --balance --niter 2
    """

    logger = get_logger(__name__)
    logger.info(f"Running snipping for: segmentation file {segmentation}, heatmap {map} in {format} format ...")

    logger.info(f"Reading {map} file with balance={balance} in {format} format ...")
    if format=='cool':
        dataset, chrms, resolution = read_cooler(map, balance=balance)
    elif format=='hiclib_heatmap':
        dataset, chrms, resolution = read_hiclib_heatmap(map, balance=balance)
    elif format=='hiclib_bychr':
        dataset, chrms, resolution = read_hiclib_bychr(map, balance=balance)
    else:
        raise Exception(f'Map format {format} is not supported ...')

    logger.info(f"Reading segmentation file: {segmentation}")
    df_segmentation = pd.read_csv(segmentation, sep='\s', header=None, engine='python')
    add_columns  = list(df_segmentation.columns[3:]) if len(df_segmentation.columns)>3 else []
    df_segmentation.columns = ['ch', 'bgn', 'end'] + add_columns

    df_segmentation.loc[:, 'bgn_bin'] = df_segmentation.bgn // resolution
    df_segmentation.loc[:, 'end_bin'] = df_segmentation.end // resolution
    df_segmentation.loc[:, 'TAD_size'] = df_segmentation.end_bin - df_segmentation.bgn_bin

    df_segmentation = df_segmentation.drop_duplicates().sort_values(['ch', 'bgn_bin']).reset_index(drop=True)

    chrms_used = np.unique(df_segmentation.loc[:, 'ch'].values)
    chrms = [ch for ch in chrms if ch in chrms_used]

    logger.info(f"Chromosomes in the dataset: {dataset.keys()}")
    logger.info(f"Lengths of chromosomes in bins of {resolution} bp: \n{[(ch, len(dataset[ch])) for ch in chrms]}")
    logger.info(f"Selected chromosomes are: {chrms}")

    # Creating shuffled segmentations
    def shuffle_segmentation_dataframe(x):
        segmentation= x[['bgn_bin', 'end_bin']].values
        shuf, order = shuffle_segmentation(segmentation)
        #print(segmentation[0:5], shuf[0:5], order[0:5])
        ret = pd.DataFrame(shuf, columns=['bgn_bin', 'end_bin']).astype(int)
        ret.loc[:, 'index'] = order
        ret.loc[:, 'TAD_size'] = ret.end_bin - ret.bgn_bin
        ret = ret.sort_values('index').reset_index(drop=True)
        assert np.all(x['TAD_size'].values==ret['TAD_size'].values)
        return ret

    for i in range(niter):
        df_segmentation_shuffled = df_segmentation.groupby('ch').apply(shuffle_segmentation_dataframe)\
            .reset_index().drop(['level_1', 'index', 'ch', 'TAD_size'], axis=1)

        df_segmentation_shuffled.columns = [f'bgn_bin_shuf{i}', f'end_bin_shuf{i}']

        df_segmentation = pd.merge(df_segmentation, df_segmentation_shuffled, left_index=True, right_index=True)

    # Computing observed over expected
    dataset_obsexp = {}
    for ch in chrms:
        mtx = numutils.observedOverExpected(dataset[ch])
        #inx_lower_triangle = np.tril_indices(len(mtx))
        #mtx[inx_lower_triangle] = np.nan
        for i in range(1, diagonals_to_remove):
            np.fill_diagonal(mtx[i:, :-i], np.nan)
            np.fill_diagonal(mtx[:-i, i:], np.nan)
        if diagonals_to_remove:
            np.fill_diagonal(mtx, np.nan)
        dataset_obsexp.update({ch:mtx})

    for mod in ['']+[f'_shuf{i}' for i in range(niter)]:
        enrichments = []
        for i, r in df_segmentation.iterrows():
            mtx = np.log2(dataset_obsexp[r.ch][r[f'bgn_bin{mod}']:r[f'end_bin{mod}'], r[f'bgn_bin{mod}']:r[f'end_bin{mod}']])
            mtx[np.isinf(mtx)] = np.nan
            enrichment = (np.nansum(mtx), np.nanmean(mtx), np.nanmedian(mtx), np.sum(np.isfinite(mtx)))
            enrichments.append(np.array(enrichment))
        enrichments = np.array(enrichments)

        df_segmentation.loc[:, f"sum{mod}"]       = enrichments[:, 0]
        df_segmentation.loc[:, f"mean{mod}"]      = enrichments[:, 1]
        df_segmentation.loc[:, f"median{mod}"]    = enrichments[:, 2]
        df_segmentation.loc[:, f"nelements{mod}"] = enrichments[:, 3]

    # Save enrichment dataframe to a file:

    cols = ['bgn_bin', 'end_bin', 'sum', 'mean', 'median', 'nelements']
    columns = ['ch', 'bgn', 'end', 'TAD_size']+add_columns+cols+[f'{x}_shuf{i}' for i in range(niter) for x in cols]
    df_segmentation[columns].to_csv(f"{output_prefix}.TADmetadata.tsv", sep='\t', index=True, header=True)

    if not enrichment_only:
        # Retrieval of snippets, log2 and filling inf with nans included:
        snips = snipper(segmentations=df_segmentation,
                        dataset=dataset_obsexp,
                        window=window)

        # Save snippets to file:
        pickle.dump(snips, open(f"{output_prefix}.TADsnips.pickle", 'wb'))

        for i in range(niter):
            # Retrieval of snippets:
            snips = snipper(segmentations=df_segmentation,
                            dataset=dataset_obsexp,
                            window=window,
                            key_bgn = f'bgn_bin_shuf{i}',
                            key_end = f'end_bin_shuf{i}')

            # Save snippets to file:
            pickle.dump(snips, open(f"{output_prefix}.TADsnips_shuf{i}.pickle", 'wb'))