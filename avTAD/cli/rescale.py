# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli, get_logger
import click

import pickle
import sys
import os
import datetime
now = datetime.datetime.now()

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

from ..tools import *

@cli.command()
@click.argument(
    "infile_pickle",
    metavar="INFILE_PICKLE")
@click.argument(
    "infile_table",
    metavar="INFILE_TABLE")
@click.argument(
    "output_prefix",
    metavar="OUTPUT_PREFIX")
# INFILE_TABLE reading/splitting parameters
@click.option(
    "--split-by",
    help="Split by groups of input columns. Split by ch is always available. Use bgn to save individual TADs.",
    metavar="SPLIT_BY",
    is_flag=False,
    default=None,
    show_default=True)
@click.option(
    "--table-is-indexed/--table-is-not-indexed",
    help="INFILE_TABLE has first index column. If not, then default sequential indexing is used. "
         "Snips are taken from pickle file according to the index. It should correspond to the step of pickle creation by avTAD snip",
    is_flag=True,
    default=True,
    show_default=True)
@click.option(
    "--table-has-header/--table-has-no-header",
    help="INFILE_TABLE has header in the first column. If not, then first columns are pre-defined as: ch, bgn, end.",
    is_flag=True,
    default=True,
    show_default=True)
@click.option(
    "--query",
    help="Query to pass to INFILE_TABLE while reading pandas dataframe, should be quoted in the terminal. "
         "Example: \"TAD_size>20 and TAD_size<30\" or \"ch=='chrX'\"",
    is_flag=False,
    default=None,
    show_default=True)
# INFILE_PICKLE snips rescaling parameters
@click.option(
    "--rescaled-size",
    help="The resulting size of the average TAD plot after rescaling. Larger size than he max TAD size is recommended.",
    is_flag=False,
    default=200,
    type=int,
    show_default=True)
@click.option(
    "--save-sum/--no-save-sum",
    help="Save sum by the zoom operation. --save-sum mode is recommended.",
    is_flag=True,
    default=True,
    show_default=True)
@click.option(
    "--smooth-order",
    help="The order of polynome to interpolate image.",
    is_flag=False,
    default=1,
    type=int,
    show_default=True)
@click.option(
    "--operation",
    help="Operation to perform over each pixel of rescaled image (mean, median, sum, count).",
    is_flag=False,
    default='mean',
    show_default=True)

def rescale(infile_pickle, infile_table, output_prefix, table_is_indexed, table_has_header, query, rescaled_size, save_sum, smooth_order, operation, split_by):
    """
    Rescale snips to the same size and average them based on index. Outputs average TAD matrix in tsv format with header with run metadata.

    Output files to be created:
    if no --split-by provided:
        {OUTPUT_PREFIX}.avTAD.tsv
    if --split-by SPLIT_BY , then a set of files will be created:
        {OUTPUT_PREFIX}.avTAD.{SPLIT_BY}:{values}.tsv

    Example usage:
       avTAD rescale OSC.TADsnips.pickle OSC.TADmetadata.tsv OSC --rescaled-size 200
    """

    logger = get_logger(__name__)
    logger.info(f"Loading pickle file {infile_pickle} ...")

    # Setting parameter for further averaging operation on snips:
    if operation=='mean':
        func = np.nanmean
    elif operation=='median':
        func = np.nanmedian
    elif operation=='sum':
        func = np.nansum
    elif operation == 'count':
        func = lambda x: np.nansum(np.isfinite(x))
    else:
        raise Exception(f'Operation {operation} is not implemented... Exiting.')

    snips = pickle.load(open(infile_pickle, 'rb'))

    logger.info(f"Performing filtering based on dataframe passed in {infile_table}")

    if infile_table:
        df_segmentation = pd.read_csv(infile_table, sep='\s',
                                  header=0 if table_has_header else None,
                                  index_col=0 if table_is_indexed else None,
                                  engine='python'
                                  )
    else:
        infile_table = 'stdout'
        inf = StringIO(sys.stdin)
        df_segmentation = pd.read_csv(inf, sep='\s',
                                  header=0 if table_has_header else None,
                                  index_col=0 if table_is_indexed else None,
                                  engine='python'
                                  ).query(query)

    if not table_has_header:
        df_segmentation = df_segmentation.loc[:, 0:3]
        df_segmentation.columns = ['ch', 'bgn', 'end']

    if query:
        df_segmentation = df_segmentation.query(query)

    assert df_segmentation.index.max() <= len(snips) # segmentation index is probably not aligned with snips if False

    if split_by is None:
        index = df_segmentation.index
        logger.info(f"Selecting index from {infile_pickle} ({len(index)} elements) ...")

        logger.info(f"Zooming snip arrays to {rescaled_size}x{rescaled_size} ...")
        snips_rescaled = np.array(
            zoom(snips[index], finalShape=(rescaled_size, rescaled_size), saveSum=save_sum, order=smooth_order)
        )

        logger.info(f"Creating average plot with function {operation} ...")
        averaged_matrix = func(snips_rescaled, axis=0)

        if os.path.isfile(f"{output_prefix}.avTAD.tsv"):
            logger.warning(f"File {output_prefix}.avTAD.tsv exists, it will be overwritten!")

        logger.info(f"Saving output to {output_prefix}.avTAD.tsv ...")
        np.savetxt(f'{output_prefix}.avTAD.tsv', averaged_matrix,
                   header=f'{len(index)} snips from {infile_pickle} as indexed in {infile_table} (query: {query}) averaged by {operation} at {now.strftime("%Y-%m-%d %H:%M")}',
                   fmt='%.6e', delimiter='\t')

    else:
        if not split_by in df_segmentation.columns:
            raise Exception(
                f"{split_by} in not in df_segmentation columns. Available choices are: {df_segmentation.columns}")
        for name, group in df_segmentation.groupby(split_by):
            index = group.index
            logger.info(f"Selecting index from {infile_pickle}, group {name} of {split_by} ({len(index)} elements) ...")

            logger.info(f"Zooming snip arrays to {rescaled_size}x{rescaled_size} ...")
            snips_rescaled = np.array(
                zoom(snips[index], finalShape=(rescaled_size, rescaled_size), saveSum=save_sum, order=smooth_order)
            )

            logger.info(f"Creating average plot with function {operation} ...")
            averaged_matrix = func(snips_rescaled, axis=0)

            fname = f"{output_prefix}.avTAD.{split_by}:{name}.tsv"
            if os.path.isfile(fname):
                logger.warning(f"File {fname} exists, it will be overwritten!")

            logger.info(f"Saving output to {fname} ...")
            np.savetxt(fname, averaged_matrix,
                       header=f'{len(index)} snips from {infile_pickle} as indexed in {infile_table} (query: {query}) averaged by {operation} at {now.strftime("%Y-%m-%d %H:%M")} group {name} of {split_by}',
                       fmt='%.6e', delimiter='\t')