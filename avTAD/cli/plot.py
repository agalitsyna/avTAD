# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli, get_logger
import click

import glob
import numpy as np
import scipy
from ..tools import plot_heatmap

@cli.command()
@click.argument(
    "infile_prefix",
    metavar="INFILE_PREFIX")
@click.argument(
    "output_prefix",
    metavar="OUTPUT_PREFIX")
@click.option(
    "--cmap",
    help="Name of colormap for plotting. Use one of seaborn/matplotlib names (Reds, RdBu_r).",
    is_flag=False,
    default='RdBu_r',
    show_default=True)
@click.option(
    "--vmax",
    help="Max value for the color scale.",
    is_flag=False,
    default=1,
    type=float,
    show_default=True)
@click.option(
    "--vmin",
    help="Min value for the color scale.",
    is_flag=False,
    default=-1,
    type=float,
    show_default=True)
@click.option(
    "--center",
    help="Center value for the color scale.",
    is_flag=False,
    default=0,
    type=float,
    show_default=True)
@click.option(
    "--autoscale",
    help="Automatic vmin and vmax based on values distribution (99-percentile). Rewrites vmax and vmin.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--figsize",
    help="Figure size passed to matplotlib.figure().",
    is_flag=False,
    default=10,
    type=float,
    show_default=True)
@click.option(
    "--cbar/--no-cbar",
    help="Whether to plot colorbar.",
    is_flag=True,
    default=True,
    show_default=True)
@click.option(
    "--window", "-w",
    help="Size of window in TAD units. Default is +-1 TAD.",
    is_flag=False,
    default=1,
    type=float,
    show_default=True)
def plot(infile_prefix, output_prefix, cmap, vmax, vmin, center, figsize, cbar, window, autoscale):
    """Plotting average TAD rescaled to the same size.
    Input is INFILE_PREFIX , output is written to OUTFILE_PREFIX

    File to be created:
        {OUTFILE_PREFIX}.avTAD.png
        {OUTFILE_PREFIX}.avTAD.ch:chrX.png etc. if multiple files are avaiable with INFILE_PREFIX

    Example usage:
        avTAD plot --cmap RdBu_r OSC OSC
    """

    logger = get_logger(__name__)
    logger.info(f'Reading infile prefixes: {infile_prefix}')

    file_list = glob.glob(f'{infile_prefix}.avTAD*.tsv')
    modes_list = ['.'+x.split('.')[-2] if not x.split('.')[-2]=='avTAD' else '' for x in file_list ]

    for f, mode in zip(file_list, modes_list):
        logger.info(f"Reading {f}")
        table = np.loadtxt(f, delimiter='\t')
        comment_line = open(f, 'r').readline()
        if comment_line.startswith('#'):
            title = comment_line[1:].split()
            title = '\n'.join([' '.join(title[(6 * i):(6 * i) + 6]) for i in range(1 + len(title) // 6)])
        else:
            title = f"Plotting matrix from {f}"

        pngname = f'{output_prefix}.avTAD{mode}.png'
        if autoscale:
            mx = np.nanpercentile(np.abs(table.flatten()), 99)
            vmax = mx
            vmin = -mx
        if not vmin:
            vmin = -vmax
        logger.info(f"Writing heatmap to {pngname} with {cmap} color map, vmax={vmax}, vmin={vmin}, center={center}, figsize={figsize}, cbar={cbar}")
        plot_heatmap(table, pngname, cmap=cmap, center=center, vmax=vmax, title=title, cbar=cbar, figsize=[figsize, figsize], window=window)