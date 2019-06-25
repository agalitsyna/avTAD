# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli, get_logger
import click

import sys
import os
import datetime
now = datetime.datetime.now()

from ..tools import *

@cli.command()
@click.argument(
    "a_prefix",
    metavar="A_PREFIX")
@click.argument(
    "b_prefix",
    metavar="B_PREFIX")
@click.argument(
    "output_prefix",
    metavar="OUTPUT_PREFIX")
@click.argument(
    "expression",
    metavar="EXPRESSION")

def evaluate(a_prefix, b_prefix, output_prefix, expression):
    """
    Evaluate matrix expression for a pair of prefixes, A_PREFIX and B_PREFIX.
    If multiple avTAD tsv files available, will be applied to corresponding files.

    Example usage:
       avTAD evaluate OSC OSC_shuf0 OSC_enrichment "a-b"
    """

    logger = get_logger(__name__)
    logger.info(f"Loading matrices from {a_prefix} and {b_prefix} for evaluation {expression} ...")

    files_list_a = glob.glob(f'{a_prefix}.avTAD*.tsv')
    files_list_b = glob.glob(f'{b_prefix}.avTAD*.tsv')
    modes = ['.'+x.split('.')[-2] if not x.split('.')[-2]=='avTAD' else '' for x in files_list_a]
    print(files_list_a, modes)

    for mode in modes:
        a_name = f'{a_prefix}.avTAD{mode}.tsv'
        b_name = f'{b_prefix}.avTAD{mode}.tsv'

        if not a_name in files_list_a or not b_name in files_list_b:
            logger.warning(f"Skipping {a_name} and {b_name}")
            continue

        a = np.loadtxt(a_name, delimiter='\t')
        b = np.loadtxt(b_name, delimiter='\t')

        if len(np.unique(expression))>5:
            raise Exception(f"Bad values in {expression}, only a, b and math operations are allowed")

        c = eval(expression)
        c[np.isinf(c)] = np.nan
        logger.info(f"Saving output to {output_prefix}.avTAD{mode}.tsv ...")
        np.savetxt(f'{output_prefix}.avTAD{mode}.tsv', c,
                   header=f'{expression} for a={a_prefix} and b={b_prefix} at {now.strftime("%Y-%m-%d %H:%M")}',
                   fmt='%.6e', delimiter='\t')