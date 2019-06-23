# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli
import click

@cli.command()
@click.argument(
    "infile",
    metavar="INFILE_PICKLE")

def plot(infile):
    """Plotting average TAD rescaled to the same size."""
    print(infile)