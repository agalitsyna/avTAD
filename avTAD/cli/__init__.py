from __future__ import division, print_function
import logging
import sys
import os
from .._version import __version__
from .._logging import get_logger
import click

# Monkey patch
click.core._verify_python3_env = lambda: None


CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}

class UnsortedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)


@click.version_option(__version__, '-V', '--version')
@click.group(context_settings=CONTEXT_SETTINGS, cls=UnsortedGroup)
@click.option(
    '-v', '--verbose',
    help="Verbose logging.",
    count=True)
# @click.option(
#     '-d', '--debug',
#     help="On error, drop into the post-mortem debugger shell.",
#     is_flag=True,
#     default=False)
# def cli(verbose, debug):
def cli(verbose):

    """
    Type -h or --help after any subcommand for more information.
    """
    # Initialize logging to stderr
    logging.basicConfig(stream=sys.stderr)
    logging.captureWarnings(True)
    root_logger = get_logger()

    # Set verbosity level
    if verbose > 0:
        root_logger.setLevel(logging.DEBUG)
    else:
        root_logger.setLevel(logging.INFO)

from . import (
    plot,
    snip
)