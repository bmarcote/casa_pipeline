"""Import Tsys information from an antab file to a Measurement Set (MS) file.

Usage: antab2ms.py  <msfile>  <antabfile>

Options:
    msdata : str
        MS data set to incorporate the Tsys information.
    antabfile : str
        ANTAB file to be read.


Version: 0.0
Date: October 2018
Written by Benito Marcote (marcote@jive.eu)
"""

import sys
import time
import argparse
import datetime as dt
from collections import defaultdict, namedtuple
import numpy as np
from pyrap import tables as pt



usage = "%(prog)s [-h] [-v]  <measurement set>  <antab file>"
description = "Import Tsys information from an antab file to a Measurement Set (MS) file."
help_msdata = "Measurement Set to be read."
help_antabfile = "ANTAB file to be read."

parser = argparse.ArgumentParser(description=description, prog='antab2ms.py', usage=usage)
parser.add_argument('msdata', type=str, help=help_msdata)
parser.add_argument('antabfile', type=str, help=help_antabfile)

args = parser.parse_args()


args.msdata = args.msdata[:-1] if args.msdata[-1]=='/' else args.msdata



class Antab():
    self.antname = None
    self.control = None
    self.baseline = None
    self.gain = {}
    self.tsys = None


def read_antabfile(antab=args.antabfile):
    """Returns the
    """
    antennas = defaultdict(dict)
    header = {}
    values = None
    with open(antab) as antabfile:
        antablines = antabfile.readlines()
        for aline in antablines:
            # Ignore if it is a comment line
            if aline[0].strip() != '!':







