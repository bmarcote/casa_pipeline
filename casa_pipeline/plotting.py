#!/usr/bin/env python3

from . import obsdata
from . import tools

class Plotting(object):
    """Functions to produce plots from a Ms object.
    """
    def __init__(self, ms: obsdata.Ms):
        self._ms = ms

    def jplotter(self, **kwargs):
        pass
