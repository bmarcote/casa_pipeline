#!/usr/bin/env python3

from . import obsdata
from . import tools




class Jplotter(object):
    """Functions to produce quick plots with jiveplot (jplotter) on the given data.
    """
    def __init__(self, ms: obsdata.Ms):
        self._ms = ms

    def open_ms(self):
        yield f"ms {str(self._ms.msfile)}"
        yield "indexr"

    def autocorr(self, sources: list = None):
        """It will create auto-correlation plots for the given sources.

        sources : list
            The list of sources to be included to be plotted.
            NOTE: It can actually contain a different set of values:
                - list of str:  the names of the sources to be plot would then be expected.
                - list of Sources: the Source-type objects for the sources to plot.
                - list of SourceType:  the type of the sources to plot.
        """
        print("Generating auto-correlation plots (amp VS frequency)")
        yield "bl auto"
        yield "fq */p;ch none"
        yield "avt scalar;avc none"
        yield "time none"
        yield "pt ampchan"
        yield "y 0 2"
        # select the scan
        yield "scan mid-30s to mid+30s where {0}".format( scansel(settings.calsrc) )
        yield "new all false bl true sb false time true"
        yield "multi true"
        yield "sort bl"
        yield "nxy 2 4"
        yield "ckey p[rr]=2 p[ll]=3 p[rl]=4 p[lr]=5 p[none]=1"
        yield "refile {0}-auto-{1}.ps/cps".format( settings.myBasename, num )
        yield "pl"
        print("done auto plots on calibrator scan")


class Plotting(object):
    """Functions to produce plots from a Ms object.
    """
    def __init__(self, ms: obsdata.Ms):
        self._ms = ms

    def jplotter(self, **kwargs):
        pass
