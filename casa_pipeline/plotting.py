#!/usr/bin/env python3
import datetime as dt
from typing import Optional, Union #, Iterable #, Tuple, NoReturn, List, Tuple
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .. import casa_pipeline as capi

# TODO: I will need to add the headermatplotlib style



class Jplotter(object):
    """Functions to produce quick plots with jiveplot (jplotter) on the given data.
    """
    def __init__(self, ms: capi.Project):
        self._ms = ms
        # TODO: It will need to run jplotter...
        if self._ms.msfile.exists():
            self.open_ms()


    def open_ms(self):
        yield f"ms {str(self._ms.msfile)}"
        yield "indexr"


    def autocorr(self, sources: Optional[list] = None, scans: Optional[list] = None):
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
        if scans is not None:
            yield f"scan {' '.join(scans)}"

        if sources is not None:
            yield f"src {' '.join(sources)}"

        # yield "scan mid-30s to mid+30s where src ~ {0}".format( scansel(settings.calsrc) )
        yield "new all false bl true sb false time true"
        yield "multi true"
        yield "sort bl"
        yield "nxy 2 4"
        yield "ckey p[rr]=2 p[ll]=3 p[rl]=4 p[lr]=5 p[none]=1"
        # yield "refile {0}-auto-{1}.ps/cps".format( settings.myBasename, num )
        yield "pl"
        print("done auto plots on calibrator scan")


    def uvplot(self, sources: Optional[list] = None, per_source: bool = True):
        """Creates a (u,v) coverage plot.
        """
        print("Generating uv plot.")
        yield "pt uv"
        yield "bl cross"
        yield "fq */p;ch none"
        yield "time none"
        yield "avt none"
        yield "avc none"
        yield "multi true"
        yield "ckey src src[none]=1"
        yield "new all false"
        if per_source:
            yield "new sb"

        if sources is not None:
            yield f"src ' '.join(sources)"

        yield "pl"
        print("(u,v) plot done.")


class Casaplot(object):
    """Functions to produce plots with the CASA-integrated tools
    """

    def casa_bandpass(self):
        """Plots the bandpass calibration
        """
        # https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.visualization.plotbandpass.html
        pass

    def casa_tsys(self):
        pass



class Plotting(object):
    """Functions to produce plots from a Ms object.
    """
    def __init__(self, ms: capi.Project):
        self._ms = ms
        self._jplotter = Jplotter(ms)

    def autocorr(self, outfile: Optional[str] = None, **kwargs):
        """Creates a plot
        """
        pass

    def crosscorr(self, outfile: Optional[str] = None, **kwargs):
        pass

    def anptime(self, outfile: Optional[str] = None, **kwargs):
        pass

    def amptime(self, outfile: Optional[str] = None, **kwargs):
        pass

    def phasetime(self, outfile: Optional[str] = None, **kargs):
        pass

    def caltable(self, calname: str, outfile: Optional[str] = None, **kargs):
        pass

    def radplot(self, mode='a&p', withmodel=False, outfile: Optional[str] = None,**kargs):
        pass

    def tplot(self, outfile: Optional[str] = None):
        """Creates a plot with the antennas in the _y_ axis and time in the _x_ axis,
        showing which antennas observed at each scan.
        """
        fig, ax = plt.subplots()
        ants_scans = self._ms.scans_in_antennas()
        # convert scans into time range
        cmap = plt.get_cmap("tab10")
        for s,a_src in enumerate(self._ms.sources.names):
            src_scans = self._ms.scans_with_source(a_src)
            if len(src_scans) > 0:
                for a,ant in enumerate(self._ms.antennas.names):
                    scans_with_src = np.intersect1d(ants_scans[ant], src_scans)
                    if len(scans_with_src) > 0:
                        x = self._ms.times_for_scans(scans_with_src)
                        ax.scatter(x=x, y=np.ones_like(x)*a, s=1, color=cmap(s), rasterized=True,
                                   label=a_src if a == 0 else None)


        ax.set_xlabel(r'Time (UTC)')
        ax.set_ylabel(r'Antenna Name')
        ax.set_yticks(range(len(self._ms.antennas)), self._ms.antennas.names)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        ax.legend()
        if outfile is None:
            fig.savefig(f"{self._ms.outdir}/plot-{self._ms.projectname}-tplot.pdf", dpi=330,
                        bbox_inches='tight', pad_inches=0.01)
        else:
            fig.savefig(outfile, dpi=330, bbox_inches='tight', pad_inches=0.01)

        plt.close()


    def uvplot(self, src: Optional[Union[list, str]] = None, outfile: Optional[str] = None):
        """Produces a (u,v) coverage plot for the given source.
        """
        pass


    def elevation(self, source: Optional[str] = None, outfile: Optional[str] = None):
        """Creates a plot with the source elevation along the time.
        If source specified, it only plots the elevation for the specified source.
        """
        pass

    def dynamic_range(self,  xaxis: str = 'time', yaxis: str = 'freq', zaxis: str = 'amp',
                      outfile: Optional[str] = None):
        """Creates a plot with the specified (x,y) axis, and with color the z axis.
        This

        Inputs
        ------
            - xaxis : str  (default = 'time')
                Variable to plot in the _x_ axis. Allowed values are 'time', 'uvdist', 'antenna'.
            - yaxis : str  (default = 'freq')
                Variable to plot in the _y_ axis. Allowed values are 'freq' or 'time'.
            - zaxis : str  (default = 'amp')
                Variable to use for the color scale of the plotted data. Allowed values are 'amp' or 'phase'.

        """
        pass




    # def phasetime(self, **kargs):
    #     pass

    # def phasetime(self, **kargs):
    #     pass
