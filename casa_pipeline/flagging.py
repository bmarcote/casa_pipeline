import casatasks
from pathlib import Path
from typing import Optional, Iterable, NoReturn, List, Union, Tuple
from rich import print as rprint
import casa_pipeline as capi

class Flagging(object):
    """Class that contains all tasks concerning flagging to the data that can be applied to
    a Ms object.
    """
    def __init__(self, ms: capi.Project):
        self._ms = ms


    def apply_flags_from_file(self, flagfile: Optional[Union[Path, str]] = None):
        """Applies the flags from a flagging file.
        If the file name is not specified, then it will assumed to be the associated
        a-priori flagging (projectname.flag).
        """
        if isinstance(flagfile, str):
            flagfile = Path(flagfile)

        flagfile = (self._ms.cwd / f"{self._ms.projectname.lower()}.flag") if flagfile is None \
                   else flagfile

        if not flagfile.exists():
            rprint(f"[bold red]The flag file {flagfile} should exist but cannot be found.[/bold red]")
            raise FileNotFoundError(f"The flagfile {flagfile} cannot be found.")

        results = casatasks.flagdata(vis=str(self._ms.msfile), mode='list', inpfile=str(flagfile),
                           reason='any', action='apply', flagbackup=False, savepars=False)
        casatasks.flagmanager(vis=str(self._ms.msfile), mode='save',
                              versionname=f"flags from {flagfile}",
                              comment='A-priori flags prior to any calibration.')
        return results


    def flagdata(self, **kwargs):
        """Runs the CASA tasks flagdata with all given parameters.
        It just saves you the option of writting the MS name.
        """
        return casatasks.flagdata(vis=str(self._ms.msfile), **kwargs)


    def autocorrelations(self):
        """Flags the auto-correlations from the data.
        """
        return casatasks.flagdata(vis=str(self._ms.msfile), mode='manual', autocorr=True)


    def edge_channels(self, edge_fraction: float = 0.1):
        """Flags the edge channels of each subband according to the specified edge_fraction.
        For example, 0.1 (default) will imply to flag the 10% of the channels next to the edge
        of each subband.

        Inputs
            edge_fraction : float [default = 0.1]
                Fraction of channels at the edge of each subband that will be flagged to all data.
        """
        edge_channel = int(self._ms.freqsetup.channels/(100*edge_fraction))
        start = str(edge_channel - 1)
        end = str(self._ms.freqsetup.channels - edge_channel)
        return casatasks.flagdata(vis=str(self._ms.msfile), mode='manual', flagbackup=False,
                           spw=f"*:0~{start};{end}~{str(self._ms.freqsetup.channels - 1)}")


    def quack(self, antenna: str, quack_interval_s: float):
        """Flags the first seconds of data for each scan, defined by quack_interval (in seconds),
        for the given antenna(s).

        Inputs
            antenna : str
                Antenna(s) to which apply the quack flagging. If more than one, they can be
                comma or space separated (within the same string).
            quack_interval_s : float

        It May Raise
            AssertionError : if the specified antenna is not present in the observation.
        """
        antennas = [ant.strip() for ant in antenna.split(',' if ',' in antenna else ' ')]
        for ant in antennas:
            assert ant in self._ms.antennas.names, \
                   f"The antenna {ant} did not participate in this observation."

        return casatasks.flagdata(vis=str(self._ms.msfile), mode='quack',
                           antenna=','.join(antennas), quackinterval=quack_interval_s,
                           flagbackup=False)


    def tfcrop(self, timecutoff=6.0, freqcutoff=5.0, timedevscale=5.0, freqdevscale=5.0):
        """Runs flagdata with the tfcrop option.
        In principle it should work fine for EVN data as long as the data are already calibrated.
        """
        return casatasks.flagdata(vis=str(self._ms.msfile), mode='tfcrop', datacolumn='corrected',
                           field='', ntime='scan', timecutoff=timecutoff, freqcutoff=freqcutoff,
                           timefit='line', freqfit='line', flagdimension='freqtime',
                           extendflags=True, timedevscale=timedevscale, freqdevscale=freqdevscale,
                           extendpols=False, growaround=False,
                           action='apply', flagbackup=True, overwrite=True, writeflags=True)


    def aoflagger(self, strategy_file: str = None):
        """Runs the AOFlagger on the associated MS.
        If you have a costumized AOflagger strategy file, you can use it.
        """
        if strategy_file is None:
            capi.tools.shell_command("aoflagger", [str(self._ms.msfile)])
        else:
            capi.tools.shell_command("aoflagger", ["-strategy", strategy_file, str(self._ms.msfile)])


    def check_unflagged_data(self):
        """Checks if there are still unflagged data.
        It will raise an error if
        """
        raise NotImplementedError

    def lsflagging(self, timeint_min: float = 1):
        """Low-Statistics Flagging scans the visibilities in the MS for all calibrator sources
        (or target if no phase-referencing is set), and analyzes the statistics along the
        different polarizations. It expects that the cross-polarizations should always exhibit
        a lower level of amplitudes than the parallel-polarizations.
        Otherwise it assumes that the antenna did not show fringes at those times.
        """
        raise NotImplementedError



