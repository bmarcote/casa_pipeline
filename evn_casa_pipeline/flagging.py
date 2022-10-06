import casatasks
from . import tools


class Flagging(object):
    def __init__(self, ms_obj):
        self._ms = ms_obj

    def apply_flags_from_file(self, flagfile: str = None):
        """Applies the flags from a flagging file.
        If the file name is not specified, then it will assumed to be the associated
        a-priori flagging (projectname.flag).
        """
        flagfile = str(self._ms.cwd / f"{self._ms.project_name}.flag") if flagfile is None else flagfile
        casatasks.flagdata(vis=str(self._ms.msfile), inpfile=flagfile,
                           reason='any', action='apply', flagbackup=False, savepars=False)
        casatasks.flagmanager(vis=str(self._ms.msfile), mode='save', versionname=f"flags from {flagfile}",
                              comment='A-priori flags prior to any calibration.')

    def edge_channels(self, edge_fraction: float = 0.1):
        """Flangs the edge channels of each subband according to the specified edge_fraction.
        For example, 0.1 (default) will imply to flag the 10% of the channels next to the edge of each subband.

        Inputs
            edge_fraction : float [default = 0.1]
                Fraction of channels at the edge of each subband that will be flagged to all data.
        """
        edge_channel = int(self._ms.freqsetup.channels/(100*edge_fraction))
        start = str(edge_channel - 1)
        end = str(self._ms.freqsetup.channels - edge_channel)
        casatasks.flagdata(vis=str(self._ms.msfile), mode='manual', flagbackup=False,
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
            assert ant in self._ms.antennas.names, f"The antenna {ant} did not participate in this observation."

        casatasks.flagdata(vis=str(self._ms.msfile), mode='quack',
                           antennas=','.join(antennas), quackinterval=quack_interval_s, flagbackup=False)

    def flag_tfcrop(self, timecutoff=6.0, freqcutoff=5.0, timedevscale=5.0, freqdevscale=5.0):
        """Runs flagdata with the tfcrop option.
        In principle it should work fine for EVN data as long as the data are already calibrated.
        """
        casatasks.flagdata(vis=str(self._ms.msfile), mode='tfcrop', datacolumn='corrected',
                           field='', ntime='scan', timecutoff=timecutoff, freqcutoff=freqcutoff,
                           timefit='line', freqfit='line', flagdimension='freqtime', extendflags=True,
                           timedevscale=timedevscale, freqdevscale=freqdevscale, extendpols=False, growaround=False,
                           action='apply', flagbackup=True, overwrite=True, writeflags=True)

    def aoflagger(self, strategy_file: str = None):
        """Runs the AOFlagger on the associated MS.
        If you have a costumized AOflagger strategy file, you can use it.
        """
        if strategy_file is None:
            tools.shell_command("aoflagger", str(self._ms.msfile))
        else:
            tools.shell_command("aoflagger", ["-strategy", strategy_file, str(self._ms.msfile)])


