import os
from astropy.io import fits
import casatasks
from . import casavlbitools
from . import project



def apply_apriori_flagging(project: project.Project):
    """Applies the flags stored in the file associated with the project.
    """
    return
    casatasks.flagdata(vis=str(project.msfile), inpfile=str(project.caldir/f"{project.project_name}.flag"),
                       reason='any', action='apply', flagbackup=False, savepars=False)
    casatasks.flagmanager(vis=str(project.msfile), mode='save', versionname='apriori_flags',
                          comment='A-priori flags prior to any calibration.')



def edge_channels(project: project.Project, edge_fraction: float=0.1):
    """Flangs the edge channels of each subband according to the specified edge_fraction.
    For example, 0.1 (default) will imply to flag the 10% of the channels next to the edge of each subband.
    """
    edge_channel = int(project.nchannels/(100*edge_fraction))
    start = str(edge_channel-1)
    end = str(project.nchannels - edge_channel)
    spw2flag = f"*:0~{start};{end}~{str(project.nchannels - 1)}"
    casatasks.flagdata(vis=str(project.msfile), mode='manual', spw=spw2flag, flagbackup=False)


def quack(project: project.Project, quack_interval_s: float=5):
    """Flags the first seconds of each scan, defined by quack_interval (in seconds).
    """
    casatasks.flagdata(vis=str(project.msfile), mode='quack', quackinterval=quack_interval_s, flagbackup=False)



def flag_tfcrop(project: project.Project):
    """Runs flagdata with the tfcrop option.
    In principle it should work fine for EVN data as long as the data are already calibrated.
    """
    casatasks.flagdata(vis=str(project.Project), mode='tfcrop', datacolumn='corrected',
                       field=','.join(project.calibrators), ntime='scan', timecutoff=6.0, freqcutoff=5.0,
                       timefit='line', freqfit='line', flagdimension='freqtime', extendflags=True,
                       timedevscale=5.0, freqdevscale=5.0, extendpols=False, growaround=False,
                       action='apply', flagbackup=True, overwrite=True, writeflags=True)

# flagdata(mode='tfcrop') can be used when the data are calibrated (also in EVN in pprinciple).
# Good results at that stage (it just looks for outliers).

