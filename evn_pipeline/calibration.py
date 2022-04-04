import shutil
from pathlib import Path
import casatasks
import casatools
from casatasks.private import tec_maps
from . import project
from .project import SourceType


def write_callib(filename, entry):
    """Appends a new entry into the callib file that is used to apply the different calibration steps.
    filename : str
        The name of the callib file.
    entry : str
        The line of text including all parameters and values that should be stored into the callib.
    """
    with open(filename, 'a') as callib_file:
        callib_file.write(f"{entry}\n")


def remove_if_exists(a_table: Path):
    """Removes the table if exists.
    """
    # a_table should always be a Path, but just in case this converts it if it is only a string
    if not isinstance(a_table, Path):
        Path(a_table)
        print(f"WARNING: {a_table} is a str, not a Path (it should, but I am converting it now).")

    if a_table.exists():
        print(f"{a_table} exists and will be removed and overwriten.")
        shutil.rmtree(a_table)


def get_spw_global_fringe(caltable: str) -> int:
    """Returns the subband (spw) where the solutions from a global fringe run have been stored.
    Note that when combining spw in fringe, it will write the solutions in the first subband with
    data (so in general spw = 0; but not always).
    This function reads the generated calibration table and returns the number of the spw where
    the solutions where actually stored.
    """
    tb = casatools.table()
    tb.open(caltable)
    tb.open(tb.getkeyword('SPECTRAL_WINDOW').replace('Table: ', ''))
    the_spw = int(tb.getcol('MEAS_FREQ_REF')[0])
    tb.close()
    tb.close()
    return the_spw


def a_priori_calibration(project: project.Project):
    """Generates the calibration tables for gain and system temperature.
    Existing versions of the calibration will be removed first.
    """
    caltsys = project.caldir/f"{project.project_name.lower()}.tsys"
    calgc = project.caldir/f"{project.project_name.lower()}.gc"
    gcfile = project.caldir/"EVN.gc"
    remove_if_exists(caltsys)
    if not caltsys.exists():
        casatasks.gencal(vis=str(project.msfile), caltable=str(caltsys), caltype='tsys', uniform=False)
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{caltsys}' tinterp='nearest'")
    else:
        print(f"{caltsys} (Tsys corrections) already exists. It will not be generated again.")

    if not calgc.exists():
        casatasks.gencal(vis=str(project.msfile), caltable=str(calgc), caltype='gc', infile=str(gcfile))
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{calgc}' tinterp='nearest' finterp='nearest'")
    else:
        print(f"{calgc} (Gain Curve corrections) already exists. It will not be generated again.")


def ionospheric_corrections(project: project.Project):
    """Runs the ionospheric calibration.
    """
    # remove_if_exists(project.caldir/f"cal-iono.IGS_TEC.im")
    # remove_if_exists(project.caldir/f"cal-iono.IGS_RMS_TEC.im")
    return
    caltec = project.caldir / f"{project.project_name.lower()}.tec"
    imtec  = project.caldir / f"{project.project_name.lower()}"
    if not caltec.exists():
        tec_maps.create(vis=str(project.msfile), doplot=False, imname=str(imtec))
        casatasks.gencal(vis=str(project.msfile), caltable=str(caltec),
                         infile=f"{imtec}.IGS_TEC.im", caltype='tecim')
        write_callib(filename=str(project.caldir / "cal_library.txt"), entry=f"caltable='{caltec}'")
    else:
        print(f"{caltec} (ionospheric corrections) already exists. It will not be generated again.")



def main_calibration(project: project.Project):
    """Runs the main calibration of the data:
    - instrumental delay corrections: in the specified time range.
    - global fringe fit: on all calibrators (and target if no phase-referencing is used).
    - bandpass calibration: using the bandpass calibrators.
    """
    cals = {'sbd': (project.caldir/f"{project.project_name.lower()}.sbd"),
            'mbd': (project.caldir/f"{project.project_name.lower()}.mbd"),
            'bp': (project.caldir/f"{project.project_name.lower()}.bp")}
    do_all_exists = []
    for a_cal in cals:
        do_all_exists.append(cals[a_cal].exists())

    if not all(do_all_exists):
        casatasks.fringefit(vis=str(project.msfile), caltable=str(cals['sbd']),
                            timerange=project.instr_delay_timerange, solint='inf', zerorates=True,
                            refant=','.join(project.refants), minsnr=50, docallib=True,
                            callib=str(project.caldir/"cal_library.txt"), parang=True)
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{cals['sbd']}' tinterp='nearest'")
        casatasks.fringefit(vis=str(project.msfile), caltable=str(cals['mbd']),
                            field=','.join(project.calibrators), solint='inf', zerorates=False,
                            refant=','.join(project.refants), combine='spw', minsnr=5, docallib=True,
                            callib=str(project.caldir/"cal_library.txt"), parang=True)
        spw_with_solutions = get_spw_global_fringe(caltable=str(cals['mbd']))
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{cals['mbd']}' tinterp='linear' " \
                           f"spwmap={project.nsubbands}*[{spw_with_solutions}]")
        casatasks.bandpass(vis=str(project.msfile), field=','.join(project.bandpass_calibrators), combine='scan',
                           caltable=str(cals['bp']), docallib=True, callib=str(project.caldir/"cal_library.txt"),
                           solnorm=True, solint='inf', refant=','.join(project.refants), bandtype='B',
                           parang=True, fillgaps=16)
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{cals['bp']}' tinterp='linear' finterp='linear'")
    else:
        print(f"The tables related to the main calibration already exists. They will not be generated again.")



def recalibration(project: project.Project):
    """Runs a second calibration of the data:
    - instrumental delay corrections: in the specified time range.
    - global fringe fit: on all calibrators (and target if no phase-referencing is used).
    """
    cals = {'sbd': (project.caldir/f"{project.project_name.lower()}.sbd2"),
            'mbd': (project.caldir/f"{project.project_name.lower()}.mbd2"),
            'bp': (project.caldir/f"{project.project_name.lower()}.bp2")}

    do_all_exists = []
    for a_cal in cals:
        do_all_exists.append(cals[a_cal].exists())

    if not all(do_all_exists):
        casatasks.fringefit(vis=str(project.msfile), caltable=str(cals['sbd']),
                            timerange=project.instr_delay_timerange, solint='inf', zerorates=True,
                            refant=','.join(project.refants), minsnr=50, docallib=True,
                            callib=str(project.caldir/"cal_library.txt"), parang=True)
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{cals['sbd']}' tinterp='nearest'")
        casatasks.fringefit(vis=str(project.msfile), caltable=str(cals['mbd']), # gaintype='T',  # combine pols
                            field=','.join(project.calibrators), solint='inf', zerorates=False,
                            refant=','.join(project.refants), combine='spw', minsnr=5, docallib=True,
                            callib=str(project.caldir/"cal_library.txt"), parang=True)
        spw_with_solutions = get_spw_global_fringe(caltable=str(cals['mbd']))
        write_callib(filename=str(project.caldir/"cal_library.txt"),
                     entry=f"caltable='{cals['mbd']}' tinterp='linear' " \
                           f"spwmap={project.nsubbands}*[{spw_with_solutions}]")


def apply_calibration(project: project.Project):
    """Applies the current calibration by using the tables written in the callib.
    """
    casatasks.applycal(vis=str(project.msfile), docallib=True,
                       callib=str(project.caldir/"cal_library.txt"), parang=True)


def split(project: project.Project):
    """Splits all the data from all calibrated sources.
    """
    for a_source in project.sources:
        if project.sources[a_source].type != SourceType.other:
            casatasks.split(vis=str(project.msfile), outputvis=str(project.msfile).replace(".ms", f".{a_source}.ms"),
                            field=a_source, datacolumn='corrected', width=project.nchannels)



def export_uvfits(project: project.Project):
    """Export the available SPLIT files into a UVFITS file.
    """
    for a_source in project.sources:
        if project.sources[a_source].type != SourceType.other:
            casatasks.exportuvfits(vis=str(project.msfile).replace(".ms", f".{a_source}.ms"),
                                   fitsfile=str(project.msfile).replace(".ms", f".{a_source}.uvfits"),
                                   multisource=False, combinespw=True, padwithflags=True)


def globalfringe_with_model(project: project.Project, calibrator_models: list):
    """Once there is a good model for the phase calibrator, this function runs the global fringe
    using the model of the phase calibrators and extrapolates the solutions to the targets.
    """
    raise NotImplementedError
    calmbd = (project.caldir/f"{project.project_name.lower()}.mbd3")
    casatasks.fringefit(vis=str(project.msfile), caltable=calmbd, # TODO: add model here
                        field=','.join(project.calibrators), solint='inf', zerorates=False,
                        refant=','.join(project.refants), combine='spw', minsnr=5, docallib=True,
                        callib=str(project.caldir/"cal_library.txt"), parang=True)
    spw_with_solutions = get_spw_global_fringe(caltable=str(calmbd))
    write_callib(filename=str(project.caldir/"cal_library.txt"),
                 entry=f"caltable='{calmbd}' tinterp='linear' spwmap={project.nsubbands}*[{spw_with_solutions}]")


def self_calibration(project: project.Project, iteraction: int, solint: int):
    raise NotImplementedError
    # Maybe just a function that calls in loop all expected iteractions: p, p, p, a&p, p, a&p, p
    # In the amplitude selfcal (not the gscale with combine scan), use solnorm=True so it only corrects for
    # time-dependent gain residuals, not the flux scale.

def self_calibration_from_final_pcal_model(project: project.Project, src_model: Path):
    raise NotImplementedError
    # Assumes the model is a CASA image
    for a_cal in project.phase_calibrators:
        calsc = project.caldir / f"{project.project_name.lower()}.{a_cal}.sc_final"
        # TODO: check if ft is in casatools or casatasks
        casatools.ft(vis=str(project.msfile), field=a_cal, model=str(src_model), usescratch=True)
        casatasks.gaincal(vis=str(project.msfile), caltable=str(calsc), field=a_cal, solint='3min',
                          refant=','.join(project.refants), gaintype='G', calmode='ap', interp=['nearest'],
                          parang=True)
        casatasks.applycal()















