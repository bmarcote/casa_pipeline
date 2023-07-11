import glob
import argparse
from pathlib import Path
from typing import Optional, Union

from AIPS import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData, AIPSImage

"""Part of this code uses or is based on the EVN Pypeline code written in ParselTongue
from JIVE.
"""



def get_antenna_number(uvdata, refant: str):
    """Returns the AIPS number of the given antenna name.
    """
    assert uvdata.exists()
    return [ant.upper() for ant in uvdata.antennas].index(refant.upper()) + 1


def max_table_no(uvdata, inext: str):
    """Returns the maximum number of the given AIPS table.
    """
    max_no = 0
    for table in uvdata.tables:
        if (table[1][-2:] == inext.upper()) and (table[0] > max_no):
            max_no = table[0]

    return max_no


def fitld(projectname: str, fitsidifiles: Optional[str] = None,
          uvfits: Optional[str] = None, aipsclass="UVDATA", digicor=-1):
    """Imports the associated FITS-IDI or UVFITS from a given project into AIPS.
    Note that one and only one of the 'fitsidifiles' or 'uvfits; parameters needs to be provided.

    Inputs
        - projectname : str  (max len <= 6)
            Name that will be associated with the created data.
        - fitsidifiles : str (optional)
            Preffix common to all the FITS-IDI files to be imported.
        - uvfits : str (optional)
            UVFITS file to be imported.

    Returns
        - uvdata : AIPSUVData
            AIPSUVData that will be created from the imported data. If exists, it will be replaced.
    """
    uvdata = AIPSUVData(projectname, aipsclass, 1, 1)
    if uvdata.exists():
        uvdata.zap(force=True)

    assert (fitsidifiles, uvfits).count(None) == 1, "One and only one of the 'fitsidifiles' or" \
            " 'uvfits' parameters need to be provided."
    fitld = AIPSTask('fitld')
    fitld.outdata = uvdata
    fitld.douvcomp = 1
    fitld.clint = 0.25
    fitld.digicor = digicor
    if fitsidifiles is not None:
        fitld.ncount = len(glob.glob(f"{fitsidifiles}*"))
        if fitld.ncount == 0:
            raise FileNotFoundError(f"Could not find files with the name {fitsidifiles}")

        fitld.infile = fitsidifiles
        fitld.doconcat = 1
    else:
        fitld.ncount = 1
        fitld.infile = uvfits
        fitld.doconcat = 1

    fitld()
    return uvdata


def indxr(uvdata, cparm=(0, 22, 0.25, 0)):
    """Runs INDXR in the associated UVDATA.
    """
    indxr = AIPSTask('indxr')
    indxr.indata = uvdata
    indxr.cparm[1:] = cparm
    indxr()


def tacop(uvdata: AIPSUVData, fromuvdata: AIPSUVData, inext: str, inver: int, ncount: int,
          outver: int):
    """Copies a table version from 'fromuvdata' into 'uvdata'.
    """
    tacop = AIPSTask('tacop')
    tacop.indata = fromuvdata
    tacop.outdata = uvdata
    tacop.inext = inext
    tacop.inver = inver
    tacop.outver = outver
    tacop.ncount = ncount
    tacop()


def ionos(uvdata: AIPSUVData, aparm=(1, 0, 1, 0.85, 56.7, 0.9782)):
    # Run tecor
    # Run clcal
    raise NotImplementedError


def fring(uvdata: AIPSUVData, calsour: list, solint: float, refant: list, docalib: int = 2,
          parmode: Optional[str] = None, aparm: list = None, dparm: list = None, snr: int = 0,
          gainuse: Optional[int] = None, flagver: Optional[int] = 0, snver: Optional[int] = None,
          doband: int = -1, bpver: int = -1, weightit: int = 1, model : Optional[AIPSImage] = None,
          timer: list = (0, 0, 0, 0, 0, 0, 0, 0), cmethod: str = 'dft', bchan: int = 0,
          echan: int = 0, antenna: list = [0]):
    """Runs FRING on the given data.
    Special functions:
    parmode : str 'mbd' or 'sbd' (default None)
        Automatically selects the aparm/dparm appropate for the instrumental delay fring ('sbd')
        or the multi-band fring ('mbd'). If None, it will pick the provided values for aparm/dparm.
    """
    fring = AIPSTask('fring')
    fring.indata = uvdata
    fring.calsour[1:] = calsour
    fring.timer[1:] = timer
    fring.antenna[1:] = antenna
    fring.docalib = docalib
    if gainuse is None:
        fring.gainuse = max_table_no(uvdata, 'CL')
    else:
        fring.gainuse = gainuse

    fring.doband = doband
    fring.bpver = bpver
    fring.bchan = bchan
    fring.echan = echan
    fring.cmethod = cmethod
    fring.weightit = weightit
    fring.solint = solint
    fring.refant = refant[0]
    fring.search[1:] = refant[1:]
    if snver is None:
        fring.snver = max_table_no(uvdata, 'SN') + 1
    else:
        fring.snver = snver

    if flagver is None:
        fring.flagver = max_table_no(uvdata, 'FG')
    else:
        fring.flagver = flagver

    if parmode == 'sbd':
        fring.aparm[1:] = (2, 0, 0, 0, 0, 1, 10 if snr == 0 else snr)
        fring.dparm[1:] = (1, 200, 50, 2, 0, 0, 1, 0, 1)
    elif parmode == 'mbd':
        fring.aparm[1:] = (2, 0, 0, 0, 1, 1, 3 if snr == 0 else snr, 0, 1)
        fring.dparm[1:] = (1, 200, 50, 2, 0, 0, 1, 0, 0)
    else:
        raise ValueError(f"'parmode' can only be 'sbd' or 'mbd' (or None). But is {parmode}.")

    if model is not None:
        fring.in2data = model

    fring()


def clcal(uvdata: AIPSUVData, calsour: list, refant: list, snver: int, gainver: int, gainuse: int,
          opcode: str = 'CALI', interpol: str = '2PT', sources: list = (''), samptype: str = '',
          doblank: int = 0, dobtween: int = 0):
    clcal = AIPSTask('clcal')
    clcal.indata = uvdata
    clcal.calsour[1:] = calsour
    clcal.sources[1:] = sources
    clcal.opcode = opcode
    clcal.interpol = interpol
    #clcal.intparm = 0.00001
    clcal.samptype = samptype
    clcal.doblank = doblank
    clcal.dobtween = dobtween
    clcal.refant = refant[0]
    clcal.snver = snver
    clcal.gainver = gainver
    clcal.gainuse = gainuse
    clcal()


def bpass(uvdata: AIPSUVData, calsour: list, refant: list, gainuse: int , bif: int = 0,
          eif: int = 0, flagver: int = 0, docalib: int = 1, solint: int = -1,
          soltype: str = 'L1R', bpver: int = 0, smooth: list = (1,),
          bpassprm: list = (0, 0, 0, 0, 1, 0, 0, 0, 1), weightit: int = 1):
    bpass = AIPSTask('bpass')
    bpass.indata = uvdata
    bpass.calsour[1:] = calsour
    bpass.bif = bif
    bpass.eif = eif
    bpass.flagver = flagver
    bpass.docalib = docalib
    bpass.gainuse = gainuse
    bpass.solint = solint
    bpass.refant = refant[0]
    bpass.soltype = soltype
    bpass.smooth[1:] = smooth
    bpass.bpassprm[1:] = bpassprm
    bpass()


def split(uvdata: AIPSUVData, sources: list = [''], bif: int = 0, eif: int = 0, bchan: int = 0,
          echan: int = 0, docal: int = 1, gainuse: int = 0, flagver: int = 0, bpver: int = 0,
          doband: int = 1, aparm: list = (2, 0, 0, 1), nchav: int = 0):
    split = AIPSTask('split')
    split.indata = uvdata
    split.sources[1:] = sources
    split.bif = bif
    split.eif = eif
    split.bchan = bchan
    split.echan = echan
    split.docal = docal
    split.gainuse = gainuse
    split.flagver = flagver
    split.doband = doband
    split.bpver = bpver
    split.aparm[1:] = aparm
    split.nchav = nchav
    split()


def fittp(uvdata: Union[AIPSUVData, AIPSImage], dataout: str):
    fittp = AIPSTask('fittp')
    fittp.indata = uvdata
    fittp.dataout = dataout
    fittp()


def apriori_calibration(aipsid: int, projectname: str,
                        fitsidifiles: Optional[str] = None,
                        uvfits: Optional[str] = None):
    """Runs the a-priori calibration in AIPS for VLBI data.
    """
    assert len(projectname) <= 6
    AIPS.userno = aipsid
    uvdata = fitld(projectname, fitsidifiles=fitsidifiles, uvfits=uvfits, aipsclass="UVDATA")
    indxr(uvdata)
    uvtasav = fitld(projectname, fitsidifiles=fitsidifiles, uvfits=uvfits, aipsclass="TASAV")
    tacop(uvdata, fromuvdata=uvtasav, inext='CL', inver=2, ncount=1, outver=2)
    tacop(uvdata, fromuvdata=uvtasav, inext='FG', inver=1, ncount=1, outver=1)
    # tecor(uvdata)
    return uvdata



def main_calibration(aipsid: int, projectname: str,
                     sbd_timerange: list, refant: list, calsour: list, solint: float,
                     target: list, bpsour : list, phaseref: Optional[list] = None,
                     import_uvfits: Optional[str, Path] = None, model: Optional[str, Path] = None,
                     bchan: int = 0, echan: int = 0):
    """bpsour are the fringe finders and calsour should contain all calibrators
    """
    assert len(projectname) <= 6
    assert (aipsid >= 100) and (aipsid < 100000)
    AIPS.userno = aipsid
    if import_uvfits is not None:
        uvdata = fitld(projectname, uvfits=import_uvfits, aipsclass="UVDATA")
    else:
        uvdata = AIPSUVData(projectname, "UVDATA", 1, 1)

    if model is not None:
        # model = AIPSImage()
        raise NotImplementedError("Using a model in fring is not implemented yet in main_calibration")

    for i in range(2):
        cl_last = max_table_no(uvdata, "CL")
        sn_last = max_table_no(uvdata, "SN")
        bp_last = max_table_no(uvdata, "BP")
        if len(sbd_timerange) > 4:
            fring(uvdata, calsour=[""], solint=10, refant=refant, parmode='sbd', gainuse=cl_last,
                  flagver=0, snver=sn_last+1, bpver = -1 if bp_last == 0 else bp_last, #model=model,
                  timer=sbd_timerange[:4])
            clcal(uvdata, calsour=calsour, refant=refant, snver=sn_last+1, gainver=cl_last,
                  gainuse=cl_last+1, opcode='CALP')
            # TODO Double-check this was the right way to do it
            fring(uvdata, calsour=[""], solint=10, refant=refant, parmode='sbd', gainuse=cl_last+1,
                  flagver=0, snver=sn_last+2, bpver = -1 if bp_last == 0 else bp_last, #model=model,
                  timer=sbd_timerange[4:])
            clcal(uvdata, calsour=calsour, refant=refant, snver=sn_last+2, gainver=cl_last+1,
                  gainuse=cl_last+2, opcode='CALP')
        else:
            fring(uvdata, calsour=[""], solint=10, refant=refant, parmode='sbd', gainuse=cl_last,
                  flagver=0, snver=sn_last+1, bpver = -1 if bp_last == 0 else bp_last, #model=model,
                  timer=sbd_timerange)
            clcal(uvdata, calsour=calsour, refant=refant, snver=sn_last+1, gainver=cl_last,
                  gainuse=cl_last+1, opcode='CALP')

        cl_last = max_table_no(uvdata, "CL")
        sn_last = max_table_no(uvdata, "SN")
        fring(uvdata, calsour=calsour, solint=solint, refant=refant, parmode='mbd', gainuse=cl_last,
              flagver=0, snver=sn_last+1, bpver = -1 if bp_last == 0 else bp_last) # model=model
        for a_source in calsour:
            clcal(uvdata, calsour=[a_source,], sources=[a_source,], refant=refant, snver=sn_last+1,
                  gainver=cl_last, gainuse=cl_last+1, opcode='CALI', interpol='self')

        if phaseref is None:
            # no phase reference!
            for a_target in target:
                assert a_target in calsour, \
                       f"There is no phaseref source. Target {a_target} should be in calsour."
        else:
            assert len(phaseref) in (1, len(target)), "The number of provided phase calibrators" \
                    " must match the number of target sources, or being only one."
            if len(phaseref) == 1:
                clcal(uvdata, calsour=phaseref, sources=target, refant=refant, snver=sn_last+1,
                      gainver=cl_last, gainuse=cl_last+1, opcode='CALI', interpol='ambg')
            else:
                for a_target, a_phaseref in zip(target, phaseref):
                    clcal(uvdata, calsour=[a_phaseref], sources=[a_target], refant=refant,
                          snver=sn_last+1, gainver=cl_last, gainuse=cl_last+1, opcode='CALI',
                          interpol='ambg')

        if i == 0:
            bpass(uvdata, gainuse=cl_last+1, calsour=bpsour, refant=refant)

    split(uvdata, sources=[''], bchan=bchan, echan=echan, gainuse=cl_last+1, doband=1,
          bpver=bp_last, aparm=(2, 0, 0, 1))
    uvsplits = []
    for a_source in list(set(target + bpsour + calsour)):
        uv_src = AIPSUVData(a_source, "SPLIT", 1, 1)
        if uv_src.exists():
            fittp(uv_src, f"PWD:{projectname}.{a_source}.SPLIT.UVFITS")
            uvsplits.append(uv_src)

    return uvsplits


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Runs a ParselTongue AIPS calibration.",
                                     prog="parseltongue.py", usage="%(prog)s [-h] [options]")
    parser.add_argument('aipsid', type=int, help="AIPS ID to use for the processing.")
    parser.add_argument('projectname', type=str, help="Codename of the project.")
    parser.add_argument('-i', '--initial', default=False, action="store_true",
                        help="Runs the a-priori calibration.")
    parser.add_argument('-c', '--calib', default=False, action="store_true",
                        help="Runs the main calibration.")
    parser.add_argument('--fits', type=str, default=None, help="FITS-IDI root names to be " \
                        "imported (following AIPS naming as used in FITLD).")
    parser.add_argument('--uvfits', type=str, default=None, help="UVFITS to be imported.")
    parser.add_argument('--target', type=str, default=None,
                        help="comma-separated list of target sources")
    parser.add_argument('--phaseref', type=str, default=None,
                        help="comma-separated list of phase-referencing sources")
    parser.add_argument('--fringefinder', type=str, default=None,
                        help="comma-separated list of fringe-finder sources")
    parser.add_argument('--sbdtime', type=str, default=None, help="Time in AIPS format " \
                        "(comma separated, no spaces) for the time range to use for the " \
                        "intrumental delay correction.")
    parser.add_argument('--refant', type=str, default=None, help="Comma-separated list " \
                        "for reference antennas")
    parser.add_argument('--solint', type=float, default=5, help="solint for the global fringe.")
    parser.add_argument('--bchan', type=int, default=0, help="BCHAN parameter for split.")
    parser.add_argument('--echan', type=int, default=0, help="ECHAN parameter for split.")
    parser.add_argument()

    args = parser.parse_args()

    assert (args.aipsid >= 100) and (args.aipsid < 100000)

    if args.initial:
        uvfiles = apriori_calibration(args.aipsid, args.projectname, fitsidifiles=args.fits,
                            uvfits=args.uvfits)

    if args.calib:
        calsources = []
        for a_src in (args.target, args.phaseref, args.fringefinder):
            if a_src is not None:
                calsources += a_src.split(',')

        uvfiles = main_calibration(args.aipsid, args.projectname, sbd_timerange=args.sbdtime,
                                   refant=args.refant, calsour=calsources, solint=args.solint,
                                   target=args.target.split(','), bpsour=args.fringefinder.split(','),
                                   phaseref=args.phaseref.split(','),
                                   import_uvfits=args.uvfits if args.initial is False else None,
                                   model=None, bchan=args.bchan, echan=args.echan)
        for uvfile in uvfiles:
            outuvfile = Path(f"{args.projectname}.{uvfile.name}.SPLIT.UVFITS")
            if outuvfile.exists():
                outuvfile.unlink()

            fittp(uvfile, f"PWD:{str(outuvfile)}")
    else:
        fittp(uvfiles, f"PWD:{args.projectname}.UVFITS")






