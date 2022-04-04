import os
from astropy.io import fits
import casatasks
from . casavlbitools import casa
from .casavlbitools import fitsidi
from .project import Project


def append_tsys(project: Project):
    """Appends the system temperature (Tsys) information to the FITS-IDI files.
    If the information is already present in the files, it will do nothing.
    Assumes that the Tsys information is available in a {project_name.lower()}.antab file in the same directory.
    """
    hdulist = fits.open(project.idi_files[0])
    if 'SYSTEM_TEMPERATURE' in hdulist:
        print('TSYS table already present in FITS-IDI files, skipping the append step.')
    else:
        print('Appending TSYS, this takes some time.')
        fitsidi.append_tsys(f"{project.project_name.lower()}.antab", project.idi_files)
        print('TSYS appended to the FITS-IDI files.')


def generate_gain_curve(project: Project):
    """Converts the gain curve for each antenna into a CASA readable gain curve table.
    Assumes that the Tsys information is available in a {project_name.lower()}.antab file in the same directory.
    """
    gctable = project.caldir / 'EVN.gc'
    if gctable.exists():
        print('GC table already created, skipping the generation of a new gain curve.')
    else:
        casa.convert_gaincurve(f"{project.project_name.lower()}.antab", str(gctable),
                               min_elevation=5.0, max_elevation=90.0)

def convert_uvflg(project: Project):
    """Converts the AIPS uvflag table into a CASA-compatible format.
    Assumes that the uvflg information is available in a {project_name.lower()}.uvflg file in the same directory.
    """
    flagfile = project.caldir / f"{project.project_name.lower()}.flag"
    if flagfile.exists():
        print('Flag file already created, skipping the conversion of the uvflg file.')
    else:
        fitsidi.convert_flags(f"{project.project_name.lower()}.uvflg", project.idi_files,
                              outfile=str(flagfile))


def importfitsidi2ms(project: Project):
    """Converts the existing IDI-FITS files into a MS file.
    """
    if project.msfile.exists():
        print(f"The MS file associated to {project.project_name} already exists. Will not be replaced.")
    else:
        casatasks.importfitsidi(fitsidifile=project.idi_files, vis=str(project.msfile), constobsid=True,
                                scanreindexgap_s=8.0, specframe='GEO')


def get_ms_info(project: Project):
    project.get_ms_info()


def listobs(project: Project):
    project.get_listobs()

