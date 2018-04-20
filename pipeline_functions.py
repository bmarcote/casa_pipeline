







def generate_gain_curve(antabfile):
    # It needs to execute:
    # casa --nologger -c ~/scripts/gc.py antabfile EVN.gc
    pass




def import_idi_files(idi_files, output_ms_file):
    """Convert the FITS-IDI files to a CASA visibility data set (MS format)
    """
    importfitsidi(fitsidifile=idi_files, vis=output_ms_file, constobsid=True,
                  scanreindexgap_s=15.0, specframe='GEO')


