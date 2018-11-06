# Add the directory containing the Mark's Python scripts to PYCAPATH and to PYTHONPATH
import glob

msfile = 'n17c3.ms'
idi_files = ['n17c3_2_1.IDI{}'.format(i) for i in range(1, 5)] # or glob.glob('n17_2_1.IDI*')
import_idi_files = True
refant = 'EF'
bpsource = 'J1848+3219'
target = bpsource
# SB 0-7, 32 MHz, 64 ch.
# Ef ( 0)  Mc ( 3)  Nt ( 4)  O8 ( 5)  T6 ( 6)  Ur ( 7)  Ys ( 9)
# Sv (10)  Zc (11)  Bd (12)  Hh (13)  Ir (14)  Cm (16)  Da (17)
# J1848+3219
# 25-Oct-2017/13:52:00.500 -> 25-Oct-2017/14:59:59.500  dT: 1.000
search_refants = ''
integration_time = '1s'
central_channels = '*:6~58'


input_uvflg_file = 'n17c3.flag'
caltable_tsys = 'cal.tsys'
caltable_gc = 'cal.gc'
caltables = {'tsys': caltable_tsys, 'gc': caltable_gc, 'fringe1': 'cal.fringe1',
             'fringe2': 'cal.fringe2', 'bpass': 'cal.bp'}


####################################################################################
# Things to do before entering CASA
####################################################################################

## Converting the uvflg to a CASA-like flag table
# python scripts/flag.py n17c3_2.uvflg n17c3_2_1.IDI1 > n17c3.flag

## Append the Tsys measurements to the FITS-IDI
# python scripts/append_tsys.py n17c3_2.antab n17c3_2_1.IDI*

## Generate the gain curve
# casa_vlbi -c scripts/gc.py n17c3_2.antab EVN.gc -l 20 -u 75



####################################################################################
# Things to run in CASA
####################################################################################


if import_idi_files:
    importfitsidi(fitsidifile=idi_files, vis=msfile, constobsid=True, scanreindexgap_s=15.0)

if input_uvflg_file is not None:
    flagdata(vis=msfile, mode='list', inpfile=input_uvflg_file, reason='any', tbuff=0.0,
         action='apply', display='', flagbackup=True, savepars=False)

# Flag autocorrelations
flagdata(vis=msfile, mode='manual', reason='autocorrelations', autocorr=True)


# Generate calibration Tsys table
gencal(vis=msfile, caltable=caltables['tsys'], caltype='tsys', uniform=False)
# Generate calibration Gain Curve table
gencal(vis=msfile, caltable=caltables['gc'], caltype='gc', infile='EVN.gc')


# Apply these corrections to be able to plot the a-priori corrected data
# applycal(vis=msfile, gaintable=[caltables['tsys'], caltables['gc']], parang=True)



## Inspection of the data
listobs(msfile)

# plotants(msfile, figfile='figs/ants.pdf')

# Probably all this part will be better done in jplotter. msplot is REALLY SLOW!!!!!!!
# plotms(msfile, xaxis='u', yaxis='v', coloraxis='field', plotfile='figs/uv-coverage.png', dpi=300, showgui=False)
# plotms(msfile, xaxis='time', yaxis='amp', coloraxis='spw', plotfile='figs/0-time-amp.png', dpi=300, showgui=False)



# Flag the edge channels. (~10% on either side)
# spw='*:0~5;58~63'
edge_channels = ''
# flagdata(vis=msfile, spw=edge_channels, autocorr=False, display='', flagbackup=True)

# Instrumental delay correction
scan_fringe = '' #'4' # ~14:20UT
fringe1_times = '14:07:00~14:08:00'
number_of_ifs = 8
# In AIPS: fring (dparm)
# no. bl comob. = 1
# delay win = 200 ns
# rate win = 50 mHz

fringefit(vis=msfile, caltable=caltables['fringe1'], scan=scan_fringe, timerange=fringe1_times, refant=refant,
          solint='inf', minsnr=3, zerorates=True, weightfactor=2, globalsolve=True,
          append=False, docallib=False, delaywindow=[], ratewindow=[],
          gaintable=[caltables['tsys'], caltables['gc']], gainfield=[], interp=[],
          spwmap=[], parang=True)


# applycal(vis=msfile, gaintable=[caltables['tsys'], caltables['gc'], caltables['fringe1']],
         # gainfield=[], interp=[], spwmap=[], calwt=[True], parang=True, flagbackup=True,
         # applymode='')

# Global fringe
fringefit(vis=msfile, caltable=caltables['fringe2'], scan='', refant=refant,
          solint='60s', minsnr=6, zerorates=False, globalsolve=True, weightfactor=2,
          delaywindow=[], ratewindow=[], append=False, docallib=False, combine='spw',
          gaintable=[caltables['tsys'], caltables['gc'], caltables['fringe1']], gainfield=[], interp=[],
          spwmap=[[], [], []], parang=True)


# applycal(vis=msfile, gaintable=[caltables['tsys'], caltables['gc'], caltables['fringe1'], caltables['fringe2']],
         # gainfield=[], interp=[], spwmap=[[],[],[],[0]*number_of_ifs], calwt=[True], parang=True, flagbackup=True,
         # applymode='')



# Band pass
bandpass(vis=msfile, caltable=caltables['bpass'], field=bpsource, solint='inf', refant=refant,
         gaintable=[caltables[i] for i in ('tsys', 'gc', 'fringe1', 'fringe2')], bandtype='B',
         spwmap=[[],[],[],[0]*number_of_ifs], parang=True)


# applycal(vis=msfile, gaintable=[caltables[i] for i in ('tsys', 'gc', 'fringe1', 'fringe2', 'bpass')],
#         gainfield=[], interp=[], spwmap=[[],[],[],[0]*number_of_ifs, []], calwt=[True], parang=True, flagbackup=True)

applycal(vis=msfile, gaintable=[caltables[i] for i in ('tsys', 'gc', 'fringe1', 'fringe2', 'bpass')],
        gainfield=[], interp=[], spwmap=[[],[],[],[0]*number_of_ifs, []], calwt=[True], parang=True, flagbackup=True)


plotms(vis=msfile)



splitfile = msfile[:-2] + bpsource + '.ms'
splitfileavg = msfile[:-2] + bpsource + '.avgchan.ms'

split(vis=msfile, outputvis=splitfile, field=bpsource, spw='*:5~58', datacolumn='corrected')
# Averaging all channels

split(vis=msfile, outputvis=splitfileavg, field=bpsource, spw='*:5~58', datacolumn='corrected', width=54)



exportuvfits(vis=splitfileavg, fitsfile=splitfile[:-2]+'uvfits', datacolumn='corrected', overwrite=True,
             multisource=False)



flagdata(vis=msfile, mode='manual', antenna='NT', reason='Nt no fringes')

tclean(vis=splitfile, imagename=splitfile[:-2]+'.clean', field=bpsource, spw='', imsize=2048, cell='0.05mas',
        weighting='briggs', robust=-2.0, interactive=False, niter=1000, interpolation='cubic', stokes='I',
        gain=0.1, specmode='mfs', nterms=1, deconvolver='hogbom', gridder='standard', savemodel='modelcolumn')



# Self-calibration
selfcalibration()











