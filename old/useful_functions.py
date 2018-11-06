
from collections import defaultdict

# NOTE I may need to modify eval if for lists (with []) it does not recognize if there are
# strings without quote marks explicitily set

def evaluate(entry):
    """Evaluate the type of the string entry and convert it to its expected one.
    e.g. entry is a string containing a float. This function would return a float
    with that value (equivalent to eval(entry)).
    The enhancement is that eval_param_type can also evaluate if a list is present
    (eval does not work if no brackets are literally written, as in 1, 2, 3.)

    Useful function to evaluate all values stores by the ConfigParser (as it assumes
    that everything is a strings.
    """
    # If it is a tentative list
    if ',' in entry:
        if ('[' in entry) and (']' in entry):
            # If has brackets, eval works fine. NO, this sentence was a lie
            #return eval(entry)
            entry = entry.replace('[', '')
            entry = entry.replace(']', '')
            try:
                return [eval(i) for i in entry.split(',')]

            except NameError:
                # it would fail if i is a str.
                return [i.strip() for i in entry.split(',')]

        elif ('[' not in entry) == (']' in entry):
            # If only one is there, something is wrong
            raise ValueError('Format of {} unclear.'.format(entry))

        else:
            # it is a list without brackets
            try:
                return [eval(i) for i in entry.split(',')]

            except NameError:
                # it would fail if i is a str.
                return [i.strip() for i in entry.split(',')]

    else:
        return eval(entry)



def build_required_file_names(inputs, extension):
    """Returns a dictionary with all the files that are going to be created within the
    pipeline. That is:
        - msdata : file with all the raw uvdata (and/or calibration applied).
        Per each specified source (in bpass, phaseref, target, sources
            - uvdata : for the calibrated uvdata
            - dirty  : dirty image of the source
            - clean  : clean image of the source
            - selfcal{i} : i-iteraction of selfcalibrated image of the source.
    """
    files = defaultdict(dict)
    files['msdata'] = '{}/{}.{}'.format(inputs['outdir'], inputs['experiment'], extension)
    files['split'] = '{}/{}.split.{}'.format(inputs['outdir'],
                      inputs['experiment'], extension)
    for a_source in inputs['bpass']:
        files[a_source]['uvdata'] = '{}/{}_{}.uvdata.{}'.format(inputs['outdir'], \
                                    inputs['experiment'], a_source, extension)
        files[a_source]['clean'] = '{}/{}_{}.clean.{}'.format(inputs['outdir'], \
                                    inputs['experiment'], a_source, extension)
        files[a_source]['dirty'] = '{}/{}_{}.dirty.{}'.format(inputs['outdir'], \
                                    inputs['experiment'], a_source, extension)

        for sciter in range(inputs['sciter']):
            files[a_source]['selfcal'+str(sciter+1)] = \
                                    '{}/{}_{}.selfcal{}.{}'.format(inputs['outdir'],
                                    inputs['experiment'], a_source, sciter+1, extension)

    return files


def build_required_msfile_names(inputs):
    """Returns a dictionary with all the MS files that are going to be created within the
    pipeline.
    """
    return build_required_file_names(inputs, 'ms')

def build_required_fitsfile_names(inputs):
    """Returns a dictionary with all the MS files that are going to be created within the
    pipeline.
    """
    return build_required_file_names(inputs, 'fits')

def build_required_calfile_names(inputs):
    """Returns a dictionary with all the calibration files that are going to be created within
    the pipeline. That is:
        tsys = %(outdir)s/%(experiment)s_caltable.tsys
        gc = %(outdir)s/%(experiment)s_caltable.gc
        tec = $(outdir)s/%(experiment)s_caltable.tec
        fringe_instr = %(outdir)s/%(experiment)s_caltable.fringe_instrumental
        fringe = %(outdir)s/%(experiment)s_caltable.fringe
        bpass = %(outdir)s/%(experiment)s_caltable.bpass
        selfcali = %(outdir)s/%(experiment)s_caltable.selfcal{i}
    """
    files = {}
    files['tsys'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.tsys'
    files['gc'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.gc'
    files['tec'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.tec'
    files['fringe_instr'] = inputs['outdir']+'/'+inputs['experiment'] + \
                            '_caltable.fringe_instrumental'
    files['fringe'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.fringe'
    files['bpass'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.bpass'
    for sciter in range(inputs['sciter']):
        files['selfcal'+str(sciter+1)] = inputs['outdir']+'/'+inputs['experiment']+ \
                                         '_caltable.selfcal'+str(sciter+1)

    return files


def build_required_plotfile_names(inputs):
    """Returns a dictionary with all the calibration files that are going to be created within
    the pipeline. That is:
        tsys = %(outdir)s/%(experiment)s_caltable.tsys
        gc = %(outdir)s/%(experiment)s_caltable.gc
        fringe_instr = %(outdir)s/%(experiment)s_caltable.fringe_instrumental
        fringe = %(outdir)s/%(experiment)s_caltable.fringe
        bpass = %(outdir)s/%(experiment)s_caltable.bpass
        selfcali = %(outdir)s/%(experiment)s_caltable.selfcal{i}
    """
    files = defaultdict(dict)
    files['tsys'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.tsys.pdf'
    files['gc'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.gc.pdf'
    files['fringe_instr']['phase'] = inputs['outdir']+'/'+inputs['experiment'] + \
                            '_caltable.fringe_instrumental.phase.pdf'
    files['fringe_instr']['delay'] = inputs['outdir']+'/'+inputs['experiment'] + \
                            '_caltable.fringe_instrumental.delay.pdf'
    files['fringe']['phase'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.fringe.phase.pdf'
    files['fringe']['delay'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.fringe.delay.pdf'
    files['fringe']['rate'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.fringe.rate.pdf'
    files['bpass'] = inputs['outdir']+'/'+inputs['experiment']+'_caltable.bpass.pdf'
    for sciter in range(inputs['sciter']):
        files['selfcal'+str(sciter+1)] = inputs['outdir']+'/'+inputs['experiment']+ \
                                         '_caltable.selfcal'+str(sciter+1)+'.pdf'

    sources = inputs['bpass']
    if inputs['phaseref'] is not None:
        sources += inputs['phaseref'] + inputs['target']

    for a_source in sources:
        files[a_source]['clean'] = '{}/{}_1_{}.clean.{}'.format(inputs['outdir'], \
                                    inputs['experiment'], a_source, 'pdf')
        files[a_source]['dirty'] = '{}/{}_1_{}.dirty.{}'.format(inputs['outdir'], \
                                    inputs['experiment'], a_source, 'pdf')

    # Plots from the uncalibrated and calibrated data
    # files['uncal']['autocorr'] = inputs['outdir']+'/'+inputs['experiment']+'_possm_autocorr.pdf'
    # files['uncal']['autocorr'] = inputs['outdir']+'/'+inputs['experiment']+'_possm_autocorr.pdf'
    return files




def print_log_header(logger, title):
    """Prints a entry in the logging formatted as a heading with title.
    """
    s = '#'*80 + '\n'
    s += '###  '+title+'\n'
    s += '#'*80 + '\n'
    logger.info(s)



def get_info_from_ms(msfile):
    """Returns a dictionary with useful information from the observations in the MS,
    as the number of channels, number of subbands, rest frequency, integration time,
    max. baseline and max. resolution of the data.
    """
    msinfo = {}
    ms.open(msfile)
    msinfo['channels'] = ms.metadata().nchan(0) #ms.getspectralwindowinfo()['0']['NumChan']
    msinfo['subbands'] = ms.metadata().nspw() #len(ms.getspectralwindowinfo())
    # CHECK SPEED WITH ONE METHOD AND THE OTHER
    # msinfo['subbands'] = len(vishead(msfiles['msdata'], mode='get', hdkey='spw_name')[0])
    msinfo['reffreq'] = 5e9#ms.getspectralwindowinfo()[str(msinfo['subbands']//2)]['RefFreq'] # In Hz
    keys = ms.getscansummary().keys()
    a_key = 1
    while str(a_key) not in keys:
        a_key += 1

    msinfo['inttime'] = ms.getscansummary()[str(a_key)]['0']['IntegrationTime']
    del a_key

    # Getting the longest baseline
    msinfo['max_baseline'] = au.getBaselineExtrema(msfiles['msdata'])[0] # In meters
    msinfo['resolution'] = (2.44*(3e8/msinfo['reffreq'])/msinfo['max_baseline'])*180*3600*1e3/np.pi # in mas

    msinfo['sources'] = ms.metadata().fieldnames()
    ms.done()
    return msinfo




