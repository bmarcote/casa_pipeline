import os
import glob
import shutil
import string
from pathlib import Path
import subprocess
import datetime as dt
from astropy.coordinates import SkyCoord
import astropy.units as u
from dataclasses import dataclass
from typing import Optional, Union
from rich import print as rprint
import casatasks
from casatasks.private import tec_maps
from casatools import table as tb
from casatools import componentlist as cl
# from casatools import msmetadata as msmd
import casa_pipeline as capi

_FILE_DIR = os.path.dirname(os.path.realpath(__file__))
"""Diferent functions that can be used directly from this module in order to calibrate the data.
Later on a class Calibration is defined, which will use these functions but by using the
obsdata/Ms objects.
"""
@dataclass
class CallibEntry:
    """Defines an entry in a calibration tables library.
    It contains the following parameters:
        - name : str
          Label that refers to which step did this calibration entry.
        - caltable : str
          Filename of the associated calibration table.
        - parameters : dict
          Parameters associated to the calibration table to used during applycal().
          It, in general, may contain the keys 'tinterp', 'finterp', and/or 'spwmap'.
    """
    name: str
    caltable: str
    parameters: dict


class CallibEntries():
    def __init__(self):
        self._entries = []

    def add(self, new_entry: CallibEntry):
        assert isinstance(new_entry, CallibEntry), f"The 'new_entry' (currently {new_entry}) " \
                                                   "should be a CallibEntry object."
        if new_entry.name in self.names:
            self.__setitem__(new_entry.name, new_entry)
        else:
            self._entries.append(new_entry)

    def delete(self, entry_name: str):
        assert entry_name in self.names, \
            f"The specified entry name ({entry_name}) is not in listed in {self.__str__}."
        self.__delitem__(entry_name)

    @property
    def names(self):
        return [entry.name for entry in self._entries]

    @property
    def caltables(self):
        return [entry.caltable for entry in self._entries]

    @property
    def parameters(self):
        return [entry.parameters for entry in self._entries]

    def get_dict(self):
        """Returns a dictionary where the keys are the names of the CallibEntry,
        and the values are the associated CallibEntry object.
        """
        return {entry.name: entry for entry in self._entries}

    def __len__(self):
        return len(self._entries)

    def __getitem__(self, entry_name: str):
        return self._entries[self.names.index(entry_name)]

    def __setitem__(self, name: str, value: CallibEntry):
        self._entries[self.names.index(name)] = value

    def __delitem__(self, entry_name: str):
        return self._entries.pop(self.names.index(entry_name))

    def __iter__(self):
        return self._entries.__iter__()

    def __reversed__(self):
        return self._entries[::-1]

    def __contains__(self, entry_name: str):
        return entry_name in self.names

    def __str__(self):
        return f"CallibEntries([{', '.join(self.names)}])"

    def __repr__(self):
        return f"CallibEntries([{', '.join(self.names)}])"



class Callib(object):
    @property
    def filename(self) -> Path:
        return self._filename

    @filename.setter
    def filename(self, new_filename: Union[Path, str]):
        if isinstance(new_filename, str):
            new_filename = Path(new_filename)

        assert isinstance(new_filename, Path)
        self._filename = new_filename

    @property
    def entries(self) -> CallibEntries:
        return self._entries

    def __init__(self, filename: Union[Path, str]):
        if isinstance(filename, str):
            filename = Path(filename)

        assert isinstance(filename, Path), \
               f"The file name {filename} must be a Path- or str-type object."
        self.filename = filename
        self._entries = CallibEntries()


    def associated_caltable(self, name: str):
        return self._entries[name].caltable

    def associated_parameters(self, name: str):
        return self._entries[name].parameters

    def names(self) -> list:
        return self.entries.names

    def gaintables(self) -> list:
        return self.entries.caltables

    def parameters(self) -> list:
        return self.entries.parameters

    def interps(self) -> list:
        """Returns a list with the interpolations to apply for each calibration table.
        This will have the format expected to run a calibration step in CASA when calling
        the 'interp' parameter.
        """
        interp_params = []
        for acal in self.entries.parameters:
            interp_params.append(f"{acal['tinterp'] if 'tinterp' in acal else 'linear'}," \
                                 f"{acal['finterp'] if 'finterp' in acal else 'linear'}")

        return interp_params


    def spwmaps(self) -> list:
        """Returns a list with the spwmap to apply for each calibration table.
        """
        spwmap_pars = []
        for acal in self.entries.parameters:
            spwmap_pars.append([] if 'spwmap' not in acal else acal['spwmap'])

        return spwmap_pars


    def new_entry(self, name: str, caltable: str, parameters: dict):
        """Adds a new callib entry to the list of associated entries.
        If there is already an entry with the same name, it will replace it.
        """
        self.entries.add(CallibEntry(name, caltable, parameters))
        self.update_file()

    def remove_entry(self, name: str, files=True):
        """Removes the entry in the callib collection with the given name.

        Inputs:
            name : str
                Name of the associated CallibEntry to be removed.
                If 'all', it will empty the list.
            files : bool [default = False]
                If True, it will also remove the calibration files associated to the callib entry.
                Otherwise, it only removes them from the list.
        """
        if name == 'all':
            if files:
                for caltable in self._entries.caltables:
                    remove_if_exists(caltable)

            self._entries = CallibEntries()
        else:
            if name in self._entries:
                if files:
                    remove_if_exists(self._entries[name].caltable)

                self._entries.__delitem__(name)

        self.update_file()


    def backup(self, filename: Union[str, Path, None] = None, replace: bool = False):
        """Makes a backup of the current callib file, copying its content into another file
        called filename.

        Inputs:
            filename : str/Path (default = None)
                Name of the associated file that will be created as a copy of the current
                callib file. If None, it will append (an incremental) number to the current
                filename of the callib.
                That is, if no previous backup exists, it will be 'callib_file_name.1'.
                Later backups will be called '*.2', '*.3', etc.
            replace : bool (default = False)
                If 'filename' is provided, and exists, then this option defines if the file
                should be overwriten. If 'filename' is None, then this option has no effect.
        """
        i = 1
        if filename is not None:
            if isinstance(filename, str):
                filename = Path(filename)

            if filename.exists():
                if replace:
                    filename.unlink()
                else:
                    raise FileExistsError(f"The file {filename} exists and the replacing option " \
                                          "has not been set.")
        else:
            filename = Path(str(self.filename) + f".{i}")
            while filename.exists():
                i += 1
                filename = Path(str(self.filename) + f".{i}")

        with open(self.filename, 'r') as fin, open(filename, 'w') as fout:
            fout.write(fin.read())

        return filename


    def update_file(self):
        with open(self.filename, 'w') as callib_file:
            for a_step in self.entries:
                callib_file.write(f"# {a_step.name}\n")
                callib_file.write(f"caltable='{a_step.caltable}' ")
                callib_file.write(' '.join([f"{k}='{a_step.parameters[k]}'" if k != 'spwmap' \
                      else f"{k}=[{','.join(str(kk) for kk in a_step.parameters[k])}]" \
                           for k in a_step.parameters]) + '\n')
                # This guarantees no spaces in the spwmap print list. Easier when importing


    def import_from_file(self, filename: Optional[Union[Path, str]] = None):
        if filename is None:
            filename = self.filename

        self._entries = CallibEntries()
        with open(self.filename, 'r') as callib_file:
            ongoing = False
            for aline in callib_file.readlines():
                if len(aline) > 0 and aline.strip()[0] == '#':
                    entry_name, ongoing = aline.strip()[1:].strip(), True
                elif ongoing:
                    entry_caltable, *entry_params = aline.split(' ')
                    entry_caltable = entry_caltable.strip().replace("'", "").split('=')[1]
                    params = {}
                    for an_entry in entry_params:
                        key, val = an_entry.strip().replace("'", "").split('=')
                        if key == 'spwmap':
                            params[key] = \
                                [int(v) for v in val.replace('[', '').replace(']', '').split(',')]
                        else:
                            params[key] = val

                    self.entries.add(CallibEntry(entry_name, entry_caltable, params))
                    ongoing = False


def remove_if_exists(a_table: Union[Path, str]):
    """Removes the table if exists.
    """
    # a_table should always be a Path, but just in case this converts it if it is only a string
    if not isinstance(a_table, Path):
        a_table = Path(a_table)

    if a_table.exists():
        print(f"{a_table} exists and will be removed.")
        shutil.rmtree(a_table)


def get_spw_global_fringe(caltable: Path) -> int:
    """Returns the subband (spw) where the solutions from a global fringe run have been stored.
    Note that when combining spw in fringe, it will write the solutions in the first subband with
    data (so in general spw = 0; but not always).
    This function reads the generated calibration table and returns the number of the spw where
    the solutions where actually stored.
    """
    a_table = tb(str(caltable))
    try:
        a_table.open(str(caltable))
        return a_table.getcol('SPECTRAL_WINDOW_ID')[0]
    finally:
        a_table.close()



class Calibration(object):
    """Class containing the calibration tools required for VLBI data reduction.
    """

    @property
    def caldir(self) -> Path:
        """Returns the Path to the directory where all calibration tables are stored.
        """
        return self._caldir

    @property
    def callib(self) -> Callib:
        """Returns the Callib file associated to the current calibration.
        """
        return self._callib


    @property
    def aips(self):
        """Access the calibration steps to be performed in AIPS through ParselTongue.
        """
        return self._aips


    def _verify(self, files: list):
        """Verifies that the given list of files (must contain a list of Path objects), exists.
        If they cannot be found, it will raise FileNotFoundError.
        """
        for a_file in files:
            if not a_file.exists():
                rprint("[bold red] -- Error --[/bold red]")
                rprint(f"[bold red]{a_file} was not created during " \
                       "a_priori_calibration.[/bold red]")
                self._ms.logger.error(f"{a_file} was not created during a_priori_calibration.")
                raise FileNotFoundError(f"The calibration table {a_file} should have been "
                                        "generated but does not exist.")

        return True


    def __init__(self, ms: capi.Project, caldir: Union[str, Path]):
        self._ms = ms
        self._caldir = caldir
        self._callib = Callib(caldir / f"callib-{self._ms.projectname}.txt")
        self._sbd_timerange = None
        self._aips = Aips(self._ms)


    # def copy_pols(self, antenna: Union[str, obsdata.Antenna], bad_pol: str, new_pol: str):
    #     """Copies the information from 'new_pol' into 'bad_pol'
    #     """
    #     raise NotImplementedError


    def a_priori_calibration(self, replace=False):
        """Generates the calibration tables for gain and system temperature.
        Existing versions of the calibration will be removed first if replace is True.
        """
        caltsys = self.caldir / f"{str(self._ms.projectname).lower()}.tsys"
        calgc = self.caldir / f"{str(self._ms.projectname).lower()}.gc"

        self._ms.logger.info(f"Running calibration.a_priori_calibration (replace={replace})")
        if replace:
            self.callib.remove_entry('tsys', files=True)
            self.callib.remove_entry('gc', files=True)

        if not caltsys.exists():
            print(f"Generating the Tsys calibration table {caltsys}.")
            casatasks.gencal(vis=str(self._ms.msfile), caltable=str(caltsys),
                             caltype='tsys', uniform=False)
            self.callib.new_entry(name='tsyscal', caltable=str(caltsys), \
                                  parameters={"tinterp": 'nearest', "finterp": 'nearest'})
            print(f"Tsys calibration table {caltsys} properly generated.")
            self._ms.logger.debug(f"New Calibration Table: {caltsys}.")
        else:
            print(f"{caltsys} (Tsys corrections) already exists. It will not be generated again.")

        if not calgc.exists():
            print(f"Generating GC calibration table {calgc}.")
            casatasks.gencal(vis=str(self._ms.msfile), caltable=str(calgc), caltype='gc')
            self.callib.new_entry(name='gccal', caltable=str(calgc),
                                  parameters={"tinterp": 'nearest'})
            print(f"GC calibration table {calgc} properly generated.")
            self._ms.logger.debug(f"New Calibration Table: {calgc}.")
        else:
            print(f"{calgc} (Gain Curve corrections) already exists. "
                  "It will not be generated again.")

        self._verify([caltsys, calgc])


    def accor(self, replace=False, **kwargs):
        """NOTE: This should never run on EVN data!!!!!!!!!!!!!!!!!!!!!!!!!

        It is written here for compatibility and internal tests.
        """
        accor_caltable = self.caldir / f"{str(self._ms.projectname).lower()}.accor"
        accor_caltable_smooth = self.caldir / f"{str(self._ms.projectname).lower()}.accor_smooth"

        if replace:
            remove_if_exists(accor_caltable)
            self.callib.remove_entry('accor', files=True)

        if not accor_caltable_smooth.exists():
            print(f"Generating ACCOR calibration table {accor_caltable_smooth}.")
            casatasks.accor(vis=str(self._ms.msfile), caltable=str(accor_caltable), solint='30s',
                            docallib=True, callib=str(self.callib.filename),
                            **self._ms.params['calibration']['accor'])
            casatasks.smoothcal(vis=str(self._ms.msfile), tablein=accor_caltable,
                                caltable=accor_caltable_smooth, smoothtype='median',
                                smoothtime=1800.0)
            self.callib.new_entry(name='accor', caltable=str(accor_caltable_smooth),
                                  parameters={"tinterp": 'nearest', "finterp": 'linear'})
            print(f"ACCOR calibration table {accor_caltable_smooth} properly generated.")
        else:
            print(f"{accor_caltable_smooth} (ACCOR corrections) already exists. " \
                  "It will not be generated again.")

        self._verify([accor_caltable_smooth, ])


    def ionospheric_corrections(self, replace=False):
        """Runs the ionospheric calibration if required (i.e. if the observation was conducted
        at less than 6 GHz). This can also be turned off in the input parameter files.
        """
        caltec = self.caldir / f"{str(self._ms.projectname).lower()}.tec"
        imtec  = self.caldir / f"{str(self._ms.projectname).lower()}"
        if replace:
            remove_if_exists(caltec)
            self.callib.remove_entry('teccor', files=True)

        if not caltec.exists():
            print("Calibrating the ionospheric TEC influence.")
            tec_maps.create(vis=str(self._ms.msfile), doplot=False, imname=str(imtec))
            casatasks.gencal(vis=str(self._ms.msfile), caltable=str(caltec),
                             infile=f"{imtec}.IGS_TEC.im", caltype='tecim')

            # spw_with_solutions = get_spw_global_fringe(caltable=str(caltec))
            self.callib.new_entry(name='teccor', caltable=str(caltec), parameters={})
                                  # parameters={"spwmap": \
                                  #       self._ms.freqsetup.n_subbands*[spw_with_solutions]})
            print(f"Ionospheric TEC correction {caltec} properly generated.")
        else:
            print(f"{caltec} (ionospheric corrections) already exists. "
                  "It will not be generated again.")

        self._verify([caltec])


    def get_sbd_timerange(self, interval='2min') -> str:
        """Localizes a good time to compute the ionospheric corrections and returns it.
        It goes through the fringe-finder sources, and then through the scans from these
        ones selecting the time interval right in the middle of the scan or at the end of the same.
        Then it will do a basic statistics you select the best interval, by the following checks:
          1) the one containing all antennas
          2) the one with a smaller rms of the rms on each baseline.
          3) the one with higher amplitudes in all baselines (checking the min amplitude).
        """
        #TODO: implement or pass it by the user
        if 'calibration' in self._ms.params and 'sbd_timerange' in self._ms.params['calibration']:
            self._sbd_timerange = self._ms.params['calibration']['sbd_timerange']
            self._ms.logger.debug(f"sbd time range (from inputs): {self._sbd_timerange}.")
            assert self._sbd_timerange is not None, "'calibration' > 'sbd_timerange' in the input" \
                                                " params file should contain a time range."
            return self._sbd_timerange
        elif self._sbd_timerange is not None:
            self._ms.logger.debug(f"sbd time range (stored): {self._sbd_timerange}.")
            return self._sbd_timerange
        else:
            raise NotImplementedError
            # TODO: check scans with all antennas on the fringe finder
            # try:
            #     m = msmd(str(self._ms.msfile))
            #     if not m.open(str(self._ms.msfile)):
            #         return ValueError(f"The MS file {self._ms.msfile} could not be openned.")
            # finally:
            #     m.close()
            #
            # self._ms.logger.debug(f"sbd time range (from MS): {self._sbd_timerange}.")



    def prioritize_ref_antennas(self) -> list:
        """Given that some calibration steps (like fringe-fit) will break if there are no data
        for the given reference antenna, this will give the full list of available antennas,
        sorted with the reference antennas first, and later a manually-prioritized list.
        """
        full_antenna_list = ['EF', 'YS', 'GB', 'O8', 'JB', 'FD', 'PT', 'LA', 'OV', 'NL', 'SR', 'MC']
        sorted_antenna_list = list(self._ms.refant)

        # First priotizes the antennas with full bandwidth
        best_ants = [ant.name for ant in self._ms.antennas if \
                (len(self._ms.antennas[ant.name].subbands) == self._ms.freqsetup.n_subbands) \
                and ant.observed and (ant.name not in sorted_antenna_list)]

        rest_ants = [ant for ant in self._ms.antennas.observed if (ant not in best_ants) and \
                (ant not in sorted_antenna_list)]

        for ant in full_antenna_list:
            if ant in best_ants:
                sorted_antenna_list.append(ant)
                best_ants.remove(ant)

        sorted_antenna_list += best_ants
        sorted_antenna_list += rest_ants
        self._ms.logger.debug(f"Prioritized reference antennas: {', '.join(sorted_antenna_list)}.")
        return sorted_antenna_list


    def main_calibration(self, replace=False, fixed_fringe_callib=False, dispersive=False,
                         parang=True):
        """Runs the main calibration of the data:
        - instrumental delay corrections: in the specified time range.
        - global fringe fit: on all calibrators (and target if no phase-referencing is used).
        - bandpass calibration: using the bandpass calibrators.
        """
        #TODO: add the possibility of using a source model for fringe.
        # from the Source.model param if not None
        cals = {'sbd': self.caldir / f"{str(self._ms.projectname).lower()}.sbd",
                'mbd': self.caldir / f"{str(self._ms.projectname).lower()}.mbd",
                'bp': self.caldir / f"{str(self._ms.projectname).lower()}.bp"}

        if replace:
            for cal in cals:
                self.callib.remove_entry(cal, files=True)

        if cals['sbd'].exists():
            print("No instrumental delay calibration performed as the table "
                  f"{cals['sbd']} already exists.")
        else:
            rprint("[bold]Running fringefit for instrumental delay correction.[/bold]")
            # TODO: this is a patch because currently fringe-fit breaks with callib.... :-(
            if fixed_fringe_callib:
                casatasks.fringefit(vis=str(self._ms.msfile), caltable=str(cals['sbd']),
                                    timerange=self.get_sbd_timerange(), solint='inf',
                                    zerorates=True, paramactive=[True, True, dispersive],
                                    refant=','.join(self.prioritize_ref_antennas()), minsnr=50,
                                    docallib=True, delaywindow=[-200, 200],
                                    ratewindow=[-5e-8, 5e-8], callib=str(self.callib.filename),
                                    corrdepflags=True, parang=parang)
            else:
                casatasks.fringefit(vis=str(self._ms.msfile), caltable=str(cals['sbd']),
                                    timerange=self.get_sbd_timerange(),
                                    solint='inf', zerorates=True,
                                    paramactive=[True, True, dispersive],
                                    refant=','.join(self.prioritize_ref_antennas()), minsnr=5,
                                    gaintable=self.callib.gaintables(),
                                    delaywindow=[-200, 200], ratewindow=[-5e-8, 5e-8],
                                    interp=self.callib.interps(), #weightfactor=1,
                                    spwmap=self.callib.spwmaps(),
                                    # TODO: weightit not implemented yet in the stable release
                                    corrdepflags=True, parang=parang)

            self.callib.new_entry(name='sbd', caltable=str(cals['sbd']),
                                  parameters={"tinterp": 'nearest'})

        self._verify([cals['sbd']])

        if cals['mbd'].exists():
            print("No multi-band delay calibration performed as the table "
                  f"{cals['mbd']} already exists.")
        else:
            rprint("[bold]Running fringefit for multi-band delay correction.[/bold]")
            if fixed_fringe_callib:
                casatasks.fringefit(vis=str(self._ms.msfile), caltable=str(cals['mbd']),
                                    field=','.join(self._ms.sources.all_calibrators.names),
                                    solint='inf', zerorates=False,
                                    paramactive=[True, True, dispersive],
                                    delaywindow=[-200, 200], ratewindow=[-5e-8, 5e-8],
                                    refant=','.join(self.prioritize_ref_antennas()), combine='spw',
                                    minsnr=3, gaintable=self.callib.gaintables(),
                                    #weightfactor=1, TODO
                                    docallib=True, callib=str(self.callib.filename),
                                    corrdepflags=True, parang=parang)
            else:
                casatasks.fringefit(vis=str(self._ms.msfile), caltable=str(cals['mbd']),
                                    field=','.join(self._ms.sources.all_calibrators.names),
                                    solint='inf', zerorates=False,
                                    paramactive=[True, True, dispersive],
                                    delaywindow=[-200, 200], ratewindow=[-5e-8, 5e-8],
                                    refant=','.join(self.prioritize_ref_antennas()), combine='spw',
                                    minsnr=3, gaintable=self.callib.gaintables(),
                                    #weightfactor=1, TODO
                                    interp=self.callib.interps(), corrdepflags=True,
                                    spwmap=self.callib.spwmaps(), parang=parang)

            spw_with_solutions = get_spw_global_fringe(caltable=str(cals['mbd']))
            self.callib.new_entry(name='mbd', caltable=str(cals['mbd']),
                                  parameters={"tinterp": 'linear', "spwmap": \
                                        self._ms.freqsetup.n_subbands*[spw_with_solutions]})

        self._verify([cals['mbd']])

        if cals['bp'].exists():
            print(f"No bandpass calibration performed as the table {cals['bp']} already exists.")
        else:
            rprint("[bold]Running bandpass calibration.[/bold]")
            casatasks.bandpass(vis=str(self._ms.msfile),
                               field=','.join(self._ms.sources.fringe_finders.names), fillgaps=16,
                               combine='scan', caltable=str(cals['bp']), docallib=True,
                               corrdepflags=True, minsnr=3, minblperant=2,
                               callib=str(self.callib.filename), solnorm=True, solint='inf',
                               refant=','.join(self._ms.refant), bandtype='B', parang=parang)
            self.callib.new_entry(name='bandpass', caltable=str(cals['bp']),
                                  parameters={"tinterp": 'linear', "finterp": 'linear'})

        self._verify([cals['bp']])


    def recalibration(self, replace=False, dispersive=False, parang=True):
        """Runs the main calibration of the data:
        - instrumental delay corrections: in the specified time range.
        - global fringe fit: on all calibrators (and target if no phase-referencing is used).
        """
        #TODO: add the possibility of using a source model for fringe.
        cals = {'sbd2': self.caldir/f"{str(self._ms.projectname).lower()}.sbd2",
                'mbd2': self.caldir/f"{str(self._ms.projectname).lower()}.mbd2"}

        if replace:
            for cal in cals:
                remove_if_exists(cals[cal])

        if cals['sbd2'].exists():
            print("No second instrumental delay calibration performed as the table "
                  f"{cals['sbd2']} already exists.")
        else:
            rprint("[bold]Running fringefit for a second instrumental delay correction.[/bold]")
            casatasks.fringefit(vis=str(self._ms.msfile), caltable=str(cals['sbd2']),
                                timerange=self.get_sbd_timerange(), solint='inf', zerorates=True,
                                refant=','.join(self.prioritize_ref_antennas()), minsnr=5,
                                gaintable=self.callib.gaintables(), #weightfactor=1,
                                delaywindow=[-200, 200], ratewindow=[-5e-8, 5e-8],
                                interp=self.callib.interps(), corrdepflags=True,
                                spwmap=self.callib.spwmaps(), parang=parang)
            self.callib.new_entry(name='sbd2', caltable=str(cals['sbd2']),
                                  parameters={"tinterp" :'nearest'})

        self._verify([cals['sbd2']])

        if cals['mbd2'].exists():
            print("No second multi-band delay calibration performed as the table "
                  f"{cals['mbd2']} already exists.")
        else:
            rprint("[bold]Running fringefit for a second multi-band delay correction.[/bold]")
            casatasks.fringefit(vis=str(self._ms.msfile), caltable=str(cals['mbd2']),
                                field=','.join(self._ms.sources.all_calibrators.names),
                                solint='inf', zerorates=False, minsnr=3,
                                refant=','.join(self.prioritize_ref_antennas()), combine='spw',
                                gaintable=self.callib.gaintables(), #weightfactor=1,
                                delaywindow=[-200, 200], ratewindow=[-5e-8, 5e-8],
                                interp=self.callib.interps(), corrdepflags=True,
                                spwmap=self.callib.spwmaps(), parang=parang)
        spw_with_solutions = get_spw_global_fringe(caltable=str(cals['mbd2']))
        self.callib.new_entry(name='mbd2', caltable=str(cals['mbd2']),
                              parameters={"tinterp": 'linear', "spwmap": \
                                    self._ms.freqsetup.n_subbands*[spw_with_solutions]})

        self._verify([cals['mbd2']])

    #TODO: implement a fringefit instance here standalone
    def fringefit(self, caltable: Union[str, Path], field='', spw='', intent='', selectdata=True,
                  timerange='', uvrange='', antenna='', scan='', observation='', msselect='',
                  solint='inf', combine='', refant=None, minsnr=3.0, zerorates=False,
                  globalsolve=True, niter=100, delaywindow=[-200, 200], ratewindow=[-5e-8, 5e-8],
                  append=False, corrdepflags=True, corrcomb='none', docallib=False,
                  callib='', gaintable='', gainfield='', interp='', spwmap='', paramactive='',
                  concatspws=True, parang=False):
        """Fringe fit CASA task.
        """
        if refant is None:
            #TODO: select the prioritized antennas
            pass

        casatasks.fringefit(vis=str(self._ms.msfile), caltable=caltable, field=field,
                            spw=spw, intent=intent, selectdata=selectdata, timerange=timerange,
                            uvrange=uvrange, antenna=antenna, scan=scan, observation=observation,
                            msselect=msselect, solint=solint, combine=combine, refant=refant,
                            minsnr=minsnr, zerorates=zerorates, globalsolve=globalsolve,
                            niter=niter, delaywindow=delaywindow, ratewindow=ratewindow,
                            append=append, corrdepflags=corrdepflags, corrcomb=corrcomb,
                            docallib=docallib, callib=callib, gaintable=gaintable,
                            gainfield=gainfield, interp=interp, spwmap=spwmap,
                            paramactive=paramactive, concatspws=concatspws, parang=parang)
        # TODO: add the caltable to the callib
        return caltable


    def bandpass(self, caltable: Union[str, Path] = None, field: Union[str, None] = None,
                 spw='', intent='', selectdata=True, timerange='', uvrange='', antenna='',
                 scan='', observation='', msselect='', solint='inf', combine='scan',
                 refant: Union[str, None] = None, minblperant=2, minsnr=3.0, solnorm=False,
                 bandtype='B', smodel='', corrdepflags=True, append=False, fillgaps=16, degamp=3,
                 degphase=3, visnorm=False, maskcenter=0, maskedge=5, docallib=False, callib='',
                 gaintable='', gainfield='', interp='', spwmap='', parang=False, replace=True):
        """
        Modifed parameters:
        - caltable : str/Path = None
        - field : str/None = None
            If None, it will pick up the list of fringe finders from the current Project.
        - refant : str/None = None
            If None, it will pick the defined referece antennas from the project.
        """
        if caltable is None:
            caltable = self.caldir / f"{str(self._ms.projectname).lower()}.bp"

        if (not replace) and caltable.exists():
            raise FileExistsError(f"The calibration table {caltable} exists and 'replace' is" \
                                  "set to False. Cannot be overwritten.")
        elif replace and caltable.exists():
            shutil.rmtree(caltable)

        if field is None:
            field =  ','.join(self._ms.sources.fringe_finders.names)

        if refant is None:
            refant = ','.join(self._ms.refant)

        rprint("\n[bold]Running bandpass calibration.[/bold]")
        casatasks.bandpass(vis=str(self._ms.msfile), caltable=caltable, field=field,
                           spw=spw, intent=intent, selectdata=selectdata, timerange=timerange,
                           uvrange=uvrange, antenna=antenna, scan=scan, observation=observation,
                           msselect=msselect, solint=solint, combine=combine, refant=refant,
                           minblperant=minblperant, minsnr=minsnr, solnorm=solnorm,
                           bandtype=bandtype, smodel=smodel, corrdepflags=corrdepflags,
                           append=append, fillgaps=fillgaps, degamp=degamp, degphase=degphase,
                           visnorm=visnorm, maskcenter=maskcenter, maskedge=maskedge,
                           docallib=docallib, callib=callib, gaintable=gaintable,
                           gainfield=gainfield, interp=interp, spwmap=spwmap, parang=parang)
        # TODO: add to the callib
        # self.callib.new_entry(name=f"bandpass{str(bpver) if bpver > 1 else ''}",
        #                       caltable=str(calbp),
        #                       parameters={"tinterp": 'linear', "finterp": 'linear'})
        #
        # self._verify([calbp])

    def scalar_bandpass(self):
        """Pure scalar bandpass (amplitude only) calibration table based on the autocorrelations.
        It is based on the Michael Janssen's task from rPicard.
        """
        pass

    def apply_calibration(self, callib: Union[str, Path, None] = None, parang=True):
        """Applies the current calibration by using the tables written in the callib.
        """
        if callib is None:
            callib = str(self.callib.filename)
        elif isinstance(callib, Path):
            callib = str(callib)

        return casatasks.applycal(vis=str(self._ms.msfile), docallib=True,
                                  callib=callib, parang=parang)


    def clearcal(self, **kwargs):
        """Re-initializes the calibration for a visibility data set.
        """
        casatasks.clearcal(vis=str(self._ms.msfile), **kwargs)


    def create_cl_file(self, difmap_mod_file: str, outfile: str = 'component_list.cl'):
        """Creates a CASA-compatible file containing the component list of
        a given model created from Difmap.
        Code firsly done by Joe Bright during its visit to JIVE in 2022.

        Parameters
        ----------

        difmap_mod_file : str
            File containing the model created by Difmap.
        outfile : str  (optional)
            Output file to write. By default 'component_list.cl'.
        """

        with open(difmap_mod_file, 'r') as mymodel:
            mylines = mymodel.readlines()

            delta_components = []
            gaussian_components = []

            for line in mylines:
                if line.startswith('! Center'):
                    ra, dec = [a_line.split(':')[1] for a_line in line.split(',')]
                    ra = "{}h{}m{}s".format(*ra.split())
                    dec = "{}d{}m{}s".format(*(dec.split('(')[0][:-1]).split())
                    map_center = SkyCoord(ra, dec, frame='icrs')

                elif '!' not in line:
                    simpl_line = " ".join(line.split())
                    component = simpl_line.replace('v', '').split()
                    if len(component) == 3:
                        delta_components.append(component)
                    if len(component) == 9:
                        gaussian_components.append(component)

        # cl.open(outfile)
        for component in delta_components:
            offset = float(component[1]) / 1000. / 60. / 60. * u.deg
            angle = float(component[2]) * u.deg
            print('Component offset ' + str(offset))
            print('Component angle ' + str(angle))
            component_coordinate = map_center.directional_offset_by(angle, offset)
            print('Component coordinate + ' + component_coordinate.to_string('hmsdms'))
            component_flux = float(component[0])
            cl.addcomponent(dir='J2000 ' + map_center.to_string('hmsdms'), flux=component_flux,
                            fluxunit='Jy', shape='point')

        for component in gaussian_components:
            offset = float(component[1]) / 1000. / 60. / 60. * u.deg
            print(offset)
            angle = float(component[2]) * u.deg
            bmaj = float(component[3]) / 1000.
            bmin = bmaj / float(component[4])
            ang = float(component[5])
            component_coordinate = map_center.directional_offset_by(angle, offset)
            print('Component coordinate ' +  component_coordinate.to_string('hmsdms'))
            component_flux = float(component[0])
            cl.addcomponent(dir='J2000 ' + map_center.to_string('hmsdms'), flux=component_flux,
                            fluxunit='Jy', shape='Gaussian', majoraxis=f"{bmaj}arcsec",
                            minoraxis=f"{bmin}arcsec", positionangle=f"{ang}deg")

        cl.rename(outfile)
        cl.close()
        cl.done()


    # def self_calibration(project: obsdata.Project, source_model: str, calmode: str, solint: int):
    #     raise NotImplementedError
    #     # Maybe just a function that calls in loop all expected iteractions:
    #     #  p, p, p, a&p, p, a&p, p
    #     # In the amplitude selfcal (not the gscale with combine scan), use solnorm=True so it
    #     # only corrects for time-dependent gain residuals, not the flux scale.
    #
    #     # First, let's refine the delays
    #     gaincal(vis='tl016b_cal1.ms', field='J1154+6022', caltable='tl016b_cal1.dcal',
    #     solint='inf', refant=project.refants, minblperant=3, gaintype='K', calmode=calmode,
    #     parang=False)


# NOTE: for selfcal phase-only, use solnorm=False. but for A&P, solnorm=True. And parang=False?

# gaincal(vis='obj_selfcal.ms', caltable='selfcal_amplitude.tb', solint='48s', refant='ea24',
# calmode='ap', solnorm=True, normtype='median', gaintype='T', minsnr=6)
#
# applycal(vis='obj_selfcal.ms', gaintable='selfcal_amplitude.tb')


    # def self_calibration_from_difmap_image(project: obsdata.Project):
    #     """
    #     Project needs to contain the difmap associated images!!!!!!!!!!!
    #
    #     This task automatically corrects for it when importing the model in the MS file.
    #     """
    #     raise NotImplementedError
    #     # Assumes the model is a FITS file image generated by Difmap
    #     for a_cal in project.phase_calibrators:
    #         calsc = project.caldir / f"{project.project_name.lower()}.{a_cal}.sc_final"
    #         # TODO: check if ft is in casatools or casatasks
    #         # NOT IN MSFILE BUT IN SPLIT FILES
    #         casatools.ft(vis=str(project.msfile), field=a_cal, model=str(src_model),
    #                      usescratch=True)
    #         # Here open the MS, getcol MODEL  divide all complex numbers by 1000 and load it again
    #         casatasks.gaincal(vis=str(project.msfile), caltable=str(calsc), field=a_cal,
    #                           solint='3min', parang=True)
    #         casatasks.applycal()
    #


# casa task uvsubtraction()



class Aips(object):
    """Runs the calibration of a dataset in AIPS through ParselTongue.
    """
    def aipsno_from_project(self):
        """Follows a personal way to get an AIPS no from EVN experiments which is:
        for XXAAAB (e.g. EM123C), it takes the numeric part (123) and the epoch is
        converted to a number following the alphabet, so the resulting AIPS no is 1233.

        Otherwise it returns 101.
        """
        alphabet = list(string.ascii_lowercase)
        if self._ms.projectname[2:5].isnumeric():
            aipsno = self._ms.projectname[2:5]
            if len(self._ms.projectname) == 6:
                aipsno += str(alphabet.index(self._ms.projectname[5].lower()) + 1)
        else:
            aipsno = '101'

        return int(aipsno)


    def converttime2aips(self, casa_time: str):
        """Covnerts a time or timerange given in CASA format to AIPS format.
        """
        dt_time = []
        aips_time = []
        for a_time in casa_time.split('~'):
            if '/' in a_time:
                assert a_time.count('/') in (2, 3)
                # it should then be YYYY/MM/DDHH:MM[:SS]
                count_hhmm = a_time.count(':')
                if count_hhmm == 0:
                    dt_time.append(dt.datetime.strptime(a_time, "%Y/%m/%d"))
                elif count_hhmm == 1:
                    dt_time.append(dt.datetime.strptime(a_time, "%Y/%m/%d/%H:%M"))
                elif count_hhmm == 2:
                    dt_time.append(dt.datetime.strptime(a_time, "%Y/%m/%d/%H:%M:%S"))
                else:
                    raise ValueError("Unrecognized sbd_timerange ({casa_time}) to parse to AIPS.")

        for a_time in dt_time:
            aips_time += [(a_time.date() - self._ms.importdata.get_obsdate_from_fitsidi()).days,
                          a_time.hour, a_time.minute, a_time.second]

        return aips_time


    def a_priori_calibration(self, aipsno: Optional[int] = None,
                             aips_fitsidifile: Optional[str] = None):
        aipsno = self.aipsno_from_project() if aipsno is None else aipsno
        cmd = ["ParselTongue", _FILE_DIR + "/parseltongue.py", str(aipsno), self._ms.projectname,
               "-i", "--replace", "--fits",
               f"{self._ms.projectname}_1_1.IDI" if aips_fitsidifile is None else aips_fitsidifile
               ]
        if ('calibration' in self._ms.params) and \
                                ('do_ionospheric' in self._ms.params['calibration']):
            if self._ms.params['calibration']['do_ionospheric'] == "default":
                if self._ms.freqsetup.frequency is not None:
                    if self._ms.freqsetup.frequency < 6*u.GHz:
                        cmd += ["--iono"]
                else:
                    if self._ms.importdata.get_freq_from_fitsidi() < 6*u.GHz:
                        cmd += ["--iono"]
            elif self._ms.params['calibration']['do_ionospheric']:
                cmd += ["--iono"]
        else:
            if self._ms.freqsetup is not None:
                if self._ms.freqsetup.frequency < 6*u.GHz:
                    cmd += ["--iono"]
            else:
                if self._ms.importdata.get_freq_from_fitsidi() < 6*u.GHz:
                    cmd += ["--iono"]

        rprint(f"\n[bold]{' '.join(cmd)}[/bold]")
        result = subprocess.run(cmd, shell=False, stdout=None, stderr=subprocess.STDOUT)
        result.check_returncode()
        return f"{self._ms.projectname}.UVFITS"


    def main_calibration(self, aipsno: Optional[int] = None, uvfits: Optional[str] = None,
                         avgchan: bool = True):
        aipsno = self.aipsno_from_project() if aipsno is None else aipsno
        sbd_timerange = self.converttime2aips(self._ms.calibrate.get_sbd_timerange())

        cmd = ["ParselTongue", _FILE_DIR + "/parseltongue.py", str(aipsno), self._ms.projectname,
               "--uvfits", f"{self._ms.projectname}.UVFITS" if uvfits is None else uvfits,
               "--target", ','.join(self._ms.sources.targets.names),
               "--fringefinder", ','.join(self._ms.sources.fringe_finders.names),
               "--refant", ','.join(self._ms.refant), "--calib",
               "--sbdtime", ','.join([str(t) for t in sbd_timerange]), "--replace"]
        if len(self._ms.sources.phase_calibrators) > 0:
            cmd += ["--phaseref", ','.join(self._ms.sources.phase_calibrators.names)]

        if avgchan:
            cmd += ["--avgchan"]

        rprint(f"\n[bold]{' '.join(cmd)}[/bold]")
        result = subprocess.run(cmd, shell=False, stdout=None, stderr=subprocess.STDOUT)
        result.check_returncode()

        splits = {}
        for a_src in self._ms.sources.names:
            if Path(f"{self._ms.projectname}.{a_src}.SPLIT.UVFITS").exists():
                splits[a_src] = capi.Project(f"{self._ms.projectname}.{a_src}", cwd=self._ms.cwd,
                                             params=self._ms.params,
                                             logging_level=self._ms._logging_level)
                splits[a_src].importdata.import_uvfits(
                        uvfitsfile=f"{self._ms.projectname}.{a_src}.SPLIT.UVFITS", delete=True)
                self._ms.splits[a_src].append(splits[a_src])

        return splits


    def selfcal_with_difmapimage(self, fitsimage: str, aipsno: Optional[int] = None):
        """Imports a FITS image into AIPS.
        """
        fitsimage = Path(fitsimage)
        aipsno = self.aipsno_from_project() if aipsno is None else aipsno
        if not fitsimage.exists():
            raise FileNotFoundError(f"The file {fitsimage} cannot be found.")

        capi.tools.fix_difmap_image(fitsimage)
        cmd = ["ParselTongue", _FILE_DIR + "/parseltongue.py", str(aipsno), self._ms.projectname,
               "--selfcal", str(fitsimage), "--target", ','.join(self._ms.sources.targets.names),
               "--fringefinder", ','.join(self._ms.sources.fringe_finders.names),
               "--refant", ','.join(self._ms.refant)]
        if len(self._ms.sources.phase_calibrators) > 0:
            cmd += ["--phaseref", ','.join(self._ms.sources.phase_calibrators.names)]

        rprint(f"\n[bold]{' '.join(cmd)}[/bold]")
        result = subprocess.run(cmd, shell=False, stdout=None, stderr=subprocess.STDOUT)
        result.check_returncode()


    def transfer_calibration(self, fitsididata: str, aipsno: Optional[int] = None,
                             source: Optional[Union[list,str]] = None, avgchan: bool = True):
        """Transfers the calibration conducted in the current dataset (within AIPS) to
        a new set of FITS-IDI files that should belong to the same observation
        (but e.g. a different correlator pass).

        It will write out a SPLIT FITS file for the given sources with the latest calibration
        applied.

        Inputs
        ------

        - fitsididata : str
            Name of the FITS-IDI/UVFITS files to import into AIPS. In case of multiple files,
            it should be the root name of the files as typically used in AIPS when importing.
            For example, if the files to import are exp_2_1.IDI1...9, this param should be
            'exp_2_1.IDI'.
        - aipsno: int   (default None)
            Number to use for the AIPS user. By default it will retrieve it from the project
            name.
        - source : list or str  (default None)
            The list of sources to which the calibration will be applied and will be split'ed
            in UVFITS files. If str, it needs to be a comma-separated succession of sources.
        - avgchan : bool  (default True)
            If True, the final SPLIT files will be single-channel multiple-IF FITS files.
            This should be consistent with the options used during selfcalibration.
        """
        fitsfiles = glob.glob(f"{fitsididata}*")
        aipsno = self.aipsno_from_project() if aipsno is None else aipsno
        if len(fitsfiles) == 0:
            raise FileNotFoundError(f"The file {fitsididata} cannot be found.")

        if source is None:
            source = ','.join(self._ms.sources.targets.names)
        elif type(source) == str:
            source = ','.join([src.strip() for src in source.split(',')])
        else:
            source = ','.join(source)

        cmd = ["ParselTongue", _FILE_DIR + "/parseltongue.py", str(aipsno), self._ms.projectname,
               "--transfer", str(fitsididata), "--target", source,
               "--refant", ','.join(self._ms.refant)]
        if avgchan:
            cmd += ["--avgchan"]

        rprint(f"\n[bold]{' '.join(cmd)}[/bold]")
        result = subprocess.run(cmd, shell=False, stdout=None, stderr=subprocess.STDOUT)
        result.check_returncode()


    def __init__(self, ms: capi.Project):
        self._ms = ms
