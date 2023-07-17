
import casa_pipeline as capi

class Imaging(object):
    """Defines all imaging tasks that can run over a Ms object
    """
    def __init__(self, ms: capi.Project):
        self._ms = ms

    def tclean(self):
        pass

    def wsclean(self, parameters: str = None):
        if parameters is None:
            capi.tools.shell_command("wsclean", [self._ms.params['imaging']['wsclean'].split(' '),
                                                 str(self._ms.msfile)])
        else:
            capi.tools.shell_command("wsclean", [parameters.split(' '), str(self._ms.msfile)])

    def import_difmap(self, fitsfile = str):
        """Imports a Difmap image to CASA format
        It will first fix the FITS image if the metadata is unexpected (as for stokes)
        """
        pass

