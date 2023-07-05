
import casa_pipeline as capi

class Imaging(object):
    """Defines all imaging tasks that can run over a Ms object
    """
    def __init__(self, ms: capi.Project):
        self._ms = ms

    def tclean(self):
        pass

    def wsclean(self, params=None):
        if params is None:
            capi.tools.shell_command("wsclean", [str(self._ms.msfile)])
        else:
            if isinstance(params, str):
                capi.tools.shell_command("wsclean", [*params.split(' '), str(self._ms.msfile)])
            elif isinstance(params, list):
                capi.tools.shell_command("wsclean", [*params, str(self._ms.msfile)])
            else:
                raise TypeError(f"The parameters for WSClean ({params}) are neither a str or a list.")

