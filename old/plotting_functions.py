import sys
import os
import jplotter




def open_ms(msdata):
    yield "ms {}".format(msdata)
    yield "indexr"


def weights():
    """Plot of weights vs time."""
    yield "pt wt"
    yield "bl auto"
    yield "fq */p"
    yield "sort bl"
    yield "ckey sb"
    yield "ptsz 4"
    yield "pl"
    print('Done weight-time plot.')



def autocorr(source, outfile=None):
    yield "pt ampchan"
    yield "bl auto"
    yield "fq */p"
    yield "src {}".format(source)
    yield "ch none"
    yield "avt vector"
    yield "ckey p"
    yield "multi true"
    yield "new sb false"
    yield "new time true"
    yield "sort time bl"
    yield "y 0 1.6"
    yield "nxy 1 4"
    yield "pl"
    print('Done amp-freq autocorrelations plot.')


def crosscorr(source, refstation, outfile=None):
    yield "pt ampchan"
    yield "bl {}* -auto".format(refstation)
    yield "fq *"
    yield "src {}".format(source)
    yield "ch none"
    yield "avt vector"
    yield "ckey p"
    yield "multi true"
    yield "new sb false"
    yield "new time true"
    yield "sort time bl"
    yield "y 0 1.6"
    yield "nxy 1 4"
    yield "pl"
    print('Done amp-freq crosscorrelations plot.')









