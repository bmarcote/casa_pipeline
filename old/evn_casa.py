"""EVN CASA Pipeline


Usage: evn_casa.py  <input file>

"""


import sys
import argparse


# input file
## at least two things to read from the input file first: showgui?  log? (outdir)


exec("casa --nogui -c scripts/gc.py <antabfile>  <outdir>/EVN.gc")
# It will probably need -l <value> -u <value> for elevation ranges.




exec("casa -c main.py  <inputfile> --logfile <logfile>")
exec("casa -c main.py --nogui <inputfile> --logfile <logfile>")

# Typically you call it like
#   !casa -c thisfile.py -i input_file
#   !mpicasa -n <num_cores> casa -c thisfile.py -i input_file
# --nogui --logfile <logfile>

