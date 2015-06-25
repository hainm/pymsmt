#!/usr/bin/env python
# Filename: espgen.py
"""
This is a module for generate the esp points from the Gaussian/GAMESS-US
output file.
"""
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

from pymsmtmol.gauio import get_esp_from_gau
from pymsmtmol.gmsio import get_esp_from_gms
from optparse import OptionParser

parser = OptionParser("usage: -i input_file -o output_file [-v software]")

parser.set_defaults(softversion='gau')

parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-o", dest="outputfile", type='string',
                  help="Output file name")
parser.add_option("-v", dest="softversion", type='string',
                  help="Software version [Default is gau (means Gaussian), \n"
                       "           other option is gms (means GAMESS-US)]")
(options, args) = parser.parse_args()

if options.softversion == 'gau':
    get_esp_from_gau(options.inputfile, options.outputfile)
elif options.softversion == 'gms':
    get_esp_from_gms(options.inputfile, options.outputfile)

quit()

