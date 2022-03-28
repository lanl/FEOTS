#!/usr/bin/python3

DOC="""feots_compare

feots_compare is use to compare two FEOTS NetCDF output files and report simple statistics.
Currently feots_compare will generate a histogram of log_{10}( |f_1 - f_2| ) where  f_1 and
f_2 are tracer fields from two FEOTS output files.

Usage:
  feots_compare absdiff <file1> <file2> [--field=<tracerfield>]

Commands:
  absdiff              Compute statistics using absolute differences between two FEOTS files

Options:
  -h --help            Display this help screen
  --field=<string>     Specification of the field in the NetCDF file to compare [default: DyeTracer_01]
"""
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
from docopt import docopt

def parse_cli():

  args = docopt(DOC,version='feots_compare 0.0.0')
  return args

#END parse_cli

def load_netcdf(filename, field):

  rootgrp = nc.Dataset(filename,"r",format="NETCDF4")
  ncdata = rootgrp.variables[field][:]
  return ncdata

def main():

  args = parse_cli()

  if args['absdiff'] :
    print('Comparing {FIELD} in {FILE1} and {FILE2}'.format(FIELD=args['--field'],
                                                            FILE1=args['<file1>'],
                                                            FILE2=args['<file2>']))
    file1_data = load_netcdf(args['<file1>'],args['--field'])
    file2_data = load_netcdf(args['<file2>'],args['--field'])

    absdiff =  np.log10(np.absolute( file1_data - file2_data ))
    absdiffHist, absdiffBins = np.histogram(absdiff, bins=50, range=(-16, 0))

    print(absdiffBins)
    print(absdiffHist)
    plt.hist(absdiff.flatten(), absdiffBins)
    plt.title("histogram")
    plt.show()

#END main

if __name__ == '__main__':

  main()
