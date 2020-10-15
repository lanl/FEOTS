#!/usr/bin/python3
DOC="""zProfile

zProfile is used to create vertical profiles of different lateral field statistics of FEOTS output.

Usage:
  zProfile plot <file> [--field=<field>] [--out=<out>] [--stats=<stats>] [--opts=<opts>]

Commands:
  plot              Create a vertical profile plot of the chosen statistics for the given FEOTS output field.
  report            Create a curve file with vertical profiles of the absolute max volume correction field.

Options:
  -h --help                 Display this help screen
  --field=<field>           The field in the FEOTS output you want to plot [default: DyeTracer_00]
  --out=<out>               The path to place the output files [default: ./]
  --stats=<stats>           Comma separated list of statistics to show in the vertical profile plot. [default: absmax]
  --opts=<opts>             Comma separated list of plot options. [default: none]
"""
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
from docopt import docopt
import feotsPostProcess as feots

def parse_cli():

  args = docopt(DOC,version='zProfile 0.0.0')
  return args

#END parse_cli

def load_netcdf(filename, field):

  rootgrp = nc.Dataset(filename,"r",format="NETCDF4")
  ncdata = rootgrp.variables[field][:]
  return ncdata

#END load_netcdf

def plotZstats(zstats, stats, z, field, opts, plotfile):

    f, ax = plt.subplots()
    for key in stats:
      if key == 'stdev':
        ax.fillbetweenx(-z/100.0,zstats['mean']-zstats['std'], zstats['mean']+zstats['std'], color=(0.8,0.8,0.8,0.8))
    else:
        ax.plot(zstats[key], -z/100.0, marker='', color='black', linewidth=2, label=key)


    if 'logx' in opts:
      ax.set(xscale='log')

    ax.grid(color='gray', linestyle='-', linewidth=1)

    ax.set(xlabel=field, ylabel='Depth (m)')

    f.savefig(plotfile)

    plt.close('all')

#END plotZstats

def main():

  args = parse_cli()

  if args['plot'] :
    field = args['--field']
    stats = args['--stats'].split(',')
    opts = args['--opts'].split(',')
    print('Plotting {FIELD} from {FILE}'.format(FIELD=field, FILE=args['<file>']))

    feotsData = feots.load_netcdf(args['<file>'],field)
    z = feots.load_netcdf(args['<file>'],'z_t')

    zProfiles = feots.VerticalProfileStats(feotsData)

    iterate = args['<file>'].split('/')[-1].split('.')[2]
    outFile = args['--out']+'/'+field+'_'+iterate+'.png'
    plotZstats(zProfiles, stats, z, field, opts, outFile)


#END main

if __name__ == '__main__':
    main()
