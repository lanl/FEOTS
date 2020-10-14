#!/usr/bin/python3
DOC="""VolumeZProfile

VolumeZProfile is used to create vertical profiles of the absolute max of the volume correction field.

Usage:
  VolumeZProfile plot <file> [--tracerid=<tracerid>] [--out=<out>]

Commands:
  plot              Create a vertical profile plot of the mean, min, and max volume correction field.
  report            Create a curve file with vertical profiles of the absolute max volume correction field.

Options:
  -h --help                 Display this help screen
  --tracerid=<tracerid>     The tracer field associated with the volume correction field you want to plot [default: 0]
  --out=<out>               The path to place the output files [default: ./]
"""
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
from docopt import docopt

def parse_cli():

  args = docopt(DOC,version='VolumeZProfile 0.0.0')
  return args

#END parse_cli

def load_netcdf(filename, field):

  rootgrp = nc.Dataset(filename,"r",format="NETCDF4")
  ncdata = rootgrp.variables[field][:]
  return ncdata

#END load_netcdf

def genVerticalProfileStats(field):

  vertProfStats = np.zeros((field.shape[1]))
  for k in range(0,field.shape[1]):
    subfield = np.squeeze(field[0,k,:,:])
    vertProfStats[k] = np.max(np.abs(subfield),(0,1))

  return vertProfStats

#END genVerticalProfileStats

def plotZstats(zstats, z, plotfile):

    f, ax = plt.subplots()
    #ax.set_title('Simple plot')
    ax.semilogx(zstats, -z/100.0, marker='', color='black', linewidth=2)
    ax.grid(color='gray', linestyle='-', linewidth=1)

    ax.set(xlabel='max(|V|)', ylabel='Depth (m)')

    f.savefig(plotfile)

    plt.close('all')

#END plotZstats

def main():

  args = parse_cli()

  if args['plot'] :
    field = 'Volume_'+str(args['--tracerid']).zfill(2)
    print('Plotting {FIELD} from {FILE}'.format(FIELD=field, FILE=args['<file>']))

    volume = load_netcdf(args['<file>'],field)
    z = load_netcdf(args['<file>'],'z_t')

    volumeZProfiles = genVerticalProfileStats(volume)

    ncFile = args['<file>'].split('/')[-1]
    outFile = args['--out']+'/'+ncFile+'.png'
    plotZstats(volumeZProfiles, z, outFile)


#END main

if __name__ == '__main__':
    main()
