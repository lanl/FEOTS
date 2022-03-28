#!/usr/bin/python3
DOC="""plotCurve

plotCurve is used to create vertical profiles of different lateral ylabel statistics of FEOTS output.

Usage:
  plotCurve plot <file>  [--out=<out>] [--opts=<opts>] [--scalex=<scalex>] [--xlabel=<xlabel>] [--ylabel=<ylabel>]

Commands:
  plot              Create a vertical profile plot of the chosen statistics for the given FEOTS output ylabel.

Options:
  -h --help                 Display this help screen
  --out=<out>               The path to place the output files [default: ./]
  --opts=<opts>             Comma separated list of plot options. [default: none]
  --scalex=<scalex>         Amount to scale the x dimension by for the plot (multiplicative). [default: 1.0]
  --xlabel=<xlabel>         Label for the x-dimension in the plot. [default: x]
  --ylabel=<ylabel>         Label for the y-dimension in the plot. [default: y]
"""
import numpy as np
from matplotlib import pyplot as plt
from docopt import docopt
import feotsPostProcess as feots

def parse_cli():

  args = docopt(DOC,version='plotCurve 0.0.0')
  return args

#END parse_cli

def loadCurve(filename):

  curveData = np.loadtxt(filename, delimiter=",", skiprows=1)
  return curveData

#END loadCurve

def plotCurve(curveData, opts, scalex, xlabel, ylabel, plotfile):

    f, ax = plt.subplots()
    ax.fillbetween(curveData[:,0]*scalex,curveData[:,1], color=(0.8,0.8,0.8,0.8))
    ax.plot(curveData[:,0]*scalex, curveData[:,1], marker='', color='black', linewidth=2)

    if 'logx' in opts:
      ax.set(xscale='log')

    if 'logy' in opts:
      ax.set(yscale='log')

    ax.grid(color='gray', linestyle='-', linewidth=1)

    ax.set(xlabel=xlabel, ylabel=ylabel)

    f.savefig(plotfile)

    plt.close('all')

#END plotCurve

def main():

  args = parse_cli()

  if args['plot'] :
    xlabel = args['--xlabel']
    scalex = args['--scalex']
    ylabel = args['--ylabel']
    opts = args['--opts'].split(',')

    curveData = loadCurve(args['<file>'])

    outFile = args['--out']
    plotCurve(curveData, opts, scalex, xlabel, ylabel, outFile)


#END main

if __name__ == '__main__':
    main()
