#!/usr/bin/python3
# Copyright 2020 Los Alamos National Laboratories
#
# Author : Dr. Joseph Schoonover, Fluid Numerics LLC
#
# About :
#   Basic support routines for working with FEOTS output in Python.
#
# //////////////////////////////////////////////////////////////////// #

import netCDF4 as nc
import numpy as np

def load_netcdf(filename, field):

  rootgrp = nc.Dataset(filename,"r",format="NETCDF4")
  ncdata = rootgrp.variables[field][:]
  return ncdata

#END load_netcdf

def VerticalProfileStats(field):

  stats = {'min': np.zeros((field.shape[1])),
           'mean': np.zeros((field.shape[1])),
           'max': np.zeros((field.shape[1])),
           'std': np.zeros((field.shape[1])),
           'absmax': np.zeros((field.shape[1]))}
#  vertProfStats = np.zeros((field.shape[1],5))
  for k in range(0,field.shape[1]):
    subfield = np.squeeze(field[0,k,:,:])
    stats['min'][k] = np.min(subfield,(0,1))
    stats['mean'][k] = np.mean(subfield,(0,1))
    stats['max'][k] = np.max(subfield,(0,1))
    stats['std'][k] = np.std(subfield,(0,1))
    stats['absmax'][k] = np.max(np.abs(subfield),(0,1))

  return stats

#END genVerticalProfileStats

#if __name__ == '__main__':
#    main()
