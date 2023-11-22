import os
import sys
import time
import subprocess
import numpy
import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
from function_io import *
from k_registration import *
import shutil
import math
from numba import jit
import copy

RAD_MAX =16
NPTS =16

CMPLX 	=0
SCMPLX  =1
FLOAT	=2
INT	=3
SHORT   =4
BYTE	=5
RASTER  =6

class SRCHTAB:
    def __init__(self,count):
        self.wgt =np.zeros(count,dtype="float32")
        self.ir = np.zeros(count, dtype=int)
        self.ix = np.zeros(count, dtype=int)
        self.iy = np.zeros(count, dtype=int)

def spiral(rmax, stab, w_mode):

    # Generate spiral search table containing locations, distance,
    # and weights to perform interpolation of data with missing values.
    #
    # spiral returns the number of entries in stab
    # rmax   maximum radius of the interpolation window
    # stab   interpolation window locations, radii, and weights
    # wmode  interpolation weighting function flag with values:
    #
    #  0: 1				  uniform
    #  1: 1 - r/(rmax+1)		  linear distance
    #  2: 1 - r*r/((rmax+1)*(rmax+1)) distance squared
    #  3: exp(-2r*r/(rmax*rmax)	  Gaussian
  	# 			21-Nov-2003 clw
    w1 = 2 * rmax + 1
    nps = 0
    mpb = np.zeros([w1 * w1],dtype=int)
    mp= np.zeros([w1 , w1],dtype=int)
    for k in range(1,rmax+1):
        r2m = k * k
        for i in range(-rmax,rmax+1):
            for j in range(-rmax, rmax + 1):
                r2 = i * i + j * j
                if (r2 <= r2m + 1):
                    i1 = i + rmax
                    j1 = j + rmax

                    if (mp[i1][j1] == 0):
                        mp[i1][j1] = k
                        stab.ix[nps] = i
                        stab.iy[nps] = j
                        stab.ir[nps] = k
                        if w_mode==0:
                            stab.wgt[nps] = 1.0
                        elif w_mode==1:
                            stab.wgt[nps] = 1.0 - sqrt(r2) / (rmax + 1)
                        elif w_mode == 2:
                            stab.wgt[nps] = 1.0 - r2 / ((rmax + 1) * (rmax + 1))
                        elif w_mode == 3:
                            stab.wgt[nps] = exp(-2. * r2 / ((rmax + 1) * (rmax)))
                        else:
                            print("ERROR: unsupported interpolation method: %d"%w_mode)
                            exit(-1)
                        nps += 1
    return (nps)

def interp_ad(data_in,data_out,height,r_max1,np_min1,np_max1,w_mode1,type,cp_data1):
    irmax = RAD_MAX
    np_min = NPTS
    np_max = NPTS
    w_mode = 2
    dtype = FLOAT
    cp_data = 1

    print("*** Weighted interpolation of gaps in 2D data using an adaptive smoothing window ***\n")
    print("*** Copyright 2008, Gamma Remote Sensing, v2.0 25-Aug-2008 clw/uw ***\n")

    print("usage: %s <data_in> <data_out> <width> [r_max] [np_min] [np_max] [w_mode] [type] [cp_data]")
    print("input parameters:")
    print("  data_in   (input) data with gaps")
    print("  data_out  (output) data with gaps filled by interpolation")
    print("  width -> height     number of samples/row -> row")    #这里我改成height(方位向行数)
    print("  r_max     maximum interpolation window radius (default(-): %d)"%irmax)
    print("  np_min    minimum number of points used for the interpolation (default(-): %d)"%np_min)
    print("  np_max    maximum number of points used for the interpolation (default(-): %d)"%np_max)
    print("  w_mode    data weighting mode (enter - for default):\n")
    print("              0: constant\n")
    print("              1: 1 - (r/r_max)\n")
    print("              2: 1 - (r/r_max)**2  (default)\n")
    print("              3: exp(-2.*(r**2/r_max**2))\n")
    print("  type      input and output data type:\n")
    print("              0: FCOMPLEX\n")
    print("              1: SCOMPLEX\n")
    print("              2: FLOAT (default)\n")
    print("              3: INT\n")
    print("              4: SHORT\n")
    print("  cp_data   copy data flag:\n")
    print("              0: do not copy input data values to output\n")
    print("              1: copy input data values to output (default)\n\n")

    print("input data file: %s\n"%data_in)
    print("output data file: %s\n"%data_out)

    if type!="-":
        dtype=type
    if ((dtype < 0) or (dtype > 4)):
        print("ERROR: invalid data type, must be in the range 0->4: %d"%dtype)
        exit(-1)
    if cp_data1 != "-":
        cp_data = cp_data1

    if cp_data==0:
        print("input data values not copied to output data file")
    elif cp_data==1:
        print("input data values copied to output data file")
    else:
        print("ERROR: invalid value for cp_data input parameter: %d\n"%cp_data)
        exit(-1)

    if dtype==0:
        da1=freadbk(data_in,">f",height,1)
        width=da1.shape[1]
        d1cpx=np.zeros([height,width],dtype=complex)
        d2cpx = np.zeros([height, width], dtype=complex)
        d1cpx=copy.copy(da1)
        print("data type: float complex\n")
    elif dtype==1:
        da1 = freadbk(data_in, ">f", height, 1)
        width = da1.shape[1]
        d1scpx = np.zeros([height, width], dtype=complex)
        d2scpx = np.zeros([height, width], dtype=complex)
        d1scpx = copy.copy(da1)
        print("data type: short complex\n")
    elif dtype == 2:
        da1 = freadbk(data_in, ">f", height, 0)
        width = da1.shape[1]
        d1flt = np.zeros([height, width], dtype="float32")
        d2flt = np.zeros([height, width], dtype="float32")
        d1flt = copy.copy(da1)
        print("data type: float\n")
    elif dtype == 3:
        da1 = freadbk(data_in, ">f", height, 0)
        width = da1.shape[1]
        d1int = np.zeros([height, width], dtype=int)
        d2int = np.zeros([height, width], dtype=int)
        d1int = copy.copy(da1)
        print("data type: integer\n")
    elif  dtype == 4:
        da1 = freadbk(data_in, ">f", height, 0)
        width = da1.shape[1]
        d1short = np.zeros([height, width], dtype=int)
        d2short = np.zeros([height, width], dtype=int)
        d1short = copy.copy(da1)
        print("data type: integer\n")
    else:
        print("\nERROR: invalid data type flag: %d\n"%dtype)
        exit(-1)

    print("input/output width: %d   lines: %d"%(width, height))
    if r_max1!="-":
        irmax=r_max1
        if (irmax < 1):
            print("ERROR: invalid value for r_max, must be > 0")
            exit(-1)
    print("interpolation window radius: %d"%irmax)

    if np_min1!="-":
        np_min=np_min1
        if (np_min < 1):
            print("ERROR: invalid value for np_min, must be > 0")
            exit(-1)
    print("minimum number of interpolation points: %d"%np_min)
    np_max = NPTS
    if np_max1!="-":
        np_max=np_max1
        if (np_max < 1):
            print("ERROR: invalid value for np_max, must be > 0")
            exit(-1)

    print("maximum number of interpolation points: %d"%np_max)

    if w_mode1!="-":
        w_mode=w_mode1

    if w_mode==0:
        print("interpolation weights: 1.0")
    elif w_mode == 1:
        print("interpolation weights: 1 - r/(rmax+1)")
    elif w_mode == 2:
        print("interpolation weights: 1 - r**2/((rmax+1)**2)")
    elif w_mode == 3:
        print("interpolation weights: Gaussian exp(-2*(r**2/(rmax**2))")
    else:
        print("ERROR: invalid interpolation mode")
        exit(-1)

    spir=SRCHTAB((2 * irmax + 1) * (2 * irmax + 1))
    nps = spiral(irmax, spir, w_mode)
    print("number of points in the interpolation window: %d"%nps)

    if dtype==0:
        for i in range(0,height):
            if (i % 400 == 0):
                print("line: %d"%i)
            for j in range(0, width):
                if ((cp_data == 1) and ((d1cpx[i][j].real != 0.0) or (d1cpx[i][j].imag != 0.0))) :
                    d2cpx[i][j]=d1cpx[i][j]
                    continue
                np1 = 0
                swv1 = 0.0
                swv2 = 0.0
                swgt = 0.0
                for k in range(0,nps):
                    if (spir.ir[k] > irmax):break
                    if ((i + spir.iy[k] < 0) or (i + spir.iy[k] >= height)):continue
                    if ((j + spir.ix[k] < 0) or (j + spir.ix[k] >= width)):continue
                    v1 = d1cpx[i + spir.iy[k]][j + spir.ix[k]].real
                    v2 = d1cpx[i + spir.iy[k]][j + spir.ix[k]].imag
                    if ((v1 != 0.0) or (v2 != 0.0)):
                        np1 +=1
                        wgt = spir.wgt[k]
                        swgt += wgt
                        swv1 += wgt * v1
                        swv2 += wgt * v2
                        if (np1 >= np_max):break
                if ((np1 >= np_min) and (swgt > 0.0)):
                    d2cpx[i][j]=complex(swv1 / swgt,swv2 / swgt)
        f = open(data_out, 'wb')
        f.write(d2cpx)
        f.close()

    elif dtype==1:
        for i in range(0,height):
            if (i % 400 == 0):
                print("line: %d"%i)
            for j in range(0, width):
                if ((cp_data == 1) and ((d1scpx[i][j].real != 0.0) or (d1scpx[i][j].imag != 0.0))) :
                    d2scpx[i][j]=d1scpx[i][j]
                    continue
                np1 = 0
                swv1 = 0.0
                swv2 = 0.0
                swgt = 0.0
                for k in range(0,nps):
                    if (spir.ir[k] > irmax): break
                    if ((i + spir.iy[k] < 0) or (i + spir.iy[k] >= height)): continue
                    if ((j + spir.ix[k] < 0) or (j + spir.ix[k] >= width)): continue
                    v1 = d1scpx[i + spir.iy[k]][j + spir.ix[k]].real
                    v2 = d1scpx[i + spir.iy[k]][j + spir.ix[k]].imag
                    if ((v1 != 0.0) or (v2 != 0.0)):
                        np1 +=1
                        wgt = spir.wgt[k]
                        swgt += wgt
                        swv1 += wgt * v1
                        swv2 += wgt * v2
                        if (np1 >= np_max):break
                if ((np1 >= np_min) and (swgt > 0.0)):
                    d2scpx[i][j]=complex(round(swv1 / swgt),round(swv2 / swgt))
        f = open(data_out, 'wb')
        f.write(d2scpx)
        f.close()

    elif dtype==2:
        print("float interpolation")
        for i in range(0,height):
            if (i % 400 == 0):
                print("line: %d"%i)
            for j in range(0, width):
                if ((cp_data == 1) and (d1flt[i][j] != 0.0)) :
                    d2flt[i][j] = d1flt[i][j]
                    continue
                np1 = 0
                swv1 = 0.0
                swgt = 0.0
                for k in range(0,nps):
                    if (spir.ir[k] > irmax): break
                    if ((i + spir.iy[k] < 0) or (i + spir.iy[k] >= height)): continue
                    if ((j + spir.ix[k] < 0) or (j + spir.ix[k] >= width)): continue
                    v1 = d1flt[i+spir.iy[k]][j+spir.ix[k]]
                    if (v1 != 0.0):
                        np1 +=1
                        wgt = spir.wgt[k]
                        swgt += wgt
                        swv1 += wgt * v1
                        if (np1 >= np_max):break
                if ((np1 >= np_min) and (swgt > 0.0)):
                    d2flt[i][j] = swv1 / swgt
        f = open(data_out, 'wb')
        f.write(d2flt)
        f.close()

    elif dtype == 3:
        for i in range(0, height):
            if (i % 400 == 0):
                print("line: %d" % i)
            for j in range(0, width):
                if ((cp_data == 1) and (d1int[i][j] != 0)):
                    d2int[i][j] = d1int[i][j]
                    continue
                np1 = 0
                swv1 = 0.0
                swgt = 0.0
                for k in range(0, nps):
                    if (spir.ir[k] > irmax): break
                    if ((i + spir.iy[k] < 0) or (i + spir.iy[k] >= height)): continue
                    if ((j + spir.ix[k] < 0) or (j + spir.ix[k] >= width)): continue
                    v1 = d1int[i + spir.iy[k]][j + spir.ix[k]]
                    if (v1 != 0.0):
                        np1 += 1
                        wgt = spir.wgt[k]
                        swgt += wgt
                        swv1 += wgt * v1
                        if (np1 >= np_max): break
                if ((np1 >= np_min) and (swgt > 0.0)):
                    d2int[i][j] = round(swv1 / swgt)

        f = open(data_out, 'wb')
        f.write(d2int)
        f.close()

    elif dtype == 4:
        for i in range(0, height):
            if (i % 400 == 0):
                print("line: %d" % i)
            for j in range(0, width):
                if ((cp_data == 1) and (d1short[i][j] != 0)):
                    d2short[i][j] = d1short[i][j]
                    continue
                np1 = 0
                swv1 = 0.0
                swgt = 0.0
                for k in range(0, nps):
                    if (spir.ir[k] > irmax): break
                    if ((i + spir.iy[k] < 0) or (i + spir.iy[k] >= height)): continue
                    if ((j + spir.ix[k] < 0) or (j + spir.ix[k] >= width)): continue
                    v1 = d1short[i + spir.iy[k]][j + spir.ix[k]]
                    if (v1 != 0.0):
                        np1 += 1
                        wgt = spir.wgt[k]
                        swgt += wgt
                        swv1 += wgt * v1
                        if (np1 >= np_max): break
                if ((np1 >= np_min) and (swgt > 0.0)):
                    s1 = round(swv1 / swgt)
                    if ((s1 == 0) and (swv1 != 0.0)):
                        if (swv1 < 0):d2short[i][j] = -1
                        else:d2short[i][j] = 1
                    else:
                        d2short[i][j] = s1
        f = open(data_out, 'wb')
        f.write(d2short)
        f.close()
    else:
        print("\nERROR: invalid data type flag: %d"%dtype)
        exit(-1)

if __name__ == '__main__':
    interp_ad("20121228.unw","20121228.unw_inter1",7051,32,8,16,2,"-","-")
    #m = freadbk('20121228.unw', ">f", 7051, 0)
    s1 = freadbk('20121228.unw_inter', ">f", 7051, 0)
    s2 = freadbk('20121228.unw_inter1', "<f", 7051, 0)
    s=s2-s1

    # ms = m - s
    plt.figure(dpi=600, figsize=(6, 10))
    plt.imshow(s2-s1, cmap=plt.cm.jet, aspect='auto',vmin=-0.0001,vmax=0.0001)
    plt.colorbar()
    # plt.savefig('jingpeizhun4.png')
    plt.show()
    plt.clf()

































































