#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 14:37:37 2024

@author: alex

Rewrite results after SExtractor for region plot on DS9, etc.

Input information format:
    X_IMAGE         Object position along x                     [pixel]
    Y_IMAGE         Object position along y                     [pixel]
    A_IMAGE         Profile RMS along major axis                [pixel]
    B_IMAGE         Profile RMS along minor axis                [pixel]
    XMIN_IMAGE      Minimum x-coordinate among detected pixels  [pixel]
    YMIN_IMAGE      Minimum y-coordinate among detected pixels  [pixel]
    XMAX_IMAGE      Maximum x-coordinate among detected pixels  [pixel]
    THETA_IMAGE     Position angle (CCW/x)                      [deg]
    FLAGS           Extraction flags                                         
    FLUX_ISOCOR     Corrected isophotal flux                   [count]

output format:
    ellipse(   269.4581,   397.0928,   46.4340,   1.7650,  -69.43)
    
Запускать под оболочкой astroconda (есть конфликт в numpy версиях  
1.21.5 в astropy vs. 1.26.4 в seiferts etc.

"""

import numpy as np
from os import system

#%%
DIR = 'TMP/'
f_res = f'{DIR}k1-impTEST.fts.sx'
f_out = f'{DIR}elik1-imp-field-TEST.reg'

X,Y,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX= np.genfromtxt(f_res, unpack=True)

output = []

for x,y,a,b,th in zip(X,Y,A,B,TH):
    output.append(f'ellipse( {x:.4f} , {y:.4f} , {a:.4f} , {b:.4f},{th:.2f}')

np.savetxt(f_out, output, fmt='%s')

#%%

fn = 'XY.txt'
arr = np.vstack((X, Y, FLUX))
np.savetxt(f'{DIR}{fn}', arr.T, fmt=('%s %s %s'), header='X_CENTER  Y_CENTER  AREA')

## use  text2fits  utility file from astrometry.net

cmd = f'/usr/local/astrometry/bin/text2fits {DIR}{fn}  {DIR}{fn.split(".")[0]}.fits'
system(cmd)