"""Holds all the technological parameters needed for wire delay extraction, both basic and derived.
"""

from device_geometry import *

nodes = [16, 7, 5, 4, 3.0, 3.1]
#3.0 is 3a in the paper, and 3.1 3b

#Entries from Table 3
MxR = [31.61, 128.69, 151.61, 392.89, 666.38, 396.65]
MxC = [0.22, 0.22, 0.22, 0.22, 0.22, 0.28]
MyR = [18.67, 21.55, 25.12, 75.72, 86.41, 18.67]
MyC = [0.24, 0.24, 0.24, 0.24, 0.24, 0.24]

#Entries from Table 3
MxMx_via = [10.88, 34.75, 39.90, 58.91, 92.92, 92.92]

#Entries from Table 4
stacked_via = [19.18, 30.54, 34.68, 69.78, 88.05, 44.87]

MyP = [80, 76, 72, 50, 48, 80]

max_V_span = [\
              {6 : {1 : 32,  2 : 16, 4 : 8,  8 : 4,  16 : 2, 32 : 1, 64 : 1}},\
              {6 : {1 : 32,  2 : 16, 4 : 8,  8 : 4,  16 : 2, 32 : 1, 64 : 1}},\
              {6 : {1 : 32,  2 : 16, 4 : 8,  8 : 4,  16 : 2, 32 : 1, 64 : 1}},\
              {6 : {1 : 16,  2 : 8,  4 : 4,  8 : 2,  16 : 1, 32 : 1, 64 : 1}},\
              {6 : {1 : 16,  2 : 8,  4 : 4,  8 : 2,  16 : 1, 32 : 1, 64 : 1}},\
              {6 : {1 : 32,  2 : 16, 4 : 8,  8 : 4,  16 : 2, 32 : 1, 64 : 1}},\
             ]

#Taken from k6_N10_mem32K_40nm and scaled using logic_delays/scale_logic_delays.py
td_lut = {\
          40 : {6 : 261e-12},\
          16 : {6 : 94e-12},\
          7  : {6 : 68e-12},\
          5  : {6 : 64e-12},\
          4  : {6 : 64e-12},\
          3.0  : {6 : 61e-12},\
          3.1  : {6 : 61e-12}\
         }
td_ff_su = {\
            40 : 66e-12,\
            16 : 24e-12,\
             7 : 17e-12,\
             5 : 16e-12,\
             4 : 16e-12,\
           3.0 : 15e-12,\
           3.1 : 15e-12\
           }
td_ff_clkq = {\
              40 : 124e-12,\
              16 : 45e-12,\
               7 : 32e-12,\
               5 : 30e-12,\
               4 : 30e-12,\
             3.0 : 29e-12,\
             3.1 : 29e-12\
             }

#I/O is difficult to model correctly and should have no influence on results.
td_inpad = 0e-12
td_outpad = 0e-12

#Mux delays will be overwritten by measured data.
td_ble_mux_lut = 0e-12
td_ble_mux_ff = 0e-12
