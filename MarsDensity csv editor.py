from math import exp
import sys
import numpy as np
import csv

alt=np.arange(0,500000,1) # In km with step of 1 m so every 0.001 km
altkm=alt/1000
den=[]

tlist=[]

def density(h):
   if h <= 7000:
      tk = -31-0.000998*h+273.15
      pkilo=0.699*exp(-0.00009*h)
   else:
      tk = -23.4 - 0.00222 * h +237.15
      if tk<=3:
         tk=3
      pkilo=0.699*exp(-0.00009*h)
      
   d = pkilo / (0.1921 * tk)
   return d


for i, h in enumerate(alt):
   den.append(density(h))


with open("MarsDensity.csv", "w", newline='') as f:
   fnames=['Altitude', 'Density']
   writer=csv.DictWriter(f, fieldnames=fnames)
   writer.writeheader()
   
   for i,h in enumerate(alt):
      writer.writerow({'Altitude': altkm[i], 'Density':den[i]})

