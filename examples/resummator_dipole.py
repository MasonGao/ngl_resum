#! /usr/bin/env python

########################################################################
#                                                                      #
# Resums the non-global logarithms, needs ngl_resum.py                 #
#                                                                      #
# If using ngl_resum, please cite                                      #
#               doi:10.1007/JHEP09(2020)029                            #
#               https://inspirehep.net/literature/1798660              #
#                                                                      #
########################################################################

__author__ = 'Marcel Balsiger'
__email__ = 'marcel.balsiger@hotmail.com'
__date__ = 'October 19, 2020'

import time
import numpy as np
import argparse
import ngl_resum as ngl

parser = argparse.ArgumentParser(description='This code shows how to '\
        'use ngl_resum to shower a single dipole aligned with the '\
        'z-axis, both legs with velocity b. The outside region is '\
        'defined by the symmetric rapidity gap from -y to y. '\
        'This code was used to produce some of the results in '\
        'Section 4 of arXiv:2006.00014')
parser.add_argument('-b','--beta', help='beta of dipole legs', \
        default=1, type=float)
parser.add_argument('-y','--ymax', help='ymax of outside region', \
        default=0.8, type=float)
parser.add_argument('-n','--nsh', help='number of showerings', \
        default=100, type=int)
parser.add_argument('-t','--tmax', help='maximal shower time tmax', \
        default=0.1, type=float)
parser.add_argument('-m','--nbins', help='number of bins in hists', \
        default=100, type=int)
parser.add_argument('-c','--cutoff', help='cutoff of shower', \
        default=6, type=float)
parser.add_argument('-s','--seed', help='random seed', \
        default=None, type=int)
args = vars(parser.parse_args())

nbins=int(args['nbins'])
tmax=float(args['tmax'])
nsh=int(args['nsh'])
showerCutoff=float(args['cutoff'])
b=float(args['beta'])


if not(args['seed'] is None) : np.random.seed(args['seed'])
   
dipole=[ngl.FourVector(1,0,0,b),ngl.FourVector(1,0,0,-b)]

ev=ngl.Event(feedDipole=dipole)

def _outside(self,vec):
    rapRangeMax=float(args['ymax'])
    rapRangeMin=0.0
    return (abs(vec.rap)<rapRangeMax) and (abs(vec.rap)>=rapRangeMin)

outsideRegion=ngl.OutsideRegion()
outsideRegion.outside = _outside.__get__(outsideRegion,\
                            ngl.OutsideRegion)

shower=ngl.Shower(ev,outsideRegion,nsh,nbins,tmax,showerCutoff)

timeStart = time.time()

shower.shower()

    
print('runtime=', time.time()-timeStart,' sec')    

print('*************************************')
print('*  t       LL(t)          dS(t)     * ')
print('*************************************\n')
print('*** Binned Result ***\n\n')


for i in range(0,shower.resLL.nbins):    
    print( round(shower.resLL.centerBinValue[i],4),' ', \
            shower.resLL.entries[i],' ', \
            np.sqrt(shower.resLL.squaredError[i]))

print('\n\n'  )  
snlo=shower.ngl1Loop
snloError=np.sqrt((shower.ngl1LoopSq-shower.ngl1Loop**2)/(nsh))
print('snlo=',snlo)
print('snloError=',snloError)


print('\n')
snnlo=shower.ngl2Loop+0.5*snlo**2
#Error(snnlo)=|d(snnlo)/d(fullNGL2Loop)*Error(fullNGL2Loop)|
#               + |d(snnlo)/d(snlo)*Error(snlo)|
snnloError=abs(np.sqrt((shower.ngl2LoopSq-shower.ngl2Loop**2)/(nsh)))\
                    +abs(snlo*snloError)
print('snnlo=',snnlo)
print('snnloError=',snnloError)
print('\n')
print('\n')
    
