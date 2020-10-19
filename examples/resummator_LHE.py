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
import pylhe
import ngl_resum as ngl

parser = argparse.ArgumentParser(description='This code shows how to '\
        'use ngl_resum in combination with LHE-files, considering '\
        'top-pair production. First, each event gets tested whether '\
        'it fulfills the conditions of Table 1 from the ATLAS paper '\
        'arXiv:1203.5015 [hep-ex] and then showers the dipoles with '\
        'the outside region defined by the symmetric rapidity gap '\
        'from -y to y with areas around the bottom quarks cut away. '\
        'Similar code was used to resum the non-global logarithms in '\
        'Section 5 of arXiv:2006.00014')
parser.add_argument('-f','--file', help='lhe event file to shower',\
        default=None,required=True)
parser.add_argument('-y','--ymax', help='ymax of outside region', \
        default=0.8, type=float)
parser.add_argument('-n','--nsh', help='number of showers per dipole', \
        default=100, type=int)
parser.add_argument('-t','--tmax', help='maximal shower time tmax', \
        default=0.1, type=float)
parser.add_argument('-m','--nbins', help='number of bins in hists', \
        default=100, type=int)
parser.add_argument('-c','--cutoff', help='cutoff of shower', \
        default=5, type=float)
parser.add_argument('-s','--seed', help='random seed', \
        default=None, type=int)
parser.add_argument('-b','--break', help='stop after so many events', \
        default=100000, type=int)
args = vars(parser.parse_args())
eventFile=args['file']

if not(args['seed'] is None) : np.random.seed(args['seed'])


showerCutoff=float(args['cutoff'])
nbins=int(args['nbins'])
tmax=float(args['tmax'])
nsh=int(args['nsh'])


def _outside(self,v):
    jetaxis1=self.event.outgoingBottom[0]/self.event.outgoingBottom[0].e
    jetaxis2=self.event.outgoingBottom[1]/self.event.outgoingBottom[1].e
    jetRadius=0.4
    rapRangeMax=float(args['ymax'])
    rapRangeMin=0.0
    return (v.R2(jetaxis1)>jetRadius**2) and \
           (v.R2(jetaxis2)>jetRadius**2) and \
           (abs(v.rap)<rapRangeMax) and (abs(v.rap)>=rapRangeMin)

def validEvent(ev): # ev is the ngl.Event we want to test
    
    # check whether we have the necessary particles
    if ev.intermediateTop == None : return False
    if ev.outgoingBottom == None : return False
    if (ev.outgoingElectron == None) and (ev.outgoingMuon == None): \
            return False
    if len(ev.intermediateTop) != 2 :  return False
    if len(ev.outgoingBottom) != 2 : return False

    momentaLeptonsOut=[]
    momentaNeutrinoOut=[]
    
    electronmuonevent=True
    if not ev.outgoingElectron==None:
        for i in ev.outgoingElectron:
            momentaLeptonsOut.append(i)
            # checks on electron(s)
            if i.eT< 25: return False
            if abs(i.rap)>2.47: return False
        for i in ev.outgoingENeutrino:
            momentaNeutrinoOut.append(i)
    else:
        electronmuonevent=False
            

    if not ev.outgoingMuon==None:
        for i in ev.outgoingMuon:
            momentaLeptonsOut.append(i)
            # checks on muon(s)
            if i.pT< 20: return False
            if abs(i.rap)>2.5: return False
        for i in ev.outgoingMNeutrino:
            momentaNeutrinoOut.append(i)
    else:
        electronmuonevent=False
    
    # check number of leptons ans neutrinos
    if len(momentaLeptonsOut) != 2 : return False
    if len(momentaNeutrinoOut) != 2 : return False
    
    dileptonmass=np.sqrt((momentaLeptonsOut[0]+momentaLeptonsOut[1])*\
                            (momentaLeptonsOut[0]+momentaLeptonsOut[1]))
    missingMomentum=(momentaNeutrinoOut[0]+momentaNeutrinoOut[1])
        
    if not electronmuonevent:
        # checks on "missing momenta" (neutrinos) and dilepton mass
        if missingMomentum.eT<40 : return False
        if (dileptonmass<15 or abs(dileptonmass-91)<10) : return False
    else:
        # check on visible transverse momentum
        if (momentaLeptonsOut[0].pT+momentaLeptonsOut[1].pT+\
            ev.outgoingBottom[0].pT+ev.outgoingBottom[1].pT)<130:
                return False

    # checks on bottom quarks
    for i in ev.outgoingBottom:
        if i.pT<25: return False
        if abs(i.rap)>2.4: return False
        for j in momentaLeptonsOut:
                    if i.R2(j)<0.4**2: return False

    return True # only gets reached, if no check failed.
    


evtFile = pylhe.readLHE(eventFile)


fullResultLL=ngl.Hist(nbins,tmax,errorHistCalc=True)
fullNGL1Loop=0.
fullNGL1LoopSq=0.
fullNGL2Loop=0.
fullNGL2LoopSq=0.

eventWeight=0.

numberEvents=0
numberValidEvents=0

timeStart = time.time()

for event in evtFile:
    numberEvents+=1
    
    
    ev=ngl.Event(eventFromFile=event,productionDipoles='intermediate',\
                    decayDipoles=False)
    
    if not eventWeight > 0:
        eventWeight=ev.weight
    if not eventWeight==ev.weight:
        print("Warning: events not of equal weight!")

    if validEvent(ev):
        numberValidEvents+=1
        
        outsideRegion=ngl.OutsideRegion(ev)
        outsideRegion.outside = _outside.__get__(outsideRegion,\
            ngl.OutsideRegion)
        shower=ngl.Shower(ev,outsideRegion,nsh,nbins,tmax,showerCutoff)
        shower.shower()
        fullResultLL+=shower.resLL
        fullNGL1Loop+=shower.ngl1Loop
        fullNGL1LoopSq+=shower.ngl1LoopSq
        fullNGL2Loop+=shower.ngl2Loop
        fullNGL2LoopSq+=shower.ngl2LoopSq
            
    if numberEvents >= int(args['break']):break
    
print('runtime=', time.time()-timeStart,' sec')    
print("of ", numberEvents," events, ", numberValidEvents," were valid.")
print("Weight of each event:", eventWeight)
print('\n\n'  )

print('*************************************')
print('*  t       LL(t)          dS(t)     * ')
print('*************************************\n')
print('*** Binned Result ***\n\n')


for i in range(0,fullResultLL.nbins):    
    print( round(fullResultLL.centerBinValue[i],4),' ', \
                fullResultLL.entries[i]/numberValidEvents,' ', \
                np.sqrt(fullResultLL.squaredError[i])/numberValidEvents)

print('\n'  )  
snlo=fullNGL1Loop/numberValidEvents
snloError=np.sqrt((fullNGL1LoopSq/numberValidEvents-\
                        (fullNGL1Loop/numberValidEvents)**2)\
                        /(nsh*numberValidEvents))
print('snlo=',snlo)
print('snloError=',snloError)

print('\n')

snnlo=fullNGL2Loop/numberValidEvents+0.5*snlo**2
#Error(snnlo)=|d(snnlo)/d(fullNGL2Loop)*Error(fullNGL2Loop)|
#               + |d(snnlo)/d(snlo)*Error(snlo)|
snnloError=abs(np.sqrt((fullNGL2LoopSq/numberValidEvents-\
                (fullNGL2Loop/numberValidEvents)**2)/\
                (nsh*numberValidEvents)))\
                +abs(snlo*snloError)
print('snnlo=',snnlo)
print('snnloError=',snnloError)
print('\n')
