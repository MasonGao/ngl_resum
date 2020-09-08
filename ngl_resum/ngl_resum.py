#! /usr/bin/env python

########################################################################
#                                                                      #
# Backbone of the Resummation of non-global logarithms                 #
#                                                                      #
# If using ngl_resum, please cite                                      #
#               doi:10.1007/JHEP09(2020)029                            #
#               https://inspirehep.net/literature/1798660              #
#                                                                      #
########################################################################

__author__ = 'Marcel Balsiger'
__email__ = 'marcel.balsiger@hotmail.com'
__date__ = 'August 1, 2020'

import numpy as np
import math
import physt

# ==================================================================== #
# Various Constants                                                    #
# ==================================================================== #

# Nc = number of colors
Nc=3

# assume p to be a massive FourVector, if p*p>cutMassive
cutMassive=1e-7

# assume p to have infinite rapidity, if (1-|cos(theta)|)> cutInfRap
cutInfRap=1e-8

# do not boost p to q, if (p-q)*(p-q)<cutBoost 
cutBoost=1e-10

# assume p = q if (p-q).e**2+...+(p-q).pz**2<toleranceSameVector
toleranceSameVector=1e-10
                                                                        
# ==================================================================== #
# FourVectors                                                          #
# ==================================================================== #
class FourVector:  
    """
    A class which represents four-vectors


    Attributes
    ----------
    e : float
        energy component of the four vector
    px : float
        px-component of the four vector
    py : float
        py-component of the four vector
    pz : float
        pz-component of the four vector
    m : float
        mass of the four vector
    vec : np.array
        fourvector as np.array([e,px,py,pz]), mainly used for boosting
    absSpace : float
        absolute value of spatial components np.sqrt(px**2+py**2+pz**2)
    beta : float
        velocity of four vector (e,e*beta*\vec{n})
    theta : float
        angle between beam axis (z-axis) and four vector
    phi : float
        angle of the vector in the plane transversal to beam (z-axis)
    pT : float
        transverse momentum of four vector (beam = z-axis)
    eT : float
        transverse energy of four vector (beam = z-axis)
    pseudorap : float
        pseudorapidity of vector in relation to beam along z-axis
        is float("inf") or float("-inf") for vectors along beam
    rap : float
        rapidity of vector in relation to beam along z-axis
        is float("inf") or float("-inf") for vectors along beam
    
    Methods
    -------
    costheta(other) : float
        returns the cosine of the angle between instance and parameter
    isMassive() : bool
        returns True if squared four vector is larger than cutMassive
    isMassless() : bool
        returns True if four vector is not massive   
    isSame(other): bool
        returns True, if the sum of the squared difference of all 
        components (e,px,py,pz) is smaller than a tolerance
    tensorProd(other): np.array(4x4)
        returns the tensor product of the instance with parameter
    R2(other) : float
        returns the squared angular distance R^2 between inst and param 
    """ 
    
    def __init__(self, e: float, px: float, py: float, pz: float):
        """
        Initialization of the FourVector with all its attributes
        
        Parameters
        ----------
        e : float
            The energy of the FourVector about to be created
        px : float
            The px-component of the FourVector about to be created
        py : float
            The py-component of the FourVector about to be created
        pz : float
            The pz-component of the FourVector about to be created
        """
        
        self.e=1.*e
        self.px=1.*px
        self.py=1.*py
        self.pz=1.*pz
        
        self.vec=np.array([e,px,py,pz])
        
        self.absSpace= np.sqrt(px**2+py**2+pz**2)
        
        if self.isMassive(): self.m=np.sqrt(abs(e*e-self.absSpace**2))
        else: self.m=0.
        
        if e>0: self.beta=self.absSpace/e
        else: self.beta=float("inf")
    
    
        if self.absSpace>0: 
            self.theta=np.arccos(pz/self.absSpace)
            self.eT=e*np.sqrt(1.-pz**2/self.absSpace**2)
        else: 
            self.theta=np.arccos(1.)
            self.eT=0.
    

    
        self.phi=math.atan2(py,px) 
        
        self.pT=np.sqrt(px**2+py**2) 
        
        if (1-abs(np.cos(self.theta)))> cutInfRap and e>0:
            self.pseudorap=-np.log(np.tan(self.theta/2))
            self.rap=1./2.*np.log((e+pz)/(e-pz))
        else:
            if pz>0:
                self.pseudorap=float("inf")
                self.rap=float("inf")
            else:
                self.pseudorap=float("-inf")
                self.rap=float("-inf")
            
    def __repr__(self):
        """
        Representation of FourVector is [e,px,py,pz]
        
        Returns
        -------
        str 
            each component of the FourVectors in brackets
        """
        
        s = f"[{self.e},{self.px},{self.py},{self.pz}]"
        return s
        
    def __add__(self,other):
        """
        Adding two FourVectors returns new FourVector with added comps 
        
        Returns
        -------
        FourVector 
            each component of other added to components of self
        """
        
        e=self.e+other.e
        px=self.px+other.px
        py=self.py+other.py
        pz=self.pz+other.pz
        r=FourVector(e,px,py,pz)
        return r
        
    def __sub__(self,other):
        """
        Subtracting two FourVectors returns new FourVector with 
        subtracted components 
        
        Returns
        -------
        FourVector 
            each component of other subtracted from components of self
        """
        
        e=self.e-other.e
        px=self.px-other.px
        py=self.py-other.py
        pz=self.pz-other.pz
        r=FourVector(e,px,py,pz)
        return r
        
    def __mul__(self,other):
        """
        Multiplying FourVector by
            - Fourvector: returns scalarproduct with metric [1,-1,-1,-1]
            - Scalar: returns FourVector with each component multiplied
        
        Returns
        -------
        FourVector or float
            if other is float each component of self multiplied by other 
            if other is FourVector scalarproduct with self     
        """

        
        if isinstance(other,FourVector):
            r=self.e*other.e-self.px*other.px-self.py*other.py\
                                                    -self.pz*other.pz
        else:
            e=self.e*other
            px=self.px*other
            py=self.py*other
            pz=self.pz*other
            r=FourVector(e,px,py,pz)
        return r
        
    def __rmul__(self,other):
        """
        Multiplying Scalar by Fourvector returns FourVector with each 
        component multiplied
        
        Returns
        -------
        FourVector
            each component of self multiplied by other 
        """
        e=self.e*other
        px=self.px*other
        py=self.py*other
        pz=self.pz*other
        r=FourVector(e,px,py,pz)
        return r
    
    def __truediv__(self,other):
        """
        Division of FourVector by Scalar returns FourVector with each 
        component divided
        
        Returns
        -------
        FourVector
            each component of self divided by other 
        """
        
        e=self.e/other
        px=self.px/other
        py=self.py/other
        pz=self.pz/other
        r=FourVector(e,px,py,pz)
        return r
    
    def tensorProd(self,other):
        """
        Tensorial product of two FourVectors, instance and parameter
        
        The FourVector given as a parameter is treated as a covector
        such that instance.tensorProd(other) returns self.metric.other
        with metric = diag(1,-1,-1,-1)
        
        Parameters
        ----------
        other: FourVector
            the FourVector you plan to make the tensorial product with
            
        Returns
        -------
        np.array(4x4)
            tensor product of instance with parameter incl metric
        """
        
        r= np.array([\
            [self.e*other.e,-self.e*other.px,\
                                    -self.e*other.py,-self.e*other.pz],\
            [self.px*other.e,-self.px*other.px,\
                                    -self.px*other.py,-self.px*other.pz],\
            [self.py*other.e,-self.py*other.px,\
                                    -self.py*other.py,-self.py*other.pz],\
            [self.pz*other.e,-self.pz*other.px,\
                                    -self.pz*other.py,-self.pz*other.pz]])
        return r
    
    def R2(self,other):
        """
        Returns the squared angular distance R^2 between the instance
        and the parameter
        
        R^2 as defined by Atlas between instance and parameter:
        R2=|self.phi-other.phi|^2+|self.pseudorap-other.pseudorap|^2
        
        Parameters
        ----------
        other : FourVector
            the FourVector you want to measure the distance to
            
        Returns
        -------
        float
            squared angular distance R^2 between instance and parameter
        """
            
        dphi=abs(self.phi-other.phi)    
        # delta phi = 0 ... Pi !
        if dphi>np.pi:        
            dphi=2*np.pi-dphi       
        return  (self.pseudorap-other.pseudorap)**2+dphi**2

    def isMassive(self):
        """
        Returns the answer to whether the instance is massive
        
        Returns
        -------
        bool
            whether instance is massive 
        """
        
        return abs(self*self)>cutMassive

    def isMassless(self):
        """
        Returns the answer to whether the instance is massless
        
        Returns
        -------
        bool
            whether instance is massless 
        """
        
        return not self.isMassive()
    
    def isSame(self,other):
        """
        Returns the answer to whether the instance and the parameter
        are the same FourVectors (up to some precision)
        
        Parameters
        ----------
        other: FourVector
            the FourVector you want to test
        
        
        Returns
        -------
        bool
            whether instance is the same as the parameter
        """
        
        dv=self-other
        sdv=dv.e*dv.e+dv.px*dv.px+dv.py*dv.py+dv.pz*dv.pz
        return sdv<toleranceSameVector
    
    def costheta(self,other):
        """
        Returns the cosine of the spatial angle between the instance 
        and the parameter
        
        Parameters
        ----------
        other: FourVector
            the FourVector you want to measure the distance to
        
        
        Returns
        -------
        float
            cosine of spatial angle between instance and parameter
        """
        
        r= self.px/self.absSpace * other.px/other.absSpace\
            + self.py/self.absSpace * other.py/other.absSpace\
            + self.pz/self.absSpace * other.pz/other.absSpace
        return r

# ==================================================================== #
# Boosting                                                             #
# ==================================================================== #
class Boost:
    """
    A class which contains all information on boosts to and from 
    the cms frame of two FourVectors (given in the labframe)


    Attributes
    ----------
    X : np.array(4x4)
        Rotation to bring the sum of the FourVectors along the x-axis
    B : np.array(4x4)
        Boost to bring the FourVectors back-to-back
    Z : np.array(4x4)
        Rotation to bring the back-to-back FourVectors along the z-axis
    LABtoCMS : np.array(4x4)
        Contains the full boost Z.B.X
    CMStoLAB : np.array(4x4)
        Contains the full inverse boost (Z.B.X)^-1
                  
    Methods
    -------
    boost(v , vp) : np.array(4x4)
        returns the Householder transformation to get from v to vp
    boostCMStoLAB(v) : np.array(4x4)
        returns the FourVector v given in cms frame boosted to lab frame
    boostLABtoCMS(v) : np.array(4x4)
        returns the FourVector v given in lab frame boosted to cms frame
    """ 
    
    def __init__(self, p1: FourVector, p2: FourVector):
        """
        Initialization of the Boost including constructing the boost
        from p1 and p2 in the lab frame to back-to-back along z-axis
        
        Parameters
        ----------
        p1 : FourVector
            The FourVector in lab frame that will be boosted to the 
            positive z-axis in the cms frame [a,0,0,a]
        p2 : FourVector
            The FourVector in lab frame that will be boosted to the 
            negative z-axis in the cms frame [a,0,0,-a]
        """
        
        P=p1+p2
        self.X=np.identity(4)
        self.B=np.identity(4)
        self.Z=np.identity(4)
        if abs(P.absSpace)>0. and abs(P.beta)<1:
            n=FourVector(1,P.px/P.absSpace,P.py/P.absSpace,P.pz/P.absSpace)
            nTar=FourVector(1,1,0,0)
            self.X= self.boost(n,nTar)
            
            gamma=1./np.sqrt(1.-P.beta**2)
            self.B= np.array([np.array([gamma,-P.beta*gamma,0,0]),\
                        np.array([-P.beta*gamma,gamma,0,0]),\
                        np.array([0,0,1,0]),np.array([0,0,0,1])])
       
        p1Tilde=np.dot(np.dot(self.B,self.X),p1.vec)
        p=FourVector(p1Tilde[0],p1Tilde[1],p1Tilde[2],p1Tilde[3])
        n=FourVector(1,p.px/p.absSpace,p.py/p.absSpace,p.pz/p.absSpace)
        nTar=FourVector(1,0,0,1)
        self.Z=self.boost(n,nTar)
        
        self.LABtoCMS=np.dot(self.Z,np.dot(self.B,self.X))
        self.CMStoLAB=np.linalg.inv(self.LABtoCMS)
    
    def __repr__(self):
        """
        Representation of the Boost is just the matrice to boost from 
        lab to cms 
        
        Returns
        -------
        str 
            matrice of the boost from lab to cms frame 
        """
        
        return str(self.LABtoCMS)
    
    def boostLABtoCMS(self,v):
        """
        Boosts the FourVector v given in the lab frame to the cms frame 
        of p1 and p2 given in the initialization
        
        Parameters
        ----------
        v: FourVector
            the FourVector in the lab frame to boost to cms frame
        
        Returns
        -------
        FourVector
            v, but boosted from lab to cms frame
        """
                
        if isinstance(v,FourVector):
            r=np.dot(self.LABtoCMS,v.vec)
        return FourVector(r[0],r[1],r[2],r[3])
    
    def boostCMStoLAB(self,v):
        """
        Boosts the FourVector v given in the cms frame of p1 and p2 
        given in the initialization to the lab frame 
        
        Parameters
        ----------
        v: FourVector
            the FourVector in the cms frame to boost to lab frame
        
        Returns
        -------
        FourVector
            v, but boosted from cms to lab frame
        """

        if isinstance(v,FourVector):
            r=np.dot(self.CMStoLAB,v.vec)
        return FourVector(r[0],r[1],r[2],r[3])
    
    def boost(self, v, vp):
        """
        Constructs the Householder transformation that boost v to vp   
        
        Parameters
        ----------
        v : FourVector
            the initial FourVector
        vp : FourVector
            the target FourVector
        
        Returns
        -------
        np.array(4x4)
            Householder transformation that boosts v to vp
        """
        
        diff = vp-v
        diffSdiff= diff*diff
        diffTdiff= diff.tensorProd(diff)
        bst =  np.identity(4)
        if abs(diffSdiff) > cutBoost :
            bst += -(2./diffSdiff)*diffTdiff
        return bst

# ==================================================================== #
# Event                                                                #
# ==================================================================== #

class Event:
    """
    A class which reads in all the information on an event from a 
    pylhe.event and generates the color dipoles. May also be fed with 
    single dipoles of the form [FourVector, FourVector] to shower.

    Attributes
    ----------
    dipoles : [[FourVector,...],[FourVector,...],...]
        color dipoles of entire event
    productionDipoles : [[FourVector,...],[FourVector,...],...]
        color dipoles associated to production - not available if Event
        gets created by feeding Dipole
    decayDipoles : [[FourVector,...],[FourVector,...],...]
        color dipoles associated to decay of produced particles - not 
        available if Event gets created by feeding Dipole or 
        disabled decayDipoles
    incomingDown : [FourVector, ...]
        list of all incoming down quarks and down antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingUp : [FourVector, ...]
        list of all incoming up quarks and up antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingStrange : [FourVector, ...]
        list of all incoming strange quarks and strange antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingCharm : [FourVector, ...]
        list of all incoming charm quarks and charm antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingBottom : [FourVector, ...]
        list of all incoming bottom quarks and bottom antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingTop : [FourVector, ...]
        list of all incoming top quarks and top antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingElectron : [FourVector, ...]
        list of all incoming electrons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingENeutrino : [FourVector, ...]
        list of all incoming electron neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingMuon : [FourVector, ...]
        list of all incoming muons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingMNeutrino : [FourVector, ...]
        list of all incoming muon neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingTau : [FourVector, ...]
        list of all incoming taus 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingTNeutrino : [FourVector, ...]
        list of all incoming tau neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingGluon : [FourVector, ...]
        list of all incoming gluons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingPhoton : [FourVector, ...]
        list of all incoming photons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingZBoson : [FourVector, ...]
        list of all incoming Z Bosons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingWBoson : [FourVector, ...]
        list of all incoming W+ Bosons and W- Bosons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    incomingHiggs : [FourVector, ...]
        list of all incoming Higgs 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateDown : [FourVector, ...]
        list of all intermediate down quarks and down antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateUp : [FourVector, ...]
        list of all intermediate up quarks and up antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateStrange : [FourVector, ...]
        list of all intermediate strange quarks and strange antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateCharm : [FourVector, ...]
        list of all intermediate charm quarks and charm antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateBottom : [FourVector, ...]
        list of all intermediate bottom quarks and bottom antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateTop : [FourVector, ...]
        list of all intermediate top quarks and top antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateElectron : [FourVector, ...]
        list of all intermediate electrons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateENeutrino : [FourVector, ...]
        list of all intermediate electron neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateMuon : [FourVector, ...]
        list of all intermediate muons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateMNeutrino : [FourVector, ...]
        list of all intermediate muon neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateTau : [FourVector, ...]
        list of all intermediate taus 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateTNeutrino : [FourVector, ...]
        list of all intermediate tau neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateGluon : [FourVector, ...]
        list of all intermediate gluons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediatePhoton : [FourVector, ...]
        list of all intermediate photons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateZBoson : [FourVector, ...]
        list of all intermediate Z Bosons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateWBoson : [FourVector, ...]
        list of all intermediate W+ Bosons and W- Bosons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    intermediateHiggs : [FourVector, ...]
        list of all intermediate Higgs 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingDown : [FourVector, ...]
        list of all outgoing down quarks and down antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingUp : [FourVector, ...]
        list of all outgoing up quarks and up antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingStrange : [FourVector, ...]
        list of all outgoing strange quarks and strange antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingCharm : [FourVector, ...]
        list of all outgoing charm quarks and charm antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingBottom : [FourVector, ...]
        list of all outgoing bottom quarks and bottom antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingTop : [FourVector, ...]
        list of all outgoing top quarks and top antiquarks 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingElectron : [FourVector, ...]
        list of all outgoing electrons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingENeutrino : [FourVector, ...]
        list of all outgoing electron neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingMuon : [FourVector, ...]
        list of all outgoing muons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingMNeutrino : [FourVector, ...]
        list of all outgoing muon neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingTau : [FourVector, ...]
        list of all outgoing taus 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingTNeutrino : [FourVector, ...]
        list of all outgoing tau neutrinos 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingGluon : [FourVector, ...]
        list of all outgoing gluons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingPhoton : [FourVector, ...]
        list of all outgoing photons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingZBoson : [FourVector, ...]
        list of all outgoing Z Bosons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingWBoson : [FourVector, ...]
        list of all outgoing W+ Bosons and W- Bosons 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    outgoingHiggs : [FourVector, ...]
        list of all outgoing Higgs 
        - not availabe if Event created by feeding dipole 
        - None if no particles of this type in this pylhe.event
    """ 
        
    def __init__(self, 
        eventFromFile=False, 
        productionDipoles: str='outgoing',
        decayDipoles: bool=False,
        feedDipole=False):
        """
        Initialization of the Event created either by feeding a 
        pylhe.event or a dipole
        
        Parameters
        ----------
        eventFromFile : pylhe.event, optional 
            pylhe.event to read in
        productionDipoles : str='outgoing', optional
            either 'outgoing' or 'intermediate', tells the Event with 
            which particles to build the dipoles
        decayDipoles: bool = False, optional
            if True, does add the decay dipoles from the intermediate 
            and outgoing particles of the event
        feedDipole: [FourVector,FourVector], optional
            dipole you want to shower
        """
        
        if eventFromFile:
            
            incoming=[None]*26
            intermediate=[None]*26
            outgoing=[None]*26
            toColorSortProd=[]
            if decayDipoles: toColorSortDec=[]
            
            
            for p in eventFromFile.particles:
                stat=int(p.status)
                color1=int(p.color1)
                color2=int(p.color2)
                energy=float(p.e)
                px=float(p.px)
                py=float(p.py)
                pz=float(p.pz)
                pid=int(p.id)
                           
                
                if pid<0:pid=-pid
                
                if stat==-1:
                    self.addToList(incoming,pid,\
                                            FourVector(energy,px,py,pz))
                    toColorSortProd.append(np.array([color1,color2,\
                          FourVector(1,px/energy,py/energy,pz/energy)]))
                
                if stat==2:
                    self.addToList(intermediate,pid,\
                                            FourVector(energy,px,py,pz))
                    if (productionDipoles=='intermediate' and not \
                                         (color1 == 0 and color2 == 0)):
                        toColorSortProd.append([color1,color2,\
                           FourVector(1,px/energy,py/energy,pz/energy)])
                    if (decayDipoles and not \
                                         (color1 == 0 and color2 == 0)): 
                        toColorSortDec.append([color1,color2,\
                           FourVector(1,px/energy,py/energy,pz/energy)])
                
                if stat==1:
                    self.addToList(outgoing,pid,\
                                            FourVector(energy,px,py,pz))
                    if (productionDipoles=='outgoing' and not \
                                         (color1 == 0 and color2 == 0)):
                        toColorSortProd.append([color1,color2,\
                           FourVector(1,px/energy,py/energy,pz/energy)])
                    if (decayDipoles and not \
                                         (color1 == 0 and color2 == 0)): 
                        toColorSortDec.append([color1,color2,\
                           FourVector(1,px/energy,py/energy,pz/energy)])
            
            self.writeAttributes(incoming,intermediate,outgoing)
            self.productionDipoles=self.colorSort(toColorSortProd)
            self.dipoles=self.productionDipoles
            if decayDipoles: 
                self.decayDipoles=self.colorSort(toColorSortDec)
                for i in range(0,len(self.decayDipoles)):
                    self.dipoles.append(self.decayDipoles[i])
        
        if feedDipole:
            self.dipoles=[feedDipole]
                      
    def __repr__(self):
        """
        Representation of Event is an array of all its dipoles each 
        stored in an array [[FourVector,...],[FourVector,...],...]
        
        Returns
        -------
        str 
            each dipole in the collection of dipoles of Event 
        """
        
        s=""
        for i in self.dipoles:
            s += f"{i}\n"
        return s
        
    def addToList(self,lst,i,toAdd):
        """
        Add particle toAdd to array lst at i th entry
        
        
        Parameters
        ----------
        lst : []
            array to which to add the particle
        i : int
            id of the particle
        toAdd : FourVector
            FourVector of the particle to add        
        """
    
        if lst[i]==None:
            lst[i]=[toAdd]
        else:
            lst[i].append(toAdd)
    
    def colorSort(self,toSort):
        """
        sorts a list of dipoles according to their colors
        
        
        Parameters
        ----------
        toSort : [np.array([int,int,FourVector]),...]
            array of FourVectors with color information to sort
            
        Returns
        -------
        [[FourVector,...],[FourVector,...],...]
            a list of all the dipoles after color sorting them
        """
        
        result = []
        while len(toSort)>0:
            isPureGluonic=True
            for i in range(0,len(toSort)):
                if toSort[i][0]==0 or toSort[i][1]==0:
                    isPureGluonic=False
                    break
            tmp=[]
            nextcolor=max(toSort[i][0],toSort[i][1])
            
            nextMom=toSort.pop(i)[2]
            tmp.append(nextMom)
            
            while nextcolor>0 and len(toSort)>0:
                for i in range(0,len(toSort)):
                    if toSort[i][0]==nextcolor:
                        nextcolor=toSort[i][1]
                        break
                    if toSort[i][1]==nextcolor:
                        nextcolor=toSort[i][0]
                        break
                nextMom=toSort.pop(i)[2]
                tmp.append(nextMom)
                
                
            if isPureGluonic: tmp.append(tmp[0])
            
            result.append(tmp)
            
        return result

    
    def writeAttributes(self,incLst,intLst,finLst):
        """
        writes all the lists into the attributes of the instance        
        
        Parameters
        ----------
        incLst : []
            list of incoming FourVectors
        intLst : []
            list of intermediate FourVectors
        finLst : []
            list of outgoing FourVectors
        """
        
        self.incomingDown=incLst[1]
        self.incomingUp=incLst[2]
        self.incomingStrange=incLst[3]
        self.incomingCharm=incLst[4]
        self.incomingBottom=incLst[5]
        self.incomingTop=incLst[6]
        self.incomingElectron=incLst[11]
        self.incomingENeutrino=incLst[12]
        self.incomingMuon=incLst[13]
        self.incomingMNeutrino=incLst[14]
        self.incomingTau=incLst[15]
        self.incomingTNeutrino=incLst[16]
        self.incomingGluon=incLst[21]
        self.incomingPhoton=incLst[22]
        self.incomingZBoson=incLst[23]
        self.incomingWBoson=incLst[24]
        self.incomingHiggs=incLst[25]
        
        self.intermediateDown=intLst[1]
        self.intermediateUp=intLst[2]
        self.intermediateStrange=intLst[3]
        self.intermediateCharm=intLst[4]
        self.intermediateBottom=intLst[5]
        self.intermediateTop=intLst[6]
        self.intermediateElectron=intLst[11]
        self.intermediateENeutrino=intLst[12]
        self.intermediateMuon=intLst[13]
        self.intermediateMNeutrino=intLst[14]
        self.intermediateTau=intLst[15]
        self.intermediateTNeutrino=intLst[16]
        self.intermediateGluon=intLst[21]
        self.intermediatePhoton=intLst[22]
        self.intermediateZBoson=intLst[23]
        self.intermediateWBoson=intLst[24]
        self.intermediateHiggs=intLst[25]
        
        self.outgoingDown=finLst[1]
        self.outgoingUp=finLst[2]
        self.outgoingStrange=finLst[3]
        self.outgoingCharm=finLst[4]
        self.outgoingBottom=finLst[5]
        self.outgoingTop=finLst[6]
        self.outgoingElectron=finLst[11]
        self.outgoingENeutrino=finLst[12]
        self.outgoingMuon=finLst[13]
        self.outgoingMNeutrino=finLst[14]
        self.outgoingTau=finLst[15]
        self.outgoingTNeutrino=finLst[16]
        self.outgoingGluon=finLst[21]
        self.outgoingPhoton=finLst[22]
        self.outgoingZBoson=finLst[23]
        self.outgoingWBoson=finLst[24]
        self.outgoingHiggs=finLst[25]
        
# ==================================================================== #
# Histograms                                                           #
# ==================================================================== #

class Hist:
    """
    A class which is used to manage the histograms of resummation
    
    Also allows to sum the error of a Histogram, according to 
    (dr/r)^2=(ds1/s1)^2+(ds2/s2)^2
    
    Attributes
    ----------
    nbins : int
        number of bins from 0 to tmax
    tmax : float
        maximal value of t of the histograms
    hist : physt.histogram
        underlying physt histogram
    entries[] : float
        entries of each bin
    lowerBinBoundary[] : float
        lower t value of the bin
    upperBinBoundary[] : float
        upper t value of the bin
    centerBinValue[] : float
        central t value of the bin 
    errorHistCalc : bool
        whether or not to also include an error histogram
    squaredHist : physt.histogram
        underlying physt histogram of error - not available if 
        errorHistCalc == False
    squaredError[] : float
        squared error of each bin - not available 
        if errorHistCalc == False
                  
    Methods
    -------
    addToBin(t,w) : 
        adds weight w to bin containing t
    setZero() :
        sets all the entries to 0
    setOne() :
        sets all the entries to 1
    """ 
    
    def __init__(self,nbins: int,tmax: float, errorHistCalc: bool=False):
        """
        Initialization of the Histogram including setting to zero
        
        Parameters
        ----------
        nbins : int
            number of bins between 0 and tmax
        tmax : float
            maximal value of t of the histograms
        errorHistCalc : bool, optional
            whether or not to also have an error histogram
        """
        self.errorHistCalc=errorHistCalc
        self.nbins=nbins
        self.tmax=tmax
        tstep=tmax/nbins
        binning=[i*tstep for i in range(0,nbins+1)] 
        self.hist=physt.histogram([],binning,dtype="float64")
        self.entries=self.hist.frequencies
        self.lowerBinBoundary=[None]*nbins
        self.upperBinBoundary=[None]*nbins
        self.centerBinValue=[None]*nbins
        for i in range(0,self.nbins):
            self.lowerBinBoundary[i]=self.hist.bins[i][0]
            self.upperBinBoundary[i]=self.hist.bins[i][1]
            self.centerBinValue[i]=\
                     (self.hist.bins[i][0]+self.hist.bins[i][1])/2.
        if errorHistCalc:
            self.squaredHist=physt.histogram([],binning,dtype="float64")
            self.squaredError=self.squaredHist.frequencies
        self.setZero()
    
    def __repr__(self):
        """
        Representation of Hist is 
        
        Returns
        -------
        str 
            the histogram with colums t and entries, one entry per row 
        """
        if self.errorHistCalc:
            
            s= "    t   | entries  |  error  \n"
            s+="--------|----------|---------\n"
            for i in range(0,self.nbins):
                s +=" " + "{:6.4f}".format(self.centerBinValue[i])\
                        +" | "+"{:8.6f}".format(self.entries[i])\
                        +" | "+"{:8.6f}".format(np.sqrt(\
                            self.squaredError[i]))+"\n" 
                                
        else:
            s= "    t   | entries  \n"
            s+="--------|----------\n"
            for i in range(0,self.nbins):
                s +=" " + "{:6.4f}".format(self.centerBinValue[i])\
                        +" | "+"{:8.6f}".format(self.entries[i])+"\n"
        return s
    
    def __add__(self,other):
        """
        Adding two Hist gives a new Hist, each entry as the sum of the 
        two entries of the previous Hists  
        
        If the error is also calculated, we also add the errors
        
        Returns
        -------
        Hist 
            histogram with entries of self and other added up
        """
        
        r=Hist(self.nbins,self.tmax,self.errorHistCalc)
        for i in range(0,r.nbins):
            r.entries[i]=self.entries[i]+other.entries[i]
        
        if self.errorHistCalc : 
            for i in range(0,r.nbins):
                r.squaredError[i]=self.squaredError[i]\
                                        +other.squaredError[i]
        
        return r
        
        
    def __sub__(self,other):
        """
        Subtracting one Hist from another gives a new Hist, each entry 
        of other subtracted from self 
        
        If the error is also calculated, we add the errors
        
        Returns
        -------
        Hist 
            histogram with entries of self minus other
        """
        
        r=Hist(self.nbins,self.tmax,self.errorHistCalc)
        for i in range(0,r.nbins):
            r.entries[i]=self.entries[i]-other.entries[i]
        
        if self.errorHistCalc : 
            for i in range(0,r.nbins):
                r.squaredError[i]=self.squaredError[i]\
                                        +other.squaredError[i]
        return r
    
    def __mul__(self,other):
        """
        Multiplying two Hists or one Hist with a scalar and gives a new 
        Hist, each entry the product of the two entries of the previous 
        Hists or each entry multiplied by the scalar.
        
        If the error is also calculated, we divide each entry of the 
        error (the squared values) by their squared entries in the 
        histogram before adding them and multiplying the sum by the 
        squared multiplied entries. 
        Note that 
        (dr/r)^2=((ds1/s1)^2+(ds2/s2)^2)
        with s1 and s2 the two observables, r = s1*s2 and dx the error
        of the observable x
        
        Returns
        -------
        Hist 
            histogram with entries of self and other multiplied
        """
        
        r=Hist(self.nbins,self.tmax,self.errorHistCalc)
        
        if isinstance(other,Hist):
            for i in range(0,r.nbins):
                r.entries[i]=self.entries[i]*other.entries[i]
            if self.errorHistCalc : 
                for i in range(0,r.nbins):
                    tmpErr=0.
                    if self.entries[i]>0:
                        tmpErr+=self.squaredError[i]/self.entries[i]**2
                    if other.entries[i]>0:
                        tmpErr+=other.squaredError[i]/other.entries[i]**2
                    r.squaredError[i]=tmpErr*r.entries[i]**2
                    
        else:
            for i in range(0,r.nbins):
                r.entries[i]=self.entries[i]*other
                
        return r
        
    def __rmul__(self,other):
        """
        Multiplying scalar with a Hist.
        
        Returns
        -------
        Hist 
            histogram with entries of self and other multiplied
        """
        
        r=Hist(self.nbins,self.tmax,self.errorHistCalc)
        
        for i in range(0,r.nbins):
            r.entries[i]=self.entries[i]*other
                
        return r
        
    def __truediv__(self,other):
        """
        Dividing Hist by a scalar.
        
        Returns
        -------
        Hist 
            hist with entries of self divided by other, error unchanged
        """
        
        r=Hist(self.nbins,self.tmax,self.errorHistCalc)
        
        for i in range(0,r.nbins):
            r.entries[i]=self.entries[i]/other
                
        return r    
      
    def addToBin(self,t,w): 
        """
        Adds the weight w to the entry of the histogram which 
        corresponds to the time t
        
        If the error is also calculated, we add to its histogram the 
        squared weight
        
        Parameters
        ----------
        t : float
            time at which the Histogram needs to be filled
        w : float
            weight to add at the time t
        """
           
        self.hist.fill(t,weight=w)
        if self.errorHistCalc:
            self.squaredHist.fill(t,weight=w*w)
    
    def setZero(self):
        """
        Sets all entries of the Hist to zero
        
        If the error is also calculated we also set its entries to zero   
        """
        
        for i in range(0,self.nbins):self.entries[i]=0
        if self.errorHistCalc:
            for i in range(0,self.nbins): self.squaredError[i]=0
        
    
    def setOne(self):
        """
        Sets all entries of the Hist to one
        
        If the error is also calculated we set its entries to zero as
        the setOne()-function is used to allow multiplication of the 
        histograms of each dipole in an event which is not needed for
        the error (see __mul__) 
        """
        
        for i in range(0,self.nbins):self.entries[i]=1
        if self.errorHistCalc:
            for i in range(0,self.nbins): self.squaredError[i]=0
    
       
# ==================================================================== #
# OutsideRegion                                                        #
# ==================================================================== #
class OutsideRegion:
    """
    A class which contains all information on the outside region that
    needs to be defined in the code.

    OutsideRegion needs to be instanciated for each event seperately, if 
    it depends on the event. If the outside region is not depending on 
    the event, one can also instanciate it once. 
    
    Attributes
    ----------
    event : Event
        contains the event with all its particles if available
                      
    Methods
    -------
    outside(v ): bool
        returns true, if FourVector v belongs to the outside region - 
        needs to be implemented in the code
    """ 
    
    def __init__(self,event: Event =  None):
        """
        Initialization of the class OutsideRegion
        
        Parameters
        ----------
        event : Event, optional
            the Event for which the outside region gets calculated (if 
            outside region depends on the event configuration)
        """ 
        self.event=event
        
    def outside(self, v ):
        """
        returns True, if FourVector v belongs to the outside region
        
        instead of changing this method, you can also implement your own
        outsideregion as follows in the code and then replace the native 
        method by yours. 
        
        First define your function as follows:
        
        def _outside(self,v):
            #
            # do something
            #
            # you can access the particles in your event in your 
            # OutsideRegion by using self.event.outgoingBottom[0] etc
            #
            return (conditions v needs to fulfill if outside)
        
        In your code, instanciate the OutsideRegion with your Event:
        
        outsideRegion=ngl_resum.OutsideRegion(Event)
        
        Then replace this method by your function using
        outsideRegion.outside = \
                 _outside.__get__(outsideRegion,ngl_resum.OutsideRegion)
        
        Parameters
        ----------
        v : FourVector
            the FourVector you want to test whether it lies outside
        """ 
        return False

# ==================================================================== #
# Shower                                                               #
# ==================================================================== #

class Shower:
    """
    the class which contains all information to shower an event and 
    then does the showering by Shower.shower()
    
    Attributes
    ----------
    event : Event
        the Event this Shower will shower
    cutoff : float 
        the collinear cutoff - new vectors can not come closer to 
        massless dipole leg p than that (p*newVec>cutoff)
    outsideRegion : OutsideRegion
        OutsideRegion of the Shower
    nsh : int
        number of showerings per dipole
    nbins : int
        number of bins to have in the final histogram
    tmax : float
        maximal value of t in the final histogram
    fixedOrderExpansion : bool
        whether or not the fixed-order coefficient of the first and 
        second order should be calculated
    virtualSubtracted : bool
        whether or not to subtract one half of the original virtual 
        dipole  corrections - see arXiv:1901.09038, equation (3.5)
    resLL : Hist
        result of the LL Resummation of all the dipoles in the given 
        event with the error
    dipsToShower : [[FourVector,...],[FourVector,...],...]
        array of all dipoles to shower, each given in an array 
    ngl1Loop : float
        first-order coefficient of the fixed-order expansion
        for this showering of the event
        - not available if fixedOrderExpansion=False
    ngl2Loop : float
        second-order non-global coefficient of the fixed-order expansion
        for this showering of the event
        - not available if fixedOrderExpansion=False
    ngl1LoopCounter : int
        number of times written into ngl1Loop
    ngl2LoopCounter : int
        number of times written into ngl2Loop
    ngl1LoopSq : float
        squared first-order coefficient of the fixed-order expansion
        for this showering of the event
        - not available if fixedOrderExpansion=False
    ngl2LoopSq : float
        squared second-order non-global coefficient of the fixed-order 
        expansion for this showering of the event
        - not available if fixedOrderExpansion=False
                      
    Methods
    -------
    shower() : 
        showers each dipole one-by-one
    showerDipole(dip) : [Hist, Hist]
        showers one dipole and returns its LL result and the error
    Wijk(i,j,k): float
        returns the dipole radiator of dipole legs i and j emitting k
    Wijkmassive(i,j,k) : float
        returns the modified dipole radiator with subtracted massive 
        monopoles of dipole i and j emitting k 
    getRapidityBoundaries(u1,u2,p1,p2) : [float,float]
        returns the rapidity boundaries yMin and yMax of the integration 
        for dipole u1 and u2 in Lab frame corresponding to p1 and p2 in
        cms frame
    realEmission(u1,u2) : [FourVector, float]
        generates a new FourVector, emitted by the dipole [u1,u2] with 
        its weight (Rijk/Vij)
    virtualCorrection(u1,u2) : float
        calculates the virtual correction of one dipole    
    """ 
    
    def __init__(self,
        event: Event,
        outsideRegion: OutsideRegion,
        nsh: int=50, 
        nbins: int=100, 
        tmax: float=0.1, 
        cut: float=5.0, 
        fixedOrderExpansion: bool=True,
        virtualSubtracted: bool=False):
        """
        Initialization of the Shower
        
        Parameters
        ----------
        event : Event
            the Event you want to shower
        outsideRegion : OutsideRegion
            the outside region of the showering
        nsh : int, optional
            number of showerings per dipole (default is 50)
        nbins : int, optional
            number of bins in the histogram (default is 100)
        tmax : float, optional
            maximal value of t of the histograms (default is 0.1)
        cut : float, optional
            collinear cutoff given in terms of rapidity of a massless 
            back-to-back dipole (default is 5.0)
        fixedOrderExpansion : bool, optional
            if True, will also calculate the first two coefficients of 
            the fixed-order expansion - to save computing time, turn off
            (default = True)
        virtualSubtracted : bool, optional
            if True, will subtract one half of the virtual correction of
            the original dipole from the virtual corrections - see 
            arXiv:1901.09038, equation (3.5) (default = False)
        """ 
        
        self.event=event
        self.cutoff=abs(1-np.tanh(cut))
        self.outsideRegion=outsideRegion
        self.nsh=nsh
        self.nbins=nbins
        self.tmax=tmax
        self.fixedOrderExpansion=fixedOrderExpansion
        self.virtualSubtracted=virtualSubtracted
        self.resLL=Hist(nbins,tmax, errorHistCalc=True)
        self.resLL.setOne()
            
    def shower(self):
        """
        Showers the event the Shower is initialized with
        """ 
        
        if not self.fixedOrderExpansion and not self.virtualSubtracted:
            self.dipsToShower=self.event.dipoles
        else:    
            tmpDip=[]
            for dip in self.event.dipoles:
                for j in range(0,len(dip)-1):
                    tmpDip.append([dip[j],dip[j+1]])
            
            self.dipsToShower=tmpDip
            self.ngl1Loop=0.
            self.ngl1LoopSq=0.
            self.ngl1LoopCounter=0
            self.ngl2Loop=0.
            self.ngl2LoopSq=0.
            self.ngl2LoopCounter=0
            
        for i in range(0,len(self.dipsToShower)):
            dipole=self.dipsToShower[i]
            showered=self.showerDipole(dipole)
            
            self.resLL*=showered[0]
                    
    def showerDipole(self,dip):
        """
        Showers one dipole or multiple dipoles that can treated as one,
        such as production dipoles from gg > tt, where we have three 
        dipoles at t=0
        
        Parameters
        ----------
        dip : [FourVector,FourVector(,...)]
            The dipole(s) to shower
        
        Returns
        -------
        [Hist,Hist] 
            Hist of the LL and the Error of this showering of the 
            dipole(s)
            dipole(s)
        """
        
        binLL  =Hist(self.nbins,self.tmax,errorHistCalc=True)
        binLLSq=Hist(self.nbins,self.tmax)
        startDipoles=dip
        for ish in range(0,self.nsh):  
            t=0.0
            uDipoles=list(startDipoles)
            vDipoles=[self.virtualCorrection(uDipoles[i],uDipoles[i+1])\
                                     for i in range(0, len(uDipoles)-1)]
            virtualtot=sum(vDipoles)
            
            if self.virtualSubtracted: 
                vsub=virtualtot/2.
                weight=1./((virtualtot-vsub)*self.nsh*\
                                                (self.tmax/self.nbins))
            else:
                weight=1./(virtualtot*self.nsh*(self.tmax/self.nbins))
            
            emissionCounter=1
            while t<self.tmax or emissionCounter<3:  
                
                if self.virtualSubtracted: 
                    dt=-np.log(np.random.random_sample())\
                                                /((virtualtot-vsub))
                else:
                    dt=-np.log(np.random.random_sample())/(virtualtot)
                
                t=t+dt
                
                binLL.addToBin(t,weight)
                binLLSq.addToBin(t,weight*weight)
                virtsum=[]
                tmpsum=0
                
                virtsum=np.cumsum(vDipoles)/virtualtot
                randomdip=np.random.random_sample()
                
                pos=len([ x for x in virtsum if x<randomdip])
                           
                newV=self.realEmission(uDipoles[pos],uDipoles[pos+1])
                newVec=newV[0]
                RijkVij=newV[1]
                if self.fixedOrderExpansion:
                    if emissionCounter==1:
                        if self.outsideRegion.outside(newVec):
                            W123=virtualtot*RijkVij
                            self.ngl1Loop-=W123/self.nsh
                            self.ngl1LoopSq+=(W123)**2/self.nsh
                            self.ngl1LoopCounter+=1
                        else:
                            W123=RijkVij
                            n1=uDipoles[pos]
                            n2=uDipoles[pos+1]
                            n3=newVec
                            
                    if emissionCounter==2:
                        if self.outsideRegion.outside(newVec):
                            #needs to be V3^2, as we have only 
                            #probability V2/V3 to get here, and need 
                            #to multiply by V2 and V3
                            Wn34=virtualtot*virtualtot*RijkVij
                            if n1.isSame(uDipoles[pos]):
                                
                                #if Theta23<Theta13:
                                if (n2*n3<n1*n3):
                                    toSub=\
                                     0.5*W123*Wn34*(1.\
                                     -self.Wijkmassive(n1,n2,newVec)\
                                     /self.Wijkmassive(n1,n3,newVec))\
                                     /self.nsh
                                else:
                                    toSub=0.5*W123*Wn34/self.nsh
                            
                            if n2.isSame(uDipoles[pos+1]):
                              
                                #if Theta13<Theta23:
                                if (n1*n3<n2*n3):
                                    toSub=\
                                     0.5*W123*Wn34*(1.\
                                     -self.Wijkmassive(n1,n2,newVec)\
                                     /self.Wijkmassive(n2,n3,newVec))\
                                     /self.nsh
                                else:
                                    toSub=0.5*W123*Wn34/self.nsh
                                    
                            self.ngl2Loop-=toSub
                            self.ngl2LoopSq+=toSub**2*self.nsh
                            self.ngl2LoopCounter+=1
                                
                if self.outsideRegion.outside(newVec):
                    break
                    
                weight=weight*RijkVij 
                
                uDipoles.insert(pos+1,newVec)
                del vDipoles[pos]
                vDipoles.insert(pos,\
                  self.virtualCorrection(uDipoles[pos],uDipoles[pos+1]))
                vDipoles.insert(pos+1,\
                self.virtualCorrection(uDipoles[pos+1],uDipoles[pos+2]))
                        
                newvirtualtot=sum(vDipoles)
                
                if self.virtualSubtracted:
                    VtotVtotNew=virtualtot/(newvirtualtot-vsub)
                else:
                    VtotVtotNew=virtualtot/newvirtualtot
                
                if VtotVtotNew > 1: 
                    weight*=VtotVtotNew
                if np.random.random_sample()>(VtotVtotNew):
                    break
                
                
                
                virtualtot=newvirtualtot
                emissionCounter+=1
                    
        return [binLL,binLLSq]

    def Wijk(self,i,j,k):
        """
        Dipole radiator
        
        Parameters
        ----------
        i : FourVector
            one of the dipole legs of the dipole radiator
        j : FourVector
            one of the dipole legs of the dipole radiator
        k : FourVector
            radiated gluon from the dipole radiator
        
        Returns
        -------
        float 
            dipole radiator i*j/((i*k)*(k*j))
        """
        
        return i*j/((i*k)*(k*j))

    def Wijkmassive(self,i,j,k):
        """
        Modified dipole radiator with subtracted massive monopoles
        
        Parameters
        ----------
        i : FourVector
            one of the dipole legs of the dipole radiator
        j : FourVector
            one of the dipole legs of the dipole radiator
        k : FourVector
            radiated gluon from the dipole radiator
        
        Returns
        -------
        float 
            dipole radiator with subtracted monopoles
            i*j/((i*k)*(k*j))-0.5(i*i/(i*k)^2+(j*j/(j*k)^2))
        """
        
        r=self.Wijk(i,j,k)
        if i.isMassive(): r -= 0.5*self.Wijk(i,i,k)
        if j.isMassive(): r -= 0.5*self.Wijk(j,j,k)
        return r
 
  
    def getRapidityBoundaries(self,u1,u2,p1,p2):
        """
        Calculates the boundaries of the integration
        
        For the massless-massless dipoles, the cutoff techniques of 
        Appendix A of arXiv 1803.07045 were applied. For 
        massive-massless dipoles, we treat the massive leg like a
        massless one to calculate the rapidity cutoffson the massless
        leg. Massive-massive dipoles do not have a cutoff. 
        Such that we do not have to boost u1 and u2 into the cms frame, 
        they are required to be giben in lab and cms frame (p1, p2).
        
        Parameters
        ----------
        u1 : FourVector
            one of the dipole legs in lab frame
        u2 : FourVector
            the other dipole leg in lab frame
        p1 : FourVector
            one of the dipole legs in cms frame
        p2 : FourVector
            the other dipole leg in cms frame
        
        Returns
        -------
        [float,float] 
            yMin and yMax as the integration boundaries on the rapidity
            of the angular integral in cms frame
        """
        
        M2=2.*(1.-u1.costheta(u2))
        beta=np.sqrt(abs(1-M2/4))
        alpha=(M2-2*self.cutoff)/(2*self.cutoff)  
        etaMax= abs(np.log(beta+np.sqrt(alpha+(beta)**2)))  
        yMax=etaMax
        yMin=-etaMax
        if p1.isMassive(): beta1=p1.beta
        else: beta1=1.
        if p2.isMassive(): beta2=p2.beta
        else: beta2=1.
        if p1.isMassive(): yMax=1./2.*np.log((1.+beta2)/(1.-beta1))
        if p2.isMassive(): yMin=1./2.*np.log((1.-beta2)/(1.+beta1))  
        
        return [yMin,yMax]
  
    def realEmission(self,u1,u2):
        """
        Generates a real emission from the dipole constructed by the 
        legs u1 and u2.
        
        Parameters
        ----------
        u1 : FourVector
            one of the dipole legs in lab frame
        u2 : FourVector
            the other dipole leg in lab frame
        
        
        Returns
        -------
        [FourVector,float] 
            the newly generated vector in the lab frame normed to have 
            energy component e=1. and the weight needed in the parton
            shower resummation Rijk/Vij
        """
    
        bst=Boost(u1,u2)
        
        p1=bst.boostLABtoCMS(u1)
        p2=bst.boostLABtoCMS(u2)
               
        if p1.isMassive(): beta1=p1.beta
        else: beta1=1.
        if p2.isMassive(): beta2=p2.beta
        else: beta2=1.
        ycuts=self.getRapidityBoundaries(u1,u2,p1,p2)
        yMin=ycuts[0]
        yMax=ycuts[1]
        
        phiMin= 0.
        phiMax= 2.*np.pi
        
        phiVal=phiMin+np.random.random_sample()*(phiMax-phiMin)
        yVal=yMin+np.random.random_sample()*(yMax-yMin)
        k3overk0=(np.exp(2.*yVal) -1.)/(beta1*np.exp(2.*yVal) + beta2)
        kToverk0=np.sqrt(1.-k3overk0*k3overk0)
        cop=np.cos(phiVal)
        sip=np.sin(phiVal)
        newvec=FourVector(1,  kToverk0*cop, kToverk0*sip, k3overk0)
        newvec=bst.boostCMStoLAB(newvec)
        newvec=newvec/newvec.e
               
        if newvec*newvec>cutMassive: 
            print ("new vector appears massive:" )
            print ("newvec=",newvec)
        
        weight=4.*Nc*(beta1*beta2+1.)/(beta1+beta2)*(yMax-yMin)
               
        W11kFac=0.
        W22kFac=0. 
        # make sure to only subtract the monopolefactor, 
        # if it really is a monopole
        if u1.isMassive(): 
            W11kFac=0.5*self.Wijk(u1,u1,newvec)/self.Wijk(u1,u2,newvec)
        if u2.isMassive(): 
            W22kFac=0.5*self.Wijk(u2,u2,newvec)/self.Wijk(u1,u2,newvec)
        
        weight=weight*(1.- W11kFac-W22kFac)
        weight=weight/(self.virtualCorrection(u1,u2))
           
        return [newvec,weight]


    def virtualCorrection(self,u1,u2):
        """
        Calculates the virtual corrections of a dipole
        
        Parameters
        ----------
        u1 : FourVector
            one of the dipole legs in lab frame
        u2 : FourVector
            the other dipole leg in lab frame
        
        
        Returns
        -------
        float 
            the virtual corrections to the dipole, which is basically 
            an emission radiated off the dipole and integrated over the 
            entire phase space
        """
    
        bst=Boost(u1,u2)
        
        p1=bst.boostLABtoCMS(u1)
        p2=bst.boostLABtoCMS(u2)
        
        
        if p1.isMassive(): beta1=p1.beta
        else: beta1=1.
        if p2.isMassive(): beta2=p2.beta
        else: beta2=1.
        
        ycuts=self.getRapidityBoundaries(u1,u2,p1,p2)
        yMin=ycuts[0]
        yMax=ycuts[1]

        weight= 4.*Nc*(beta1*beta2+1.)/(beta1+beta2)*(yMax-yMin)
        
        if p1.isMassive(): weight-=2.*Nc
        if p2.isMassive(): weight-=2.*Nc

        
        if weight<0: 
            print("negative virtual corrections: ",weight)
            print("u1=",u1)
            print("u1.u1=",u1*u1)
            print("u2=",u2)
            print("u2.u2=",u2*u2)
        
        return weight   
