#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for MainAdmittance

@author: pdesma
"""
import numpy as np


from pyfilon.pyfilon import filon_tab_sin, filon_tab_cos

def WKACF(x,time):
    """
    ==========================================================================================
    Computes the autocorrelation function of a time series using the Wiener-Khinchin theorem
    ==========================================================================================
    From Giovanni Pireddu Q2Z https://github.com/gpireddu/Q2Z/blob/main/Q2Z.py# Reference: 
    Frequency-Dependent Impedance of Nanocapacitors from Electrode Charge Fluctuations as a Probe of Electrolyte Dynamics
    Giovanni Pireddu and Benjamin Rotenberg
    Phys. Rev. Lett. 130, 098001, 2023 
    DOI: https://doi.org/10.1103/PhysRevLett.130.098001

    Parameters
    ----------
    x : float np.array
        quantity to plot the ACF of
    time : float np.array
        time array.

    Returns
    -------
    ACF :  float np.array
        ACF
    time :  float np.array
        time array.

    """
    Ctt = np.fft.fftn(x)
    CC = Ctt[:] * np.conjugate(Ctt[:])
    CC[:] = np.fft.ifftn(CC[:])
    ACF = (CC[:np.shape(CC)[0]//2]).real /np.shape(CC)[0] 
    time= time[:len(time)//2]
    return ACF, time



def FilonLagrange(DFT,freq,Time,functab):
    """
    
    Perform filon lagrange from library
    Parameters
    ----------
    DFT : numpy complex array
        array to be populated with the result
    freq : numpy float array
        list of frequency
    Time : numpy float array
        Time array (must be evenly spaces)
    functab : numpy float array
        arrazy to compute the Lagrange transform of 
        

    Returns
    -------
    DFT : numpy complex array
        array to be populated with the result

    """
    for i,f in enumerate(freq):
      DFT[i]=( filon_tab_cos(functab, Time[0], Time[-1], f) -
           1j*filon_tab_sin(functab, Time[0], Time[-1], f))
    return DFT
def AdmFromQ(time,QACF,freq):
    """
    ==========================================================================================
    Computes the admittance from the total charge autocorrelation function
    ==========================================================================================
    
    Inspired from Giovanni Pireddu Q2Z https://github.com/gpireddu/Q2Z/blob/main/Q2Z.py# Reference: 
    Frequency-Dependent Impedance of Nanocapacitors from Electrode Charge Fluctuations as a Probe of Electrolyte Dynamics
    Giovanni Pireddu and Benjamin Rotenberg
    Phys. Rev. Lett. 130, 098001, 2023 
    DOI: https://doi.org/10.1103/PhysRevLett.130.098001
    
    Parameters
    ----------
   time :  float np.array
       time array.
    QACF : float np.array
        Charge autocorrelation fonction
    freq : float np.array
        list of frequencies to compute the admittance at

    Returns
    -------
    Adm : complexe np array
        Admittance

    """
    #Forward filtering
    DFT= FilonLagrange(np.zeros(len(freq),dtype=complex),freq,time,QACF)
        
    # DFT=  FilonLagrange2(np.zeros(len(freq),dtype=complex),freq,time,QACF)
         
    Adm= np.zeros(len(DFT),dtype=complex)
    for i in range(len(freq)):
        Adm[i]= ((freq[i]**2)*DFT[i] + 1j * freq[i] * QACF[0] )
    return Adm

def VarRedTot(freq,CAdm,Adm):
    """
    Reduce Variate globally and separtion of real and imaganinary ADMITTANCE
    https://arxiv.org/abs/2607.00932see article in prep

    Parameters
    ----------
    freq : numpy array of float 1D
       Frequency list
    CAdm : numpy array complex 2D 
        Admittance from charge (N realisation of)
    Adm : numpy array 2D
        Admittance from Ion forces (N realisation of)
    Returns
    -------
    optAdm : numpy array complex 2D 
        optimal admittance based on inputs
    optVar : numpy array complex 2D 
        optimal variance of real and imaginary admittance based on inputs
    lambdat :  numpy array of float 1D
        proportion of each contribution
    
    """
    N=min(len(Adm[0,:]),len(CAdm[0,:]))
    Adm=Adm[:,:N]    
    CAdm=CAdm[:,:N]    
    delta=Adm-CAdm

    optAdm=np.zeros(len(freq),dtype='complex128')
    optVar=np.zeros(len(freq),dtype='complex128')
    lambdat=np.zeros(len(freq))
    for i in range(0,len(freq)):
        temp=np.sum(((CAdm[i,:]-np.mean(CAdm[i,:])).real*(delta[i,:]-np.mean(delta[i,:])).real)+
              ((CAdm[i,:]-np.mean(CAdm[i,:])).imag*(delta[i,:]-np.mean(delta[i,:])).imag))/(N-1)
        lbd=-temp/np.var(delta[i,:],ddof=1)
     
        lambdat[i]=lbd
        
        optAdm[i]=np.mean(CAdm[i,:])+np.mean(lbd*delta[i,:])
        optVar[i]=(np.var(CAdm[i,:].real,ddof=1)+2*lbd*np.sum((CAdm[i,:]-np.mean(CAdm[i,:])).real*(delta[i,:]-np.mean(delta[i,:])).real)/(N-1)
                   +lbd**2*np.var(delta[i,:].real,ddof=1)
                  ) + 1j*(np.var(CAdm[i,:].imag,ddof=1)+2*lbd*np.sum((CAdm[i,:]-np.mean(CAdm[i,:])).imag*(delta[i,:]-np.mean(delta[i,:])).imag)/(N-1)
                             +lbd**2*np.var(delta[i,:].imag,ddof=1))
    return optAdm,optVar,lambdat
