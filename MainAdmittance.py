#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 08:45:13 2026

Computes and combines admittance based on the charge fluctuaction and forces acting on the ions
See https://arxiv.org/abs/2607.00932 for details

Based on an orginal development of Giovanni Pireddu 

https://github.com/gpireddu/Q2Z

@author: pdesma
"""

import numpy as np

from funcAdmittance import WKACF,FilonLagrange,AdmFromQ,VarRedTot

# Constants
k = 1.3806485279e-23 #J.K-1
e = 1.602176620898e-19 #C
eps0=8.8541878128e-12
eV=1.602177e-19 # Joule to eV
au2A=5.2917721092e-1
epse_e2_pereV_perang=eps0/(e**2)*eV*1e-10
eps_r=78
hartree2J=4.3597482e-18
temperature= 298
beta= 1/(k*temperature)

# Geometry 
area=1.288037300000000E+02*6.971707000000001E+01*(5.2917721092e-1)**2
L=39.72

V=L*area*1e30

# System and simulation characteristics

tstep=5e-15 
rate=50

name_global="1M0ltf"
taufit=5e-10
Leff=((L)*1e-10)     
D1=1.1262314150958323e-09
D2=1.1262314150958323e-09
q1=1
q2=-1
ltf=np.array([0,1,2,0,0,0])
N=102

lb=e**2/(4*np.pi*eps0*eps_r)*beta #bjerrum
# Folder names
file="./"

ionic_file="/ionic_current.out.gz"
charge_file="/total_charges.out.gz"



# Caracteristics of the block averaging
corr_len= 1.1e-9#length of each average in SI

eqtime=int(1e-10//(tstep*rate))# equilibration time
block=int(corr_len//(tstep*rate))
if block%4 ==0 or block%4 ==1 : #insure oddity of the half
    block+=2
# frequency choice
nfreq=101
hif= (2*np.pi)/((tstep*rate)*3)
lof= (2*np.pi)/(block*(tstep*rate))
#shifting freqs for bloc average



#Computes the admittance with both methods for all extract and save it

freq= np.logspace(np.log10(lof),np.log10(hif),nfreq)

nameapp="len%.1enfreq%d"%(corr_len,nfreq)


# Read the output files, here the sum of the forces for each ion species
Data=np.loadtxt(file+'/'+ionic_file,comments=['*',"#"]) 
ind= (np.arange(0,len(Data)))
if len(ind)%2==1:
    ind=ind[:-1]
    Data=Data[:-1,:]
#unit conversion + conversion in flux
flux=(Data[ind%2==0,2]*q1*D1+Data[ind%2==1,2]*q2*D2)*hartree2J/(au2A*1e-10)*beta*e
time= (np.arange(0,len(flux)))*rate*tstep
flux=flux[eqtime:]
time=time[eqtime:]
time-=time[0]
time=np.abs(time)


nsplit=int(len(time)/block)
   
AdmI=np.zeros((nfreq,nsplit),dtype="complex128")
FFTI=np.zeros((nfreq,nsplit),dtype="complex128")

print("block %d,len %d, split %d"%(block,len(time),nsplit))
   
NEpart=beta * e**2* (51) * (D1 + D1)/(Leff)**2 
for isplit in range(nsplit):
    start= isplit * block
    stop= (isplit+1) * block
    
    # Autocorrelation 
    wQACF, redTime = WKACF(flux[start:stop],time[start:stop])
    time=time[:block//2]
    print(start,stop,len(wQACF),block)
   
    toIntegrate=(wQACF*beta)/(Leff**2)
    FFTI[:,isplit]=FilonLagrange(np.zeros(len(freq),dtype=complex),freq,
                                  time,toIntegrate)
    
    AdmI[:,isplit]=NEpart-FFTI[:,isplit]


# # Charge based
# https://github.com/gpireddu/Q2Z/tree/main
# using the method of Giovnni Pireddu with a tweaked Laplace transform
# using pyfilon 1.3.0 , https://pypi.org/project/pyfilon/ of Alex Room
# https://github.com/alexhroom/pyfilon for more information
Data= np.loadtxt(file+'/'+charge_file,comments=['*',"#"])            
Charges= Data[eqtime:,1]*e
time= (np.arange(0,len(Charges)))*rate*tstep
del Data

print("block %d,len %d, split %d"%(block,len(time),nsplit))

Adm=np.zeros((nfreq,nsplit),dtype="complex128")
for isplit in range(nsplit):

    start= isplit * block
    stop= (isplit+1) * block
    
    # Autocorrelation
    wQACF, redtime = WKACF(Charges[start:stop],time[start:stop])
    

    Adm[:,isplit]=  beta *AdmFromQ(time[:block//2]-time[0],wQACF, freq)


#%%  Perform the averages for each configuration and performs reduce variate
ltf=np.array([0])
NEpart=beta * e**2* (51) * (D1 + D1)/(Leff)**2 

# read all files and store them in memory

nrep=len(FFTI[0,:])

#Compute the variance
varIimFFTI=np.sqrt(np.var(-FFTI.imag,axis=1,ddof=1))/np.sqrt(nrep)*1j*2
varIrealFFTI=np.sqrt(np.var(-FFTI.real,axis=1,ddof=1))/np.sqrt(nrep)*2
varIimAdm=np.sqrt(np.var(-Adm.imag,axis=1,ddof=1))/np.sqrt(nrep)*1j*2
varIrealAdm=np.sqrt(np.var(-Adm.real,axis=1,ddof=1))/np.sqrt(nrep)*2
AdmMean=np.mean(Adm,axis=1)
FFTIMean=np.mean(FFTI,axis=1)

name_global= ""

nrep=min(len(Adm[0,:]),len(AdmI[0,:]))

optAdmTot,optVarTot,lbdTot=VarRedTot(freq,Adm,NEpart-FFTI)
optVar1=2*(np.sqrt(optVarTot.real)+1j*np.sqrt(optVarTot.imag))/np.sqrt(nrep)
   

np.savetxt(file+"/Adm%s%s.data"%(name_global,nameapp),optAdmTot)
np.savetxt(file+"/VarAdm%s%s.data"%(name_global,nameapp),optVar1)


 
