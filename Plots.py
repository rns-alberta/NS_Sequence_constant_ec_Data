#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Code to compare various plots given files with different quantities 
of rotating neutron stars, considering ten EOS
"""

import future 
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import future 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import numpy
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from scipy.optimize import curve_fit
# Size of the figure
import numpy, scipy, scipy.optimize
import matplotlib
from mpl_toolkits.mplot3d import  Axes3D
from matplotlib import cm # to colormap 3D surfaces from blue to red
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from scipy.interpolate import make_interp_spline, BSpline
from plotnine import *
from plotnine.data import mpg
graphWidth = 800 # units are pixels
graphHeight = 600 # units are pixels

# 3D contour plot lines
numberOfContourLines = 16



#Mass Best fit

def func(data,C0, C1,C2,C3,C4,C5):
    x = data[0]
    y = data[1]
    return 1+(np.exp(C0*x**2)-1)*(C1+C2*y+C3*y**2+C4*y**3+C5*y**4)
   
def zfromfitM(x,y,C):
    return 1+(np.exp(C[0]*x**2)-1)*(C[1]+C[2]*y+C[3]*y**2+C[4]*y**3+C[5]*y**4) 

def DivSepM(x,y,z,fitpar):
    div = np.zeros(len(x))

    for i in range(len(z)):
      if ((zfromfitM(x[i],y[i],fitpar)!=0)):
       #if np.abs((z[i] - zfromfitM(x,y,fitpar,i))/zfromfitM(x,y,fitpar,i))<10:
        div[i]=(z[i] - zfromfitM(x[i],y[i],fitpar))/z[i]

    return div











#Compactness Best fit


def funcComp(data,C0, C1,C2,C3):
    x = data[0]
    y = data[1]
    return y+(np.log(1-(x/1.1)**3))*(C0*y+C1*y**2+C2*y**3+C3*y**4)

def zfromfitC(x,y,C):
    return y+(np.log(1-(x/1.1)**3))*(C[0]*y+C[1]*y**2+C[2]*y**3+C[3]*y**4)


def DivSepC(x,y,z,fitpar):
    div = np.zeros(len(x))

    for i in range(len(z)):
      if ((zfromfitC(x[i],y[i],fitpar)!=0)):
       #if np.abs((z[i] - zfromfitM(x,y,fitpar,i))/zfromfitM(x,y,fitpar,i))<10:
        div[i]=(z[i] - zfromfitC(x[i],y[i],fitpar))/z[i]

    return div
    
    
    
    
    
    
    
    
    
#Radius Best fi 
    
    
def funcR(data,C0, C1,C2,C3,C4,C5,C6):
    x = data[0]
    y = data[1]
    #return 1+((np.exp(C0*x**2+C7*x**20)-1)-C1*np.log(1-(x/1.1)**4)**2)*(1+C2*y+C3*y**2+C4*y**3+C5*y**4+C6*y**5)
    return 1+(np.exp(C0*x**2)-1-C1*np.log(1-(x/1.1)**4)**2)*(1+C2*y+C3*y**2+C4*y**3+C5*y**4+C6*y**5)
def zfromfitR(x,y,C):
    #return 1+((np.exp(C[0]*x**2+C[7]*x**20)-1)-C[1]*np.log(1-(x/1.1)**4)**2)*(1+C[2]*y+C[3]*y**2+C[4]*y**3+C[5]*y**4+C[6]*y**5) 
    return 1+(np.exp(C[0]*x**2)-1-C[1]*np.log(1-(x/1.1)**4)**2)*(1+C[2]*y+C[3]*y**2+C[4]*y**3+C[5]*y**4+C[6]*y**5) 
def DivSepR(x,y,z,fitpar):
    div = np.zeros(len(x))

    for i in range(len(z)):
      if ((zfromfitR(x[i],y[i],fitpar)!=0)):
       #if np.abs((z[i] - zfromfitR(x,y,fitpar,i))/zfromfitR(x,y,fitpar,i))<10:
        div[i]=(z[i] - zfromfitR(x[i],y[i],fitpar))/z[i]

    return div

Rcoef=7



#General
    
    
def funcGen(data,C0, C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,C34,C35,C36,C37,C38,C39,C40,C41,C42,C43,C44,C45,C46, C47,C48,C49,C50,C51,C52,C53,C54,C55):
    x = data[0]
    y = data[1]
    return C0+C1*x+C2*y+C3*x*y+C4*x**2+C5*y**2+C6*x**3+C7*y**3+C8*x**2*y+C9*x*y**2+C10*x**4+C11*y**4+C12*x*y**3+C13*x**2*y**2+C14*x**3*y+C15*x**5+C16*y**5+C17*x**4*y+C18*x**3*y**2+C19*x**2*y**3+C20*x*y**4 +C21*x**6+C22*y**6+C23*x**5*y+C24*x**4*y**2+C25*x**3*y**3+C26*x**2*y**4+C27*x*y**5+C28*x**7+C29*y**7+C30*x**6*y+C31*x**5*y**2+C32*x**4*y**3+C33*x**3*y**4+C34*x**2*y**5+C35*x*y**6+ C36*x**8+C37*y**8+C38*x**7*y+C39*x**6*y**2+C40*x**5*y**3+C42*x**4*y**4+C43*x**3*y**5+C44*x**2*y**6+C45*x*y**7 +C46*x**9+C47*y**9+C48*x**8*y+C49*x**7*y**2+C50*x**6*y**3+C51*x**5*y**4+C52*x**4*y**5+C53*x**3*y**6+C54*x**2*y**7+C55*x*y**8
    
def zfromfitGen(x,y,C):
    return C[0]+C[1]*x+C[2]*y+C[3]*x*y+C[4]*x**2+C[5]*y**2+C[6]*x**3+C[7]*y**3+C[8]*x**2*y+C[9]*x*y**2+C[10]*x**4+C[11]*y**4+C[12]*x*y**3+C[13]*x**2*y**2+C[14]*x**3*y+C[15]*x**5+C[16]*y**5+C[17]*x**4*y+C[18]*x**3*y**2+C[19]*x**2*y**3+C[20]*x*y**4 +C[21]*x**6+C[22]*y**6+C[23]*x**5*y+C[24]*x**4*y**2+C[25]*x**3*y**3+C[26]*x**2*y**4+C[27]*x*y**5+C[28]*x**7+C[29]*y**7+C[30]*x**6*y+C[31]*x**5*y**2+C[32]*x**4*y**3+C[33]*x**3*y**4+C[34]*x**2*y**5 +C[35]*x*y**6 +C[36]*x**8+C[37]*y**8+C[38]*x**7*y+C[39]*x**6*y**2+C[40]*x**5*y**3+C[42]*x**4*y**4+C[43]*x**3*y**5+C[44]*x**2*y**6+C[45]*x*y**7 +C[46]*x**9+C[47]*y**9+C[48]*x**8*y+C[49]*x**7*y**2+C[50]*x**6*y**3+C[51]*x**5*y**4+C[52]*x**4*y**5+C[53]*x**3*y**6+C[54]*x**2*y**7+C[55]*x*y**8

def DivSepGen(x,y,z,fitpar):
    div = np.zeros(len(x))

    for i in range(len(z)):
      if (z[i]!=0):
       #if np.abs((z[i] - zfromfitGen(x[i],y[i],fitpar))/z[i])<1:
        div[i]=(z[i] - zfromfitGen(x[i],y[i],fitpar))/z[i]

    return div










#Other Best fit relations 
def funcPol(data, C1,C2,C3,C4,C5,C6,C7,C8,C9,C10):
    x = data[0]
    y = data[1]
    #Moment of Inertia
    #return (C1/(C2+y)+C3/(C4+y)**2+C5/(C6+y)**3)#*np.exp(x*C7)#+C7/(C8+y)**4#(C1/y+C2/y**2+C3/y**3+C4/y**4+C5/y**5+C6/y**6+C7/y**7+C8/y**8+C9/y**9+C10/y**10)#*(1+C11*x)#+C12*x**2)#C1*y/(1-C2*y+C3*y**2+C4*y**3)+C5*x#+C5*y*np.log(1-(x/1.1)**2)
    #Binding
    #return C1+C2*y+C3*y**2#/(1+C2*y)#+C3*y**2+C4*y**3+C5*y**4+C6*y**5)#+C5*x
    #Rratio
    return 1+(C1*x**2+C2*np.log(1-(x/1.1)**2))*(1+C3*y+C4*y**2)
    #Vp
    #return C1+C2*y+C3*y**2#(C1*y+C5+C6*y**2)/(1+C2*y+C3*y**2+C4*y**3)
    #Virial
    #return C3*x**2+C2*x +C1

def zfromfitPol(x,y,C,i):
    #return (C[0]/(C[1]+y[i])+C[2]/(C[3]+y[i])**2+C[4]/(C[5]+y[i])**3)#*np.exp(x[i]*C[6])#+C[6]/(C[7]+y[i])**4#(C[0]/y[i]+C[1]/y[i]**2+C[2]/y[i]**3+C[3]/y[i]**4+C[4]/y[i]**5+C[5]/y[i]**6+C[6]/y[i]**7+C[7]/y[i]**8+C[8]/y[i]**9+C[9]/y[i]**10)#*(1+C[10]*x[i])#+C[11]*x[i]**2)#C[0]*y[i]/(1-C[1]*y[i]+C[2]*y[i]**2+C[3]*y[i]**3)+C[4]*x[i]#+
    #return C[0]+C[1]*y[i]+C[2]*y[i]**2#/(1+C[1]*y[i])#+C[2]*y[i]**2+C[3]*y[i]**3+C[4]*y[i]**4+C[5]*y[i]**5)#+C[4]*x[i]
    return 1+(C[0]*x[i]**2+C[1]*np.log(1-(x[i]/1.1)**2))*(1+C[2]*y[i]+C[3]*y[i]**2)
    #return C[0]+C[1]*y[i]+C[2]*y[i]**2#(C[0]*y[i]+C[4]+C[5]*y[i]**2)/(1+C[1]*y[i]+C[2]*y[i]**2+C[3]*y[i]**3)
    #return C[1]*x[i]**2+C[0]*x[i]+C[2]

def DivSepPol(x,y,z,fitpar):
    div = np.zeros(len(x))
    for i in range(len(z)):
      if ((z[i]!=0)):
       #if np.abs((z[i] - zfromfitM(x,y,fitpar,i))/zfromfitM(x,y,fitpar,i))<10:
        div[i]=(z[i] - zfromfitPol(x,y,fitpar,i))/z[i]
    return div









#Mass Best fit Reverse

def func1(data,C0, C1,C2,C3,C4,C5,C6):
    x = data[0]
    y = data[1]
    return 1+(x+C0*x**2+C1*x**3+C2*x**4)*(C3*y+C4*y**2+C5*y**3+C6*y**4)
   
def zfromfitM1(x,y,C):
    return 1+(x+C[0]*x**2+C[1]*x**3+C[2]*x**4)*(C[3]*y+C[4]*y**2+C[5]*y**3+C[6]*y**4) 

def DivSepM1(x,y,z,fitpar):
    div = np.zeros(len(x))

    for i in range(len(z)):
      if ((zfromfitM1(x[i],y[i],fitpar)!=0)):
       #if np.abs((z[i] - zfromfitM1(x,y,fitpar,i))/zfromfitM(x,y,fitpar,i))<10:
        div[i]=(z[i] - zfromfitM1(x[i],y[i],fitpar))/z[i]

    return div












# Size of the figure
fig = plt.figure(figsize = (12,12))

# Files to read the data from
dataGreif0 = np.loadtxt('NS_data_eosGreif0.txt', unpack=True)
dataGreif1 = np.loadtxt('NS_data_eosGreif1.txt', unpack=True)
dataGreif2 = np.loadtxt('NS_data_eosGreif2.txt', unpack=True)
dataGreif3 = np.loadtxt('NS_data_eosGreif3.txt', unpack=True)
dataGreif4 = np.loadtxt('NS_data_eosGreif4.txt', unpack=True)
dataGreif5 = np.loadtxt('NS_data_eosGreif5.txt', unpack=True)
dataGreif6 = np.loadtxt('NS_data_eosGreif6.txt', unpack=True)
dataGreif7 = np.loadtxt('NS_data_eosGreif7.txt', unpack=True)
dataGreif8 = np.loadtxt('NS_data_eosGreif8.txt', unpack=True)
dataGreif9 = np.loadtxt('NS_data_eosGreif9.txt', unpack=True)
dataGreif10 = np.loadtxt('NS_data_eosGreif10.txt', unpack=True)
dataGreif11 = np.loadtxt('NS_data_eosGreif11.txt', unpack=True)
dataGreif12 = np.loadtxt('NS_data_eosGreif12.txt', unpack=True)



#dataPol0nr = np.loadtxt('NS_data_eosPol0nonrot.txt', unpack=True)
#dataPol0K = np.loadtxt('NS_data_eosPol0Kepler.txt', unpack=True)
#dataPol0ins = np.loadtxt('NS_data_eosPol0ins.txt', unpack=True)
#dataPol0J = np.loadtxt('NS_data_eosPol0J.txt', unpack=True)
#dataP0 = np.loadtxt('eosPol0table', unpack=True)

dataPol0 = np.loadtxt('NS_data_eosPol0.txt', unpack=True)
dataPol1 = np.loadtxt('NS_data_eosPol1.txt', unpack=True)
dataPol2 = np.loadtxt('NS_data_eosPol2.txt', unpack=True)
dataPol3 = np.loadtxt('NS_data_eosPol3.txt', unpack=True)
dataPol4 = np.loadtxt('NS_data_eosPol4.txt', unpack=True)
dataPol5 = np.loadtxt('NS_data_eosPol5.txt', unpack=True)
dataPol6 = np.loadtxt('NS_data_eosPol6.txt', unpack=True)
dataPol7 = np.loadtxt('NS_data_eosPol7.txt', unpack=True)
dataPol8 = np.loadtxt('NS_data_eosPol8.txt', unpack=True)
dataPol9 = np.loadtxt('NS_data_eosPol9.txt', unpack=True)
dataPol10 = np.loadtxt('NS_data_eosPol10.txt', unpack=True)
dataPol11 = np.loadtxt('NS_data_eosPol11.txt', unpack=True)
dataPol12 = np.loadtxt('NS_data_eosPol12.txt', unpack=True)
dataPol13 = np.loadtxt('NS_data_eosPol13.txt', unpack=True)
dataPol14 = np.loadtxt('NS_data_eosPol14.txt', unpack=True)
dataPol15 = np.loadtxt('NS_data_eosPol15.txt', unpack=True)
dataPol16 = np.loadtxt('NS_data_eosPol16.txt', unpack=True)
dataPol17 = np.loadtxt('NS_data_eosPol17.txt', unpack=True)
dataPol18 = np.loadtxt('NS_data_eosPol18.txt', unpack=True)
dataAPR = np.loadtxt('NS_data_eosAPR.txt', unpack=True)
dataBBB1 = np.loadtxt('NS_data_eosBBB1.txt', unpack=True)
dataBBB2 = np.loadtxt('NS_data_eosBBB2.txt', unpack=True)
dataL = np.loadtxt('NS_data_eosL.txt', unpack=True)
dataHLPS1 = np.loadtxt('NS_data_eosHLPS1.txt', unpack=True)
dataHLPS2 = np.loadtxt('NS_data_eosHLPS2.txt', unpack=True)
#dataHLPS3 = np.loadtxt('NS_data_eosHLPS3.txt', unpack=True)
dataABPR1 = np.loadtxt('NS_data_eosABPR1.txt', unpack=True)
dataABPR2 = np.loadtxt('NS_data_eosABPR2.txt', unpack=True)
dataQHCD = np.loadtxt('NS_data_eosQHC_D.txt', unpack=True)
dataH0 = np.loadtxt('NS_data_eosH0.txt', unpack=True)
#dataH0 = np.loadtxt('NS_data_eosQ160.txt', unpack=True)
#dataQ1 = np.loadtxt('NS_data_eosQ160.txt', unpack=True)
#dataQ1 = np.loadtxt('NS_data_QMF18.txt', unpack=True)
#dataPol0rho1 = np.loadtxt('NS_data_eosPol0_rho1.txt', unpack=True)
#dataPol0rho2 = np.loadtxt('NS_data_eosPol0_rho2.txt', unpack=True)



# Some constants
G = 6.67e-8
c = 2.99792458e10
Msun = 1



'''

ec = data[0,:]*10**15   
M = data[1,:]           
M0 = data[2,:]          
Mstat = data[3,:]       
#Mmax = data[4,:]        
R = data[5,:]           
Rratio = data[6,:]      
Rstat = data[7,:]       
freq = data[8,:]        
kfreq = data[9,:]       
J = data[10,:]          
T = data[11,:]          
W = data[12,:]          
#Rmax = data[13,:]       
#Qmax = data[14,:]       
Vp = data[15,:] 
'''
data=numpy.concatenate((dataGreif0, dataGreif1, dataGreif2, dataGreif3, dataGreif4, dataGreif5, dataGreif6, dataGreif7, dataGreif8, dataGreif9, dataGreif10,dataGreif11,dataGreif12,dataPol0 ,dataPol1 ,dataPol2 ,dataPol3 ,dataPol4 ,dataPol5 ,dataPol6 ,dataPol7  ,dataPol8 ,dataPol9 ,dataPol10 ,dataPol11 ,dataPol12 ,dataPol13 ,dataPol14 ,dataPol15 ,dataPol16 ,dataPol17 ,dataPol18),axis=1)

dataPol=numpy.concatenate((dataPol0 ,dataPol1 ,dataPol2 ,dataPol3 ,dataPol4 ,dataPol5 ,dataPol6 ,dataPol7  ,dataPol8 ,dataPol9 ,dataPol10 ,dataPol11 ,dataPol12 ,dataPol13 ,dataPol14 ,dataPol15 ,dataPol16 ,dataPol17 ,dataPol18),axis=1)

#dataPol=dataPol14
#12
dataGr=numpy.concatenate((dataGreif0, dataGreif1, dataGreif2, dataGreif3, dataGreif4, dataGreif5, dataGreif6, dataGreif7, dataGreif8, dataGreif9, dataGreif10,dataGreif11,dataGreif12),axis=1)

#dataP0 = np.loadtxt('eosPol0table', unpack=True)
#dataPol0 = np.loadtxt('NS_data_eosPol4.txt', unpack=True)

ec = data[0,:]
M = data[1,:]
M0= data[2,:]
Mstat = data[3,:]
R = data[4,:]
Rratio = data[5,:]
Rstat = data[6,:]
freq =data[7,:]
kfreq = data[8,:]
J = data[9,:]
T = data[10,:]
W = data[11,:]
RratioS = data[12,:]
Mp=M+W-T

Mmax0=0;
for i in M[:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in data[4,:]:
##  data[4,i]=Mmax0;

Rmax0=0;
for i in R[:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in data[13,:]:
##  Rmax[i]=Mmax0;
Mmax=Mmax0;
Rmax=Rmax0;
Qmax=Mmax0/Rmax0;

    

#Wstat


Wstat=np.zeros(len(W))
Wstat[0]=W[0]
Wdum=W[0]
edum=ec[0]
for i in range(len(W)):
  if(i>=1):
    if(ec[i]==edum):
      Wstat[i]=Wdum
    else:
      edum=ec[i]
      Wdum=W[i]
      Wstat[i]=W[i]


#Mpstat


Mpstat=np.zeros(len(Mp))
Mpstat[0]=Mp[0]
Mpdum=Mp[0]
edum=ec[0]
for i in range(len(Mp)):
  if(i>=1):
    if(ec[i]==edum):
      Mpstat[i]=Mpdum
    else:
      edum=ec[i]
      Mpdum=Mp[i]
      Mpstat[i]=Mp[i]




#M0stat
size0=0
for i in freq[:]:
 size0=size0+1

M0_stat=np.zeros(size0)
M0_stat0=M0[0];
count0=0;

for i in freq[:]:
 if i<1:
  M0_stat0=M0[count0]
 M0_stat[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# EOS 
m = M * 1.9884e33   
r = R * 1.0e5       
mmax = Mmax * 1.9884e33   
rstat = Rstat * 1.0e5   
mstat = Mstat * 1.9884e33   
delta = freq * (2*np.pi) * ((rstat**3)/(G*mstat))**(0.5)  
m0 = M0 * 1.9884e33     
normJ = (c*J)/(G*m0**2)    
a = (c*J)/(G*m**2)         










#dame
NorK=np.zeros(len(delta))
ComK=np.zeros(len(delta))
countk=0
for i in range(len(delta)):
  if i>1:
    if (Rratio[i]==1) and (Rratio[i-1]!=1) and (kfreq[i-1]-freq[i-1]<100):
      NorK[i-1] = delta[i-1]#freq[i-1] * (2*np.pi) * ((r[i-1]**3)/(G*m[i-1]))**(0.5)
      ComK[i-1] = Mstat[i-1]/Rstat[i-1]
      countk=countk+1

Nor2Kepler=np.zeros(countk)
ComKepler=np.zeros(countk)
countk=0
for i in range(len(delta)):
  if NorK[i]!=0:
    Nor2Kepler[countk]=NorK[i]
    ComKepler[countk]=ComK[i]# * G* 1.9884e33/ (1.0e5*c**2) 
    countk=countk+1


SlopK4,SlopK3,SlopK2,SlopK,InterK=np.polyfit(ComKepler,Nor2Kepler,4)

print(np.max(np.abs(((SlopK4*(ComKepler)**4+SlopK3*(ComKepler)**3+SlopK2*(ComKepler)**2+SlopK*ComKepler+InterK)-Nor2Kepler)/Nor2Kepler)))

NormD2=(SlopK4*(Mstat/Rstat)**4+SlopK3*(Mstat/Rstat)**3+SlopK2*(Mstat/Rstat)**2+SlopK*Mstat/Rstat+InterK)

delta=delta/NormD2#**0.5




absError = (SlopK4*(ComKepler)**4+SlopK3*(ComKepler)**3+SlopK2*(ComKepler)**2+SlopK*ComKepler+InterK) - Nor2Kepler

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
RsquaredK = 1.0 - (numpy.var(absError) / numpy.var(Nor2Kepler))




ec_Pol0 = dataPol0[0,:]*10**15 
M_Pol0 = dataPol0[1,:]
M0_Pol0 = dataPol0[2,:]
Mstat_Pol0 = dataPol0[3,:]
R_Pol0 = dataPol0[4,:]
Rratio_Pol0 = dataPol0[5,:]
Rstat_Pol0 = dataPol0[6,:]
freq_Pol0 =dataPol0[7,:]
kfreq_Pol0 = dataPol0[8,:]
J_Pol0 = dataPol0[9,:]
T_Pol0 = dataPol0[10,:]
W_Pol0 = dataPol0[11,:]
RratioS_Pol0 = dataPol0[12,:]
Mmax0=0;
for i in dataPol0[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataPol0[4,:]:
##  dataPol0[4,i]=Mmax0;

Rmax0=0;
for i in dataPol0[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataPol0[13,:]:
##  Rmax_Pol0[i]=Mmax0;
Mmax_Pol0=Mmax0;
Rmax_Pol0=Rmax0;
Qmax_Pol0=Mmax0/Rmax0;

size0=0
for i in dataPol0[8,:]:
 size0=size0+1

M0_stat_Pol0=np.zeros(size0)
M0_stat0=dataPol0[2,0];
count0=0;

for i in dataPol0[8,:]:
 if i<1:
  M0_stat0=M0_Pol0[count0]
 M0_stat_Pol0[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_Pol0 = M_Pol0 * 1.9884e33   
r_Pol0 = R_Pol0 * 1.0e5       
mmax_Pol0 = Mmax_Pol0 * 1.9884e33   
rstat_Pol0 = Rstat_Pol0 * 1.0e5   
mstat_Pol0 = Mstat_Pol0 * 1.9884e33   
delta_Pol0 = freq_Pol0 * (2*np.pi) * ((rstat_Pol0**3)/(G*mstat_Pol0))**(0.5)/((SlopK4*(Mstat_Pol0/Rstat_Pol0)**4+SlopK3*(Mstat_Pol0/Rstat_Pol0)**3+SlopK2*(Mstat_Pol0/Rstat_Pol0)**2+SlopK*Mstat_Pol0/Rstat_Pol0+InterK))#**0.5  
m0_Pol0 = M0_Pol0 * 1.9884e33     
normJ_Pol0 = (c*J_Pol0)/(G*m0_Pol0**2)    
a_Pol0 = (c*J_Pol0)/(G*m_Pol0**2)         



'''
ec_Pol0J = dataPol0J[0,:]*10**15   
M_Pol0J = dataPol0J[1,:]           
M0_Pol0J = dataPol0J[2,:]          
Mstat_Pol0J = dataPol0J[3,:]       
#Mmax_Pol0J = dataPol0J[4,:]        
R_Pol0J = dataPol0J[5,:]           
Rratio_Pol0J = dataPol0J[6,:]      
Rstat_Pol0J = dataPol0J[7,:]       
freq_Pol0J = dataPol0J[8,:]        
kfreq_Pol0J = dataPol0J[9,:]       
J_Pol0J = dataPol0J[10,:]          
T_Pol0J = dataPol0J[11,:]          
W_Pol0J = dataPol0J[12,:]          
#Rmax_Pol0J = dataPol0J[13,:]       
#Qmax_Pol0J = dataPol0J[14,:]       

Mmax0=0;
for i in dataPol0J[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataPol0J[4,:]:
##  dataPol0J[4,i]=Mmax0;

Rmax0=0;
for i in dataPol0J[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataPol0J[13,:]:
##  Rmax_Pol0J[i]=Mmax0;
Mmax_Pol0J=Mmax0;
Rmax_Pol0J=Rmax0;
Qmax_Pol0J=Mmax0/Rmax0;

size0=0
for i in dataPol0J[8,:]:
 size0=size0+1

M0_stat_Pol0J=np.zeros(size0)
M0_stat0=dataPol0J[2,0];
count0=0;

for i in dataPol0J[8,:]:
 if i<1:
  M0_stat0=M0_Pol0J[count0]
 M0_stat_Pol0J[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# EOS Pol0J
m_Pol0J = M_Pol0J * 1.9884e33   
r_Pol0J = R_Pol0J * 1.0e5       
mmax_Pol0J = Mmax_Pol0J * 1.9884e33   
rstat_Pol0J = Rstat_Pol0J * 1.0e5   
mstat_Pol0J = Mstat_Pol0J * 1.9884e33   
delta_Pol0J = freq_Pol0J /(1223.47) 
m0_Pol0J = M0_Pol0J * 1.9884e33     
normJ_Pol0J = (c*J_Pol0J)/(G*m0_Pol0J**2)    
a_Pol0J = (c*J_Pol0J)/(G*m_Pol0J**2)         


ec_Pol0nr = dataPol0nr[0,:]*10**15   
M_Pol0nr = dataPol0nr[1,:]           
M0_Pol0nr = dataPol0nr[2,:]          
Mstat_Pol0nr = dataPol0nr[3,:]       
R_Pol0nr = dataPol0nr[4,:]           
Rratio_Pol0nr = dataPol0nr[5,:]      
Rstat_Pol0nr = dataPol0nr[6,:]       
freq_Pol0nr = dataPol0nr[7,:]        
kfreq_Pol0nr = dataPol0nr[8,:]       
J_Pol0nr = dataPol0nr[9,:]          
T_Pol0nr = dataPol0nr[10,:]          
W_Pol0nr = dataPol0nr[11,:]          
#Rmax_Pol0nr = dataPol0nr[13,:]       
#Qmax_Pol0nr = dataPol0nr[14,:]       
Vp_Pol0nr = dataPol0nr[12,:] 

Mmax0=0;
for i in dataPol0nr[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataPol0nr[4,:]:
##  dataPol0nr[4,i]=Mmax0;

Rmax0=0;
for i in dataPol0nr[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataPol0nr[13,:]:
##  Rmax_Pol0nr[i]=Mmax0;
Mmax_Pol0nr=Mmax0;
Rmax_Pol0nr=Rmax0;
Qmax_Pol0nr=Mmax0/Rmax0;

size0=0
for i in dataPol0nr[8,:]:
 size0=size0+1

M0_stat_Pol0nr=np.zeros(size0)
M0_stat0=dataPol0nr[2,0];
count0=0;

for i in dataPol0nr[8,:]:
 if i<1:
  M0_stat0=M0_Pol0nr[count0]
 M0_stat_Pol0nr[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# EOS Pol0nr
m_Pol0nr = M_Pol0nr * 1.9884e33   
r_Pol0nr = R_Pol0nr * 1.0e5       
mmax_Pol0nr = Mmax_Pol0nr * 1.9884e33   
rstat_Pol0nr = Rstat_Pol0nr * 1.0e5   
mstat_Pol0nr = Mstat_Pol0nr * 1.9884e33   
delta_Pol0nr = freq_Pol0nr * (2*np.pi) * ((rstat_Pol0nr**3)/(G*mstat_Pol0nr))**(0.5)/((SlopK4*(Mstat_Pol0nr/Rstat_Pol0nr)**4+SlopK3*(Mstat_Pol0nr/Rstat_Pol0nr)**3+SlopK2*(Mstat_Pol0nr/Rstat_Pol0nr)**2+SlopK*Mstat_Pol0nr/Rstat_Pol0nr+InterK))  
m0_Pol0nr = M0_Pol0nr * 1.9884e33     
normJ_Pol0nr = (c*J_Pol0nr)/(G*m0_Pol0nr**2)    
a_Pol0nr = (c*J_Pol0nr)/(G*m_Pol0nr**2)         




ec_Pol0K = dataPol0K[0,:]*10**15   
M_Pol0K = dataPol0K[1,:]           
M0_Pol0K = dataPol0K[2,:]          
Mstat_Pol0K = dataPol0K[3,:]       
R_Pol0K = dataPol0K[4,:]           
Rratio_Pol0K = dataPol0K[5,:]      
Rstat_Pol0K = dataPol0K[6,:]       
freq_Pol0K = dataPol0K[7,:]        
kfreq_Pol0K = dataPol0K[8,:]       
J_Pol0K = dataPol0K[9,:]          
T_Pol0K = dataPol0K[10,:]          
W_Pol0K = dataPol0K[11,:]          
#Rmax_Pol0K = dataPol0K[13,:]       
#Qmax_Pol0K = dataPol0K[14,:]       
Vp_Pol0K = dataPol0K[12,:] 

Mmax0=0;
for i in dataPol0K[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataPol0K[4,:]:
##  dataPol0K[4,i]=Mmax0;

Rmax0=0;
for i in dataPol0K[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataPol0K[13,:]:
##  Rmax_Pol0K[i]=Mmax0;
Mmax_Pol0K=Mmax0;
Rmax_Pol0K=Rmax0;
Qmax_Pol0K=Mmax0/Rmax0;

size0=0
for i in dataPol0K[8,:]:
 size0=size0+1

M0_stat_Pol0K=np.zeros(size0)
M0_stat0=dataPol0K[2,0];
count0=0;

for i in dataPol0K[8,:]:
 if i<1:
  M0_stat0=M0_Pol0K[count0]
 M0_stat_Pol0K[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# EOS Pol0K
m_Pol0K = M_Pol0K * 1.9884e33   
r_Pol0K = R_Pol0K * 1.0e5       
mmax_Pol0K = Mmax_Pol0K * 1.9884e33   
rstat_Pol0K = Rstat_Pol0K * 1.0e5   
mstat_Pol0K = Mstat_Pol0K * 1.9884e33   
delta_Pol0K = freq_Pol0K * (2*np.pi) * ((rstat_Pol0K**3)/(G*mstat_Pol0K))**(0.5)/((SlopK4*(Mstat_Pol0K/Rstat_Pol0K)**4+SlopK3*(Mstat_Pol0K/Rstat_Pol0K)**3+SlopK2*(Mstat_Pol0K/Rstat_Pol0K)**2+SlopK*Mstat_Pol0K/Rstat_Pol0K+InterK)) 
m0_Pol0K = M0_Pol0K * 1.9884e33     
normJ_Pol0K = (c*J_Pol0K)/(G*m0_Pol0K**2)    
a_Pol0K = (c*J_Pol0K)/(G*m_Pol0K**2)         





ec_Pol0ins = dataPol0ins[0,:]*10**15   
M_Pol0ins = dataPol0ins[1,:]           
M0_Pol0ins = dataPol0ins[2,:]          
Mstat_Pol0ins = dataPol0ins[3,:]       
R_Pol0ins = dataPol0ins[4,:]           
Rratio_Pol0ins = dataPol0ins[5,:]      
Rstat_Pol0ins = dataPol0ins[6,:]       
freq_Pol0ins = dataPol0ins[7,:]        
kfreq_Pol0ins = dataPol0ins[8,:]       
J_Pol0ins = dataPol0ins[9,:]          
T_Pol0ins = dataPol0ins[10,:]          
W_Pol0ins = dataPol0ins[11,:]          
#Rmax_Pol0ins = dataPol0ins[13,:]       
#Qmax_Pol0ins = dataPol0ins[14,:]       
Vp_Pol0ins = dataPol0ins[12,:] 

Mmax0=0;
for i in dataPol0ins[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataPol0ins[4,:]:
##  dataPol0ins[4,i]=Mmax0;

Rmax0=0;
for i in dataPol0ins[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataPol0ins[13,:]:
##  Rmax_Pol0ins[i]=Mmax0;
Mmax_Pol0ins=Mmax0;
Rmax_Pol0ins=Rmax0;
Qmax_Pol0ins=Mmax0/Rmax0;

size0=0
for i in dataPol0ins[8,:]:
 size0=size0+1

M0_stat_Pol0ins=np.zeros(size0)
M0_stat0=dataPol0ins[2,0];
count0=0;

for i in dataPol0ins[8,:]:
 if i<1:
  M0_stat0=M0_Pol0ins[count0]
 M0_stat_Pol0ins[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# EOS Pol0ins
m_Pol0ins = M_Pol0ins * 1.9884e33   
r_Pol0ins = R_Pol0ins * 1.0e5       
mmax_Pol0ins = Mmax_Pol0ins * 1.9884e33   
rstat_Pol0ins = Rstat_Pol0ins * 1.0e5   
mstat_Pol0ins = Mstat_Pol0ins * 1.9884e33   
delta_Pol0ins = freq_Pol0ins * (2*np.pi) * ((rstat_Pol0ins**3)/(G*mstat_Pol0ins))**(0.5)/((SlopK4*(Mstat_Pol0ins/Rstat_Pol0ins)**4+SlopK3*(Mstat_Pol0ins/Rstat_Pol0ins)**3+SlopK2*(Mstat_Pol0ins/Rstat_Pol0ins)**2+SlopK*Mstat_Pol0ins/Rstat_Pol0ins+InterK))  
m0_Pol0ins = M0_Pol0ins * 1.9884e33     
normJ_Pol0ins = (c*J_Pol0ins)/(G*m0_Pol0ins**2)    
a_Pol0ins = (c*J_Pol0ins)/(G*m_Pol0ins**2)         



ec_Pol0_rho1 = dataPol0rho1[0,:]*10**15 
M_Pol0_rho1 = dataPol0rho1[1,:]
M0_Pol0_rho1 = dataPol0rho1[2,:]
Mstat_Pol0_rho1 = dataPol0rho1[3,:]
R_Pol0_rho1 = dataPol0rho1[4,:]
Rratio_Pol0_rho1 = dataPol0rho1[5,:]
Rstat_Pol0_rho1 = dataPol0rho1[6,:]
freq_Pol0_rho1 =dataPol0rho1[7,:]
kfreq_Pol0_rho1 = dataPol0rho1[8,:]
J_Pol0_rho1 = dataPol0rho1[9,:]
T_Pol0_rho1 = dataPol0rho1[10,:]
W_Pol0_rho1 = dataPol0rho1[11,:]
RratioS_Pol0_rho1 = dataPol0rho1[12,:]


# Converting parameters to CGS units
# EOS Pol0_rho1
m_Pol0_rho1 = M_Pol0_rho1 * 1.9884e33   
r_Pol0_rho1 = R_Pol0_rho1 * 1.0e5       
rstat_Pol0_rho1 = Rstat_Pol0_rho1 * 1.0e5   
mstat_Pol0_rho1 = Mstat_Pol0_rho1 * 1.9884e33   
delta_Pol0_rho1 = freq_Pol0_rho1 * (2*np.pi) * ((rstat_Pol0_rho1**3)/(G*mstat_Pol0_rho1))**(0.5)/((SlopK4*(Mstat_Pol0_rho1/Rstat_Pol0_rho1)**4+SlopK3*(Mstat_Pol0_rho1/Rstat_Pol0_rho1)**3+SlopK2*(Mstat_Pol0_rho1/Rstat_Pol0_rho1)**2+SlopK*Mstat_Pol0_rho1/Rstat_Pol0_rho1+InterK))  
m0_Pol0_rho1 = M0_Pol0_rho1 * 1.9884e33     
normJ_Pol0_rho1 = (c*J_Pol0_rho1)/(G*m0_Pol0_rho1**2)    
a_Pol0_rho1 = (c*J_Pol0_rho1)/(G*m_Pol0_rho1**2)         

ec_Pol0_rho2 = dataPol0rho2[0,:]*10**15 
M_Pol0_rho2 = dataPol0rho2[1,:]
M0_Pol0_rho2 = dataPol0rho2[2,:]
Mstat_Pol0_rho2 = dataPol0rho2[3,:]
R_Pol0_rho2 = dataPol0rho2[4,:]
Rratio_Pol0_rho2 = dataPol0rho2[5,:]
Rstat_Pol0_rho2 = dataPol0rho2[6,:]
freq_Pol0_rho2 =dataPol0rho2[7,:]
kfreq_Pol0_rho2 = dataPol0rho2[8,:]
J_Pol0_rho2 = dataPol0rho2[9,:]
T_Pol0_rho2 = dataPol0rho2[10,:]
W_Pol0_rho2 = dataPol0rho2[11,:]
RratioS_Pol0_rho2 = dataPol0rho2[12,:]

# Converting parameters to CGS units
# EOS Pol0_rho2
m_Pol0_rho2 = M_Pol0_rho2 * 1.9884e33   
r_Pol0_rho2 = R_Pol0_rho2 * 1.0e5       
rstat_Pol0_rho2 = Rstat_Pol0_rho2 * 1.0e5   
mstat_Pol0_rho2 = Mstat_Pol0_rho2 * 1.9884e33   
delta_Pol0_rho2 = freq_Pol0_rho2 * (2*np.pi) * ((rstat_Pol0_rho2**3)/(G*mstat_Pol0_rho2))**(0.5)/((SlopK4*(Mstat_Pol0_rho2/Rstat_Pol0_rho2)**4+SlopK3*(Mstat_Pol0_rho2/Rstat_Pol0_rho2)**3+SlopK2*(Mstat_Pol0_rho2/Rstat_Pol0_rho2)**2+SlopK*Mstat_Pol0_rho2/Rstat_Pol0_rho2+InterK))  
m0_Pol0_rho2 = M0_Pol0_rho2 * 1.9884e33     
normJ_Pol0_rho2 = (c*J_Pol0_rho2)/(G*m0_Pol0_rho2**2)    
a_Pol0_rho2 = (c*J_Pol0_rho2)/(G*m_Pol0_rho2**2)         

'''
ec_Pol = dataPol[0,:]*10**15 
M_Pol = dataPol[1,:]
M0_Pol = dataPol[2,:]
Mstat_Pol = dataPol[3,:]
R_Pol = dataPol[4,:]
Rratio_Pol = dataPol[5,:]
Rstat_Pol = dataPol[6,:]
freq_Pol =dataPol[7,:]
kfreq_Pol = dataPol[8,:]
J_Pol = dataPol[9,:]
T_Pol = dataPol[10,:]
W_Pol = dataPol[11,:]
RratioS_Pol = dataPol[12,:]
Mmax0=0;
for i in dataPol[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataPol[4,:]:
##  dataPol[4,i]=Mmax0;

Rmax0=0;
for i in dataPol[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataPol[13,:]:
##  Rmax_Pol[i]=Mmax0;
Mmax_Pol=Mmax0;
Rmax_Pol=Rmax0;
Qmax_Pol=Mmax0/Rmax0;

size0=0
for i in dataPol[8,:]:
 size0=size0+1

M0_stat_Pol=np.zeros(size0)
M0_stat0=dataPol[2,0];
count0=0;

for i in dataPol[8,:]:
 if i<1:
  M0_stat0=M0_Pol[count0]
 M0_stat_Pol[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_Pol = M_Pol * 1.9884e33   
r_Pol = R_Pol * 1.0e5       
mmax_Pol = Mmax_Pol * 1.9884e33   
rstat_Pol = Rstat_Pol * 1.0e5   
mstat_Pol = Mstat_Pol * 1.9884e33   
delta_Pol = freq_Pol * (2*np.pi) * ((rstat_Pol**3)/(G*mstat_Pol))**(0.5)#**0.5  
m0_Pol = M0_Pol * 1.9884e33     
normJ_Pol = (c*J_Pol)/(G*m0_Pol**2)    
a_Pol = (c*J_Pol)/(G*m_Pol**2)         


NorKP=np.zeros(len(delta_Pol))
ComKP=np.zeros(len(delta_Pol))
countk=0
for i in range(len(delta_Pol)):
  if i>1:
    if (Rratio_Pol[i]==1) and (Rratio_Pol[i-1]!=1) and (kfreq_Pol[i-1]-freq_Pol[i-1]<100):
      NorKP[i-1] = delta_Pol[i-1]#**2
      ComKP[i-1] = Mstat_Pol[i-1]/Rstat_Pol[i-1]
      countk=countk+1

Nor2KeplerP=np.zeros(countk)
ComKeplerP=np.zeros(countk)
countk=0
for i in range(len(delta_Pol)):
  if NorKP[i]!=0:
    Nor2KeplerP[countk]=NorKP[i]
    ComKeplerP[countk]=ComKP[i]
    countk=countk+1


delta_Pol = freq_Pol * (2*np.pi) * ((rstat_Pol**3)/(G*mstat_Pol))**(0.5)/((SlopK4*(Mstat_Pol/Rstat_Pol)**4+SlopK3*(Mstat_Pol/Rstat_Pol)**3+SlopK2*(Mstat_Pol/Rstat_Pol)**2+SlopK*Mstat_Pol/Rstat_Pol+InterK))#**0.5  


ec_Gr = dataGr[0,:]*10**15 
M_Gr = dataGr[1,:]
M0_Gr = dataGr[2,:]
Mstat_Gr = dataGr[3,:]
R_Gr = dataGr[4,:]
Rratio_Gr = dataGr[5,:]
Rstat_Gr = dataGr[6,:]
freq_Gr =dataGr[7,:]
kfreq_Gr = dataGr[8,:]
J_Gr = dataGr[9,:]
T_Gr = dataGr[10,:]
W_Gr = dataGr[11,:]
RratioS_Gr = dataGr[12,:]


Mmax0=0;
for i in dataGr[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataGr[4,:]:
##  dataGr[4,i]=Mmax0;

Rmax0=0;
for i in dataGr[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataGr[13,:]:
##  Rmax_Gr[i]=Mmax0;
Mmax_Gr=Mmax0;
Rmax_Gr=Rmax0;
Qmax_Gr=Mmax0/Rmax0;

size0=0
for i in dataGr[8,:]:
 size0=size0+1

M0_stat_Gr=np.zeros(size0)
M0_stat0=dataGr[2,0];
count0=0;

for i in dataGr[8,:]:
 if i<1:
  M0_stat0=M0_Gr[count0]
 M0_stat_Gr[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_Gr = M_Gr * 1.9884e33   
r_Gr = R_Gr * 1.0e5       
mmax_Gr = Mmax_Gr * 1.9884e33   
rstat_Gr = Rstat_Gr * 1.0e5   
mstat_Gr = Mstat_Gr * 1.9884e33   
delta_Gr = freq_Gr * (2*np.pi) * ((rstat_Gr**3)/(G*mstat_Gr))**(0.5)#**0.5  
m0_Gr = M0_Gr * 1.9884e33     
normJ_Gr = (c*J_Gr)/(G*m0_Gr**2)    
a_Gr = (c*J_Gr)/(G*m_Gr**2)         



NorKG=np.zeros(len(delta_Gr))
ComKG=np.zeros(len(delta_Gr))
countk=0
for i in range(len(delta_Gr)):
  if i>1:
    if (Rratio_Gr[i]==1) and (Rratio_Gr[i-1]!=1) and (kfreq_Gr[i-1]-freq_Gr[i-1]<100):
      NorKG[i-1] = delta_Gr[i-1]#**2
      ComKG[i-1] = Mstat_Gr[i-1]/Rstat_Gr[i-1]
      countk=countk+1

Nor2KeplerG=np.zeros(countk)
ComKeplerG=np.zeros(countk)
countk=0
for i in range(len(delta_Gr)):
  if NorKG[i]!=0:
    Nor2KeplerG[countk]=NorKG[i]
    ComKeplerG[countk]=ComKG[i]
    countk=countk+1


delta_Gr = freq_Gr * (2*np.pi) * ((rstat_Gr**3)/(G*mstat_Gr))**(0.5)/((SlopK4*(Mstat_Gr/Rstat_Gr)**4+SlopK3*(Mstat_Gr/Rstat_Gr)**3+SlopK2*(Mstat_Gr/Rstat_Gr)**2+SlopK*Mstat_Gr/Rstat_Gr+InterK))#**0.5  





ec_QHCD = dataQHCD[0,:]*10**15   
M_QHCD = dataQHCD[1,:]           
M0_QHCD = dataQHCD[2,:]          
Mstat_QHCD = dataQHCD[3,:]       
#Mmax_QHCD = dataQHCD[4,:]        
R_QHCD = dataQHCD[5,:]           
Rratio_QHCD = dataQHCD[6,:]      
Rstat_QHCD = dataQHCD[7,:]       
freq_QHCD = dataQHCD[8,:]        
kfreq_QHCD = dataQHCD[9,:]       
J_QHCD = dataQHCD[10,:]          
T_QHCD = dataQHCD[11,:]          
W_QHCD = dataQHCD[12,:]          
#Rmax_QHCD = dataQHCD[13,:]       
#Qmax_QHCD = dataQHCD[14,:]       

Mmax0=0;
for i in dataQHCD[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataQHCD[4,:]:
##  dataQHCD[4,i]=Mmax0;

Rmax0=0;
for i in dataQHCD[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataQHCD[13,:]:
##  Rmax_QHCD[i]=Mmax0;
Mmax_QHCD=Mmax0;
Rmax_QHCD=Rmax0;
Qmax_QHCD=Mmax0/Rmax0;

size0=0
for i in dataQHCD[8,:]:
 size0=size0+1

M0_stat_QHCD=np.zeros(size0)
M0_stat0=dataQHCD[2,0];
count0=0;

for i in dataQHCD[8,:]:
 if i<1:
  M0_stat0=M0_QHCD[count0]
 M0_stat_QHCD[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_QHCD = M_QHCD * 1.9884e33   
r_QHCD = R_QHCD * 1.0e5       
mmax_QHCD = Mmax_QHCD * 1.9884e33   
rstat_QHCD = Rstat_QHCD * 1.0e5   
mstat_QHCD = Mstat_QHCD * 1.9884e33   
delta_QHCD = freq_QHCD * (2*np.pi) * ((rstat_QHCD**3)/(G*mstat_QHCD))**(0.5)/((SlopK4*(Mstat_QHCD/Rstat_QHCD)**4+SlopK3*(Mstat_QHCD/Rstat_QHCD)**3+SlopK2*(Mstat_QHCD/Rstat_QHCD)**2+SlopK*Mstat_QHCD/Rstat_QHCD+InterK))#**0.5  
m0_QHCD = M0_QHCD * 1.9884e33     
normJ_QHCD = (c*J_QHCD)/(G*m0_QHCD**2)    
a_QHCD = (c*J_QHCD)/(G*m_QHCD**2)         





# Extracting columns from the data files of each EOS
# Greif EOS
ec_ABPR2 = dataABPR2[0,:]*10**15   
M_ABPR2 = dataABPR2[1,:]           
M0_ABPR2 = dataABPR2[2,:]          
Mstat_ABPR2 = dataABPR2[3,:]       
#Mmax_ABPR2 = dataABPR2[4,:]        
R_ABPR2 = dataABPR2[5,:]           
Rratio_ABPR2 = dataABPR2[6,:]      
Rstat_ABPR2 = dataABPR2[7,:]       
freq_ABPR2 = dataABPR2[8,:]        
kfreq_ABPR2 = dataABPR2[9,:]       
J_ABPR2 = dataABPR2[10,:]          
T_ABPR2 = dataABPR2[11,:]          
W_ABPR2 = dataABPR2[12,:]          
#Rmax_ABPR2 = dataABPR2[13,:]       
#Qmax_ABPR2 = dataABPR2[14,:]       

Mmax0=0;
for i in dataABPR2[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataABPR2[4,:]:
##  dataABPR2[4,i]=Mmax0;

Rmax0=0;
for i in dataABPR2[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataABPR2[13,:]:
##  Rmax_ABPR2[i]=Mmax0;
Mmax_ABPR2=Mmax0;
Rmax_ABPR2=Rmax0;
Qmax_ABPR2=Mmax0/Rmax0;

size0=0
for i in dataABPR2[8,:]:
 size0=size0+1

M0_stat_ABPR2=np.zeros(size0)
M0_stat0=dataABPR2[2,0];
count0=0;

for i in dataABPR2[8,:]:
 if i<1:
  M0_stat0=M0_ABPR2[count0]
 M0_stat_ABPR2[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_ABPR2 = M_ABPR2 * 1.9884e33   
r_ABPR2 = R_ABPR2 * 1.0e5       
mmax_ABPR2 = Mmax_ABPR2 * 1.9884e33   
rstat_ABPR2 = Rstat_ABPR2 * 1.0e5   
mstat_ABPR2 = Mstat_ABPR2 * 1.9884e33   
delta_ABPR2 = freq_ABPR2 * (2*np.pi) * ((rstat_ABPR2**3)/(G*mstat_ABPR2))**(0.5)/((SlopK4*(Mstat_ABPR2/Rstat_ABPR2)**4+SlopK3*(Mstat_ABPR2/Rstat_ABPR2)**3+SlopK2*(Mstat_ABPR2/Rstat_ABPR2)**2+SlopK*Mstat_ABPR2/Rstat_ABPR2+InterK))  
m0_ABPR2 = M0_ABPR2 * 1.9884e33     
normJ_ABPR2 = (c*J_ABPR2)/(G*m0_ABPR2**2)    
a_ABPR2 = (c*J_ABPR2)/(G*m_ABPR2**2)         



# Extracting columns from the data files of each EOS
# Greif EOS
ec_HLPS2 = dataHLPS2[0,:]*10**15   
M_HLPS2 = dataHLPS2[1,:]           
M0_HLPS2 = dataHLPS2[2,:]          
Mstat_HLPS2 = dataHLPS2[3,:]       
#Mmax_HLPS2 = dataHLPS2[4,:]        
R_HLPS2 = dataHLPS2[5,:]           
Rratio_HLPS2 = dataHLPS2[6,:]      
Rstat_HLPS2 = dataHLPS2[7,:]       
freq_HLPS2 = dataHLPS2[8,:]        
kfreq_HLPS2 = dataHLPS2[9,:]       
J_HLPS2 = dataHLPS2[10,:]          
T_HLPS2 = dataHLPS2[11,:]          
W_HLPS2 = dataHLPS2[12,:]          
#Rmax_HLPS2 = dataHLPS2[13,:]       
#Qmax_HLPS2 = dataHLPS2[14,:]       

Mmax0=0;
for i in dataHLPS2[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataHLPS2[4,:]:
##  dataHLPS2[4,i]=Mmax0;

Rmax0=0;
for i in dataHLPS2[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataHLPS2[13,:]:
##  Rmax_HLPS2[i]=Mmax0;
Mmax_HLPS2=Mmax0;
Rmax_HLPS2=Rmax0;
Qmax_HLPS2=Mmax0/Rmax0;

size0=0
for i in dataHLPS2[8,:]:
 size0=size0+1

M0_stat_HLPS2=np.zeros(size0)
M0_stat0=dataHLPS2[2,0];
count0=0;

for i in dataHLPS2[8,:]:
 if i<1:
  M0_stat0=M0_HLPS2[count0]
 M0_stat_HLPS2[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_HLPS2 = M_HLPS2 * 1.9884e33   
r_HLPS2 = R_HLPS2 * 1.0e5       
mmax_HLPS2 = Mmax_HLPS2 * 1.9884e33   
rstat_HLPS2 = Rstat_HLPS2 * 1.0e5   
mstat_HLPS2 = Mstat_HLPS2 * 1.9884e33   
delta_HLPS2 = freq_HLPS2 * (2*np.pi) * ((rstat_HLPS2**3)/(G*mstat_HLPS2))**(0.5)/((SlopK4*(Mstat_HLPS2/Rstat_HLPS2)**4+SlopK3*(Mstat_HLPS2/Rstat_HLPS2)**3+SlopK2*(Mstat_HLPS2/Rstat_HLPS2)**2+SlopK*Mstat_HLPS2/Rstat_HLPS2+InterK))  
m0_HLPS2 = M0_HLPS2 * 1.9884e33     
normJ_HLPS2 = (c*J_HLPS2)/(G*m0_HLPS2**2)    
a_HLPS2 = (c*J_HLPS2)/(G*m_HLPS2**2)         


# Extracting columns from the data files of each EOS
# Greif EOS
'''
ec_HLPS3 = dataHLPS3[0,:]*10**15   
M_HLPS3 = dataHLPS3[1,:]           
M0_HLPS3 = dataHLPS3[2,:]          
Mstat_HLPS3 = dataHLPS3[3,:]       
#Mmax_HLPS3 = dataHLPS3[4,:]        
R_HLPS3 = dataHLPS3[5,:]           
Rratio_HLPS3 = dataHLPS3[6,:]      
Rstat_HLPS3 = dataHLPS3[7,:]       
freq_HLPS3 = dataHLPS3[8,:]        
kfreq_HLPS3 = dataHLPS3[9,:]       
J_HLPS3 = dataHLPS3[10,:]          
T_HLPS3 = dataHLPS3[11,:]          
W_HLPS3 = dataHLPS3[12,:]          
#Rmax_HLPS3 = dataHLPS3[13,:]       
#Qmax_HLPS3 = dataHLPS3[14,:]       

Mmax0=0;
for i in dataHLPS3[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataHLPS3[4,:]:
##  dataHLPS3[4,i]=Mmax0;

Rmax0=0;
for i in dataHLPS3[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataHLPS3[13,:]:
##  Rmax_HLPS3[i]=Mmax0;
Mmax_HLPS3=Mmax0;
Rmax_HLPS3=Rmax0;
Qmax_HLPS3=Mmax0/Rmax0;

size0=0
for i in dataHLPS3[8,:]:
 size0=size0+1

M0_stat_HLPS3=np.zeros(size0)
M0_stat0=dataHLPS3[2,0];
count0=0;

for i in dataHLPS3[8,:]:
 if i<1:
  M0_stat0=M0_HLPS3[count0]
 M0_stat_HLPS3[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_HLPS3 = M_HLPS3 * 1.9884e33   
r_HLPS3 = R_HLPS3 * 1.0e5       
mmax_HLPS3 = Mmax_HLPS3 * 1.9884e33   
rstat_HLPS3 = Rstat_HLPS3 * 1.0e5   
mstat_HLPS3 = Mstat_HLPS3 * 1.9884e33   
delta_HLPS3 = freq_HLPS3 * (2*np.pi) * ((rstat_HLPS3**3)/(G*mstat_HLPS3))**(0.5)/((SlopK4*(Mstat_HLPS3/Rstat_HLPS3)**4+SlopK3*(Mstat_HLPS3/Rstat_HLPS3)**3+SlopK2*(Mstat_HLPS3/Rstat_HLPS3)**2+SlopK*Mstat_HLPS3/Rstat_HLPS3+InterK))   
m0_HLPS3 = M0_HLPS3 * 1.9884e33     
normJ_HLPS3 = (c*J_HLPS3)/(G*m0_HLPS3**2)    
a_HLPS3 = (c*J_HLPS3)/(G*m_HLPS3**2)         
'''


# Extracting columns from the data files of each EOS
# Greif EOS
ec_L = dataL[0,:]*10**15   
M_L = dataL[1,:]           
M0_L = dataL[2,:]          
Mstat_L = dataL[3,:]       
#Mmax_L = dataL[4,:]        
R_L = dataL[5,:]           
Rratio_L = dataL[6,:]      
Rstat_L = dataL[7,:]       
freq_L = dataL[8,:]        
kfreq_L = dataL[9,:]       
J_L = dataL[10,:]          
T_L = dataL[11,:]          
W_L = dataL[12,:]          
#Rmax_L = dataL[13,:]       
#Qmax_L = dataL[14,:]       

Mmax0=0;
for i in dataL[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataL[4,:]:
##  dataL[4,i]=Mmax0;

Rmax0=0;
for i in dataL[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataL[13,:]:
##  Rmax_L[i]=Mmax0;
Mmax_L=Mmax0;
Rmax_L=Rmax0;
Qmax_L=Mmax0/Rmax0;

size0=0
for i in dataL[8,:]:
 size0=size0+1

M0_stat_L=np.zeros(size0)
M0_stat0=dataL[2,0];
count0=0;

for i in dataL[8,:]:
 if i<1:
  M0_stat0=M0_L[count0]
 M0_stat_L[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_L = M_L * 1.9884e33   
r_L = R_L * 1.0e5       
mmax_L = Mmax_L * 1.9884e33   
rstat_L = Rstat_L * 1.0e5   
mstat_L = Mstat_L * 1.9884e33   
delta_L = freq_L * (2*np.pi) * ((rstat_L**3)/(G*mstat_L))**(0.5)/((SlopK4*(Mstat_L/Rstat_L)**4+SlopK3*(Mstat_L/Rstat_L)**3+SlopK2*(Mstat_L/Rstat_L)**2+SlopK*Mstat_L/Rstat_L+InterK))   
m0_L = M0_L * 1.9884e33     
normJ_L = (c*J_L)/(G*m0_L**2)    
a_L = (c*J_L)/(G*m_L**2)         




# Extracting columns from the data files of each EOS
# Greif EOS
ec_HLPS1 = dataHLPS1[0,:]*10**15   
M_HLPS1 = dataHLPS1[1,:]           
M0_HLPS1 = dataHLPS1[2,:]          
Mstat_HLPS1 = dataHLPS1[3,:]       
#Mmax_HLPS1 = dataHLPS1[4,:]        
R_HLPS1 = dataHLPS1[5,:]           
Rratio_HLPS1 = dataHLPS1[6,:]      
Rstat_HLPS1 = dataHLPS1[7,:]       
freq_HLPS1 = dataHLPS1[8,:]        
kfreq_HLPS1 = dataHLPS1[9,:]       
J_HLPS1 = dataHLPS1[10,:]          
T_HLPS1 = dataHLPS1[11,:]          
W_HLPS1 = dataHLPS1[12,:]          
#Rmax_HLPS1 = dataHLPS1[13,:]       
#Qmax_HLPS1 = dataHLPS1[14,:]       

Mmax0=0;
for i in dataHLPS1[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataHLPS1[4,:]:
##  dataHLPS1[4,i]=Mmax0;

Rmax0=0;
for i in dataHLPS1[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataHLPS1[13,:]:
##  Rmax_HLPS1[i]=Mmax0;
Mmax_HLPS1=Mmax0;
Rmax_HLPS1=Rmax0;
Qmax_HLPS1=Mmax0/Rmax0;

size0=0
for i in dataHLPS1[8,:]:
 size0=size0+1

M0_stat_HLPS1=np.zeros(size0)
M0_stat0=dataHLPS1[2,0];
count0=0;

for i in dataHLPS1[8,:]:
 if i<1:
  M0_stat0=M0_HLPS1[count0]
 M0_stat_HLPS1[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_HLPS1 = M_HLPS1 * 1.9884e33   
r_HLPS1 = R_HLPS1 * 1.0e5       
mmax_HLPS1 = Mmax_HLPS1 * 1.9884e33   
rstat_HLPS1 = Rstat_HLPS1 * 1.0e5   
mstat_HLPS1 = Mstat_HLPS1 * 1.9884e33   
delta_HLPS1 = freq_HLPS1 * (2*np.pi) * ((rstat_HLPS1**3)/(G*mstat_HLPS1))**(0.5)/((SlopK4*(Mstat_HLPS1/Rstat_HLPS1)**4+SlopK3*(Mstat_HLPS1/Rstat_HLPS1)**3+SlopK2*(Mstat_HLPS1/Rstat_HLPS1)**2+SlopK*Mstat_HLPS1/Rstat_HLPS1+InterK))   
m0_HLPS1 = M0_HLPS1 * 1.9884e33     
normJ_HLPS1 = (c*J_HLPS1)/(G*m0_HLPS1**2)    
a_HLPS1 = (c*J_HLPS1)/(G*m_HLPS1**2)         



# Extracting columns from the data files of each EOS
# Greif EOS
ec_ABPR1 = dataABPR1[0,:]*10**15   
M_ABPR1 = dataABPR1[1,:]           
M0_ABPR1 = dataABPR1[2,:]          
Mstat_ABPR1 = dataABPR1[3,:]       
#Mmax_ABPR1 = dataABPR1[4,:]        
R_ABPR1 = dataABPR1[5,:]           
Rratio_ABPR1 = dataABPR1[6,:]      
Rstat_ABPR1 = dataABPR1[7,:]       
freq_ABPR1 = dataABPR1[8,:]        
kfreq_ABPR1 = dataABPR1[9,:]       
J_ABPR1 = dataABPR1[10,:]          
T_ABPR1 = dataABPR1[11,:]          
W_ABPR1 = dataABPR1[12,:]          
#Rmax_ABPR1 = dataABPR1[13,:]       
#Qmax_ABPR1 = dataABPR1[14,:]       

Mmax0=0;
for i in dataABPR1[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataABPR1[4,:]:
##  dataABPR1[4,i]=Mmax0;

Rmax0=0;
for i in dataABPR1[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataABPR1[13,:]:
##  Rmax_ABPR1[i]=Mmax0;
Mmax_ABPR1=Mmax0;
Rmax_ABPR1=Rmax0;
Qmax_ABPR1=Mmax0/Rmax0;

size0=0
for i in dataABPR1[8,:]:
 size0=size0+1

M0_stat_ABPR1=np.zeros(size0)
M0_stat0=dataABPR1[2,0];
count0=0;

for i in dataABPR1[8,:]:
 if i<1:
  M0_stat0=M0_ABPR1[count0]
 M0_stat_ABPR1[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_ABPR1 = M_ABPR1 * 1.9884e33   
r_ABPR1 = R_ABPR1 * 1.0e5       
mmax_ABPR1 = Mmax_ABPR1 * 1.9884e33   
rstat_ABPR1 = Rstat_ABPR1 * 1.0e5   
mstat_ABPR1 = Mstat_ABPR1 * 1.9884e33   
delta_ABPR1 = freq_ABPR1 * (2*np.pi) * ((rstat_ABPR1**3)/(G*mstat_ABPR1))**(0.5) /((SlopK4*(Mstat_ABPR1/Rstat_ABPR1)**4+SlopK3*(Mstat_ABPR1/Rstat_ABPR1)**3+SlopK2*(Mstat_ABPR1/Rstat_ABPR1)**2+SlopK*Mstat_ABPR1/Rstat_ABPR1+InterK))  
m0_ABPR1 = M0_ABPR1 * 1.9884e33     
normJ_ABPR1 = (c*J_ABPR1)/(G*m0_ABPR1**2)    
a_ABPR1 = (c*J_ABPR1)/(G*m_ABPR1**2)         


# Extracting columns from the data files of each EOS
# Greif EOS
ec_BBB1 = dataBBB1[0,:]*10**15   
M_BBB1 = dataBBB1[1,:]           
M0_BBB1 = dataBBB1[2,:]          
Mstat_BBB1 = dataBBB1[3,:]       
#Mmax_BBB1 = dataBBB1[4,:]        
R_BBB1 = dataBBB1[5,:]           
Rratio_BBB1 = dataBBB1[6,:]      
Rstat_BBB1 = dataBBB1[7,:]       
freq_BBB1 = dataBBB1[8,:]        
kfreq_BBB1 = dataBBB1[9,:]       
J_BBB1 = dataBBB1[10,:]          
T_BBB1 = dataBBB1[11,:]          
W_BBB1 = dataBBB1[12,:]          
#Rmax_BBB1 = dataBBB1[13,:]       
#Qmax_BBB1 = dataBBB1[14,:]       

Mmax0=0;
for i in dataBBB1[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataBBB1[4,:]:
##  dataBBB1[4,i]=Mmax0;

Rmax0=0;
for i in dataBBB1[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataBBB1[13,:]:
##  Rmax_BBB1[i]=Mmax0;
Mmax_BBB1=Mmax0;
Rmax_BBB1=Rmax0;
Qmax_BBB1=Mmax0/Rmax0;

size0=0
for i in dataBBB1[8,:]:
 size0=size0+1

M0_stat_BBB1=np.zeros(size0)
M0_stat0=dataBBB1[2,0];
count0=0;

for i in dataBBB1[8,:]:
 if i<1:
  M0_stat0=M0_BBB1[count0]
 M0_stat_BBB1[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_BBB1 = M_BBB1 * 1.9884e33   
r_BBB1 = R_BBB1 * 1.0e5       
mmax_BBB1 = Mmax_BBB1 * 1.9884e33   
rstat_BBB1 = Rstat_BBB1 * 1.0e5   
mstat_BBB1 = Mstat_BBB1 * 1.9884e33   
delta_BBB1 = freq_BBB1 * (2*np.pi) * ((rstat_BBB1**3)/(G*mstat_BBB1))**(0.5) /((SlopK4*(Mstat_BBB1/Rstat_BBB1)**4+SlopK3*(Mstat_BBB1/Rstat_BBB1)**3+SlopK2*(Mstat_BBB1/Rstat_BBB1)**2+SlopK*Mstat_BBB1/Rstat_BBB1+InterK))  
m0_BBB1 = M0_BBB1 * 1.9884e33     
normJ_BBB1 = (c*J_BBB1)/(G*m0_BBB1**2)    
a_BBB1 = (c*J_BBB1)/(G*m_BBB1**2)         


# Extracting columns from the data files of each EOS
# Greif EOS
ec_BBB2 = dataBBB2[0,:]*10**15   
M_BBB2 = dataBBB2[1,:]           
M0_BBB2 = dataBBB2[2,:]          
Mstat_BBB2 = dataBBB2[3,:]       
#Mmax_BBB2 = dataBBB2[4,:]        
R_BBB2 = dataBBB2[5,:]           
Rratio_BBB2 = dataBBB2[6,:]      
Rstat_BBB2 = dataBBB2[7,:]       
freq_BBB2 = dataBBB2[8,:]        
kfreq_BBB2 = dataBBB2[9,:]       
J_BBB2 = dataBBB2[10,:]          
T_BBB2 = dataBBB2[11,:]          
W_BBB2 = dataBBB2[12,:]          
#Rmax_BBB2 = dataBBB2[13,:]       
#Qmax_BBB2 = dataBBB2[14,:]       

Mmax0=0;
for i in dataBBB2[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataBBB2[4,:]:
##  dataBBB2[4,i]=Mmax0;

Rmax0=0;
for i in dataBBB2[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataBBB2[13,:]:
##  Rmax_BBB2[i]=Mmax0;
Mmax_BBB2=Mmax0;
Rmax_BBB2=Rmax0;
Qmax_BBB2=Mmax0/Rmax0;

size0=0
for i in dataBBB2[8,:]:
 size0=size0+1

M0_stat_BBB2=np.zeros(size0)
M0_stat0=dataBBB2[2,0];
count0=0;

for i in dataBBB2[8,:]:
 if i<1:
  M0_stat0=M0_BBB2[count0]
 M0_stat_BBB2[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_BBB2 = M_BBB2 * 1.9884e33   
r_BBB2 = R_BBB2 * 1.0e5       
mmax_BBB2 = Mmax_BBB2 * 1.9884e33   
rstat_BBB2 = Rstat_BBB2 * 1.0e5   
mstat_BBB2 = Mstat_BBB2 * 1.9884e33   
delta_BBB2 = freq_BBB2 * (2*np.pi) * ((rstat_BBB2**3)/(G*mstat_BBB2))**(0.5) /((SlopK4*(Mstat_BBB2/Rstat_BBB2)**4+SlopK3*(Mstat_BBB2/Rstat_BBB2)**3+SlopK2*(Mstat_BBB2/Rstat_BBB2)**2+SlopK*Mstat_BBB2/Rstat_BBB2+InterK))    
m0_BBB2 = M0_BBB2 * 1.9884e33     
normJ_BBB2 = (c*J_BBB2)/(G*m0_BBB2**2)    
a_BBB2 = (c*J_BBB2)/(G*m_BBB2**2)         



# Extracting columns from the data files of each EOS
# EOS APR
ec_APR = dataAPR[0,:]*10**15   
M_APR = dataAPR[1,:]           
M0_APR = dataAPR[2,:]          
Mstat_APR = dataAPR[3,:]       
#Mmax_APR = dataAPR[4,:]        
R_APR = dataAPR[5,:]           
Rratio_APR = dataAPR[6,:]      
Rstat_APR = dataAPR[7,:]       
freq_APR = dataAPR[8,:]        
kfreq_APR = dataAPR[9,:]       
J_APR = dataAPR[10,:]          
T_APR = dataAPR[11,:]          
W_APR = dataAPR[12,:]          
#Rmax_APR = dataAPR[13,:]       
#Qmax_APR = dataAPR[14,:]       

Mmax2=0;
for i in dataAPR[1,:]:
 if i > Mmax2 :
   Mmax2=i;

##for i in dataAPR[4,:]:
##  dataAPR[4,i]=Mmax2;

Rmax2=0;
for i in dataAPR[5,:]:
 if i > Rmax2 :
   Rmax2=i;

##for i in dataAPR[13,:]:
##  Rmax_APR[i]=Mmax2;
Mmax_APR=Mmax2;
Rmax_APR=Rmax2;
Qmax_APR=Mmax2/Rmax2;

size2=0
for i in dataAPR[8,:]:
 size2=size2+1

M0_stat_APR=np.zeros(size2)
M0_stat2=dataAPR[2,0];
count2=0;

for i in dataAPR[8,:]:
 if i<1:
  M0_stat2=M0_APR[count2]
 M0_stat_APR[count2]=M0_stat2
 count2=count2+1;

# Converting parameters to CGS units
# EOS APR
m_APR = M_APR * 1.9884e33   
r_APR = R_APR * 1.0e5       
mmax_APR = Mmax_APR * 1.9884e33   
rstat_APR = Rstat_APR * 1.0e5   
mstat_APR = Mstat_APR * 1.9884e33   
delta_APR = freq_APR * (2*np.pi) * ((rstat_APR**3)/(G*mstat_APR))**(0.5)/((SlopK4*(Mstat_APR/Rstat_APR)**4+SlopK3*(Mstat_APR/Rstat_APR)**3+SlopK2*(Mstat_APR/Rstat_APR)**2+SlopK*Mstat_APR/Rstat_APR+InterK))     
m0_APR = M0_APR * 1.9884e33     
normJ_APR = (c*J_APR)/(G*m0_APR**2)    
a_APR = (c*J_APR)/(G*m_APR**2)         




'''
# Extracting columns from the data files of each EOS
# Greif EOS
ec_Q1 = dataQ1[0,:]*10**15   
M_Q1 = dataQ1[1,:]           
M0_Q1 = dataQ1[2,:]          
Mstat_Q1 = dataQ1[3,:]       
#Mmax_Q1 = dataQ1[4,:]        
R_Q1 = dataQ1[5,:]           
Rratio_Q1 = dataQ1[6,:]      
Rstat_Q1 = dataQ1[7,:]       
freq_Q1 = dataQ1[8,:]        
kfreq_Q1 = dataQ1[9,:]       
J_Q1 = dataQ1[10,:]          
T_Q1 = dataQ1[11,:]          
W_Q1 = dataQ1[12,:]          
#Rmax_Q1 = dataQ1[13,:]       
#Qmax_Q1 = dataQ1[14,:]       

Mmax0=0;
for i in dataQ1[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataQ1[4,:]:
##  dataQ1[4,i]=Mmax0;

Rmax0=0;
for i in dataQ1[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataQ1[13,:]:
##  Rmax_Q1[i]=Mmax0;
Mmax_Q1=Mmax0;
Rmax_Q1=Rmax0;
Qmax_Q1=Mmax0/Rmax0;

size0=0
for i in dataQ1[8,:]:
 size0=size0+1

M0_stat_Q1=np.zeros(size0)
M0_stat0=dataQ1[2,0];
count0=0;

for i in dataQ1[8,:]:
 if i<1:
  M0_stat0=M0_Q1[count0]
 M0_stat_Q1[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_Q1 = M_Q1 * 1.9884e33   
r_Q1 = R_Q1 * 1.0e5       
mmax_Q1 = Mmax_Q1 * 1.9884e33   
rstat_Q1 = Rstat_Q1 * 1.0e5   
mstat_Q1 = Mstat_Q1 * 1.9884e33   
delta_Q1 = freq_Q1 * (2*np.pi) * ((rstat_Q1**3)/(G*mstat_Q1))**(0.5)/((SlopK4*(Mstat_Q1/Rstat_Q1)**4+SlopK3*(Mstat_Q1/Rstat_Q1)**3+SlopK2*(Mstat_Q1/Rstat_Q1)**2+SlopK*Mstat_Q1/Rstat_Q1+InterK))   
m0_Q1 = M0_Q1 * 1.9884e33     
normJ_Q1 = (c*J_Q1)/(G*m0_Q1**2)    
a_Q1 = (c*J_Q1)/(G*m_Q1**2)         


ec_Q1 = dataQ1[0,:]*10**15 
M_Q1 = dataQ1[1,:]
M0_Q1 = dataQ1[2,:]
Mstat_Q1 = dataQ1[3,:]
R_Q1 = dataQ1[4,:]
Rratio_Q1 = dataQ1[5,:]
Rstat_Q1 = dataQ1[6,:]
freq_Q1 =dataQ1[7,:]
kfreq_Q1 = dataQ1[8,:]
J_Q1 = dataQ1[9,:]
T_Q1 = dataQ1[10,:]
W_Q1 = dataQ1[11,:]
RratioS_Q1 = dataQ1[12,:]
Mmax0=0;
for i in dataQ1[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataQ1[4,:]:
##  dataQ1[4,i]=Mmax0;

Rmax0=0;
for i in dataQ1[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataQ1[13,:]:
##  Rmax_Q1[i]=Mmax0;
Mmax_Q1=Mmax0;
Rmax_Q1=Rmax0;
Qmax_Q1=Mmax0/Rmax0;

size0=0
for i in dataQ1[8,:]:
 size0=size0+1

M0_stat_Q1=np.zeros(size0)
M0_stat0=dataQ1[2,0];
count0=0;

for i in dataQ1[8,:]:
 if i<1:
  M0_stat0=M0_Q1[count0]
 M0_stat_Q1[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_Q1 = M_Q1 * 1.9884e33   
r_Q1 = R_Q1 * 1.0e5       
mmax_Q1 = Mmax_Q1 * 1.9884e33   
rstat_Q1 = Rstat_Q1 * 1.0e5   
mstat_Q1 = Mstat_Q1 * 1.9884e33   
delta_Q1 = freq_Q1 * (2*np.pi) * ((rstat_Q1**3)/(G*mstat_Q1))**(0.5)/((SlopK4*(Mstat_Q1/Rstat_Q1)**4+SlopK3*(Mstat_Q1/Rstat_Q1)**3+SlopK2*(Mstat_Q1/Rstat_Q1)**2+SlopK*Mstat_Q1/Rstat_Q1+InterK))#**0.5  
m0_Q1 = M0_Q1 * 1.9884e33     
normJ_Q1 = (c*J_Q1)/(G*m0_Q1**2)    
a_Q1 = (c*J_Q1)/(G*m_Q1**2)         
'''



# Greif EOS
ec_H0 = dataH0[0,:]*10**15
M_H0 = dataH0[1,:]
M0_H0= dataH0[2,:]
Mstat_H0 = dataH0[3,:]
R_H0 = dataH0[4,:]
Rratio_H0 = dataH0[5,:]
Rstat_H0 = dataH0[6,:]
freq_H0 =dataH0[7,:]
kfreq_H0 = dataH0[8,:]
J_H0 = dataH0[9,:]
T_H0 = dataH0[10,:]
W_H0 = dataH0[11,:]
RratioS_H0 = dataH0[12,:]
Mp=M+W-T
#Rmax_H0 = dataH0[13,:]       
#Qmax_H0 = dataH0[14,:]       

Mmax0=0;
for i in dataH0[1,:]:
 if i > Mmax0 :
   Mmax0=i;

##for i in dataH0[4,:]:
##  dataH0[4,i]=Mmax0;

Rmax0=0;
for i in dataH0[5,:]:
 if i > Rmax0 :
   Rmax0=i;

##for i in dataH0[13,:]:
##  Rmax_H0[i]=Mmax0;
Mmax_H0=Mmax0;
Rmax_H0=Rmax0;
Qmax_H0=Mmax0/Rmax0;

size0=0
for i in dataH0[8,:]:
 size0=size0+1

M0_stat_H0=np.zeros(size0)
M0_stat0=dataH0[2,0];
count0=0;

for i in dataH0[8,:]:
 if i<1:
  M0_stat0=M0_H0[count0]
 M0_stat_H0[count0]=M0_stat0
 count0=count0+1;

# Converting parameters to CGS units
# Greif EOS
m_H0 = M_H0 * 1.9884e33   
r_H0 = R_H0 * 1.0e5       
mmax_H0 = Mmax_H0 * 1.9884e33   
rstat_H0 = Rstat_H0 * 1.0e5   
mstat_H0 = Mstat_H0 * 1.9884e33   
delta_H0 = freq_H0 * (2*np.pi) * ((rstat_H0**3)/(G*mstat_H0))**(0.5)/((SlopK4*(Mstat_H0/Rstat_H0)**4+SlopK3*(Mstat_H0/Rstat_H0)**3+SlopK2*(Mstat_H0/Rstat_H0)**2+SlopK*Mstat_H0/Rstat_H0+InterK))   
m0_H0 = M0_H0 * 1.9884e33     
normJ_H0 = (c*J_H0)/(G*m0_H0**2)    
a_H0 = (c*J_H0)/(G*m_H0**2)         







deltanuc= numpy.concatenate((delta_APR,delta_BBB1,delta_BBB2,delta_L,delta_HLPS1,delta_HLPS2,delta_ABPR1,delta_ABPR2,delta_QHCD,delta_H0))

Rnuc=  numpy.concatenate((R_APR,R_BBB1,R_BBB2,R_L,R_HLPS1,R_HLPS2,R_ABPR1,R_ABPR2,R_QHCD,R_H0))

Mnuc= numpy.concatenate((M_APR,M_BBB1,M_BBB2,M_L,M_HLPS1,M_HLPS2,M_ABPR1,M_ABPR2,M_QHCD,M_H0))
   

Rstatnuc = numpy.concatenate((Rstat_APR,Rstat_BBB1,Rstat_BBB2,Rstat_L,Rstat_HLPS1,Rstat_HLPS2,Rstat_ABPR1,Rstat_ABPR2,Rstat_QHCD,Rstat_H0))

Mstatnuc=  numpy.concatenate((Mstat_APR,Mstat_BBB1,Mstat_BBB2,Mstat_L,Mstat_HLPS1,Mstat_HLPS2,Mstat_ABPR1,Mstat_ABPR2,Mstat_QHCD,Mstat_H0))






#Comp
#All EOS
xData = numpy.array(delta)
yData = numpy.array(Mstat/Rstat)
zData = numpy.array(M/R)

data = [xData, yData, zData]

initialParameters = np.ones(4) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersComp, pcov = scipy.optimize.curve_fit(funcComp, [xData, yData], zData, p0 = initialParameters)

#print('Compactnes fitted prameters', fittedParametersComp)

modelPredictionsComp = funcComp(data, *fittedParametersComp) 

absError = modelPredictionsComp - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
RsquaredComp = 1.0 - (numpy.var(absError) / numpy.var(zData))
#print('RMSE:', RMSE)
#print('R-squared Comp:', Rsquared)

'''
#Polytropes

xData = numpy.array(delta_Pol)
yData = numpy.array(Mstat_Pol/Rstat_Pol)
zData = numpy.array(M_Pol/R_Pol)

data = [xData, yData, zData]

initialParameters = np.ones(7) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersCompPol, pcov = scipy.optimize.curve_fit(funcComp, [xData, yData], zData, p0 = initialParameters)

print('Pol Compactnes fitted prameters', fittedParametersCompPol)

modelPredictionsCompPol = funcComp(data, *fittedParametersCompPol) 

absError = modelPredictionsCompPol - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
print('RMSE:', RMSE)
print('R-squared Comp:', Rsquared)



#Greif

xData = numpy.array(delta_Gr)
yData = numpy.array(Mstat_Gr/Rstat_Gr)
zData = numpy.array(M_Gr/R_Gr)

data = [xData, yData, zData]

initialParameters = np.ones(5) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersCompGr, pcov = scipy.optimize.curve_fit(funcComp, [xData, yData], zData, p0 = initialParameters)

print('Greif Compactnes fitted prameters', fittedParametersCompGr)

modelPredictionsCompGr = funcComp(data, *fittedParametersCompGr) 

absError = modelPredictionsCompGr - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
print('RMSE:', RMSE)
print('R-squared Comp:', Rsquared)

'''


















#Mass
#All
xData = numpy.array(delta)
yData = numpy.array(Mstat/Rstat)
zData = numpy.array((M)/Mstat)

data = [xData, yData, zData]

initialParameters = np.ones(6) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersBoth, pcov = scipy.optimize.curve_fit(func, [xData, yData], zData, p0 = initialParameters)

#print('Mass fitted prameters', fittedParametersBoth)

modelPredictionsBoth = func(data, *fittedParametersBoth) 

absError = modelPredictionsBoth - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
RsquaredMass = 1.0 - (numpy.var(absError) / numpy.var(zData))
#print('RMSE:', RMSE)
#print('R-squared Both:', Rsquared)


'''
#politropes

xData = numpy.array(delta_Pol)
yData = numpy.array(Mstat_Pol/Rstat_Pol)
zData = numpy.array(M_Pol/Mstat_Pol)

data = [xData, yData, zData]

initialParameters = np.ones(6) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersBothPol, pcov = scipy.optimize.curve_fit(func, [xData, yData], zData, p0 = initialParameters)

print('Polytropes Mass fitted prameters', fittedParametersBothPol)

modelPredictionsBothPol = func(data, *fittedParametersBothPol) 

absError = modelPredictionsBothPol - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
print('RMSE:', RMSE)
print('R-squared Both:', Rsquared)

#Greif

xData = numpy.array(delta_Gr)
yData = numpy.array(Mstat_Gr/Rstat_Gr)
zData = numpy.array(M_Gr/Mstat_Gr)

data = [xData, yData, zData]

initialParameters = np.ones(6) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersBothGr, pcov = scipy.optimize.curve_fit(func, [xData, yData], zData, p0 = initialParameters)

print('Greif Mass fitted prameters', fittedParametersBothGr)

modelPredictionsBothGr = func(data, *fittedParametersBothGr) 

absError = modelPredictionsBothGr - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
print('RMSE:', RMSE)
print('R-squared Both:', Rsquared)
'''












#Radius Both
#All
xData = numpy.array(delta)
yData = numpy.array(Mstat/Rstat)
zData = numpy.array(R/Rstat)

data = [xData, yData, zData]

initialParameters = np.ones(Rcoef) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersBothR, pcov = scipy.optimize.curve_fit(funcR, [xData, yData], zData, p0 = initialParameters)

#print('Radius fitted prameters', fittedParametersBothR)

modelPredictionsBothR = funcR(data, *fittedParametersBothR) 

absError = modelPredictionsBothR - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
RsquaredRadi = 1.0 - (numpy.var(absError) / numpy.var(zData))
#print('RMSE:', RMSE)
#print('R-squared Both:', Rsquared)

#andreas
print('a1 a2 a3 a4 a5 a6 a7 a8 R^2')
print('Equation 1',SlopK4,SlopK3,SlopK2,SlopK,InterK,'- -', RsquaredK)
print('Equation 6',fittedParametersComp[0],fittedParametersComp[1],fittedParametersComp[2],fittedParametersComp[3],'- - -', RsquaredComp)
print('Equation 2',fittedParametersBoth[0],fittedParametersBoth[1],fittedParametersBoth[2],fittedParametersBoth[3],fittedParametersBoth[4],fittedParametersBoth[5],'-' ,RsquaredMass)
print('Equation 4',fittedParametersBothR[0],fittedParametersBothR[1],fittedParametersBothR[2],fittedParametersBothR[3],fittedParametersBothR[4],fittedParametersBothR[5],fittedParametersBothR[6], RsquaredRadi)
'''
#Polytropes
xData = numpy.array(delta_Pol)
yData = numpy.array(Mstat_Pol/Rstat_Pol)
zData = numpy.array(R_Pol/Rstat_Pol)

data = [xData, yData, zData]

initialParameters = np.ones(Rcoef) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersBothRPol, pcov = scipy.optimize.curve_fit(funcR, [xData, yData], zData, p0 = initialParameters)

print('Polytropes Radius fitted prameters', fittedParametersBothRPol)

modelPredictionsBothRPol = funcR(data, *fittedParametersBothRPol) 

absError = modelPredictionsBothRPol - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
print('RMSE:', RMSE)
print('R-squared Both:', Rsquared)

#Greif
xData = numpy.array(delta_Gr)
yData = numpy.array(Mstat_Gr/Rstat_Gr)
zData = numpy.array(R_Gr/Rstat_Gr)

data = [xData, yData, zData]

initialParameters = np.ones(Rcoef) # these are the same as scipy default values in this example

# here a non-linear surface fit is made with scipy's curve_fit()
fittedParametersBothRGr, pcov = scipy.optimize.curve_fit(funcR, [xData, yData], zData, p0 = initialParameters)

print('Polytropes Radius fitted prameters', fittedParametersBothRGr)

modelPredictionsBothRGr = funcR(data, *fittedParametersBothRGr) 

absError = modelPredictionsBothRGr - zData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
print('RMSE:', RMSE)
print('R-squared Both:', Rsquared)
'''




if(0):
  x = delta
  y = (Mstat/Rstat)
  z = np.abs(DivSepM(delta,(Mstat/Rstat),M/Mstat,fittedParametersBoth))
  xmax=np.max(x)
  xmin=np.min(x)
  ymax=np.max(y)
  ymin=np.min(y)
  xbin=10
  ybin=20
    
  x_bins = np.linspace(xmin, xmax, xbin)
  y_bins = np.linspace(ymin, ymax,ybin)

  HM, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])



  z = np.abs(DivSepC(delta,(Mstat/Rstat),M/R,fittedParametersComp))

  HC, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
   
  z = np.abs(DivSepR(delta,(Mstat/Rstat),R/Rstat,fittedParametersBothR))

  HR, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
     
  print(x_bins)
  print(y_bins)
  print(HM)
  xint=[0,1,4,7,8]
  yint=[0,4,9,13,18]
  for i in range(len(yint)):
    for j in range(len(xint)):
     a=xint[j]
     b=yint[i]
     xbin=(x_bins[a+1]+x_bins[a])/2.
     ybin=(y_bins[b+1]+y_bins[b])/2.
     errM=HM[a,b]*100.
     errR=HR[a,b]*100.
     errC=HC[a,b]*100.
     print(f'{xbin:.4}','&',f'{ybin:.4}','&',f'{zfromfitM(xbin,ybin,fittedParametersBoth):.4}','&',f'{errM:.2}','&',f'{zfromfitR(xbin,ybin,fittedParametersBothR):.4}','&',f'{errR:.2}','&',f'{zfromfitC(xbin,ybin,fittedParametersComp):.4}','&',f'{errC:.2}','\\\\')







# Menu
print('Choose the plot to generate (Choose only integer numbers)')
print('0.  |2R*-Re-Rp|/2R* vs Omega_n vs C* - 2D Histogram (max)') 
print('1.  W/M vs E/M') 
print('2.  W/M vs C* - if you modified it you can also extract E/m vs C*')
print('3.  I/M^3 vs C*')
print('    Be sure that you have the correct equations in funcPol and zfromfitPol functions')
print('4.  Eb vs GC*/c^2') 
print('    Be sure that you have the correct equations in funcPol and zfromfitPol functions')
print('5.  Mmax(rot) vs Mmax(stat)') 
print('6.  Rmax(rot) vs Rmax(stat)')
print('7.  Re/R* vs Omega_n vs C*')
print('8.  Re/R* vs Omega_n vs C* - 2D Histogram (max) If you modified it you can also extract M/M and Ce/C*')
print('9.  Dev(Re/R*) vs Omega_n vs C*')
print('10. Dev(Re/R*) vs Omega_n vs C* - 2D Histogram (max)') 
print('11. Dev(Re/R*) vs Omega_n vs C* - from other EOSs')
print('12. Dev(Re/R*) vs Omega_n vs C* - 2D Histogram (max) - from other EOSs') 
print('13. Re/R* vs Omega_n - for EOS PP0')
print('14. Dev(M/M*) vs Omega_n vs C*')
print('15. Dev(M/M*) vs Omega_n vs C* - 2D Histogram (max)') 
print('16. Dev(M/M*) vs Omega_n vs C* - from other EOSs')
print('17. Dev(M/M*) vs Omega_n vs C* - 2D Histogram (max) - from other EOSs') 
print('18. M/M* vs Omega_n - for EOS PP0')
print('19. Ce vs Omega_n vs C*')
print('20. Dev(Ce) vs Omega_n vs C*')
print('21. Dev(Ce) vs Omega_n vs C* - 2D Histogram (max)') 
print('22. Dev(Ce) vs Omega_n vs C* - from other EOSs')
print('23. Dev(Ce) vs Omega_n vs C* - 2D Histogram (max) - from other EOSs') 
print('24. Ce vs Omega_n - for EOS PP0')
print('25. Dev(Ce) vs Omega_n vs C* - 2D Histogram (max) - from R/R* and M/M* best fit equations') 
print('26. Mass - Radius curves of other EOSs')
print('27. Finds Re/R* vs Ce and Omega_n2')
print('28. Finds M/M* vs Ce and Omega_n2 and then does inverse mapping on EOS PP0')
print('29. Application of our best fit equations that depend on C* and Omega_n')
print('30. Rratio vs Omega_n vs C*')
print('31. Dev(Rratio) vs Omega_n vs C* - 2D Histogram (max)')
print('    Be sure that you have the correct equations in funcPol and zfromfitPol functions')
print('32. Dev(RratioS) vs Omega_n vs C* - 2D Histogram (max)')
print('    Be sure that you have the correct equations in funcPol and zfromfitPol functions')
print('33. T/M vs Omega_n vs C*')
print('34. W/M vs Omega_n vs C*')
print('35. M0/M vs Omega_n vs C*')
print('36. E/M vs Omega_n vs C*')
print('37. Omega_K vs C*')
print('---------------------------------------------------------------')

# Choose the number of the desired plot 
method = eval(input())
#method = 55



if method == 0:
    '''
    d1=np.array([dataGreif0[1,0],dataGreif1[1,0],dataGreif2[1,0],dataGreif3[1,0],dataGreif4[1,0],dataGreif5[1,0],dataGreif6[1,0],dataGreif7[1,0], dataGreif8[1,0],dataGreif9[1,0],dataGreif10[1,0],
    dataGreif11[1,0] ,dataPol0[1,0],dataPol1[1,0],dataPol2[1,0],dataPol3[1,0],dataPol4[1,0],dataPol5[1,0],dataPol6[1,0],dataPol7[1,0],dataPol8[1,0],dataPol9[1,0],dataPol10[1,0],
    dataPol11[1,0],dataPol12[1,0],dataPol13[1,0],dataPol14[1,0],dataPol15[1,0],dataPol16[1,0],dataPol17[1,0],dataPol18[1,0]
    ])
    d2=np.array([dataGreif0[5,0],dataGreif1[5,0],dataGreif2[5,0],dataGreif3[5,0],dataGreif4[5,0],dataGreif5[5,0],dataGreif6[5,0],dataGreif7[5,0], dataGreif8[5,0],dataGreif9[5,0],dataGreif10[5,0],
    dataGreif11[5,0],dataPol0[5,0],dataPol1[5,0],dataPol2[5,0],dataPol3[5,0],dataPol4[5,0],dataPol5[5,0],dataPol6[5,0],dataPol7[5,0],dataPol8[5,0],dataPol9[5,0],dataPol10[5,0],
    dataPol11[5,0],dataPol12[5,0],dataPol13[5,0],dataPol14[5,0],dataPol15[5,0],dataPol16[5,0],dataPol17[5,0],dataPol18[5,0]])
    dum3=d1/d2
    dum1=[2.32196,2.14741,2.23218,1.81746,1.68266,1.88388,2.6834,1.57344,2.63287,2.31666,2.74131 ,2.09868,2.05954,2.74468,2.06436,2.66948,1.79879,2.26074,1.75407,2.72281,1.36081,2.90151,1.83953 ,2.13547,2.53598,2.40429,2.52533,1.31001,1.4324,1.55035,1.51305]#*10**15
    dum2=[11.3212,6.58378,6.91275, 8.38318,6.24665,7.90123,9.58624,5.68621,8.32294,4.42874,11.6974,10.5364,9.85551,12.7102,5.96563 ,12.3574,7.61734,11.702,8.29364,9.83551,4.03766,10.4537,3.84867,4.59831,11.215,10.1417,11.7277 ,3.68243,6.52862,5.97707,2.99321]#*10**35
    '''
    '''
    ax = fig.add_subplot(projection='3d')
    
    ax.scatter(delta, Mstat/Rstat,np.log10(np.abs(2*Rstat-R-RratioS*R)/(2*Rstat)), label = 'All EOS')
     
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$(M*/R*)$', fontsize='24')
    ax.set_zlabel(r'$|2R_*-Re-Rp|/2R_*$ ', fontsize='24')
    ax.view_init(azim=-106,elev=1)
    #ax.axes.set_zlim3d(bottom=0, top=0.01) 
    #ax.axes.set_xlim3d(left=0, right=0.6) 
    plt.legend()
    '''  
    
    dum1=np.zeros(len(R))
    dum2=np.zeros(len(R))
    dum3=np.zeros(len(R))
    countk=0
    for i in range(len(R)):
     if(delta[i]<0.9 ):
          dum1[i] = delta[i]
          dum2[i] = Mstat[i]/Rstat[i]
          dum3[i] = (np.abs(2*Rstat[i]-R[i]-RratioS[i]*R[i])/(2*Rstat[i]))
          countk=countk+1

    Dc0=np.zeros(countk)
    Cc0=np.zeros(countk)
    Rc0=np.zeros(countk)
    countk=0
    for i in range(len(R)):
      if dum2[i]!=0:
        Dc0[countk]=dum1[i]
        Cc0[countk]=dum2[i]
        Rc0[countk]=dum3[i]
        countk=countk+1

        
    x=Dc0#delta
    y=Cc0#Mstat/Rstat
    z=Rc0#(np.abs(2*Rstat-R-Rratio*R)/(2*Rstat))
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31






    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin


    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels  
   

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels   
   

    
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$|2R_*-Re-Rp|/2R_*$", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.savefig('ReRpC.png', format = 'png', transparent=False)
    plt.show()
    
    '''
    xData = numpy.array(delta)
    yData = numpy.array(Mstat/Rstat)
    zData = numpy.array(T/W)

    data = [xData, yData, zData]

    initialParameters = np.ones(56) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fpar, pcov = scipy.optimize.curve_fit(funcGen, [xData, yData], zData, p0 = initialParameters)


    modelPar = funcGen(data, *fpar) 

    absError = modelPar - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
    '''

    '''
   
    x = delta
    y = ec
    z = M*Rstat/(R*Mstat)
    

    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=10
    ybin=15
    
    x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    #y1=list(range(1,ybins+1))*(ymax-ymin)/20
    #x1=list(range(1,ybins+1))*(xmax-xmin)/20
    
    #z1=[[0]*xbin]*ybin


    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='min', bins=[x_bins, y_bins])

   

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$\epsilon_c$ [$10^{15}g/cm^3$]', fontsize='24')
    plt.title("$C_e/C_*$")
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    '''
    '''
    dum1=np.zeros(len(R))
    dum2=np.zeros(len(R))
    countk=0
    for i in range(len(R)):
     if(Mstat[i]/Rstat[i]>0.10 ):
          dum1[i] = M[i]/Mstat[i]
          dum2[i] = R[i]/Rstat[i]
          countk=countk+1

    NC0=np.zeros(countk)
    ec0=np.zeros(countk)
    countk=0
    for i in range(len(R)):
      if dum2[i]!=0:
        ec0[countk]=dum1[i]
        NC0[countk]=dum2[i]
        countk=countk+1

    plt.scatter(NC0,ec0, s=3, label = 'All EOS')



    a2,a1,a0=np.polyfit(NC0,ec0,2)
    

    print(np.max(np.abs((a2*NC0**2+a1*NC0+a0-ec0)/ec0)))
    
   

    xar=np.array(range(100, 140))/100.
    plt.plot(xar,a2*xar**2+a1*xar+a0, 'k', label = ' {:.2f}($R_e/R_*$)^2+{:.2f}($R_e/R_*$){:.2f}'.format(a2,a1,a0))  
    '''
    '''
    plt.scatter(NC0,ec0, s=3, label = 'All EOS')
   
   
    plt.xlabel(r'$R_e/R_*$ ', fontsize='24')
    plt.ylabel(r'$M/M_*$', fontsize='24')
    plt.legend()
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Things.png', format = 'png', transparent=False)
    plt.show()
    '''
elif method == 1:
    

    yData = numpy.array(Mstat/Rstat)
    xData = numpy.array(-(W)/(M))
    zData = numpy.array((M-M0+W-T)/(M))


    a4,a3,a2,a1,a0=np.polyfit(xData,zData,4)
    

    print(np.max(np.abs((a4*xData**4+a3*xData**3+a2*xData**2+a1*xData+a0-zData)/zData)))
    

    plt.scatter(xData,zData, s=7,color='blue', label = 'All EOS')
    xar=np.array(range(-736, -44))/1000.
    plt.plot(xar,a4*xar**4+a3*xar**3+a2*xar**2+a1*xar+a0, 'k', label = '{:.2f}$(M_*/R_*)^4$+{:.2f}$(M_*/R_*)^3$ +{:.2f}$(M_*/R_*)^2$ {:.2f}(M_*/R_*)+{:.2f}'.format(a4,a3,a2,a1,a0))  
    
    plt.xlabel(r'$W/M$', fontsize='24')
    plt.ylabel(r'$E/M$', fontsize='24')
    plt.legend(loc='best')
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Vir2d.png', format = 'png', transparent=False)
    plt.show()
elif method == 2:
    
    xData = numpy.array(Mstat/Rstat)
    yData = numpy.array(Mstat/Rstat)
    zData = numpy.array(-W/(M))
    #zData = numpy.array((M-M0+W-T)/M)


    a4,a3,a2,a1,a0=np.polyfit(xData,zData,4)
    

    print(np.max(np.abs((a4*xData**4+a3*xData**3+a2*xData**2+a1*xData+a0-zData)/zData)))
    

    plt.scatter(xData,zData, s=7,color='blue', label = 'All EOS')
    xar=np.array(range(2, 220))/1000.
    plt.plot(xar,a4*xar**4+a3*xar**3+a2*xar**2+a1*xar+a0, 'k', label = '{:.2f}$(M_*/R_*)^4$ + {:.2f}$(M_*/R_*)^3$  {:.2f}$(M_*/R_*)^2$  {:.2f}($M_*/R_*$){:.2f}'.format(a4,a3,a2,a1,a0))  
    
    plt.xlabel(r'$C_*$ ', fontsize='24')
    plt.ylabel(r'$W/M$', fontsize='24')
    #plt.ylabel(r'$E/M$', fontsize='24')
    plt.legend(loc='upper right')
    plt.legend(prop={"size":12})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Gr2d.png', format = 'png', transparent=False)
    #plt.savefig('InEn2d.png', format = 'png', transparent=False)
    plt.show()    
elif method == 3:

    
    Idum=np.zeros(len(J))
    Cdum=np.zeros(len(J))
    Fdum=np.zeros(len(J))

    countk=0
    for i in range(len(J)):
        if (freq[i]!=0 ):
          Fdum[i]=delta[i]
          Cdum[i]=Mstat[i]/Rstat[i]
          Idum[i]=J[i]*10**(-10)/(1.9884e33*2*np.pi*freq[i]*M[i]**3)
          countk=countk+1

    Fnor=np.zeros(countk)
    Cnor=np.zeros(countk)
    Inor=np.zeros(countk)
    countk=0
    for i in range(len(T)):
      if Idum[i]!=0:
        Fnor[countk]=Fdum[i]
        Cnor[countk]=Cdum[i]
        Inor[countk]=Idum[i]
        countk=countk+1
    
    xData = numpy.array(Fnor)
    yData = numpy.array(Cnor)
    zData = numpy.array(Inor)
    data = [xData, yData, zData]

    initialParameters = np.ones(10) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersPol, pcov = scipy.optimize.curve_fit(funcPol, [xData, yData], zData, p0 = initialParameters)
    
    print('fitted prameters', fittedParametersPol)

    modelPredictionsPol = funcPol(data, *fittedParametersPol) 

    absError = modelPredictionsPol - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
    print('RMSE:', RMSE)
    print('R-squared Comp:', Rsquared)


    fit=fittedParametersPol





    print(np.max(np.abs(((fit[0]/(fit[1]+yData)+fit[2]/(fit[3]+yData)**2+fit[4]/(fit[5]+yData)**3-zData)/zData))))







    plt.scatter(Mstat_Pol/Rstat_Pol,J_Pol*10**-10/(1.9884e33*2*np.pi*freq_Pol*M_Pol**3), s=3, label = 'PP EOS')
    plt.scatter(Mstat_Gr/Rstat_Gr,J_Gr*10**-10/(1.9884e33*2*np.pi*freq_Gr*M_Gr**3), s=3, label = '$c_s$ EOS')
    #xar=G*M * 1.9884e33/(R *c*c*1e5)
    xar=np.array(range(24, 220))/(1000.)
    xar1=np.array(range(60, 220))/(1000.)
    xar2=np.array(range(2, 220))/(1000.)
    #plt.plot(xar2,fittedParametersPol[0]/xar2+fittedParametersPol[1]/xar2**2+fittedParametersPol[2]/xar2**3+fittedParametersPol[3]/xar2**4+fittedParametersPol[4]/xar2**5+fittedParametersPol[5]/xar2**6+fittedParametersPol[6]/xar2**7+fittedParametersPol[7]/xar2**8+fittedParametersPol[8]/xar2**9+fittedParametersPol[9]/xar2**10, 'k', label ='Best Fit') #)) 
    #plt.plot(xar2,fit[0]/xar2+fit[1]/xar2**2+fit[2]/xar2**3+fit[3]/xar2**4+fit[4]/xar2**5+fit[5]/xar2**6+fit[6]/xar2**7+fit[7]/xar2**8+fit[8]/xar2**9+fit[9]/xar2**10, 'k', label ='Best Fit')
    
    plt.plot(xar2,fit[0]/(fit[1]+xar2)+fit[2]/(fit[3]+xar2)**2+fit[4]/(fit[5]+xar2)**3, 'k', label ='Best Fit')
    plt.plot(xar,0.21/(xar**2*(1-2*G*xar/(c*c))), 'r', label = ' Ravenhall & Pethick (1994)') 
    plt.plot(xar1,0.8134/xar1+0.2101/xar1**2+0.003175/xar1**3-0.0002717/xar1**4, 'b', label = ' Breu & Rezzolla (2016) [slow rot.]') 
    
    plt.xlabel(r'$C_*$ ', fontsize='24')
    plt.ylabel(r'$I/M^3 [km^2/M_\odot^2]$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Im.png', format = 'png', transparent=False)
    plt.show()


elif method == 4:


    xData = numpy.array(G*Mstat * 1.9884e33/(Rstat*c*c*1e5))
    yData = numpy.array(G*Mstat * 1.9884e33/(Rstat*c*c*1e5))
    zData = numpy.array( (M0-M)/M)
    data = [xData, yData, zData]

    initialParameters = np.ones(10) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersPol, pcov = scipy.optimize.curve_fit(funcPol, [xData, yData], zData, p0 = initialParameters)
    
    print('fitted prameters', fittedParametersPol)

    modelPredictionsPol = funcPol(data, *fittedParametersPol) 

    absError = modelPredictionsPol - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
    print('RMSE:', RMSE)
    print('R-squared Comp:', Rsquared)


    fit=fittedParametersPol


    print(np.max(np.abs(((fit[0]+fit[1]*yData+fit[2]*yData**2-zData)/zData))))

    xar=np.array(range(14, 320))/1000.
    plt.scatter(G*Mstat_Pol * 1.9884e33/(Rstat_Pol*c*c*1e5), (M0_Pol-M_Pol)/M_Pol, s=3, label = 'PP EOS')
    plt.scatter(G*Mstat_Gr * 1.9884e33/(Rstat_Gr*c*c*1e5), (M0_Gr-M_Gr)/M_Gr, s=3, label = '$c_s$ EOS')
    plt.plot(xar,fit[0]+fit[1]*xar+fit[2]*xar**2, 'k', label ='Best Fit')
    #xar=G*M * 1.9884e33/(R *c*c*1e5)

    plt.plot(xar,(0.6*xar)/(1-0.5*xar), 'r', label = ' Lattimer & Prakash (2001)') 
    plt.plot(xar,(0.619*xar)+(0.1359*xar**2), 'b', label = ' Breu & Rezzolla (2016)') 
    
    plt.xlabel(r'$GC_*/(c^2)$', fontsize='24')
    plt.ylabel(r'$E_b/M$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('BindEnergyComparison.png', format = 'png', transparent=False)
    plt.show()  
        
elif method == 5:

    xar=np.array([dataPol0[1,0],dataPol1[1,0],dataPol2[1,0],dataPol3[1,0],dataPol4[1,0],dataPol5[1,0],dataPol6[1,0],dataPol7[1,0],dataPol8[1,0],dataPol9[1,0],dataPol10[1,0],
    dataPol11[1,0],dataPol12[1,0],dataPol13[1,0],dataPol14[1,0],dataPol15[1,0],dataPol16[1,0],dataPol17[1,0],dataPol18[1,0],
    dataGreif0[1,0],dataGreif1[1,0],dataGreif2[1,0],dataGreif3[1,0],dataGreif4[1,0],dataGreif5[1,0],dataGreif6[1,0],dataGreif7[1,0], dataGreif8[1,0],dataGreif9[1,0],dataGreif10[1,0],
    dataGreif11[1,0],dataGreif12[1,0] ])

    yar=np.array([np.max(dataPol0[1,:]),   np.max(dataPol1[1,:]),np.max(dataPol2[1,:]),np.max(dataPol3[1,:]),np.max(dataPol4[1,:]),np.max(dataPol5[1,:]),np.max(dataPol6[1,:]),np.max(dataPol7[1,:]),np.max(dataPol8[1,:]),np.max(dataPol9[1,:]),np.max(dataPol10[1,:]),    np.max(dataPol11[1,:]),np.max(dataPol12[1,:]),np.max(dataPol13[1,:]),np.max(dataPol14[1,:]),np.max(dataPol15[1,:]),np.max(dataPol16[1,:]),np.max(dataPol17[1,:]),np.max(dataPol18[1,:]),    np.max(dataGreif0[1,:]),np.max(dataGreif1[1,:]),np.max(dataGreif2[1,:]),np.max(dataGreif3[1,:]),np.max(dataGreif4[1,:]),np.max(dataGreif5[1,:]),np.max(dataGreif6[1,:]),np.max(dataGreif7[1,:]), np.max(dataGreif8[1,:]),np.max(dataGreif9[1,:]),np.max(dataGreif10[1,:]),    np.max(dataGreif11[1,:]),    np.max(dataGreif12[1,:]) ])

    SlopM,InterM=np.polyfit(xar,yar,1)
    print(np.max(np.abs((SlopM*xar+InterM-yar)/yar)))
    plt.scatter(xar, yar, s=7,color='black', label = 'All EOS')
    plt.plot(xar,SlopM*xar+InterM, 'k', label = 'Best Fit')  
    plt.plot(xar,1.18*xar, 'r', label = ' Lasota et al. (1996)')  
    plt.plot(xar,1.203*xar, 'b', label = ' Breu & Rezzolla (2016)')  

    plt.xlabel(r'$M_{Max}(stat) [M_\odot]$', fontsize='24')
    plt.ylabel(r'$M_{Max}(rot) [M_\odot]$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Maxmass.png', format = 'png', transparent=False)
    plt.show()
elif method == 6:


    xar=np.array([dataPol0[4,0],dataPol1[4,0],dataPol2[4,0],dataPol3[4,0],dataPol4[4,0],dataPol5[4,0],dataPol6[4,0],dataPol7[4,0],dataPol8[4,0],dataPol9[4,0],dataPol10[4,0],
    dataPol11[4,0],dataPol12[4,0],dataPol13[4,0],dataPol14[4,0],dataPol15[4,0],dataPol16[4,0],dataPol17[4,0],dataPol18[4,0],
    dataGreif0[4,0],dataGreif1[4,0],dataGreif2[4,0],dataGreif3[4,0],dataGreif4[4,0],dataGreif5[4,0],dataGreif6[4,0],dataGreif7[4,0], dataGreif8[4,0],dataGreif9[4,0],dataGreif10[4,0],
    dataGreif11[4,0],dataGreif12[4,0]])

    yar=np.array([ dataPol0[4,np.argmax(dataPol0[1,:])],   dataPol1[4,np.argmax(dataPol1[1,:])],dataPol2[4,np.argmax(dataPol2[1,:])],dataPol3[4,np.argmax(dataPol3[1,:])],dataPol4[4,np.argmax(dataPol4[1,:])],dataPol5[4,np.argmax(dataPol5[1,:])],dataPol6[4,np.argmax(dataPol6[1,:])],dataPol7[4,np.argmax(dataPol7[1,:])],dataPol8[4,np.argmax(dataPol8[1,:])],dataPol9[4,np.argmax(dataPol9[1,:])],dataPol10[4,np.argmax(dataPol10[1,:])],    dataPol11[4,np.argmax(dataPol11[1,:])],dataPol12[4,np.argmax(dataPol12[1,:])],dataPol13[4,np.argmax(dataPol13[1,:])],dataPol14[4,np.argmax(dataPol14[1,:])],dataPol15[4,np.argmax(dataPol15[1,:])],dataPol16[4,np.argmax(dataPol16[1,:])],dataPol17[4,np.argmax(dataPol17[1,:])],dataPol18[4,np.argmax(dataPol18[1,:])],    dataGreif0[4,np.argmax(dataGreif0[1,:])],dataGreif1[4,np.argmax(dataGreif1[1,:])],dataGreif2[4,np.argmax(dataGreif2[1,:])],dataGreif3[4,np.argmax(dataGreif3[1,:])],dataGreif4[4,np.argmax(dataGreif4[1,:])],dataGreif5[4,np.argmax(dataGreif5[1,:])],dataGreif6[4,np.argmax(dataGreif6[1,:])],dataGreif7[4,np.argmax(dataGreif7[1,:])], dataGreif8[4,np.argmax(dataGreif8[1,:])],dataGreif9[4,np.argmax(dataGreif9[1,:])],dataGreif10[4,np.argmax(dataGreif10[1,:])],    dataGreif11[4,np.argmax(dataGreif11[1,:])],dataPol0[4,np.argmax(dataGreif11[1,:])] ])


    SlopR,InterR=np.polyfit(xar,yar,1)
    print(np.max(np.abs((SlopR*xar+InterR-yar)/yar)))
    plt.scatter(xar, yar, s=7,color='black', label = 'All EOS')
    plt.plot(xar,SlopR*xar+InterR, 'k', label = 'Best Fit')  
    plt.plot(xar,1.34*xar, 'r', label = ' Lasota et al. (1996)') 
     
    plt.xlabel(r'$R_{Max}(stat) [km]$', fontsize='24')
    plt.ylabel(r'$R_{Max}(rot) [km]$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Maxradius.png', format = 'png', transparent=False)
    plt.show()   



 
elif method == 7:
    ax = fig.add_subplot(projection='3d')
    #ax.plot_trisurf(xData, yData, zP7_fit)
    #ax.plot_trisurf(xData, yData, zG8_fit)
    #ax.plot_surface(xData, yData, zP7_fit,rstride=1, cstride=1, linewidth=0.25)
    #ax.plot_surface(xData, yData, zG8_fit,rstride=1, cstride=1, linewidth=0.25)
    ax.scatter(delta_Pol, (Mstat_Pol/Rstat_Pol),R_Pol/Rstat_Pol, label = 'PP EOS')
    ax.scatter(delta_Gr, (Mstat_Gr/Rstat_Gr),(R_Gr)/Rstat_Gr, label = '$c_s$ EOS')

    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'$R_e/R_*$', fontsize='24')
    #
    ax.view_init(azim=147,elev=17)
    plt.legend()
    plt.savefig('RadialChangeBoth.png', format = 'png', transparent=False)
    plt.show()
    

elif method == 8:
    #ax = fig.add_subplot(projection='3d')
    #ax.plot_trisurf(xData, yData, zP7_fit)
    #ax.plot_trisurf(xData, yData, zG8_fit)
    #ax.plot_surface(xData, yData, zP7_fit,rstride=1, cstride=1, linewidth=0.25)
    #ax.plot_surface(xData, yData, zG8_fit,rstride=1, cstride=1, linewidth=0.25)
    '''
    ax.scatter(delta_Pol, (Mstat_Pol/Rstat_Pol),M_Pol/Mstat_Pol, label = 'PP EOS')
    ax.scatter(delta_Gr, (Mstat_Gr/Rstat_Gr),(M_Gr)/Mstat_Gr, label = '$c_s$ EOS')
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'$M/M_*$', fontsize='24')
    ax.view_init(azim=147,elev=17)
    plt.legend()
    plt.savefig('MassChangeBoth.png', format = 'png', transparent=False)
    plt.show()
    '''
    
    x=delta
    y=Mstat/Rstat
    #z=M/Mstat
    z=R/Rstat
    #z=M*Rstat/(R*Mstat)
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31






    x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    #x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))
    '''    
    #y1=list(range(1,ybins+1))*(ymax-ymin)/20
    #x1=list(range(1,ybins+1))*(xmax-xmin)/20
    
    #z1=[[0]*xbin]*ybin
    
    H, xedges, yedges = np.histogram2d(x, y, bins = [x_bins, y_bins], weights = z)
    H_counts, xedges, yedges = np.histogram2d(x, y, bins = [x_bins, y_bins]) 
    H = H/H_counts
    plt.xlabel(r'$r_x$ [km]', fontsize='24')
    plt.ylabel(r'$r_z$ [km]', fontsize='24')
    plt.title("$\Omega$ = 0 Hz")
    plt.legend(loc='best')
    levels = (pow(10,14.2)/(np.max(e1)*c**2/(Kappa*G)), pow(10,14.66)/(np.max(e1)*c**2/(Kappa*G)), pow(10,14.899)/(np.max(e1)*c**2/(Kappa*G)))
    cset = plt.contour(H.T, levels, origin='lower',colors=['black','green','blue'],linewidths=(1.9, 1.6,  1.4),extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]])
    plt.clabel(cset, inline=1, fontsize=10, fmt='%2.1f')
    plt.imshow(H.T, origin='lower',  cmap='jet',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], vmin = 0, vmax =1)

    '''

    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='mean', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    levels = ( 1.050,1.150,1.250)
    #levels = ( 0.95,0.98,0.99)
    cset = plt.contour( H.T, levels, colors=['black','black','black'],linewidths=(2, 2,  2),extent = [xmin, xmax,ymin, ymax])
    plt.clabel(cset, inline=1, fontsize=15, fmt='%1.3f')
    

    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='mean', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    
    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    #plt.title("$M/M_*$", fontsize='24')
    plt.title("$R_e/R_*$", fontsize='24')
    #plt.title("$C_e/C_*$", fontsize='24')
    
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    
    
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)

    #plt.savefig('Mass2D.png', format = 'png', transparent=False)
    plt.savefig('Radius2D.png', format = 'png', transparent=False)
    #plt.savefig('Compact2D.png', format = 'png', transparent=False)
    
    plt.show()

    '''
    fig = plt.figure(figsize=(8, 6))
    hb = plt.hexbin(delta, y=(Mstat/Rstat), C= M/Mstat,gridsize=50, mincnt=0, cmap='jet')
    #,reduce_C_function=np.max
    #plt.zlabel(r'$M/M_*$', fontsize='24')
    fig.colorbar(hb,cmap='jet').set_label(label='M/$M_*$',size=15)
    
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.savefig('MassChange2D.png', format = 'png', transparent=False)
    plt.show()
    
    
    fig = plt.figure(figsize=(8, 6))
    hb = plt.hexbin(delta, y=(Mstat/Rstat), C= R/Rstat,gridsize=50, mincnt=0,cmap='jet')
    #plt.zlabel(r'$M/M_*$', fontsize='24')
    fig.colorbar(hb,cmap='jet').set_label(label='$R_e/R_*$',size=15)
    
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.savefig('RadiusChange2D.png', format = 'png', transparent=False)
    plt.show()    

    fig = plt.figure(figsize=(8, 6))
    hb = plt.hexbin(delta, y=(Mstat/Rstat), C= M/R,gridsize=50, mincnt=0,cmap='jet')
    #plt.zlabel(r'$M/M_*$', fontsize='24')
    fig.colorbar(hb,cmap='jet').set_label(label='M/$R_e$ ',size=15)
    
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.savefig('CompactChange2D.png', format = 'png', transparent=False)
    plt.show()
    '''

elif method == 9:
    ax = fig.add_subplot(projection='3d')
    #ax.scatter(xData, yData,divP0, label = 'EOS Pol0')

    ax.scatter(delta, (Mstat/Rstat),DivSepR(delta,(Mstat/Rstat),(R)/Rstat,fittedParametersBothR), label = 'All EOS')

    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'Dev($R/R_*$)', fontsize='24')

    ax.view_init(azim=85,elev=3)
    plt.legend()
    plt.savefig('Div_from_best_fitR_both.png', format = 'png', transparent=False)
    plt.show()    
elif method==10:    
    x = delta
    y = (Mstat/Rstat)
    z = np.abs(DivSepR(delta,(Mstat/Rstat),(R)/Rstat,fittedParametersBothR))
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31
    
    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin



    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
        
   

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.clim(0, 0.0814)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
        

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(R_e/R_*)$", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto', vmin = 0, vmax =0.0814)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.clim(0, 0.0814)   
    plt.savefig('DiverRAll.png', format = 'png', transparent=False)
    plt.show()
elif method == 11:

    ax = fig.gca(projection='3d')
    #ax.axes.set_xlim3d(left=0, right=1) 
    #ax.axes.set_ylim3d(bottom=0, top=0.25) 
    #ax.axes.set_zlim3d(bottom=0, top=10) 

    ax.scatter(delta_APR, Mstat_APR/Rstat_APR,DivSepR(delta_APR,(Mstat_APR/Rstat_APR),(R_APR)/Rstat_APR,fittedParametersBothR),s=3, label = 'APR EOS')
    ax.scatter(delta_BBB1, Mstat_BBB1/Rstat_BBB1,DivSepR(delta_BBB1,(Mstat_BBB1/Rstat_BBB1),(R_BBB1)/Rstat_BBB1,fittedParametersBothR),s=3, label = 'BBB1 EOS')
    ax.scatter(delta_BBB2, Mstat_BBB2/Rstat_BBB2,DivSepR(delta_BBB2,(Mstat_BBB2/Rstat_BBB2),(R_BBB2)/Rstat_BBB2,fittedParametersBothR),s=3, label = 'BBB2 EOS')
    ax.scatter(delta_ABPR1, Mstat_ABPR1/Rstat_ABPR1,DivSepR(delta_ABPR1,(Mstat_ABPR1/Rstat_ABPR1),(R_ABPR1)/Rstat_ABPR1,fittedParametersBothR),s=3, label = 'ABPR1 EOS')
    ax.scatter(delta_ABPR2, Mstat_ABPR2/Rstat_ABPR2,DivSepR(delta_ABPR2,(Mstat_ABPR2/Rstat_ABPR2),(R_ABPR2)/Rstat_ABPR2,fittedParametersBothR),s=3, label = 'ABPR2 EOS')
    ax.scatter(delta_HLPS1, Mstat_HLPS1/Rstat_HLPS1,DivSepR(delta_HLPS1,(Mstat_HLPS1/Rstat_HLPS1),(R_HLPS1)/Rstat_HLPS1,fittedParametersBothR),s=3, label = 'HLPS1 EOS')
    ax.scatter(delta_HLPS2, Mstat_HLPS2/Rstat_HLPS2,DivSepR(delta_HLPS2,(Mstat_HLPS2/Rstat_HLPS2),(R_HLPS2)/Rstat_HLPS2,fittedParametersBothR),s=3, label = 'HLPS2 EOS')
    #ax.scatter(delta_HLPS3, Mstat_HLPS3/Rstat_HLPS3,DivSepR(delta_HLPS3,(Mstat_HLPS3/Rstat_HLPS3),(R_HLPS3)/Rstat_HLPS3,fittedParametersBothR),s=3, label = 'HLPS3 EOS')
    ax.scatter(delta_L, Mstat_L/Rstat_L,DivSepR(delta_L,(Mstat_L/Rstat_L),(R_L)/Rstat_L,fittedParametersBothR),s=3, label = 'L EOS')
    ax.scatter(delta_H0, Mstat_H0/Rstat_H0,DivSepR(delta_H0,(Mstat_H0/Rstat_H0),(R_H0)/Rstat_H0,fittedParametersBothR),s=3, label = 'H0 EOS')
    ax.scatter(delta_QHCD, Mstat_QHCD/Rstat_QHCD,DivSepR(delta_QHCD,(Mstat_QHCD/Rstat_QHCD),(R_QHCD)/Rstat_QHCD,fittedParametersBothR),s=3, label = 'QHCD EOS')
    #ax.scatter(delta_Q1, Mstat_Q1/Rstat_Q1,DivSepR(delta_Q1,(Mstat_Q1/Rstat_Q1),(R_Q1)/Rstat_Q1,fittedParametersBothR),s=3, label = 'QS1 EOS')

    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$(C_*$ ', fontsize='24')
    ax.set_zlabel(r'Dev($R/R_*$)', fontsize='24')
    ax.view_init(azim=85,elev=3)
    plt.legend()
    plt.savefig('Div_from_best_fitR_both_nuc3.png', format = 'png', transparent=False)
    plt.show()
elif method==12: 


   
    x = deltanuc
    y = (Mstatnuc/Rstatnuc)
    z = np.abs(DivSepR(deltanuc,(Mstatnuc/Rstatnuc),Rnuc/Rstatnuc,fittedParametersBothR))
    
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31
    
    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin



    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])

    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
       

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.clim(0, 0.0814)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(R_e/R_*)$ [other NS EOSs]", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto', vmin = 0, vmax =0.0814)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.savefig('DiverRAllnuc.png', format = 'png', transparent=False)
    plt.show()
elif method==13:
    plt.scatter(delta_Pol0, R_Pol0/Rstat_Pol0, s=3, label = 'EOS Pol0')
    #plt.plot(delta_Pol0nr, M_Pol0nr/R_Pol0nr,'k-',linewidth=3)
    #plt.plot(delta_Pol0K, M_Pol0K/R_Pol0K,'k--',linewidth=3)
    #plt.plot(delta_Pol0ins, M_Pol0ins/R_Pol0ins,'k--')
    #plt.plot(delta_Pol0J, M_Pol0J/R_Pol0J,'k-')
    
            
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$R/R_*$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})  
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('R-frac_dim-ang-vel_alleosgreif.png', format = 'png', transparent=False)
    plt.show()
elif method == 14:
    ax = fig.add_subplot(projection='3d')
    #ax.scatter(xData, yData,divP0, label = 'EOS Pol0')

    ax.scatter(delta, (Mstat/Rstat),DivSepM(delta,(Mstat/Rstat),(M)/Mstat,fittedParametersBoth), label = 'All EOS')

    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$', fontsize='24')
    ax.set_zlabel(r'Dev($M/M_*$)', fontsize='24')

    ax.view_init(azim=130,elev=0)
    plt.legend()
    plt.savefig('Div_from_best_fit_both.png', format = 'png', transparent=False)
    plt.show()
elif method==15:    

    x = delta
    y = (Mstat/Rstat)
    z = np.abs(DivSepM(delta,(Mstat/Rstat),(M)/Mstat,fittedParametersBoth))
    '''
    fig = plt.figure(figsize=(8, 6))
    hb = plt.hexbin(x, y=y, C=z,gridsize=10, mincnt=0, cmap='jet',reduce_C_function=np.max)
    #
    #plt.zlabel(r'$M/M_*$', fontsize='24')
    fig.colorbar(hb,cmap='jet').set_label(label='Dev(M/$M_*)$',size=15)
    
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.savefig('DiverMAll2.png', format = 'png', transparent=False)
    plt.show()
    

    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=10
    ybin=20
    
    x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    #y1=list(range(1,ybins+1))*(ymax-ymin)/20
    #x1=list(range(1,ybins+1))*(xmax-xmin)/20
    '''
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31






    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin


    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])
    #levels = ( 0.010,0.020,0.030)
    #cset = plt.contour(H.T, levels, colors=['black','black','black'],linewidths=(2, 2,  2),extent = [xmin, xmax,ymin, ymax])
    #plt.clabel(cset, inline=1, fontsize=15, fmt='%1.3f')
    

    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.clim(0, 0.06)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(M/M_*)$", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    
    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto', vmin = 0, vmax =0.06)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)

    plt.savefig('DiverMAll.png', format = 'png', transparent=False)
    plt.show()
elif method == 16:
    
    ax = fig.add_subplot(projection='3d')
    ax.scatter(delta_APR, Mstat_APR/Rstat_APR,DivSepM(delta_APR,(Mstat_APR/Rstat_APR),(M_APR)/Mstat_APR,fittedParametersBoth),s=3, label = 'APR EOS')
    ax.scatter(delta_BBB1, Mstat_BBB1/Rstat_BBB1,DivSepM(delta_BBB1,(Mstat_BBB1/Rstat_BBB1),(M_BBB1)/Mstat_BBB1,fittedParametersBoth),s=3, label = 'BBB1 EOS')
    ax.scatter(delta_BBB2, Mstat_BBB2/Rstat_BBB2,DivSepM(delta_BBB2,(Mstat_BBB2/Rstat_BBB2),(M_BBB2)/Mstat_BBB2,fittedParametersBoth),s=3, label = 'BBB2 EOS')
    ax.scatter(delta_ABPR1, Mstat_ABPR1/Rstat_ABPR1,DivSepM(delta_ABPR1,(Mstat_ABPR1/Rstat_ABPR1),(M_ABPR1)/Mstat_ABPR1,fittedParametersBoth),s=3, label = 'ABPR1 EOS')
    ax.scatter(delta_ABPR2, Mstat_ABPR2/Rstat_ABPR2,DivSepM(delta_ABPR2,(Mstat_ABPR2/Rstat_ABPR2),(M_ABPR2)/Mstat_ABPR2,fittedParametersBoth),s=3, label = 'ABPR2 EOS')
    ax.scatter(delta_HLPS1, Mstat_HLPS1/Rstat_HLPS1,DivSepM(delta_HLPS1,(Mstat_HLPS1/Rstat_HLPS1),(M_HLPS1)/Mstat_HLPS1,fittedParametersBoth),s=3, label = 'HLPS1 EOS')
    ax.scatter(delta_HLPS2, Mstat_HLPS2/Rstat_HLPS2,DivSepM(delta_HLPS2,(Mstat_HLPS2/Rstat_HLPS2),(M_HLPS2)/Mstat_HLPS2,fittedParametersBoth),s=3, label = 'HLPS2 EOS')
    #ax.scatter(delta_HLPS3, Mstat_HLPS3/Rstat_HLPS3,DivSepM(delta_HLPS3,(Mstat_HLPS3/Rstat_HLPS3),(M_HLPS3)/Mstat_HLPS3,fittedParametersBoth),s=3, label = 'HLPS3 EOS')
    ax.scatter(delta_L, Mstat_L/Rstat_L,DivSepM(delta_L,(Mstat_L/Rstat_L),(M_L)/Mstat_L,fittedParametersBoth),s=3, label = 'L EOS')
    ax.scatter(delta_H0, Mstat_H0/Rstat_H0,DivSepM(delta_H0,(Mstat_H0/Rstat_H0),(M_H0)/Mstat_H0,fittedParametersBoth),s=3, label = 'H0 EOS')
    #ax.scatter(delta_Q1, Mstat_Q1/Rstat_Q1,DivSepR(delta_Q1,(Mstat_Q1/Rstat_Q1),(R_Q1)/Rstat_Q1,fittedParametersBothR),s=3, label = 'QS1 EOS')
    ax.scatter(delta_QHCD, Mstat_QHCD/Rstat_QHCD,DivSepM(delta_QHCD,(Mstat_QHCD/Rstat_QHCD),(M_QHCD)/Mstat_QHCD,fittedParametersBoth),s=3, label = 'QHCD EOS')
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'Dev($M/M_*$)', fontsize='24')
    
    ax.view_init(azim=130,elev=0)
    plt.legend()
    plt.savefig('Div_from_best_fitM_both_nuc4.png', format = 'png', transparent=False)
    plt.show()
elif method==17: 


   
    x = deltanuc
    y = (Mstatnuc/Rstatnuc)
    z = np.abs(DivSepM(deltanuc,(Mstatnuc/Rstatnuc),Mnuc/Mstatnuc,fittedParametersBoth))
    


    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31
    
    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin



    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])

    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    plt.clim(0, 0.06)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(M/M_*)$ (other NS EOSs)", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto', vmin = 0, vmax =0.06)
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto', vmin = 0, vmax =0.06)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)

    plt.savefig('DiverMAllnuc.png', format = 'png', transparent=False)
    plt.show()

elif method==18:
    plt.scatter(delta_Pol0, M_Pol0/Mstat_Pol0, s=3, label = 'EOS Pol0')
    #plt.plot(delta_Pol0nr, M_Pol0nr/Mstat_Pol0nr,'k-',linewidth=3)
    #plt.plot(delta_Pol0K, M_Pol0K/Mstat_Pol0K,'k--',linewidth=3)
    #plt.plot(delta_Pol0ins, M_Pol0ins/Mstat_Pol0ins,'k--')
    #plt.plot(delta_Pol0J, M_Pol0J/Mstat_Pol0J,'k-')
    
            
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$M/M_*$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})  
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Mstat-frac_Ang-vel_alleospol.png', format = 'png', transparent=False)
    plt.show()
elif method == 19:
    ax = fig.add_subplot(projection='3d')

    #ax.plot_trisurf(delta**2, Mstat/Rstat, zC_fit)
    ax.scatter(delta_Pol, (Mstat_Pol/Rstat_Pol),M_Pol/R_Pol, label = 'PP EOS')
    ax.scatter(delta_Gr, (Mstat_Gr/Rstat_Gr),M_Gr/R_Gr, label = '$c_s$ EOS')
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'$C_e$ ', fontsize='24')
    
    ax.view_init(azim=-119,elev=15)
    plt.legend()
    plt.savefig('Compact.png', format = 'png', transparent=False)
    plt.show()

elif method == 20:
    ax = fig.add_subplot(projection='3d')
    ax.scatter(delta, (Mstat/Rstat),DivSepC(delta,(Mstat/Rstat),M/R,fittedParametersComp), label = 'All EOS')
 
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'$Dev(C_e)$', fontsize='24')
    
    ax.view_init(azim=-95,elev=0)
    plt.legend()
    plt.savefig('Diverg_for_Comp.png', format = 'png', transparent=False)
    plt.show()

elif method==21:    
    x = delta
    y = (Mstat/Rstat)
    z = np.abs(DivSepC(delta,(Mstat/Rstat),M/R,fittedParametersComp))
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31
    
    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin




    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])

    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    #plt.clim(0, 0.1431)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(C_e)$", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')#, vmin = 0, vmax =0.1431)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.savefig('DiverCAll.png', format = 'png', transparent=False)
    plt.show()


elif method == 22:
    
    ax = fig.add_subplot(projection='3d')
    ax.scatter(delta_APR, Mstat_APR/Rstat_APR,DivSepC(delta_APR,(Mstat_APR/Rstat_APR),(M_APR)/R_APR,fittedParametersComp),s=3, label = 'APR EOS')
    ax.scatter(delta_BBB1, Mstat_BBB1/Rstat_BBB1,DivSepC(delta_BBB1,(Mstat_BBB1/Rstat_BBB1),(M_BBB1)/R_BBB1,fittedParametersComp),s=3, label = 'BBB1 EOS')
    ax.scatter(delta_BBB2, Mstat_BBB2/Rstat_BBB2,DivSepC(delta_BBB2,(Mstat_BBB2/Rstat_BBB2),(M_BBB2)/R_BBB2,fittedParametersComp),s=3, label = 'BBB2 EOS')
    ax.scatter(delta_ABPR1, Mstat_ABPR1/Rstat_ABPR1,DivSepC(delta_ABPR1,(Mstat_ABPR1/Rstat_ABPR1),(M_ABPR1)/R_ABPR1,fittedParametersComp),s=3, label = 'ABPR1 EOS')
    ax.scatter(delta_ABPR2, Mstat_ABPR2/Rstat_ABPR2,DivSepC(delta_ABPR2,(Mstat_ABPR2/Rstat_ABPR2),(M_ABPR2)/R_ABPR2,fittedParametersComp),s=3, label = 'ABPR2 EOS')
    ax.scatter(delta_HLPS1, Mstat_HLPS1/Rstat_HLPS1,DivSepC(delta_HLPS1,(Mstat_HLPS1/Rstat_HLPS1),(M_HLPS1)/R_HLPS1,fittedParametersComp),s=3, label = 'HLPS1 EOS')
    ax.scatter(delta_HLPS2, Mstat_HLPS2/Rstat_HLPS2,DivSepC(delta_HLPS2,(Mstat_HLPS2/Rstat_HLPS2),(M_HLPS2)/R_HLPS2,fittedParametersComp),s=3, label = 'HLPS2 EOS')
    #ax.scatter(delta_HLPS3, Mstat_HLPS3/Rstat_HLPS3,DivSepC(delta_HLPS3,(Mstat_HLPS3/Rstat_HLPS3),(M_HLPS3)/R_HLPS3,fittedParametersComp),s=3, label = 'HLPS3 EOS')
    ax.scatter(delta_L, Mstat_L/Rstat_L,DivSepC(delta_L,(Mstat_L/Rstat_L),(M_L)/R_L,fittedParametersComp),s=3, label = 'L EOS')
    ax.scatter(delta_H0, Mstat_H0/Rstat_H0,DivSepC(delta_H0,(Mstat_H0/Rstat_H0),(M_H0)/R_H0,fittedParametersComp),s=3, label = 'H0 EOS')
    #ax.scatter(delta_Q1, Mstat_Q1/Rstat_Q1,DivSepC(delta_Q1,(Mstat_Q1/Rstat_Q1),(M_Q1)/R_Q1,fittedParametersComp),s=3, label = 'QS1 EOS')
    ax.scatter(delta_QHCD, Mstat_QHCD/Rstat_QHCD,DivSepC(delta_QHCD,(Mstat_QHCD/Rstat_QHCD),(M_QHCD)/R_QHCD,fittedParametersComp),s=3, label = 'QHCD EOS')
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'$Dev(C_e)$', fontsize='24')
    
    ax.view_init(azim=-95,elev=0)
    plt.legend()
    plt.savefig('Div_from_best_fitComp_both_nuc.png', format = 'png', transparent=False)
    plt.show()
elif method==23: 


   
    x = deltanuc
    y = (Mstatnuc/Rstatnuc)
    z = np.abs(DivSepC(deltanuc,(Mstatnuc/Rstatnuc),Mnuc/Rnuc,fittedParametersComp))
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31
    
    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin



    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    
 

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.clim(0, 0.1431)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(C_e)$ (other NS EOSs)", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto', vmin = 0, vmax =0.1431)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.savefig('DiverCAllnuc.png', format = 'png', transparent=False)
    plt.show()
elif method==24:
    plt.scatter(delta_Pol0, M_Pol0/R_Pol0, s=3, label = 'EOS Pol0')
    plt.plot(delta_Pol0nr, M_Pol0nr/R_Pol0nr,'k-',linewidth=3)
    plt.plot(delta_Pol0K, M_Pol0K/R_Pol0K,'k--',linewidth=3)
    #plt.plot(delta_Pol0ins, M_Pol0ins/R_Pol0ins,'k--')
    plt.plot(delta_Pol0J, M_Pol0J/R_Pol0J,'k-')
    plt.plot(delta_Pol0_rho1, M_Pol0_rho1/R_Pol0_rho1,'r')
    plt.plot(delta_Pol0_rho2, M_Pol0_rho2/R_Pol0_rho2,'y')
    x=delta_Pol0ins
    y=M_Pol0ins/R_Pol0ins
    xnew = np.linspace(x.min(), x.max(), 200) 

    #define spline
    spl = make_interp_spline(x, y, k=3)
    y_smooth = spl(xnew)
    plt.plot(xnew,y_smooth,'k--')
            
    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_e\ $ ', fontsize='24')
    plt.legend(loc='upper right')
    plt.legend(prop={"size":16}) 
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Comp_eosPP0.png', format = 'png', transparent=False)
    plt.show()
    
elif method==25:    
    x = delta
    y = (Mstat/Rstat)
    
    
    Memp=(1+(np.exp(1.127*x**2)-1)*(-0.0160+3.123*y-20.721*y**2+41.202*y**3-6.464 *y**4))*Mstat

    Remp=(1+(np.exp(0.203*x**2)-1+0.1611*np.log(1-(x/1.1)**4)**2)*(1-15.496*y+442.6*y**2-4945.62 *y**3+23458.06*y**4-40544.25*y**5))*Rstat 

    
    
    z = np.abs((M/R-Memp/Remp)/(M/R))
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31
    
    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin




    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])

    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    #plt.clim(0, 0.1431)   
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    plt.title("$Dev(C_e)$ from R/$R_*$ and M/$M_*$ best fit equations", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')#, vmin = 0, vmax =0.1431)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.savefig('DiverCAllBasedonMandR.png', format = 'png', transparent=False)
    plt.show()
elif method == 26:


        
    dum1=np.zeros(len(R_APR))
    dum2=np.zeros(len(R_APR))
    countk=0
    for i in range(len(R_APR)):
      if (Rratio_APR[i]==1)  and (freq_APR[i]==0):
          dum1[i] = M_APR[i]
          dum2[i] = R_APR[i]
          countk=countk+1

    M_APR0=np.zeros(countk)
    R_APR0=np.zeros(countk)
    countk=0
    for i in range(len(R_APR)):
      if dum1[i]!=0:
        M_APR0[countk]=dum1[i]
        R_APR0[countk]=dum2[i]
        countk=countk+1
    
    dum1=np.zeros(len(R_BBB1))
    dum2=np.zeros(len(R_BBB1))
    countk=0
    for i in range(len(R_BBB1)):
      if (Rratio_BBB1[i]==1)  and (freq_BBB1[i]==0):
          dum1[i] = M_BBB1[i]
          dum2[i] = R_BBB1[i]
          countk=countk+1

    M_BBB10=np.zeros(countk)
    R_BBB10=np.zeros(countk)
    countk=0
    for i in range(len(R_BBB1)):
      if dum1[i]!=0:
        M_BBB10[countk]=dum1[i]
        R_BBB10[countk]=dum2[i]
        countk=countk+1
        
    dum1=np.zeros(len(R_BBB2))
    dum2=np.zeros(len(R_BBB2))
    countk=0
    for i in range(len(R_BBB2)):
      if (Rratio_BBB2[i]==1)  and (freq_BBB2[i]==0):
          dum1[i] = M_BBB2[i]
          dum2[i] = R_BBB2[i]
          countk=countk+1

    M_BBB20=np.zeros(countk)
    R_BBB20=np.zeros(countk)
    countk=0
    for i in range(len(R_BBB2)):
      if dum1[i]!=0:
        M_BBB20[countk]=dum1[i]
        R_BBB20[countk]=dum2[i]
        countk=countk+1
        
    dum1=np.zeros(len(R_HLPS1))
    dum2=np.zeros(len(R_HLPS1))
    countk=0
    for i in range(len(R_HLPS1)):
      if (Rratio_HLPS1[i]==1)  and (freq_HLPS1[i]==0):
          dum1[i] = M_HLPS1[i]
          dum2[i] = R_HLPS1[i]
          countk=countk+1

    M_HLPS10=np.zeros(countk)
    R_HLPS10=np.zeros(countk)
    countk=0
    for i in range(len(R_HLPS1)):
      if dum1[i]!=0:
        M_HLPS10[countk]=dum1[i]
        R_HLPS10[countk]=dum2[i]
        countk=countk+1
        
    dum1=np.zeros(len(R_HLPS2))
    dum2=np.zeros(len(R_HLPS2))
    countk=0
    for i in range(len(R_HLPS2)):
      if (Rratio_HLPS2[i]==1)  and (freq_HLPS2[i]==0):
          dum1[i] = M_HLPS2[i]
          dum2[i] = R_HLPS2[i]
          countk=countk+1

    M_HLPS20=np.zeros(countk)
    R_HLPS20=np.zeros(countk)
    countk=0
    for i in range(len(R_HLPS2)):
      if dum1[i]!=0:
        M_HLPS20[countk]=dum1[i]
        R_HLPS20[countk]=dum2[i]
        countk=countk+1
        
    dum1=np.zeros(len(R_HLPS3))
    dum2=np.zeros(len(R_HLPS3))
    countk=0
    for i in range(len(R_HLPS3)):
      if (Rratio_HLPS3[i]==1)  and (freq_HLPS3[i]==0):
          dum1[i] = M_HLPS3[i]
          dum2[i] = R_HLPS3[i]
          countk=countk+1

    M_HLPS30=np.zeros(countk)
    R_HLPS30=np.zeros(countk)
    countk=0
    for i in range(len(R_HLPS3)):
      if dum1[i]!=0:
        M_HLPS30[countk]=dum1[i]
        R_HLPS30[countk]=dum2[i]
        countk=countk+1
        
    dum1=np.zeros(len(R_ABPR1))
    dum2=np.zeros(len(R_ABPR1))
    countk=0
    for i in range(len(R_ABPR1)):
      if (Rratio_ABPR1[i]==1)  and (freq_ABPR1[i]==0):
          dum1[i] = M_ABPR1[i]
          dum2[i] = R_ABPR1[i]
          countk=countk+1

    M_ABPR10=np.zeros(countk)
    R_ABPR10=np.zeros(countk)
    countk=0
    for i in range(len(R_ABPR1)):
      if dum1[i]!=0:
        M_ABPR10[countk]=dum1[i]
        R_ABPR10[countk]=dum2[i]
        countk=countk+1

    dum1=np.zeros(len(R_ABPR2))
    dum2=np.zeros(len(R_ABPR2))
    countk=0
    for i in range(len(R_ABPR2)):
      if (Rratio_ABPR2[i]==1)  and (freq_ABPR2[i]==0):
          dum1[i] = M_ABPR2[i]
          dum2[i] = R_ABPR2[i]
          countk=countk+1

    M_ABPR20=np.zeros(countk)
    R_ABPR20=np.zeros(countk)
    countk=0
    for i in range(len(R_ABPR2)):
      if dum1[i]!=0:
        M_ABPR20[countk]=dum1[i]
        R_ABPR20[countk]=dum2[i]
        countk=countk+1
        

    dum1=np.zeros(len(R_QHCD))
    dum2=np.zeros(len(R_QHCD))
    countk=0
    for i in range(len(R_QHCD)):
      if (Rratio_QHCD[i]==1)  and (freq_QHCD[i]==0):
          dum1[i] = M_QHCD[i]
          dum2[i] = R_QHCD[i]
          countk=countk+1

    M_QHCD0=np.zeros(countk)
    R_QHCD0=np.zeros(countk)
    countk=0
    for i in range(len(R_QHCD)):
      if dum1[i]!=0:
        M_QHCD0[countk]=dum1[i]
        R_QHCD0[countk]=dum2[i]
        countk=countk+1

    dum1=np.zeros(len(R_L))
    dum2=np.zeros(len(R_L))
    countk=0
    for i in range(len(R_L)):
      if (Rratio_L[i]==1)  and (freq_L[i]==0):
          dum1[i] = M_L[i]
          dum2[i] = R_L[i]
          countk=countk+1

    M_L0=np.zeros(countk)
    R_L0=np.zeros(countk)
    countk=0
    for i in range(len(R_L)):
      if dum1[i]!=0:
        M_L0[countk]=dum1[i]
        R_L0[countk]=dum2[i]
        countk=countk+1
              
    dum1=np.zeros(len(R_Q1))
    dum2=np.zeros(len(R_Q1))
    countk=0
    for i in range(len(R_Q1)):
      if (Rratio_Q1[i]==1)  and (freq_Q1[i]==0):
          dum1[i] = M_Q1[i]
          dum2[i] = R_Q1[i]
          countk=countk+1

    M_Q10=np.zeros(countk)
    R_Q10=np.zeros(countk)
    countk=0
    for i in range(len(R_Q1)):
      if dum1[i]!=0:
        M_Q10[countk]=dum1[i]
        R_Q10[countk]=dum2[i]
        countk=countk+1
                                     
    dum1=np.zeros(len(R_H0))
    dum2=np.zeros(len(R_H0))
    countk=0
    for i in range(len(R_H0)):
      if (Rratio_H0[i]==1)  and (freq_H0[i]==0):
          dum1[i] = M_H0[i]
          dum2[i] = R_H0[i]
          countk=countk+1

    M_H00=np.zeros(countk)
    R_H00=np.zeros(countk)
    countk=0
    for i in range(len(R_H0)):
      if dum1[i]!=0:
        M_H00[countk]=dum1[i]
        R_H00[countk]=dum2[i]
        countk=countk+1


    plt.plot(R_APR0, M_APR0, '-', label = 'EOS APR')
    plt.plot(R_BBB10, M_BBB10, '-', label = 'EOS BBB1')
    plt.plot(R_BBB20, M_BBB20, '-', label = 'EOS BBB2')
    plt.plot(R_L0, M_L0, '-', label = 'EOS L')
    plt.plot(R_Q10, M_Q10, '-', label = 'EOS Q160')
    plt.plot(R_QHCD0, M_QHCD0, '-', label = 'EOS QHCD')
    plt.plot(R_H00, M_H00, '-', label = 'EOS H0')
    plt.plot(R_ABPR10, M_ABPR10, '-', label = 'EOS ABPR1')
    plt.plot(R_ABPR20, M_ABPR20, '-', label = 'EOS ABPR2')
    plt.plot(R_HLPS10, M_HLPS10, '-', label = 'EOS HLPS1')
    plt.plot(R_HLPS20, M_HLPS20, '--', label = 'EOS HLPS2')
    
    plt.xlim([6, 16])    
    plt.xlabel(r'$R$ [km]', fontsize='24')
    plt.ylabel(r'$M$ [M$_\odot$]', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})  
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Mass_Radius_allothereos.png', format = 'png', transparent=False)
    plt.show()    


        

    
elif method == 27:    
    delta=freq * (2*np.pi) * ((r**3)/(G*m))**(0.5)
    #dame
    NorK=np.zeros(len(delta))
    ComK=np.zeros(len(delta))
    countk=0
    for i in range(len(delta)):
      if i>1:
        if (Rratio[i]==1) and (Rratio[i-1]!=1) and (kfreq[i-1]-freq[i-1]<100):
          NorK[i-1] = freq[i-1] * (2*np.pi) * ((r[i-1]**3)/(G*m[i-1]))**(0.5)
          ComK[i-1] = M[i-1]/R[i-1]
          countk=countk+1

    Nor2Kepler=np.zeros(countk)
    ComKepler=np.zeros(countk)
    countk=0
    for i in range(len(delta)):
      if NorK[i]!=0:
        Nor2Kepler[countk]=NorK[i]
        ComKepler[countk]=ComK[i]# * G* 1.9884e33/ (1.0e5*c**2) 
        countk=countk+1


    SlopK4,SlopK3,SlopK2,SlopK,InterK=np.polyfit(ComKepler,Nor2Kepler,4)

    print(np.max(np.abs(((SlopK4*(ComKepler)**4+SlopK3*(ComKepler)**3+SlopK2*(ComKepler)**2+SlopK*ComKepler+InterK)-Nor2Kepler)/(SlopK4*(ComKepler)**4+SlopK3*(ComKepler)**3+SlopK2*(ComKepler)**2+SlopK*ComKepler+InterK))))

    NormD2=(SlopK4*(Mstat/Rstat)**4+SlopK3*(Mstat/Rstat)**3+SlopK2*(Mstat/Rstat)**2+SlopK*Mstat/Rstat+InterK)

    delta=delta/NormD2#**0.5
    
    
    
    
    xData = numpy.array(delta)
    yData = numpy.array(M/R)
    zData = numpy.array(R/Rstat)

    data = [xData, yData, zData]

    initialParameters = np.ones(Rcoef) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersBothR, pcov = scipy.optimize.curve_fit(funcR, [xData, yData], zData, p0 = initialParameters)

    print('Radius fitted prameters', fittedParametersBothR)

    modelPredictionsBothR = funcR(data, *fittedParametersBothR) 

    absError = modelPredictionsBothR - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
    print('RMSE:', RMSE)
    print('R-squared Both:', Rsquared)




    x = delta
    y = (M/R)
    z = np.abs(DivSepR(delta,(M/R),(R)/Rstat,fittedParametersBothR))
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31


    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin



    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels     

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels    
    
    plt.xlabel(r'$\Omega_{n2}$', fontsize='24')
    plt.ylabel(r'$C_e$ ', fontsize='24')
    plt.title("$Dev(R_e/R_*)$", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)

    plt.savefig('DiverRAll2.png', format = 'png', transparent=False)
    plt.show()
    
   
    
elif method == 28:    

    delta=freq * (2*np.pi) * ((r**3)/(G*m))**(0.5)
    #dame
    NorK=np.zeros(len(delta))
    ComK=np.zeros(len(delta))
    countk=0
    for i in range(len(delta)):
      if i>1:
        if (Rratio[i]==1) and (Rratio[i-1]!=1) and (kfreq[i-1]-freq[i-1]<100):
          NorK[i-1] = freq[i-1] * (2*np.pi) * ((r[i-1]**3)/(G*m[i-1]))**(0.5)
          ComK[i-1] = M[i-1]/R[i-1]
          countk=countk+1

    Nor2Kepler=np.zeros(countk)
    ComKepler=np.zeros(countk)
    countk=0
    for i in range(len(delta)):
      if NorK[i]!=0:
        Nor2Kepler[countk]=NorK[i]
        ComKepler[countk]=ComK[i]# * G* 1.9884e33/ (1.0e5*c**2) 
        countk=countk+1


    SlopK4,SlopK3,SlopK2,SlopK,InterK=np.polyfit(ComKepler,Nor2Kepler,4)

    print(np.max(np.abs(((SlopK4*(ComKepler)**4+SlopK3*(ComKepler)**3+SlopK2*(ComKepler)**2+SlopK*ComKepler+InterK)-Nor2Kepler)/Nor2Kepler)))

 
    '''
    plt.scatter(ComKepler,Nor2Kepler,  s=3, label = 'PP EOS')
    plt.scatter(ComKeplerG,Nor2KeplerG,  s=3, label = '$c_s$ EOS')   
    
    
    x=np.array(range(0, 220))/1000.
    plt.plot(x,SlopK4*x*x*x*x+SlopK3*x*x*x+SlopK2*x*x+SlopK*x+InterK, 'k', label = 'Best Fit') 

    plt.ylabel(r'$\Omega_K \sqrt{(R^3 / GM)}$', fontsize='24')
    plt.xlabel(r'$C_e$ ', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Kepler2.png', format = 'png', transparent=False)
    plt.show()    
    '''
    NormD2=(SlopK4*(Mstat/Rstat)**4+SlopK3*(Mstat/Rstat)**3+SlopK2*(Mstat/Rstat)**2+SlopK*Mstat/Rstat+InterK)

    delta=delta/NormD2#**0.5
    
        
    absError = (SlopK4*(ComKepler)**4+SlopK3*(ComKepler)**3+SlopK2*(ComKepler)**2+SlopK*ComKepler+InterK) - Nor2Kepler

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    RsquaredK = 1.0 - (numpy.var(absError) / numpy.var(Nor2Kepler))


    
    
    
     
    xData = numpy.array(delta)
    yData = numpy.array(M/R)
    zData = numpy.array((M)/Mstat)

    data = [xData, yData, zData]

    initialParameters = np.ones(7) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersBoth, pcov = scipy.optimize.curve_fit(func1, [xData, yData], zData, p0 = initialParameters)

    #print('Mass fitted prameters', fittedParametersBoth)

    modelPredictionsBoth = func1(data, *fittedParametersBoth) 

    absError = modelPredictionsBoth - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    RsquaredMass = 1.0 - (numpy.var(absError) / numpy.var(zData))
    #print('RMSE:', RMSE)
    #print('R-squared Both:', Rsquared)





    '''
    x = delta
    y = (M/R)
    z = np.abs(DivSepM1(delta,(M/R),M/Mstat,fittedParametersBoth))
    

    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31


    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin


    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels 
    

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels     

    plt.xlabel(r'$\Omega_{n2}$', fontsize='24')
    plt.ylabel(r'$C_e$ ', fontsize='24')
    plt.title("$Dev(M/M_*)$", fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
    plt.savefig('DiverMAll2.png', format = 'png', transparent=False)
    plt.show()

    '''
    
    
    
    
        
    xData = numpy.array(delta)
    yData = numpy.array(M/R)
    zData = numpy.array(R/Rstat)

    data = [xData, yData, zData]

    initialParameters = np.ones(Rcoef) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersBothR, pcov = scipy.optimize.curve_fit(funcR, [xData, yData], zData, p0 = initialParameters)

    #print('Radius fitted prameters', fittedParametersBothR)

    modelPredictionsBothR = funcR(data, *fittedParametersBothR) 

    absError = modelPredictionsBothR - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    RsquaredRadi = 1.0 - (numpy.var(absError) / numpy.var(zData))
    #print('RMSE:', RMSE)
    #print('R-squared Both:', Rsquared)
#andreas
    print('Equation 7',SlopK4,SlopK3,SlopK2,SlopK,InterK,'- -', RsquaredK)
    print('Equation 8',fittedParametersBoth[0],fittedParametersBoth[1],fittedParametersBoth[2],fittedParametersBoth[3],fittedParametersBoth[4],fittedParametersBoth[5],fittedParametersBoth[6] ,RsquaredMass)
    print('Equation 9',fittedParametersBothR[0],fittedParametersBothR[1],fittedParametersBothR[2],fittedParametersBothR[3],fittedParametersBothR[4],fittedParametersBothR[5],fittedParametersBothR[6], RsquaredRadi)


    dum1=np.zeros(len(R_Pol0))
    dum2=np.zeros(len(R_Pol0))
    dum3=np.zeros(len(R_Pol0))
    countk=0
    for i in range(len(R_Pol0)):
     if(freq_Pol0[i]==0 and Rratio_Pol0[i]==1):
          dum1[i] = M_Pol0[i]
          dum2[i] = R_Pol0[i]
          dum3[i] = ec_Pol0[i]
          countk=countk+1

    M_Pl0=np.zeros(countk)
    R_Pl0=np.zeros(countk)
    ec_Pl0=np.zeros(countk)
    countk=0
    for i in range(len(R_Pol0)):
      if dum2[i]!=0:
        M_Pl0[countk]=dum1[i]
        R_Pl0[countk]=dum2[i]
        ec_Pl0[countk]=dum3[i]
        countk=countk+1
    
    dum1=np.zeros(len(R_Pol0))
    dum2=np.zeros(len(R_Pol0))
    dum3=np.zeros(len(R_Pol0))
    dum4=np.zeros(len(R_Pol0))
    countk=0
    for i in range(len(R_Pol0)):
     #if(delta_Pol0[i]<=0.9):
          dum1[i] = M_Pol0[i]
          dum2[i] = R_Pol0[i]
          dum3[i] = freq_Pol0[i]
          dum4[i] = ec_Pol0[i]
          countk=countk+1

    M_Pl1=np.zeros(countk)
    R_Pl1=np.zeros(countk)
    freq_Pl1=np.zeros(countk)
    ec_Pl1=np.zeros(countk)
    countk=0
    for i in range(len(R_Pol0)):
      if dum2[i]!=0:
        M_Pl1[countk]=dum1[i]
        R_Pl1[countk]=dum2[i]
        freq_Pl1[countk]=dum3[i]
        ec_Pl1[countk]=dum4[i]
        countk=countk+1
    



    Comp_Pl0=M_Pl1/R_Pl1
    
    C=fittedParametersBoth    
    
    NormD2= (SlopK4*(Comp_Pl0)**4+SlopK3*(Comp_Pl0)**3+SlopK2*(Comp_Pl0)**2+SlopK*Comp_Pl0+InterK)


    delta_Pl0=freq_Pl1 * (2*np.pi) * (( (R_Pl1* 1.0e5)**3)/(G*M_Pl1*1.9884e33))**(0.5)/NormD2

    
    Memp=M_Pl1/(1+(delta_Pl0+C[0]*delta_Pl0**2+C[1]*delta_Pl0**3+C[2]*delta_Pl0**4)*(C[3]*Comp_Pl0+C[4]*Comp_Pl0**2+C[5]*Comp_Pl0**3+C[6]*Comp_Pl0**4))
    
    
    C=fittedParametersBothR
    Remp=R_Pl1/(1+(np.exp(C[0]*delta_Pl0**2)-1-C[1]*np.log(1-(delta_Pl0/1.1)**4)**2)*(1+C[2]*Comp_Pl0+C[3]*Comp_Pl0**2+C[4]*Comp_Pl0**3+C[5]*Comp_Pl0**4+C[6]*Comp_Pl0**5)) 
    
    
    plt.scatter(R_Pl1,M_Pl1,  s=3, label = 'EOS Pol 0')
    plt.plot(R_Pl0,M_Pl0,  'k', label = 'M-R curve from TOV equations',linewidth=3)        
    #plt.plot(Remp,Memp,  'k--', label = 'Empirical equation',linewidth=3)
    plt.scatter(Remp,Memp,  s=3, label = 'Empirical equation')



    plt.xlabel(r'$R_e$ [km] ', fontsize='24')
    plt.ylabel(r'M [$M_\odot$]', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('MRcurve3.png', format = 'png', transparent=False)
    plt.show()
    
    '''
    plt.scatter(ec_Pl1,Memp,  s=3,color='orange', label = 'Empirical equation')
    plt.plot(ec_Pl0,M_Pl0,  'k', label = 'Non-rotating neutron stars',linewidth=2)        


    plt.xlabel(r'$e_c$ [Log(g/cm$^3$)] ', fontsize='24')
    plt.ylabel(r'M [$M_\odot$]', fontsize='24')
    plt.legend()
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('MRcurve4.png', format = 'png', transparent=False)
    plt.show()
    '''
elif method == 29:    
    



    dum1=np.zeros(len(R_Pol))
    dum2=np.zeros(len(R_Pol))
    dum3=np.zeros(len(R_Pol))
    countk=0
    for i in range(len(R_Pol)):
     if i>1:
      if (Rratio_Pol[i]==1) and (Rratio_Pol[i-1]!=1) and (kfreq_Pol[i-1]-freq_Pol[i-1]<100):
          dum1[i-1] = Mstat_Pol[i-1]/Rstat_Pol[i-1]
          dum2[i-1] = R_Pol[i-1]/Rstat_Pol[i-1]
          dum3[i-1] = M_Pol[i-1]/Mstat_Pol[i-1]
          countk=countk+1

    CIn_Pol=np.zeros(countk)
    RK_Pol=np.zeros(countk)
    MK_Pol=np.zeros(countk)
    countk=0
    for i in range(len(R_Pol)):
      if dum2[i]!=0:
        CIn_Pol[countk]=dum1[i]
        RK_Pol[countk]=dum2[i]
        MK_Pol[countk]=dum3[i]
        countk=countk+1


    dum1=np.zeros(len(R_Gr))
    dum2=np.zeros(len(R_Gr))
    dum3=np.zeros(len(R_Gr))
    countk=0
    for i in range(len(R_Gr)):
     if i>1:
      if (Rratio_Gr[i]==1) and (Rratio_Gr[i-1]!=1) and (kfreq_Gr[i-1]-freq_Gr[i-1]<100):
          dum1[i-1] = Mstat_Gr[i-1]/Rstat_Gr[i-1]
          dum2[i-1] = R_Gr[i-1]/Rstat_Gr[i-1]
          dum3[i-1] = M_Gr[i-1]/Mstat_Gr[i-1]
          countk=countk+1

    CIn_Gr=np.zeros(countk)
    RK_Gr=np.zeros(countk)
    MK_Gr=np.zeros(countk)
    countk=0
    for i in range(len(R_Gr)):
      if dum2[i]!=0:
        CIn_Gr[countk]=dum1[i]
        RK_Gr[countk]=dum2[i]
        MK_Gr[countk]=dum3[i]
        countk=countk+1




    dum1=np.zeros(len(R))
    dum2=np.zeros(len(R))
    dum3=np.zeros(len(R))
    countk=0
    for i in range(len(R)):
     if i>1:
      if (Rratio[i]==1) and (Rratio[i-1]!=1) and (kfreq[i-1]-freq[i-1]<100):
          dum1[i-1] = Mstat[i-1]/Rstat[i-1]
          dum2[i-1] = R[i-1]/Rstat[i-1]
          dum3[i-1] = M[i-1]/Mstat[i-1]
          countk=countk+1

    CIn=np.zeros(countk)
    RK=np.zeros(countk)
    MK=np.zeros(countk)
    countk=0
    for i in range(len(R)):
      if dum2[i]!=0:
        CIn[countk]=dum1[i]
        RK[countk]=dum2[i]
        MK[countk]=dum3[i]
        countk=countk+1
    
    SRK5,SRK4,SRK3,SRK2,SRK,IRK=np.polyfit(CIn,RK,5)

    print(np.max(np.abs(((SRK5*(CIn)**5+ SRK4*(CIn)**4+SRK3*(CIn)**3+SRK2*(CIn)**2+SRK*CIn+IRK)-RK)/RK)))

    absError = (SRK5*(CIn)**5+SRK4*(CIn)**4+SRK3*(CIn)**3+SRK2*(CIn)**2+SRK*CIn+IRK) - RK

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    RsquaredRK = 1.0 - (numpy.var(absError) / numpy.var(RK))
    print('Equation 5',SRK5,SRK4,SRK3,SRK2,SRK,IRK,'-', RsquaredRK)
    '''
    plt.scatter(CIn_Pol,RK_Pol,  s=3, label = 'PP EOS')
    plt.scatter(CIn_Gr,RK_Gr,  s=3, label = '$c_s$ EOS')
    
    x=np.array(range(0, 220))/1000.
    plt.plot(x,SRK5*x**5+ SRK4*x*x*x*x+SRK3*x*x*x+SRK2*x*x+SRK*x+IRK, 'k', label = 'Best fit') 
    
    plt.xlabel(r'$C_*$  ', fontsize='24')
    plt.ylabel(r'$R_e/R_*|_K$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('RchangeKepler.png', format = 'png', transparent=False)
    plt.show()
    '''
    
    
    
    SMK3,SMK2,SMK,IMK=np.polyfit(CIn,MK,3)

    print(np.max(np.abs(((SMK3*(CIn)**3+SMK2*(CIn)**2+SMK*CIn+IMK)-MK)/MK)))
    absError = (SMK3*(CIn)**3+SMK2*(CIn)**2+SMK*CIn+IMK) - MK

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    RsquaredMK = 1.0 - (numpy.var(absError) / numpy.var(MK))
    print('Equation 3',SMK3,SMK2,SMK,IMK,'- - -', RsquaredMK)
    '''
    plt.scatter(CIn_Pol,MK_Pol,  s=3, label = 'PP EOS')
    plt.scatter(CIn_Gr,MK_Gr,  s=3, label = '$c_s$ EOS')
    
    x=np.array(range(0, 220))/1000.
    plt.plot(x,SMK3*x*x*x+SMK2*x*x+SMK*x+IMK, 'k', label = 'Best fit') 
    
    plt.xlabel(r'$C_*$  ', fontsize='24')
    plt.ylabel(r'$M/M_*|_K$', fontsize='24')
    plt.legend(loc='best')
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('MchangeKepler.png', format = 'png', transparent=False)
    plt.show()
    '''
    dum1=np.zeros(len(R_Pol0))
    dum2=np.zeros(len(R_Pol0))
    countk=0
    for i in range(len(R_Pol0)):
     if(delta_Pol0[i]<=0.95):
          dum1[i] = M_Pol0[i]
          dum2[i] = R_Pol0[i]
          countk=countk+1

    Ms09=np.zeros(countk)
    Rs09=np.zeros(countk)
    countk=0
    for i in range(len(R_Pol0)):
      if dum2[i]!=0:
        Ms09[countk]=dum1[i]
        Rs09[countk]=dum2[i]
        countk=countk+1
    
    dum1=np.zeros(len(R_Pol0))
    dum2=np.zeros(len(R_Pol0))
    countk=0
    for i in range(len(R_Pol0)):
     if(delta_Pol0[i]>0.95):
          dum1[i] = M_Pol0[i]
          dum2[i] = R_Pol0[i]
          countk=countk+1

    Ml09=np.zeros(countk)
    Rl09=np.zeros(countk)
    countk=0
    for i in range(len(R_Pol0)):
      if dum2[i]!=0:
        Ml09[countk]=dum1[i]
        Rl09[countk]=dum2[i]
        countk=countk+1

    







    
   
    dum1=np.zeros(len(R_Pol0))
    dum2=np.zeros(len(R_Pol0))
    countk=0
    for i in range(len(R_Pol0)):
     if(freq_Pol0[i]==0 and Rratio_Pol0[i]==1):
          dum1[i] = M_Pol0[i]
          dum2[i] = R_Pol0[i]
          countk=countk+1

    M_Pl0=np.zeros(countk)
    R_Pl0=np.zeros(countk)
    countk=0
    for i in range(len(R_Pol0)):
      if dum2[i]!=0:
        M_Pl0[countk]=dum1[i]
        R_Pl0[countk]=dum2[i]
        countk=countk+1
    
    Comp_Pl0=M_Pl0/R_Pl0
    
    C=fittedParametersBoth
    #print(C)
    Memp=(1+(np.exp(C[0]*0.95**2)-1)*(C[1]+C[2]*Comp_Pl0+C[3]*Comp_Pl0**2+C[4]*Comp_Pl0**3+C[5]*Comp_Pl0**4))*M_Pl0
    C=fittedParametersBothR
    #print(C)
    Remp=(1+(np.exp(C[0]*0.95**2)-1-C[1]*np.log(1-(0.95/1.1)**4)**2)*(1+C[2]*Comp_Pl0+C[3]*Comp_Pl0**2+C[4]*Comp_Pl0**3+C[5]*Comp_Pl0**4+C[6]*Comp_Pl0**5))*R_Pl0 
    #Remp=(1+(np.exp(-0.0633397572*1**2)-1-0.201646944*np.log(1-(1/1.1)**2))*(1+8.10970064*Comp_Pl0-69.2609273*Comp_Pl0**2+130.809880*Comp_Pl0**3))*R_Pl0 
       # Rd=Rstat*( 1+(np.exp(-0.0633397572*O**2)-1-0.201646944*np.log(1-(O/1.1)**2))*(1+8.10970064*C-69.2609273*C**2+130.809880*C**3) )
       

                            
    
    plt.plot(R_Pl0,M_Pl0,  'k', label = 'M-R curve from TOV equations',linewidth=3)        
    plt.plot(Remp,Memp,  'b-', label = 'Empirical equation for $\Omega_n = 0.95$',linewidth=3)
    plt.plot((SRK5*(Comp_Pl0)**5+ SRK4*(Comp_Pl0)**4+SRK3*(Comp_Pl0)**3+SRK2*(Comp_Pl0)**2+SRK*Comp_Pl0+IRK)*R_Pl0 ,(SMK3*Comp_Pl0**3+SMK2*Comp_Pl0**2+SMK*Comp_Pl0+IMK)*M_Pl0,  'r-', label = 'Empirical M-R curve at the Kepler limit',linewidth=3)
    

    plt.scatter(Rs09,Ms09,  s=3, label = 'EOS Pol 0 ($\Omega_n<=0.95$)',color='blue')
    plt.scatter(Rl09,Ml09,  s=3, label = 'EOS Pol 0 ($\Omega_n>0.95$)',color='red')


    plt.xlabel(r'$R_e$ [km] ', fontsize='24')
    plt.ylabel(r'M [$M_\odot$]', fontsize='24')
    plt.legend(loc='upper right')
    plt.legend(prop={"size":14})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('MRcurve.png', format = 'png', transparent=False)
    plt.show() 


elif method ==30:
    ax = fig.gca(projection='3d')
    #ax.scatter(delta_Pol, (Mstat_Pol/Rstat_Pol),(M0_Pol-M_Pol)/M_Pol, label = 'PP EOS')
    #ax.scatter(delta_Gr, (Mstat_Gr/Rstat_Gr),(M0_Gr-M_Gr)/M_Gr, label = '$c_s$ EOS')
    ax.scatter(delta_Pol, (Mstat_Pol/Rstat_Pol),Rratio_Pol, label = 'PP EOS')
    ax.scatter(delta_Gr, (Mstat_Gr/Rstat_Gr),Rratio_Gr, label = '$c_s$ EOS')
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$C_*$ ', fontsize='24')
    ax.set_zlabel(r'$Rratio$', fontsize='24')
    #ax.set_zlabel(r'$E_b/M$', fontsize='24')
    ax.view_init(azim=-90,elev=16)
    #ax.view_init(azim=-120,elev=17)
    plt.legend(loc='best')
    plt.savefig('Rratio3d.png', format = 'png', transparent=False)
    #plt.savefig('BindEChangeBothM.png', format = 'png', transparent=False)
    plt.show()


elif method==31:
    xData = numpy.array(delta)
    yData = numpy.array((Mstat/Rstat))
    zData = numpy.array(Rratio)
    #zData = numpy.array((M0-M)/M)
    data = [xData, yData, zData]

    initialParameters = np.ones(10) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersPol, pcov = scipy.optimize.curve_fit(funcPol, [xData, yData], zData, p0 = initialParameters)
    
    print('fitted prameters', fittedParametersPol)

    modelPredictionsPol = funcPol(data, *fittedParametersPol) 

    absError = modelPredictionsPol - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
    print('RMSE:',   RMSE)
    print('R-squared Comp:', Rsquared)


    '''
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xData, yData,DivSepPol(xData,yData,zData,fittedParametersPol), label = 'All EOS')
   
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$(M*/R*)$', fontsize='24')
    ax.set_zlabel(r'Rratio Divergence', fontsize='24')
    
    ax.view_init(azim=134,elev=7)
    plt.legend()
    plt.savefig('Diverg_for_Rratio_altern.png', format = 'png', transparent=False)
    plt.show()
    


    '''
    x = delta
    y = Mstat/Rstat
    z = np.abs(DivSepPol(delta,Mstat/Rstat,Rratio,fittedParametersPol))
    #z = np.abs(DivSepPol(delta,Mstat/Rstat,(M0-M)/M,fittedParametersPol))    
    #z = np.abs(DivSepPol(delta,Mstat/Rstat,(M0-M)/M,fittedParametersPol))    
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31


    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin


    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])

    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels     

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels      

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    #plt.title("$E_b$/M Deviation")
    plt.title("Dev($R_{ratio}$)",fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
   
    plt.savefig('DiverRratio.png', format = 'png', transparent=False)
    #plt.savefig('DiverAllEb.png', format = 'png', transparent=False)
    plt.show()

elif method==32:
    xData = numpy.array(delta)
    yData = numpy.array((Mstat/Rstat))
    zData = numpy.array(RratioS)
    #zData = numpy.array((M0-M)/M)
    data = [xData, yData, zData]

    initialParameters = np.ones(10) # these are the same as scipy default values in this example

    # here a non-linear surface fit is made with scipy's curve_fit()
    fittedParametersPol, pcov = scipy.optimize.curve_fit(funcPol, [xData, yData], zData, p0 = initialParameters)
    
    print('fitted prameters', fittedParametersPol)

    modelPredictionsPol = funcPol(data, *fittedParametersPol) 

    absError = modelPredictionsPol - zData

    SE = numpy.square(absError) # squared errors
    MSE = numpy.mean(SE) # mean squared errors
    RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (numpy.var(absError) / numpy.var(zData))
    print('RMSE:',   RMSE)
    print('R-squared Comp:', Rsquared)


    '''
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xData, yData,DivSepPol(xData,yData,zData,fittedParametersPol), label = 'All EOS')
   
    
    ax.set_xlabel(r'$\Omega_n$', fontsize='24')
    ax.set_ylabel(r'$(M*/R*)$', fontsize='24')
    ax.set_zlabel(r'Rratio Divergence', fontsize='24')
    
    ax.view_init(azim=134,elev=7)
    plt.legend()
    plt.savefig('Diverg_for_Rratio_altern.png', format = 'png', transparent=False)
    plt.show()
    


    '''
    x = delta
    y = Mstat/Rstat
    z = np.abs(DivSepPol(delta,Mstat/Rstat,RratioS,fittedParametersPol))
    #z = np.abs(DivSepPol(delta,Mstat/Rstat,(M0-M)/M,fittedParametersPol))    
    #z = np.abs(DivSepPol(delta,Mstat/Rstat,(M0-M)/M,fittedParametersPol))    
    
    xmax=np.max(x)
    xmin=np.min(x)
    ymax=np.max(y)
    ymin=np.min(y)
    xbin=61
    ybin=31






    #x_bins = np.linspace(xmin, xmax, xbin)
    y_bins = np.linspace(ymin, ymax,ybin)
    
    x_bins=numpy.concatenate((np.linspace(xmin, 0.2, 2, endpoint=False),np.linspace(0.2, 0.4, 4, endpoint=False), np.linspace(0.4, xmax, xbin, endpoint=False)))
    #y_bins=numpy.concatenate((np.linspace(ymin, 0.2, 20, endpoint=False), np.linspace(0.2, 0.4, 40, endpoint=False),np.linspace(0.4, ymax,ybin, endpoint=False)))    
    #z1=[[0]*xbin]*ybin


    #for i in range(len(x)):
    #  for i in range(len(ybin)):
    z1=np.zeros(len(x))
    H, xedges, yedges, binnumbers = stats.binned_statistic_2d(x, y, values=z, statistic='max', bins=[x_bins, y_bins])

    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels     

    #statistic1, xedges1, yedges1, binnumber1 = stats.binned_statistic_2d(x.ravel(), y.ravel(), values=z.ravel(), statistic='mean',bins=[np.arange(xmin, xmax, xbin), np.arange(ymin, ymax,ybin)])


    #plt.figure(2)
    #plt.pcolormesh(ret.x_edge,ret.y_edge,ret.statistic, vmin = 0, vmax = 1)
    #H = np.ma.masked_where(H==0, H)
    plt.pcolormesh(xedges,yedges,H.T,  cmap='jet')
    #plt.colorbar(plot,ax=ax2, pad = .015, aspect=10)
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels      

    plt.xlabel(r'$\Omega_n$', fontsize='24')
    plt.ylabel(r'$C_*$ ', fontsize='24')
    #plt.title("$E_b$/M Deviation")
    plt.title("Dev($R_{ratio}$_S)",fontsize='24')
    plt.legend(loc='best')

    plt.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    #plt.colorbar()

    aspect = 10
    pad_fraction = 0.5
    ax = plt.gca()
    im = ax.imshow(H.T, origin='lower',  cmap='jet',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], interpolation='nearest', aspect='auto')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)

    plt.colorbar(im, cax=cax)
   
    plt.savefig('DiverRratioS.png', format = 'png', transparent=False)
    #plt.savefig('DiverAllEb.png', format = 'png', transparent=False)
    plt.show()



elif method == 33:
    ax = fig.add_subplot(projection='3d')
    #ax.plot_trisurf(xData, yData, zP7_fit)
    #ax.plot_trisurf(xData, yData, zG8_fit)
    #ax.plot_surface(xData, yData, zP7_fit,rstride=1, cstride=1, linewidth=0.25)
    #ax.plot_surface(xData, yData, zG8_fit,rstride=1, cstride=1, linewidth=0.25)
    ax.scatter(delta, (Mstat/Rstat),T/(M), label = 'All EOS')

    ax.set_xlabel(r'$\Omega_n$',labelpad=12, fontsize='18')
    ax.set_ylabel(r'$C_*$ ',labelpad=12, fontsize='18')
    ax.set_zlabel(r'$T/M$',labelpad=12, fontsize='18')


    #ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0],fontsize=24)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="z", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    #ax.set_zticks(18)
        
    ax.view_init(azim=-106,elev=10)
    plt.legend(loc='best')
    plt.legend(prop={"size":16}) 

    plt.savefig('KineticEnergyChangeBoth.png', format = 'png', transparent=False)
    plt.show()
elif method == 34:
    ax = fig.add_subplot(projection='3d')
    #ax.plot_trisurf(xData, yData, zP7_fit)
    #ax.plot_trisurf(xData, yData, zG8_fit)
    #ax.plot_surface(xData, yData, zP7_fit,rstride=1, cstride=1, linewidth=0.25)
    #ax.plot_surface(xData, yData, zG8_fit,rstride=1, cstride=1, linewidth=0.25)
    ax.scatter(delta, (Mstat/Rstat),-W/(M), label = 'All EOS')
 
 
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="z", labelsize=15)
    ax.tick_params(axis="y", labelsize=15) 
    
    ax.set_xlabel(r'$\Omega_n$',labelpad=12, fontsize='18')
    ax.set_ylabel(r'$C_*$ ',labelpad=12, fontsize='18')
    ax.set_zlabel(r'$W/M$',labelpad=12, fontsize='18')
    
    ax.view_init(azim=-153,elev=9)
    plt.legend(loc='best')
    plt.legend(prop={"size":16}) 
    plt.savefig('GrBindEnChangeBoth.png', format = 'png', transparent=False)
    plt.show()
elif method == 35:
    ax = fig.add_subplot(projection='3d')
    #ax.plot_trisurf(xData, yData, zP7_fit)
    #ax.plot_trisurf(xData, yData, zG8_fit)
    #ax.plot_surface(xData, yData, zP7_fit,rstride=1, cstride=1, linewidth=0.25)
    #ax.plot_surface(xData, yData, zG8_fit,rstride=1, cstride=1, linewidth=0.25)
    ax.scatter(delta, (Mstat/Rstat),M0/(M), label = 'All EOS')

    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="z", labelsize=15)
    ax.tick_params(axis="y", labelsize=15) 
    
    
    ax.set_xlabel(r'$\Omega_n$',labelpad=12, fontsize='18')
    ax.set_ylabel(r'$C_*$ ',labelpad=12, fontsize='18')
    ax.set_zlabel(r'$M_0/M$',labelpad=12, fontsize='18')
    
    ax.view_init(azim=-146,elev=7)
    plt.legend(loc='best')
    plt.legend(prop={"size":16}) 
    plt.savefig('BarMassChangeBoth.png', format = 'png', transparent=False)
    plt.show()
elif method == 36:
    ax = fig.add_subplot(projection='3d')
    #ax.plot_trisurf(xData, yData, zP7_fit)
    #ax.plot_trisurf(xData, yData, zG8_fit)
    #ax.plot_surface(xData, yData, zP7_fit,rstride=1, cstride=1, linewidth=0.25)
    #ax.plot_surface(xData, yData, zG8_fit,rstride=1, cstride=1, linewidth=0.25)
    ax.scatter(delta, (Mstat/Rstat),(M-M0+W-T)/(M), label = 'All EOS')


    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="z", labelsize=15)
    ax.tick_params(axis="y", labelsize=15) 
    

    ax.set_xlabel(r'$\Omega_n$',labelpad=12, fontsize='18')
    ax.set_ylabel(r'$C_*$ ',labelpad=12, fontsize='18')
    ax.set_zlabel(r'$E/M$',labelpad=12, fontsize='18')
    
    ax.view_init(azim=-153,elev=15)
    plt.legend(loc='best')
    plt.legend(prop={"size":16}) 
    plt.savefig('InEnChangeBoth.png', format = 'png', transparent=False)
    plt.show()


    
elif method == 37:


    plt.scatter(ComKepler,Nor2Kepler,  s=3, label = 'PP EOS')
    plt.scatter(ComKeplerG,Nor2KeplerG,  s=3, label = '$c_s$ EOS')
    x=np.array(range(0, 220))/1000.
    plt.plot(x,SlopK4*x*x*x*x+SlopK3*x*x*x+SlopK2*x*x+SlopK*x+InterK, 'k', label = ' Best Fit'.format(SlopK4,SlopK3,SlopK2,SlopK,InterK)) 

    plt.ylabel(r'$\Omega_K \sqrt{(R_*^3 / GM_*)}$', fontsize='24')
    plt.xlabel(r'$C_*$ ', fontsize='24')
    plt.legend()
    plt.legend(prop={"size":16})
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('Kepler.png', format = 'png', transparent=False)
    plt.show()



    

    

else:
    print('Enter a number from 1 to 38')
