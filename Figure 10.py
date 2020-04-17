# -*- coding: cp1250 -*-
import math
r1 = 5.0*pow(10, -6)
r1_max = 5.2*pow(10, -6)
r1_min = 4.8*pow(10, -6)

nFeSiB = r1
nFeSiB_max = r1_max
nFeSiB_min = r1_min

nPyrex = 8.5*pow(10, -6)
nPyrex_max = 8.7*pow(10, -6)
nPyrex_min = 8.3*pow(10, -6)

nAu = 0.03*pow(10, -6)
nAu_max = 0.033*pow(10, -6)
nAu_min = 0.027*pow(10, -6)

nCoNi = 2.0*pow(10, -6)
nCoNi_max = 2.2*pow(10, -6)
nCoNi_min = 1.8*pow(10, -6)

#raduis
r2 = r1 + nPyrex + nAu + nCoNi
r2_max = r1_max + nPyrex_max + nAu_max + nCoNi_max
r2_min = r1_min+ nPyrex_min + nAu_min + nCoNi_min

r2_av = (r2_max + r2_min)/2


r2_error = (r2_max - r2_min)/2 


print ('(r2_max-r2_min)/2 =', (r2_max- r2_min)/2)

print ('r2 = ', r2)

print ('r2_av = ', r2_av)

def average_error(max, min):
    av = (max + min)/2
    er = (max - min)/2
    return av, 'Error =', er

print ('Total radius')
print ('r2 = ', average_error(r2_max, r2_min),'[m]')

SFeSiB = math.pi * r1**2
SFeSiB_max = math.pi*r1_max**2
SFeSiB_min = math.pi*r1_min**2

SPyrex = math.pi*(r1+nPyrex)**2 - math.pi*r1**2
SPyrex_max = math.pi*(r1_max+nPyrex_max)**2 - math.pi*r1_max**2
SPyrex_min = math.pi*(r1_min+nPyrex_min)**2 - math.pi*r1_min**2

SAu = math.pi*(r1+nPyrex+nAu)**2- math.pi*(r1+nPyrex)**2
SAu_max = math.pi*(r1_max+nPyrex_max+nAu_max)**2- math.pi*(r1_max+nPyrex_max)**2
SAu_min = math.pi*(r1_min+nPyrex_min+nAu_min)**2- math.pi*(r1_min+nPyrex_min)**2

SCoNi = math.pi*(r1+nPyrex+nAu+ nCoNi)**2- math.pi*(r1+nPyrex+nAu)**2
SCoNi_max = math.pi*(r1_max+nPyrex_max+nAu_max+ nCoNi_max)**2- math.pi*(r1_max+nPyrex_max+nAu_max)**2
SCoNi_min = math.pi*(r1_min+nPyrex_min+nAu_min+ nCoNi_min)**2- math.pi*(r1_min+nPyrex_min+nAu_min)**2


print ('Cross section areas')
print ('SFeSiB = ', SFeSiB, average_error(SFeSiB_max, SFeSiB_min),'[m^2]')
print ('SAu = ', SAu, average_error(SAu_max, SAu_min), '[m^2]')
print ('SCoNi = ', SCoNi,average_error(SCoNi_max, SCoNi_min),'[m^2]')
print ('SPyrex = ', SPyrex, average_error(SPyrex_max, SPyrex_min), '[m^2]')


# density of layers

DFeSiB = 7400.0
DAu = 19300.0
DCoNi = 8900.0
l = 0.01

# mass of layers
mFeSiB = DFeSiB*SFeSiB*l
mFeSiB_min = DFeSiB*SFeSiB_min*l
mFeSiB_max = DFeSiB*SFeSiB_max*l

mAu    = DAu      *SAu      *l
mAu_min= DAu      *SAu_min      *l
mAu_max= DAu      *SAu_max      *l

mCoNi     = DCoNi  *SCoNi  *l
mCoNi_min = DCoNi  *SCoNi_min  *l
mCoNi_max = DCoNi  *SCoNi_max  *l


print ('Mass of layers')
print ('mFeSiB =', mFeSiB, average_error(mFeSiB_max, mFeSiB_min), '[kg.m^-3]' )
print ('mAu = ',   mAu, average_error(mAu_max, mAu_min), '[kg.m^-3]')
print ('mCoNi = ', mCoNi, average_error(mCoNi_max, mCoNi_min),'[kg.m^-3]')
####################################################
####################################################

#volume
VFeSiB = SFeSiB*l
VFeSiB_max = SFeSiB_max*l
VFeSiB_min = SFeSiB_min*l

VPyrex = SPyrex*l
VPyrex_max = SPyrex_max*l
VPyrex_min = SPyrex_min*l

VAu = SAu*l
VAu_max = SAu_max*l
VAu_min = SAu_min*l

VCoNi = SCoNi*l
VCoNi_max = SCoNi_max*l
VCoNi_min = SCoNi_min*l

Vtot = VFeSiB + VPyrex + VAu + VCoNi
Vtot_max = VFeSiB_max + VPyrex_max + VAu_max + VCoNi_max
Vtot_min = VFeSiB_min + VPyrex_min + VAu_min + VCoNi_min

print ('Vtot =', Vtot, average_error(Vtot_max, Vtot_min), '[m^3]')

#fraction of layers in percent
FeSiB_frac = VFeSiB*100/Vtot
FeSiB_frac_max = VFeSiB_max*100/Vtot
FeSiB_frac_min = VFeSiB_min*100/Vtot

Pyrex_frac = VPyrex*100/Vtot
Pyrex_frac_max = VPyrex_max*100/Vtot
Pyrex_frac_min = VPyrex_min*100/Vtot

Au_frac = VAu*100/Vtot
Au_frac_max = VAu_max*100/Vtot
Au_frac_min = VAu_min*100/Vtot

CoNi_frac = VCoNi*100/Vtot
CoNi_frac_max = VCoNi_max*100/Vtot
CoNi_frac_min = VCoNi_min*100/Vtot
print ('FeSiB_frac = ', FeSiB_frac/100)
print ('Pyrex_frac = ', Pyrex_frac/100)
print ('Au_frac = ', Au_frac/100)
print ('CoNi_frac =', CoNi_frac/100)
print ('tot_frac = ', FeSiB_frac/100+Pyrex_frac/100+Au_frac/100+CoNi_frac/100)
####################################################
####################################################

#resistance
roCoNi  = 0.10*pow(10, -6)
roFeSiB = 1.14*pow(10, -6)
roFeSiB_max = 1.14*pow(10, -6)
roFeSiB_min = 1.14*pow(10, -6)
roAu    = 2.44*pow(10, -8)

RCoNi = roCoNi*l/SCoNi
RCoNi_min = roCoNi*l/SCoNi_max
RCoNi_max = roCoNi*l/SCoNi_min

RFeSiB = roFeSiB*l/SFeSiB
RFeSiB_min = roFeSiB_min*l/SFeSiB_max
RFeSiB_max = roFeSiB_max*l/SFeSiB_min

RAu = roAu*l/SAu
RAu_min = roAu*l/SAu_max
RAu_max = roAu*l/SAu_min

Rtot = 1.0/(1.0/RFeSiB + 1.0/RAu + 1.0/RCoNi)
Rtot_min = 1.0/(1.0/RFeSiB_min + 1.0/RAu_min + 1.0/RCoNi_min)
Rtot_max = 1.0/(1.0/RFeSiB_max + 1.0/RAu_max + 1.0/RCoNi_max)

print ('Rtot =', Rtot)
RtotA_min = Rtot - 0.2
RtotA_max = Rtot + 0.2


while 1.0/(1.0/RFeSiB_min + 1.0/RAu_min + 1.0/RCoNi_min) < RtotA_min:
    RAu_min = RAu_min + RAu_min*0.01204
    RCoNi_min = RCoNi_min + RCoNi_min*0.01013
    RFeSiB_min = RFeSiB_min + RFeSiB_min*0.00570

while 1.0/(1.0/RFeSiB_max + 1.0/RAu_max + 1.0/RCoNi_max) > RtotA_max:
    RAu_max = RAu_max - RAu_max*0.01204
    RCoNi_max = RCoNi_max- RCoNi_max*0.01013
    RFeSiB_max = RFeSiB_max - RFeSiB_max*0.00570

print ('resistance')
print ('RFeSiB = ', RFeSiB, average_error(RFeSiB_max, RFeSiB_min), 'Ohm')
print ('RAu = ', RAu, average_error(RAu_max, RAu_min), 'Ohm')
print ('RCoNi = ', RCoNi, average_error(RCoNi_max, RCoNi_min), 'Ohm')


I = 0.005 
print ('current')

ICoNi = I/(RCoNi/RFeSiB + RCoNi/RAu + 1.0)
ICoNi_min = I/(RCoNi_max/RFeSiB_min + RCoNi_max/RAu_min + 1.0)
ICoNi_max = I/(RCoNi_min/RFeSiB_max + RCoNi_min/RAu_max + 1)

IAu = RCoNi*ICoNi/RAu
IAu_min = RCoNi_min*ICoNi_min/RAu_max
IAu_max = RCoNi_max*ICoNi_max/RAu_min

IFeSiB = RCoNi*ICoNi/RFeSiB
IFeSiB_min = RCoNi_min*ICoNi_min/RFeSiB_max
IFeSiB_max = RCoNi_max*ICoNi_max/RFeSiB_min
ro=I/SFeSiB
print ('ro =', ro)

print ('IFeSiB = ', IFeSiB, average_error(IFeSiB_max, IFeSiB_min), '[A]')
print ('IAu = ', IAu, average_error(IAu_max, IAu_min), '[A]')
print ('ICoNi = ', ICoNi, average_error(ICoNi_max, ICoNi_min), '[A]')


QFeSiB = IFeSiB**2*RFeSiB
QFeSiB_min = IFeSiB_min**2*RFeSiB_min
QFeSiB_max = IFeSiB_max**2*RFeSiB_max

QAu = IAu**2*RAu
QAu_min = IAu_min**2*RAu_min
QAu_max = IAu_max**2*RAu_max

QCoNi = ICoNi**2*RCoNi
QCoNi_min = ICoNi_min**2*RCoNi_min
QCoNi_max = ICoNi_max**2*RCoNi_max

print ('heat')
print ('QFeSiB = ', QFeSiB, average_error(QFeSiB_max, QFeSiB_min), '[J]')
print ('QAu = ', QAu, average_error(QAu_max, QAu_min), '[J]')
print ('QCoNi = ', QCoNi, average_error(QCoNi_max, QCoNi_min), '[J]')

cCoNi = 423.3 #J.kg-1.K-1
cAu = 129
cFeSiB = 2200
cFeSiB_max = 2220
cFeSiB_min = 2180

deltaTPyrex = 0.0

##deltaTFeSiB = QFeSiB/(mFeSiB*cFeSiB)
#deltaTFeSiB_min = QFeSiB_min/(mFeSiB_max*cFeSiB_max)
#deltaTFeSiB_max = QFeSiB_max/(mFeSiB_min*cFeSiB_min)

##deltaTAu = QAu/(mAu*cAu)
##deltaTAu_min = QAu_min/(mAu_max*cAu)
##deltaTAu_max = QAu_max/(mAu_min*cAu)

##deltaTCoNi = QCoNi/(mCoNi*cCoNi)
##deltaTCoNi_min = QCoNi_min/(mCoNi_max*cCoNi)
##deltaTCoNi_max = QCoNi_max/(mCoNi_min*cCoNi)

deltaTFeSiB = 7.9472477282054941
deltaTFeSiB_min = 7.9472477282054941
deltaTFeSiB_max = 7.9472477282054941

deltaTAu = 10.467156725238567
deltaTAu_min = 10.467156725238567
deltaTAu_max = 10.467156725238567

deltaTCoNi = 8.0987191321449821
deltaTCoNi_min = 8.0987191321449821
deltaTCoNi_max = 8.0987191321449821

print ('temperature')
print ('deltaTFeSiB =', deltaTFeSiB, average_error(deltaTFeSiB_max, deltaTFeSiB_min), '[K]')
print ('deltaTAu =', deltaTAu, average_error(deltaTAu_max, deltaTAu_min), '[K]')
print ('deltaTCoNi =', deltaTCoNi, average_error(deltaTCoNi_max, deltaTCoNi_min), '[K]')

####################################################
####################################################

deltaTPyrex_max = 0
deltaTPyrex_min = 0
deltaT = (deltaTFeSiB*FeSiB_frac/100 + deltaTPyrex*Pyrex_frac/100 + deltaTAu*Au_frac/100 + deltaTCoNi*CoNi_frac/100)/4

deltaT_max = (deltaTFeSiB_max*FeSiB_frac_min/100 + deltaTPyrex_max*Pyrex_frac_min/100 + deltaTAu_max*Au_frac_min/100 +
              deltaTCoNi_max*CoNi_frac_min/100)/4
deltaT_min = (deltaTFeSiB_min*FeSiB_frac_max/100 + deltaTPyrex_min*Pyrex_frac_max/100 + deltaTAu_min*Au_frac_max/100 +
              deltaTCoNi_min*CoNi_frac_max/100)/4

print ('deltaT = ', deltaT, '[K] in 1 s')
print ('deltaT = ', deltaT, average_error(deltaT_max, deltaT_min), '[K]')

alphaFeSiB = 3.3*pow(10,-6)
alphaAu = 14.1*pow(10, -6)
alphaCoNi = 13*pow(10, -6)
alphaPyrex = 3.2*pow(10, -6)

EFeSiB = 160*pow(10, 9)
EPyrex = 64 *pow(10, 9)
EAu    = 80 *pow(10, 9)
ECoNi  = 140*pow(10, 9)

Es = (EFeSiB*FeSiB_frac_min + EPyrex*Pyrex_frac_min + EAu*Au_frac_min + ECoNi*CoNi_frac_min)/4
print ('Es =', Es)


a = [0.0, alphaFeSiB, alphaPyrex, alphaAu, alphaCoNi]
v = [0.0, VFeSiB_min, VPyrex_min, VAu_min, VCoNi_min]
E = [0.0, EFeSiB, EPyrex, EAu, ECoNi]

alphas = (a[1]*v[1]*E[1] + a[2]*v[2]*E[2] + a[3]*v[3]*E[3] + a[4]*v[4]*E[4])/(a[1]*v[1] + a[2]*v[2]+ a[3]*v[3] + a[4]*v[4])

##############################################################################

import numpy as np
import pylab
import matplotlib.pyplot as plt


c13 = ECoNi*(nCoNi_max-nCoNi_min)*abs(alphaCoNi-alphaFeSiB)*deltaT/(EFeSiB*(2*r1_min)**2)
c13_max = ECoNi*(nCoNi_max-nCoNi_min)*abs(alphaCoNi-alphaFeSiB)*deltaT_max/(EFeSiB*(2*r1_min)**2)
c13_min = ECoNi*(nCoNi_max-nCoNi_min)*abs(alphaCoNi-alphaFeSiB)*deltaT_min/(EFeSiB*(2*r1_min)**2)
c23 = EAu*(nAu_max-nAu_min)*abs(alphaAu-alphaFeSiB)*deltaT/(EFeSiB*(2*r1_min)**2)
c23_max = EAu*(nAu_max-nAu_min)*abs(alphaAu-alphaFeSiB)*deltaT_max/(EFeSiB*(2*r1_min)**2)
c23_min = EAu*(nAu_max-nAu_min)*abs(alphaAu-alphaFeSiB)*deltaT_min/(EFeSiB*(2*r1_min)**2)
c33 = EPyrex*(nPyrex_max-nPyrex_min)*abs(alphaPyrex-alphaFeSiB)*deltaT/(EFeSiB*(2*r1_min)**2)
c33_max = EPyrex*(nPyrex_max-nPyrex_min)*abs(alphaPyrex-alphaFeSiB)*deltaT_max/(EFeSiB*(2*r1_min)**2)
c33_min = EPyrex*(nPyrex_max-nPyrex_min)*abs(alphaPyrex-alphaFeSiB)*deltaT_min/(EFeSiB*(2*r1_min)**2)

radius_a4 = 1/(6*(c13 + c23 + c33))
radius_a4_max = 1/(6*(c13_min + c23_min + c33_min))
radius_a4_min = 1/(6*(c13_max + c23_max + c33_max))

print ('radius_a4=', radius_a4)
print ('l/2=', l/2)
print ('radius_a4**2-(l/2)**2 =', radius_a4**2-(l/2)**2)
delta_a4 = radius_a4 - (radius_a4**2-(l/2)**2)**0.5
delta_a4_min = radius_a4_max - (radius_a4_max**2-(l/2)**2)**0.5
delta_a4_max = radius_a4_min - (radius_a4_min**2-(l/2)**2)**0.5

print ('delta_a4 =', delta_a4, average_error(delta_a4_max, delta_a4_min))
print ('deltaT =', deltaT)
x = [deltaT]
y = [delta_a4*1e6]

yerr = [0.0]

nn = 0 #number of runs of loop

while deltaT < 3.2:
    deltaT = deltaT + 0.1
    c13 = ECoNi*(nCoNi_max-nCoNi_min)*(alphaCoNi-alphaFeSiB)*deltaT/(EFeSiB*(2*r1_min)**2)
    c13_max = ECoNi*(nCoNi_max-nCoNi_min)*(alphaCoNi-alphaFeSiB)*deltaT_max/(EFeSiB*(2*r1_min)**2)
    c13_min = ECoNi*(nCoNi_max-nCoNi_min)*(alphaCoNi-alphaFeSiB)*deltaT_min/(EFeSiB*(2*r1_min)**2)
    c23 = EAu*(nAu_max-nAu_min)*(alphaAu-alphaFeSiB)*deltaT/(EFeSiB*(2*r1_min)**2)
    c23_max = EAu*(nAu_max-nAu_min)*(alphaAu-alphaFeSiB)*deltaT_max/(EFeSiB*(2*r1_min)**2)
    c23_min = EAu*(nAu_max-nAu_min)*(alphaAu-alphaFeSiB)*deltaT_min/(EFeSiB*(2*r1_min)**2)
    c33 = EPyrex*(nPyrex_max-nPyrex_min)*(alphaPyrex-alphaFeSiB)*deltaT/(EFeSiB*(2*r1_min)**2)
    c33_max = EPyrex*(nPyrex_max-nPyrex_min)*(alphaPyrex-alphaFeSiB)*deltaT_max/(EFeSiB*(2*r1_min)**2)
    c33_min = EPyrex*(nPyrex_max-nPyrex_min)*(alphaPyrex-alphaFeSiB)*deltaT_min/(EFeSiB*(2*r1_min)**2)

    radius_a4 = 1/(6*(c13 + c23 + c33))
    radius_a4_max = 1/(6*(c13_min + c23_min + c33_min))
    radius_a4_min = 1/(6*(c13_max + c23_max + c33_max))
    delta_a4 = radius_a4 - (radius_a4**2-(l/2)**2)**0.5
    delta_a4_min = radius_a4_max - (radius_a4_max**2-(l/2)**2)**0.5
    delta_a4_max = radius_a4_min - (radius_a4_min**2-(l/2)**2)**0.5
    y.append(delta_a4*1e6)#-y
    nn = nn + 1
    x.append(deltaT)
    yerr.append((delta_a4_max - delta_a4_min)*1e6)
print ('len(x)=', len(x))
  
##########################
kk = 3.5
k_max = 4.0
k_min = 3.0

A = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9]
T1 = A[0]/kk
T = [T1]
for ss in range(1,19,1):
    a1 =A[ss]/kk
    a1_max = A[ss]/k_min
    a1_min = A[ss]/k_max
    T.append(a1)
print ('average_error =', average_error(a1_max, a1_min))


###############################
'''Magnetostriction'''
##h = np.linspace(2.84885714e-07,6.1e-7,26, True)
h = np.linspace(8.985714e-08,2.471408e-07,26, True)
hx = np.linspace(0.57142857,3.202653658056983,26, True)
print ('len(h)=', len(h), 'len(y)=', len(y))
'''Magnetostriction + Thermal expansion = MT'''
MT = []
for i in range(len(h)):
    hy = h[i] + y[i]
    MT.append(hy)
print ('MT=', MT)    
print ('delta_a4 = ', delta_a4)
print ('deltaT = ',deltaT)


plt.rc('font', size=35)
plt.figure(figsize=(10, 10), dpi=60)
pylab.ylabel (r' A [$\mu$m]', fontsize = 40)
pylab.xlabel(r'$\Delta$ T [K]',fontsize = 40)

##plt.plot(x,MT,'rs',markersize=26, label='calc.')#-y
plt.plot(T,A,'ks', markersize = 15)#, label='exp')
plt.xticks([0,1,2,3])
plt.yticks([2,4,6,8,10])
ax = plt.gca() # get the current axes
for l in ax.get_xticklines() + ax.get_yticklines():
    l.set_markersize(10)
    l.set_markeredgewidth(2) 

##plt.legend(loc =2)

pylab.savefig('C:/Users/kolesar/Documents/cv/51_Schönaich bei Böblingen/01_nanuscript/Figure 10.png')
##plt.show()

