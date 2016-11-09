# -*- coding: utf-8 -*-

"""
Created 2014

@author: Lidia Bressan

DOP...


Vreal = Vmis /cos(ang_dop)

ang_dop = 45° ==> cos(ang_dop) = 1.0/sqrt(2.0)

  ==>  Vreal = Vmis /cos(ang_dop) = Vmis*sqrt(2.0)


fd = 2 fe v cos_a /c
fe - frequenza di emissione
c - velocita' del suono
fd - frequenza doppler?

Tprf - periodo di pulse repetition frequency
Fprf - periodo di pulse repetition frequency

Pmax = Tprf*c/2.0 = c/2.0/Fprf
Pmax*Vmax = c**2/8.0/fe

Vmax = Fprf*c/4.0/fe

fe = c**2/8.0/Pmax/Vmax

"""
import numpy as np

c = 1500. ## m/s

fe = np.array([1.0, 4.0])*10**6 ## MHz
##Fprf = np.array([100., 15625.]) ## MHz
Fprf = np.array([5000., 10000., 13000.]) ## MHz
Tprf = 1/Fprf


Pmax_x_Vmax = c**2/8.0/fe

def p_max(c, Fprf):
    Pmax = c/2/Fprf
    return Pmax

Pmax = p_max(c, Fprf)

#print("P max = {0} cm".format(Pmax*100))

Vmax = np.array([Fprf*c/4.0/fe[0], Fprf*c/4.0/fe[1]])

#print Vmax


def prf(Vmax, fe):
    return 4.0*fe*Vmax/c

for Vmax in [1.0, 2.0, 2.5, 5.0]: #m/s
    print ('Vmax, Fprf, fe', Vmax, prf(Vmax, fe), fe)



#Fprf, fe = 15000, 4*10**6
#print ([Fprf, fe, p_max(c, Fprf)])

#Fprf, fe = 13000, 4*10**6
#print ([Fprf, fe, p_max(c, Fprf)])

#Fprf, fe = 5000, 1*10**6
#print ([Fprf, fe, p_max(c, Fprf)])


