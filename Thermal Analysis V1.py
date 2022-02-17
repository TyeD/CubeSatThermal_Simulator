# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 14:50:47 2022

@author: Tye Dougherty 
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
#Orbital parameters:
#=================================#
#Constants 
G = 6.67430e-11


print("Orbital Parameters:") 
print("------------------------------------------------------\n")
#https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.ht

#Earth Parameters
massE = 5.972e24 #kg
d_sun = 149.6e9
blk_body_temp = 254.0 #K
solar_flux = 1361 # W/m^2 
albedo = .3
solar_emission = 237
r_Earth = 6361 #km
Umbra = 2160 #s
Penumbra = 8 #s 

#Orbital Altitude from Earth's surface: 
#LEO (Low Earth Orbit)
h = 300 #km

#Orbital Radius 
r_orb1 = r_Earth + h #km
r_orb = r_orb1*10**3 #m
print("Orbital Radius is: %.0f km \n" %(r_orb1))

#Flight Velocity 
c1 = 398600.5 #Constant 1 
c2 = 6378.14 #Constant 2
orb_vel = np.sqrt(G*massE/r_orb) #m/s
print("Orbital Velocity is: %.2f m/s\n" % (orb_vel))

#Orbital Period 
P = np.sqrt((4*(np.pi**2)*r_orb**3)/(G*massE))
print("Orbital Period: %.1f s" %(P))

#Asorbitivity, assume 6061 Aluminium 
alpha_lower = .02
alpha_upper = .2  
#=====================================================#

dt = 10000#Time steps 
t = np.linspace(0, P, dt)
phi = np.abs(np.pi - t)
#print(phi)
#print(phi)
#print(t)

k = r_Earth/(r_Earth+h)
beta = np.linspace(0, 2*np.pi, dt)

beta1 = t 
beta2 = 2*np.pi - t
#print(beta1, beta2)

'''
print(np.arccos(k))
print((np.pi + np.arccos(k)))
print(np.pi)
print(beta)
'''
print((np.pi + np.arccos(k)))
Xe1 = []
Xe2 = []
Xe3 = []
Xe4 = []
Xe5 = []
Xe = []
for i in range(len(beta)):
    #print(i)
    #print(beta[i])
    if 0 <= beta[i] <= np.arccos(k):
        #print(beta[i])
        Xe1.append((k**2)*np.cos(beta[i]))
        #print(Xe1)
    
    elif np.arccos(k) < beta[i] < (np.pi - np.arccos(k)):
        #print(beta[i])
        Xe2.append((k**2)*np.cos(beta[i]) + (1/np.pi)*((np.pi/2)-((1-k**2)*(k**2 - np.cos(beta[i])**2))**(1/2)-np.arcsin(((1-k**2))**(1/2)/np.sin(beta[i]))-((k**2)*np.cos(beta[i])*np.arccos(((1/k**2)-1)**(1/2)*(np.cos(beta[i]))/np.sin(beta[i])))))
        #print(Xe2)
    
    elif (np.pi - np.arccos(k)) <= beta[i] <= (np.pi + np.arccos(k)):
        #print(beta[i])
        Xe3.append(0*beta[i])
    
    elif (np.pi + np.arccos(k)) < beta[i] < (2*np.pi - np.arccos(k)):
        Xe4.append((k**2)*np.cos(beta[i]) + (1/np.pi)*((np.pi/2)-((1-k**2)*(k**2 - np.cos(beta[i])**2))**(1/2)-np.arcsin(((1-k**2))**(1/2)/np.sin(beta[i]))-(k**2)*np.cos(beta[i])*np.arccos((1/((k**2))-1)**(1/2)*(np.cos(beta[i]))/np.sin(beta[i]))))
    
    elif (2*np.pi - np.arccos(k)) <= beta[i] <= (2*np.pi):
        #print(beta[i])
        Xe5.append((k**2)*np.cos(beta[i]))
'''
#print(Xe1)
print(Xe2)
print(Xe3)
print(Xe4)
print(Xe5)
'''
Xe = np.append(Xe, Xe1)
Xe = np.append(Xe, Xe2)
Xe = np.append(Xe, Xe3)
Xe = np.append(Xe, Xe4)
Xe = np.append(Xe, Xe5)
#print(Xe2)
#print("Hi", Xe)

#Xa = Xe*np.cos(phi)

#qe = albedo*solar_emission*alpha_lower*Xe
plt.plot(t, Xe)
#plt.plot(t,Xa)