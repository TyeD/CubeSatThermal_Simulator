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
Boltz = 5.67e-8 #W/m^2 K^4

print("Orbital Parameters:") 
print("------------------------------------------------------\n")
#https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.ht

#Earth Parameters
#==========================#
massE = 5.972e24 #kg
d_sun = 149.6e9 #m
blk_body_temp_earth = 254.0 #K
solar_flux = 1361 # W/m^2 
solar_emission = 237
albedo_avg = .3
r_Earth = 6361 #km
Umbra = 2160 #s
Penumbra = 8 #s 
#==========================#
#Satellite Dimension










#======================================#
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

#=====================================================#
#Asorbitivity, 6061 Aluminium 
absorb = 0.031
#Emissivity, 6061 Aluminium
emit = 0.039
#=====================================================#
#Simulation Time Step
dt = 1000#Time steps 
t = np.linspace(0, P, dt)
#======================================================#
#View Factor Cacluation
phi = np.abs(np.pi - t)
Phi = t 
#NOTE: Phi is currently unknown need to ask ADCS
# about if they have any goals for the 
#attitude of the ship (orententation w.r.t. some frame) 
#over the orbit.


#Trigonometric relations for 
k = r_Earth/(r_Earth+h)
#Alpha is the angle the sat will fly past in single cycle.
alpha = np.linspace(0, 2*np.pi, dt)



#Generate empty lists for each if statment
Xe1 = []
Xe2 = []
Xe3 = []
Xe4 = []
Xe5 = []
Xe = []
#Loop through all alpha values for each condition. 
for i in range(len(alpha)):
    #Loop through conditions given in paper 
    if 0 <= alpha[i] <= np.arccos(k):
        Xe1.append((k**2)*np.cos(alpha[i]))

    elif np.arccos(k) < alpha[i] < (np.pi - np.arccos(k)):
        Xe2.append((k**2)*np.cos(alpha[i]) + (1/np.pi)*((np.pi/2)-((1-k**2)*(k**2 - np.cos(alpha[i])**2))**(1/2)-np.arcsin(((1-k**2))**(1/2)/np.sin(alpha[i]))-((k**2)*np.cos(alpha[i])*np.arccos(((1/k**2)-1)**(1/2)*(np.cos(alpha[i]))/np.sin(alpha[i])))))
    
    elif (np.pi - np.arccos(k)) <= alpha[i] <= (np.pi + np.arccos(k)):
        Xe3.append(0*alpha[i])
    
    elif (np.pi + np.arccos(k)) < alpha[i] < (2*np.pi - np.arccos(k)):
        Xe4.append((k**2)*np.cos(2*np.pi - alpha[i]) + (1/np.pi)*((np.pi/2)-((1-k**2)*(k**2 - np.cos(2*np.pi - alpha[i])**2))**(1/2) - (np.arcsin((((1-k**2))**(1/2))/np.sin(2*np.pi - alpha[i]))) - ((k**2)*np.cos(2*np.pi - alpha[i])*np.arccos(((1/((k**2))-1)**(1/2))*(np.cos(2*np.pi - alpha[i]))/np.sin(2*np.pi - alpha[i])))))
    
    elif (2*np.pi - np.arccos(k)) <= alpha[i] <= (2*np.pi):
        Xe5.append((k**2)*np.cos(alpha[i]))
#Append Values to Xe
Xe = np.append(Xe, Xe1)
Xe = np.append(Xe, Xe2)
Xe = np.append(Xe, Xe3)
Xe = np.append(Xe, Xe4)
Xe = np.append(Xe, Xe5)


Xa = Xe*np.cos(phi)

qe = solar_emission*absorb*Xe
qa = albedo_avg*solar_flux*absorb*Xa
plt.plot(t, Xe)
plt.show()


plt.plot(t,Xe, 'y')
plt.xlabel("Time (s)")
plt.ylabel(r"Sun View Factor $X_e$")
plt.title("Veiw Factor as function of time over one period of a LEO")
plt.show()


plt.plot(t, qe, 'r')
plt.show()

plt.plot(t, Xa, 'k')
plt.show()
