# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 14:50:47 2022

@author: Tye Dougherty 
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
'''
Write function to pass varity of emissivity and absorbtivity constants in and
plot the output on the same graph. This will require some research into the statistical
approch that should be made.

'''
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
earth_emission = 237
albedo_avg = .3
r_Earth = 6361 #km
Umbra = 2160 #s
Penumbra = 8 #s 
#==========================#
#Satellite Dimension DELETE 
#As faces exposed to sun Ae is faces exposed to earth A1f is a single face At is the total area (Note this should be improved)
A1f = 22919 #mm 
A1f = A1f/1e6                         #Play with this parameter
At = 4*A1f                            #Play with this parameter
As = 2 *A1f 
Ae = A1f                              #Play with this parameter

#======================================#

#=====================================================#
#Asorbitivity
absorb1 = .14          #Play with this parameter
#Emissivity
emitAl = 0.22      #Play with this parameter 
'''
Usally spacecraft remain within 126 degrees but this depends on electrical
system temperature operatating requirements. A material with a emissivity of
greater then .3 would keep the space craft below 110 degrees C. 
'''
#=====================================================#
#Simulation Time Step
dt = 1000#Time steps                   #Play with this parameter
##Orbital Altitude from Earth's surface: 
#LEO (Low Earth Orbit)
h = 300 #km
def runSimulation(G, Boltz, massE, r_Earth, albedo_avg, earth_emission, solar_flux, A1f, absorb, emitAl, dt):
    #Spacecraft dimensions again 
    A1f = A1f/1e6                         #Play with this parameter
    At = 4*A1f                            #Play with this parameter
    As = 2 *A1f 
    Ae = A1f 
    #=====================================================#
    print("Orbit height above Earth's surface': " + str(h)+ " km \n")
    
    #Orbital Radius 
    r_orb1 = r_Earth + h #km
    r_orb = r_orb1*10**3 #m
    #print("Orbital Radius is: %.0f m \n" %(r_orb1))
    
    #Flight Velocity 
    orb_vel = np.sqrt(G*massE/r_orb) #m/s
    #print("Orbital Velocity is: %.2f m/s\n" % (orb_vel))
    
    #Orbital Period 
    P = np.sqrt((4*(np.pi**2)*r_orb**3)/(G*massE))
    #print("Orbital Period: %.1f s" %(P))
    
    #=====================================================#
    t = np.linspace(0, P, dt)
    #======================================================#
    #View Factor Cacluation
    
    #Phi = t 
    #NOTE: Phi is currently unknown need to ask ADCS
    # about if they have any goals for the 
    #attitude of the ship (orententation w.r.t. some frame) 
    #over the orbit.
    
    #Trigonometric relation 
    k = r_Earth/(r_Earth+h)
    #Alpha is the angle the sat will fly past in single cycle.
    alpha = np.linspace(0, 2*np.pi, dt)
    phi = np.abs(np.pi - alpha)
    
    
    ##Cacluate sun emitted view factor
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
            Xe2.append((k**2)*np.cos(alpha[i]) + 
            (1/np.pi)*((np.pi/2)-((1-k**2)*(k**2 - np.cos(alpha[i])**2))**(1/2)-
            np.arcsin(((1-k**2))**(1/2)/np.sin(alpha[i]))-
            ((k**2)*np.cos(alpha[i])*np.arccos(((1/k**2)-1)**(1/2)*
            (np.cos(alpha[i]))/np.sin(alpha[i])))))
        
        elif (np.pi - np.arccos(k)) <= alpha[i] <= (np.pi + np.arccos(k)):
            Xe3.append(0*alpha[i])
        
        elif (np.pi + np.arccos(k)) < alpha[i] < (2*np.pi - np.arccos(k)):
            Xe4.append((k**2)*np.cos(2*np.pi - alpha[i]) + 
            (1/np.pi)*((np.pi/2)-((1-k**2)*(k**2 - np.cos(2*np.pi - alpha[i])**2))**(1/2) - 
            (np.arcsin((((1-k**2))**(1/2))/np.sin(2*np.pi - alpha[i]))) - 
            ((k**2)*np.cos(2*np.pi - alpha[i])*np.arccos(((1/((k**2))-1)**(1/2))*(np.cos(2*np.pi - 
            alpha[i]))/np.sin(2*np.pi - alpha[i])))))
        
        elif (2*np.pi - np.arccos(k)) <= alpha[i] <= (2*np.pi):
            Xe5.append((k**2)*np.cos(alpha[i]))
    
    #Append Values to Xe
    Xe = np.append(Xe, Xe1)
    Xe = np.append(Xe, Xe2)
    Xe = np.append(Xe, Xe3)
    Xe = np.append(Xe, Xe4)
    Xe = np.append(Xe, Xe5)
    #=============================================================================#
    #Albedo View Factor
    #(Note may not be required, may just be constant, depends on ADCS)
    #====================================================#
    #Required to make albedo view factor zero 
    Xa1 = Xe*np.cos(phi)
    Xa = []
    for i in range(len(t)):
        if Xa1[i] > 0:    
            Xa.append(Xe[i]*np.cos(phi[i]))
    
        elif Xa1[i] <= 0:
    
            Xa.append(0*phi[i])
    
    Xa = np.array(Xa)
    
    
    #Heat flux calcluations
    #==================================================#
    #Xe = np.array(Xe, dtype = np.int32)
    #print(np.dtype(Xe[1]))
    #print(np.dtype(absorb))
    #print(np.shape(Xe))
    qs = solar_flux*.14*Xe
    q_earth_emitted = earth_emission*absorb
    q_earth_reflected = albedo_avg*solar_flux*absorb*Xa
    #==================================================#
    
    
    #Electronics Power and survival temps 
    #=================================================#
    Qelec = 3*1.9 #W Base assmption more detailed model can be applied below in space below 1U cube sat = 1.9W 
    Tmax = 40; Tmin = 10
    #
    #
    #
    #=================================================#
    
    
    #Temperature Function (Heat balence)
    #https://www.alternatewars.com/BBOW/Space/Spacecraft_Ext_Temps.htm
    #=========================================#
    Ts = (((((qs*As*.3))+
            (q_earth_reflected*Ae)+Qelec)/(emitAl*Boltz*At*.75))**(1/4)) - 273.13 #Degree C
    #=========================================#
    
    
    #Produce Plots
    #=======================================================================#
    plt.plot(t,Xe, 'y', label = r'Sun emitted radiation $X_{emitted}$')
    plt.plot(t, Xa, 'k', label = r'Earth reflected radiation $X_{albedo}$')
    plt.xlabel("Time (s)")
    plt.ylabel(r"View Factor")
    plt.title("View Factor as function of time over one period of a LEO")
    plt.legend()
    plt.show()
    
    plt.plot(t, qs, 'b')
    plt.xlabel("Time (s)")
    plt.ylabel(r'Heat flux $(W/m^2)$')
    plt.title("Sun's heat flux as function of time over one period of a LEO")
    plt.show()
    
    '''
    plt.plot(t, q_earth_emitted , 'r')
    plt.xlabel("Time (s)")
    plt.ylabel(r'Heat flux $(W/m^2)$')
    plt.title("Earths's emitted heat flux as function of time over one full orbit")
    plt.show()
    '''
    '''
    plt.plot(t, q_earth_reflected, 'b')
    plt.xlabel("Time (s)")
    plt.ylabel(r'Heat flux $(W/m^2)$')
    plt.title("Earths's reflected heat flux as function of time over one full orbit")
    plt.show()
    '''
    
    plt.plot(t, q_earth_reflected+q_earth_emitted+qs, 'r')
    plt.xlabel("Time (s)")
    plt.ylabel(r'Heat flux $(W/m^2)$')
    plt.title("Total heat flux as function of time over one full orbit")
    plt.show()
    
    plt.plot(t, Ts)
    plt.plot(t, Tmax*(t/t), 'r--', label = r'$T_{max}$ and $T_{min}$ for Electronics')
    plt.plot(t, Tmin*(t/t), 'r--')
    plt.xlabel("Time (s)")
    plt.ylabel(r'Temperature $(^\circ C)$')
    plt.title("Temperature as function of time over one full orbit")
    plt.legend(loc = 'best')
    plt.show()
    return r_orb, P, orb_vel
#=======================================================================#

runSimulation(G, Boltz, massE, r_Earth, albedo_avg, earth_emission, solar_flux, A1f, absorb1, emitAl, dt)
