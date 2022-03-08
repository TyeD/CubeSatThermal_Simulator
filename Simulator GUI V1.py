# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 18:35:55 2022

@author: abdla
"""
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
#Satellite Dimension
#As faces exposed to sun Ae is faces exposed to earth A1f is a single face At is the total area (Note this should be improved)
#A1f = 22919 #mm 
#A1f = A1f/1e6                         #Play with this parameter
#At = 4*A1f                            #Play with this parameter
#As = 2 *A1f 
#Ae = A1f                              #Play with this parameter

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



import tkinter as tk 
from tkinter import ttk 
from ThermalAnalysisV1 import runSimulation
class ToggledFrame(tk.Frame):

    def __init__(self, parent, text="", *args, **options):
        tk.Frame.__init__(self, parent, *args, **options)

        self.show = tk.IntVar()
        self.show.set(0)

        self.title_frame = ttk.Frame(self)
        self.title_frame.pack(fill="x", expand=1)

        ttk.Label(self.title_frame, text=text).pack(side="left", fill="x", expand=1)

        self.toggle_button = ttk.Checkbutton(self.title_frame, width=2, text='+', command=self.toggle,
                                            variable=self.show, style='Toolbutton')
        self.toggle_button.pack(side="left")

        self.sub_frame = tk.Frame(self, relief="sunken", borderwidth=1)

    def toggle(self):
        if bool(self.show.get()):
            self.sub_frame.pack(fill="x", expand=1)
            self.toggle_button.configure(text='-')
        else:
            self.sub_frame.forget()
            self.toggle_button.configure(text='+')

if __name__ == "__main__":
    #Create gui Window
    gui = tk.Tk()
    # set the background colour of GUI window
    gui.configure(background="grey")
    
    # set the title of GUI window
    gui.title("3U CubeSat Thermal Simulator")
    
    
    userinputFrame = tk.Frame(master=gui, width=50, height=50)

    t = ToggledFrame(userinputFrame, text='Orbital Parameters', relief="raised", borderwidth=1)
    h_var = tk.IntVar()
    t.pack(fill="x", expand=1, pady=2, padx=2, anchor="n")

    ttk.Label(t.sub_frame, text= "Orbit Height (km):").pack(side="left", fill="x", expand=1)
    ttk.Entry(t.sub_frame).pack(side="left")
    h = ttk.Entry(t.sub_frame).get()
    
    

    
    t2 = ToggledFrame(userinputFrame, text='Thermal Matieral Parameters', relief="raised", borderwidth=1)
    t2.pack(fill="x", expand=1, pady=2, padx=2, anchor="n")
    
    ttk.Label(t2.sub_frame, text= "Absorbitivity:").pack(side="left", expand=1)
    ttk.Entry(t2.sub_frame).pack(side="left")
    absorb = ttk.Entry(t2.sub_frame).get()
    
    
    ttk.Label(t2.sub_frame, text= "Emissivity").pack(side="left", expand=1)
    ttk.Entry(t2.sub_frame).pack(side="left")
    emiss = ttk.Entry(t2.sub_frame).get()
    
    
    
    
    t3 = ToggledFrame(userinputFrame, text='Dimensions of CubeSat', relief="raised", borderwidth=1)
    t3.pack(fill="x", expand=1, pady=2, padx=2, anchor="n")
    
    Length_var = tk.DoubleVar()
    Width_var = tk.DoubleVar()
    Height_var = tk.DoubleVar()
    
    #Make a formula so this depends on the above three inputed values and changes the area in gui
    A1f_var = tk.DoubleVar()
    
    ttk.Label(t3.sub_frame, text= "Length:").pack(side="left", expand=1)
    ttk.Entry(t3.sub_frame).pack(side="left")
    Length = ttk.Entry(t3.sub_frame).get()
    
    
    ttk.Label(t3.sub_frame, text= "Width").pack(side="left", expand=1)
    ttk.Entry(t3.sub_frame).pack(side="left")
    Width= ttk.Entry(t3.sub_frame).get()
    
    
    ttk.Label(t3.sub_frame, text= "Height").pack(side="left", expand=1)
    ttk.Entry(t3.sub_frame).pack(side="left")
    Height = ttk.Entry(t3.sub_frame).get()
    
    ttk.Label(t3.sub_frame, text= "Area of 1 Face(m^2)").pack(side="left", expand=1)
    ttk.Entry(t3.sub_frame, textvariable = A1f_var).pack(side="left")
    
    
    
    
    t4 = ToggledFrame(userinputFrame, text='Electrical Parameters', relief="raised", borderwidth=1)
    t4.pack(fill="x", expand=1, pady=2, padx=2, anchor="n")
    
    ttk.Label(t4.sub_frame, text= "Electrical Heat Output (W):").pack(side="left", expand=1)
    ttk.Entry(t4.sub_frame).pack(side="left")
    Qelec = ttk.Entry(t4.sub_frame).get()
    
    
    ttk.Label(t4.sub_frame, text= "Tmax").pack(side="left", expand=1)
    ttk.Entry(t4.sub_frame).pack(side="left")
    Tmax = ttk.Entry(t4.sub_frame).get()
    
    ttk.Label(t4.sub_frame, text= "Tmin").pack(side="left", expand=1)
    ttk.Entry(t4.sub_frame).pack(side="left")
    Tmin = ttk.Entry(t4.sub_frame).get()
    
    def sim():
        A1f = A1f_var.get()
        rorb, P, orbVel = runSimulation(G, Boltz, massE, r_Earth, albedo_avg, earth_emission, solar_flux, A1f, absorb1, emitAl, dt)
        orbitalRadius = tk.Label(userinputFrame, text = "Orbital Radius (km): {0}".format(rorb/1000))
        orbitalRadius.pack()
        orbVel = tk.Label(userinputFrame, text = "Orbital Velocity(m/s): {0}".format(round(int(orbVel), 0)))
        orbVel.pack()
        period = tk.Label(userinputFrame, text = "Orbital Period (s): {0}".format(int(P)))
        period.pack()
        return
    runSim = tk.Button(userinputFrame, text = 'Run', command = sim)
    runSim.pack()
    userinputFrame.pack(fill=tk.BOTH, side=tk.LEFT)


    

gui.mainloop()



