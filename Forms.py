# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 16:54:04 2023

@author: mauri
"""

#This module, 'Forms', is used to define the formulas of the aircraft dynamics without having to include them 
#in the main script. This module imports variables from the 'constants'module so that it recognizes variables such as 
#CL0, CLa, and so on. 
import numpy as np
import constants

def Coefficient_of_Lift(a, d):
    return constants.CL0 + constants.CLa * a + constants.CLde * d

def Coefficient_of_Moment(a, d):
    return constants.CM0 + constants.CMa * a + constants.CMde * d

def Coefficient_of_Drag(a, d):
    return constants.CD0 + constants.K * (Coefficient_of_Lift(a, d))**2

def Lift(a, d, velocity, gamma):
    return (0.5 * constants.air_density * velocity**2 * constants.wing_surface *
            Coefficient_of_Lift(a, d))

def Drag(a, d, velocity, gamma):
    return (0.5 * constants.air_density * velocity**2 * constants.wing_surface *
            Coefficient_of_Drag(a, d))

def Moment(a, d, velocity, gamma):
    return (0.5 * constants.air_density * velocity**2 * constants.wing_surface *
           constants.cbar * Coefficient_of_Moment(a, d))

def Engine_Thrust(a, d, the, velocity, gamma):
    return (0.5 * constants.air_density * velocity**2 * constants.wing_surface *
            Coefficient_of_Drag(a, d) * np.cos(a) +
            constants.mass * constants.gravity * np.sin(the) -
            0.5 * constants.air_density * velocity**2 * constants.wing_surface *
            Coefficient_of_Lift(a, d) * np.sin(a))

def Equilibrium(a, velocity, gamma):
    return (-0.5 * constants.air_density * velocity**2 * constants.wing_surface *
            (constants.CL0 + constants.CLa * a - constants.CLde * (constants.CM0 + constants.CMa * a) / constants.CMde) * np.cos(a) -
            0.5 * constants.air_density * velocity**2 * constants.wing_surface *
            (constants.CD0 + constants.K * (constants.CL0 + constants.CLa * a - constants.CLde * (constants.CM0 + constants.CMa * a) / constants.CMde)**2) * np.sin(a) +
            constants.mass * constants.gravity * np.cos(a + gamma))
