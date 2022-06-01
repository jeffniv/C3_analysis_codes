#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Jeff Nivitanont, U. Wyoming 2022

import numpy as np

def gplume_fwd(x, y, z, H, Q, U, SC):
    #  Compute contaminant concentration (kg/m**3) at a given  
    #   set of receptor locations using the standard Gaussian plume 
    #   solution.  This code handles a single source (located at the 
    #   origin) and multiple receptors. Assumes complete reflectoin from the
    #   surface and does not account for reflectio from PBL.
    #
    #  Input parameters: 
    #     x - receptor locations: distance along the wind direction, with 
    #         the source at x=0 (m)  
    #     y - receptor locations: cross-wind direction (m) 
    #     z - receptor locations: vertical height (m) 
    #     H - source height (m) 
    #     Q - contaminant emission rate (kg/s) 
    #     U - wind velocity (m/s) 
    #     SC - stability class 1-5(see below)
    # 
    #Output: 
    #     C - contaminant concentration (kg/m**3) 
    # 
    # code from John Stockie except for dispersion coeffient calculations
    #
    # First, define the cut-off velocity, below which concentration = 0. 
    # Determine the sigma coefficients using rural dispersion coefficients 
    # based on stability class   
    #
    #*********Stability Classs Determination*************
    # Go to http://www.ready.noaa.gov/READYpgclass.php for descriptions
    # of each stability class
    # Parameters calculated according to Briggs (1973) functions for rural
    # conditions

    Umin = 0.0; 

    if SC==1:
        #stability class 1: A (Extremely Unstable)
        sigmay = 0.22*x*np.abs(1+0.0001*x)**(-0.5) ;
        sigmaz = 0.2*x;

    elif SC==2:
        #stability class 2: B (Moderately Unstable)
        sigmay = 0.16*x*np.abs(1+0.0001*x)**(-0.5);
        sigmaz = 0.12*x;

    elif SC==3:
        #stability class 3: C (Slightly Unstable)
        sigmay = 0.11*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.08*x*np.abs(1+0.0002*x)**(-0.5);

    elif SC==4:
        #stability class 4: D (Neutral)
        sigmay = 0.08*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.06*x*np.abs(1+0.0015*x)**(-0.5); 

    elif SC==5:
        #stability class 5: E (Slightly Stable) *Night Only
        sigmay = 0.06*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.03*x*np.abs(1+0.0003*x)**(-1); 

    else:
        #stability class 6: F (Moderately Stable) *Night Only
        sigmay = 0.04*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.016*x*np.abs(1+0.0003*x)**(-1); 
    #end


    #Calculate the contaminant concentration (kg/m**3) using Ermak's formula. 
    if U < Umin:   
        C = 0*z; 
    else:
        C  = Q/(2*np.pi*U*sigmay*sigmaz)*np.exp( -0.5*y**2./sigmay**2 )*(np.exp(-0.5*(z-H)**2./sigmaz**2)+np.exp(-0.5*(z+H)**2./sigmaz**2) );   
        ii = np.logical_or(np.isnan(C), np.isinf(C));   
        C[ii] = 0;  
        # Set all NaN or inf values to zero. 

    return C


def gplume_inv(x, y, z, H, C, U, SC):
    #  Compute contaminant concentration (kg/m**3) at a given  
    #   set of receptor locations using the standard Gaussian plume 
    #   solution.  This code handles a single source (located at the 
    #   origin) and multiple receptors. Assumes complete reflectoin from the
    #   surface and does not account for reflectio from PBL.
    #
    #  Input parameters: 
    #     x - receptor locations: distance along the wind direction, with 
    #         the source at x=0 (m)  
    #     y - receptor locations: cross-wind direction (m) 
    #     z - receptor locations: vertical height (m) 
    #     H - source height (m)
    #     C - contaminant concentration (kg/m**3)
    #     U - wind velocity (m/s) 
    #     SC - stability class 1-5(see below)
    # 
    #Output: 
    #     Q - contaminant emission rate (kg/s) 
    # 
    # code from John Stockie except for dispersion coeffient calculations
    #
    # First, define the cut-off velocity, below which concentration = 0. 
    # Determine the sigma coefficients using rural dispersion coefficients 
    # based on stability class   
    #
    #*********Stability Classs Determination*************
    # Go to http://www.ready.noaa.gov/READYpgclass.php for descriptions
    # of each stability class
    # Parameters calculated according to Briggs (1973) functions for rural
    # conditions

    Umin = 0.0; 

    if SC==1:
        #stability class 1: A (Extremely Unstable)
        sigmay = 0.22*x*np.abs(1+0.0001*x)**(-0.5) ;
        sigmaz = 0.2*x;

    elif SC==2:
        #stability class 2: B (Moderately Unstable)
        sigmay = 0.16*x*np.abs(1+0.0001*x)**(-0.5);
        sigmaz = 0.12*x;

    elif SC==3:
        #stability class 3: C (Slightly Unstable)
        sigmay = 0.11*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.08*x*np.abs(1+0.0002*x)**(-0.5);

    elif SC==4:
        #stability class 4: D (Neutral)
        sigmay = 0.08*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.06*x*np.abs(1+0.0015*x)**(-0.5); 

    elif SC==5:
        #stability class 5: E (Slightly Stable) *Night Only
        sigmay = 0.06*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.03*x*np.abs(1+0.0003*x)**(-1); 

    else:
        #stability class 6: F (Moderately Stable) *Night Only
        sigmay = 0.04*x*np.abs(1+0.0001*x)**(-0.5); 
        sigmaz = 0.016*x*np.abs(1+0.0003*x)**(-1); 
    #end


    #Calculate the contaminant concentration (kg/m**3) using Ermak's formula. 
    if U < Umin:   
        Q = 0*z; 
    else:
        Q  = C*(2*np.pi*U*sigmay*sigmaz)*np.exp(-0.5*y**2./sigmay**2)*(np.exp(-0.5*(z-H)**2./sigmaz**2)+np.exp(-0.5*(z+H)**2./sigmaz**2));   
        ii = np.logical_or(np.isnan(Q), np.isinf(Q));   
        Q[ii] = 0;  
        # Set all NaN or inf values to zero. 

    return Q






