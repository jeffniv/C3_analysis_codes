#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Jeff Nivitanont, U. Wyoming 2022

import numpy as np

re = 6371e3 # diameter of Earth in meters
dtor = np.pi/180.
def haversine_dist(lat1, lat2, lon1, lon2):
    '''
    Given two lat/lon coordinates, this function calculates the arc distance in meters on a great sphere.
    '''
    phi1 = lat1*dtor
    phi2 = lat2*dtor
    lam1 = lon1*dtor
    lam2 = lon2*dtor
    h = np.sin((phi2-phi1)/2.)**2 + np.cos(phi1)*np.cos(phi2)*np.sin((lam2-lam1)/2)**2
    d = 2*re*np.arcsin(np.sqrt(h))
    return d