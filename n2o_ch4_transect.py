#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Jeff Nivitanont, U. Wyoming 2022

import numpy as np
import pandas as pd
from scipy.integrate import simpson
from .haversine import haversine_dist

def log_interp(x, xp, yp):
    ylog = np.interp(x, xp, np.log(yp))
    return np.exp(ylog)

def retrieve_transect(df, sttime, endtime):
    year, month, day = df.index[0].year, df.index[0].month, df.index[0].day
    t0, t1 = np.datetime64(f'{year}-{month:02d}-{day:02d}T{sttime}'), np.datetime64(f'{year}-{month:02d}-{day:02d}T{endtime}')
    return df.loc[t0:t1].copy()

def calc_transect(df, sourcept, angle_correction=None, interp=True, centroid=False):
    slat, slon = sourcept
    lat = df['Lat'].copy()
    lon = df['Lon'].copy()
    # c2h2_0 = df['C2H2']
    n2o_0 = df['N2O'].copy()
    ch4_0 = df['CH4'].copy()

    # find the lowest 10% of the spike, use as enhancement reference
    # c0 = c2h2_0.min()*0.9 + c2h2_0.max()*0.1
    n0 = n2o_0.min()*0.9 + n2o_0.max()*0.1
    m0 = ch4_0.min()*0.9 + ch4_0.max()*0.1
    n2o_1 = n2o_0-n0
    n2o_1[n2o_1 < 0] = 0
    # c2h2_1 = c2h2_0-c0
    # c2h2_1[c2h2_1 < 0] = 0
    ch4_1 = ch4_0-m0
    ch4_1[ch4_1 < 0] = 0

    #calculate transect slant distances
    if centroid:
        x = np.arange(n2o_0.size)
        maxn2o = int(np.round(np.trapz(x*n2o_0, x)/np.trapz(n2o_0, x)))
    else:
        maxn2o = n2o_0.argmax()
    ctrlat = lat[maxn2o]
    ctrlon = lon[maxn2o]
    y0 = haversine_dist(ctrlat, slat, ctrlon, slon)
    dh0 = haversine_dist(ctrlat, lat, ctrlon, lon)
    dh0[0:maxn2o] = -dh0[0:maxn2o]

    if interp:
        #interpolate slant concentrations
        dh = np.arange(-250, 251)
        n2o_1 = np.interp(dh, dh0, n2o_1)
        ch4_1 = np.interp(dh, dh0, ch4_1)
    else:
        dh = dh0

    #calculate crosswind and downwind distances for interpolated points
    dy = np.zeros(dh.size)
    dx = np.zeros(dh.size)
    if angle_correction is not None:
        dtor = np.pi/180.
        dy = y0-dh*np.sin(angle_correction*dtor)
        dx = dh*np.abs(np.cos(angle_correction*dtor))
    else:
        dx = dh.copy()

    return [n2o_1, ch4_1, dh, dx, dy, y0]


if __name__ == '__main__':
    '''
    Main script is not fully operational.
    '''
    sttime = sys.argv[1]
    endtime = sys.argv[2]
    save = bool(int(sys.argv[3]))
    sourcept = (40.484295, -104.771116) #lat, lon
    """
    Find max, find distances from max, integrate.
    """
    dft = retrieve_transect(df1, sttime, endtime)
    [n2o_1, ch4_1, dh, dx, dy, y0] = calc_transect(dft, sourcept, angle_correction=21)

    n_intg = simpson(n2o_1, dh)
    # c_intg = simpson(c2h2_1)
    m_intg = simpson(ch4_1, dh)

    # tr_ratio = n_intg/c_intg
    nc_ratio = m_intg/n_intg


    """
    Plot plume.
    """
    fig = plt.figure(figsize=(8, 6), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    twin2 = ax.twinx()

    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines['right'].set_position(("axes", .91))
    p1, = ax.plot(dh, n2o_1, 'b', label="N2O")
    p3, = twin2.plot(dh, ch4_1, 'g', label="CH4")

    ax.set_xlabel("distance from center(meters)")
    ax.set_ylabel("N2O (ppbv)")
    top = n2o_1.max() + 5 - (n2o_1.max() % 10)
    ax.set_ylim([0, top])
    ax.grid(axis='y')
    twin2.set_ylabel("CH4 (ppbv)")
    top = ch4_1.max() + 20 - (ch4_1.max() % 10)
    twin2.set_ylim([0, top])

    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())

    ax.text(.05, .95, f'CH4 plume intg.(ppm m): {m_intg: 2.4f}\nN2O plume intg.(ppm m): \
    {n_intg: 2.4f}\nCH4/N2O ratio: {nc_ratio: 2.4f}',
            transform=ax.transAxes, bbox=dict(boxstyle="round", ec='k', fc='lightcyan', alpha=0.5), ha='left', va='top')
    # ax.set_title('Integration method')
    plt.show()
