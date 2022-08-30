#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Jeff Nivitanont, U. Wyoming 2022

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from .haversine import haversine_dist

verbose = False

c2h2_flux = 0
n2o_flux = 10

def enhwind(tildf,
            wxdf,
            date,
            index,
            sttime,
            timewindow=10,
            plot=False,
            satellite=False,
            verbose=False):
    if not type(sttime)==str:# missing times in .xls
        area, mmin, mmax, mq5, mq95, enh5, time, avgwind = np.full(8, np.nan)
        returndf = pd.DataFrame({'avgwind(m/s)': avgwind, 'ch4plumeintg(ppb m)': area, 'ch4min(ppb)': mmin,
            'ch4max(ppb)': mmax, 'ch4bot5(ppb)': mq5, 'ch4top5(ppb)': mq95, 'ch4enh5(ppb)': enh5}, index=[index])
        return returndf
    year, month, day = date[:4], date[4:6], date[6:]
    t0 = np.datetime64(f'{year}-{month}-{day}T{sttime}')
    t1 = t0 + np.timedelta64(timewindow, 'm')
    lats = tildf.loc[t0:t1,'Lat']
    lons = tildf.loc[t0:t1,'Lon']
    ch4 = tildf.loc[t0:t1,'CH4']
    dx = haversine_dist(lats, lats.shift(-1), lons, lons.shift(-1))
    dx.fillna(0, inplace=True)
    x = np.cumsum(dx.values)
    if verbose:
        print(f'{area: 2.3f} ppb meters')

    try:
        area = np.trapz(ch4, x)
        mq5, mq95 = np.nanquantile(ch4, [.05, .95])
        mmin = np.nanmin(ch4)
        mmax = np.nanmax(ch4)
        enhmax5 = mmax-mq5
        enh5 = mq95-mq5
        time = ch4.index
        avgwind = np.nanmean(wxdf.loc[t0:t1, 'GPSCorWindSpeed (m/s)'])
    except:
        area, mmin, mmax, mq5, mq95, enh5, enhmax5, time, avgwind = np.full(9, np.nan)

    if plot:
        plotch4 = ch4-mq5
        plotch4[plotch4<0.] = 0.
        if satellite:
            plt.style.use('seaborn-whitegrid')
            fig = plt.figure(figsize=(8, 4), dpi=150)
            ax1 = fig.add_subplot(1, 2, 1)
            terrain = cimgt.GoogleTiles(style='satellite')
            ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
            sc = ax2.scatter(lons, lats, c=plotch4, s=2,  cmap='hot_r',)
            ax2.add_image(terrain, 18)
            ax2.gridlines()
            cbar = plt.colorbar(sc, orientation = 'horizontal', shrink=.3)
            cbar.set_label('CH4 enhancement (ppbv)')
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=(5, 4), dpi=150)
        ax1.plot(time, plotch4)
        ax1.set_xlabel('Time')
        ax1.set_ylabel('CH4 enhancement (ppbv)')
        ax1.tick_params(axis='x', rotation=45)
        plt.show()
    returndf = pd.DataFrame({'avgwind(m/s)': avgwind, 'ch4plumeintg(ppb m)': area, 'ch4min(ppb)': mmin,
        'ch4max(ppb)': mmax, 'ch4bot5(ppb)': mq5, 'ch4top5(ppb)': mq95, 'ch4enh5(ppb)': enh5, 'ch4enhmax5(ppb)': enhmax5}, index=[index])
    return returndf

def intg_peak(tildf):
    '''
    Plot a slice of TILDAS data and returns integrated peak in (ppm m).User supplied left and right bounds for integration.
    '''
    ch4 = tildf.loc[:,'CH4']
    time = ch4.index
    mq5 = np.nanquantile(ch4, [.05]) #use this as background CH4

    fig, ax1 = plt.subplots(1, 1, figsize=(5, 4), dpi=150)
    ax1.plot(time, ch4-mq5)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('CH4 enhancement (ppbv)')
    ax1.tick_params(axis='x', rotation=45)
    plt.show()

    t0 = input('Enter integration start time (HH:MM:SS):')
    t1 = input('Enter integration end time (HH:MM:SS)')

    ch4 = tildf.loc[t0:t1,'CH4']
    lats = tildf.loc[t0:t1,'Lat']
    lons = tildf.loc[t0:t1,'Lon']
    dx = haversine_dist(lats, lats.shift(-1), lons, lons.shift(-1))
    dx.fillna(0, inplace=True)
    x = np.cumsum(dx.values)
    area = np.trapz(ch4-mq5, x)
    return area



if __name__ == '__main__':
    sttime = sys.argv[1]
    endtime = sys.argv[2]
    save = bool(int(sys.argv[3]))
    year, month, day = date[:4], date[4:6], date[6:]
    if save:
        outtxt = open(f'./figures/{date}/{date}_{sttime.replace(":", "")}T{endtime.replace(":", "")}_regression_output.txt', 'w')

    """
    Start integration method.

    """
    t0, t1 = np.datetime64(f'{year}-{month}-{day}T{sttime}'), np.datetime64(f'{year}-{month}-{day}T{endtime}')
    # find the lowest 10% of the spike, use as enhancement reference
    c0 = df1.loc[t0:t1,'C2H2'].min()*0.9 + df1.loc[t0:t1,'C2H2'].max()*0.1
    n0 = df1.loc[t0:t1,'N2O'].min()*0.9 + df1.loc[t0:t1,'N2O'].max()*0.1
    m0 = df1.loc[t0:t1,'CH4'].min()*0.9 + df1.loc[t0:t1,'CH4'].max()*0.1

    fig = plt.figure(figsize=(13, 9), dpi=100)
    ax = fig.add_subplot(2, 2, 1)
    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines['right'].set_position(("axes", .91))
    n2o_1 = df1.loc[t0:t1, 'N2O']-n0
    n2o_1[n2o_1 < 0] = 0
    c2h2_1 = df1.loc[t0:t1, 'C2H2']-c0
    c2h2_1[c2h2_1 < 0] = 0
    ch4_1 = df1.loc[t0:t1, 'CH4']-m0
    ch4_1[ch4_1 < 0] = 0

    time = n2o_1.index
    p1, = ax.plot(time, n2o_1, 'b', label="N2O")
    p2, = twin1.plot(time, c2h2_1, 'r',  label="C2H2")
    p3, = twin2.plot(time, ch4_1, 'g', label="CH4")

    ax.set_xlabel("index")
    ax.set_ylabel("N2O (ppbv)")
    top = n2o_1.max() + 10 - (n2o_1.max() % 10)
    ax.set_ylim([0, top])
    ax.set_xlim([t0, t1])
    ax.grid(axis='y')
    twin1.set_ylabel("C2H2 (ppbv)")
    top = c2h2_1.max() + 10 - (c2h2_1.max() % 10)
    twin1.set_ylim([0, top])
    twin2.set_ylabel("CH4 (ppbv)")
    top = ch4_1.max() + 10 - (ch4_1.max() % 10)
    twin2.set_ylim([0, top])

    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())

    nc_ratio = simps(ch4_1)/simps(n2o_1)
    cc_ratio = simps(ch4_1)/simps(c2h2_1)
    tr_ratio = simps(n2o_1)/simps(c2h2_1)

    nc_emit = nc_ratio*n2o_flux/22.4*16.043*60/1e3 # mol/L * g/mol * min/hr *kg/g
    cc_emit = cc_ratio*c2h2_flux/22.4*16.043*60/1e3

    ax.text(.05, .95, f'CH4/N2O ratio: {nc_ratio: 2.4f}\nCH4/C2H2 ratio: {cc_ratio: 2.4f}\nN2O/C2H2 ratio: \
    {tr_ratio: 2.4f}\nCH4 emissions (N2O): {nc_emit: 2.4f} kg/hr\nCH4 emissions (C2H2): {cc_emit: 2.4f} kg/hr',
            transform=ax.transAxes, bbox=dict(boxstyle="round", ec='k', fc='lightcyan', alpha=0.5), ha='left', va='top')
    ax.set_title('Integration method')

    """
    Start single tracer regression (N2O).
    """

    results = sm.OLS(ch4_1, n2o_1).fit()

    beta1 = results.predict([1])[0]
    ch4_flux = beta1*n2o_flux/22.4*16.043*60/1e3
    r2 = results.rsquared

    x = [n2o_1.min(), n2o_1.max()]
    y = results.predict(x)

    ax = fig.add_subplot(2, 2, 2)
    ax.scatter(n2o_1, ch4_1)
    ax.plot(x, y, 'r--')
    ax.set_ylabel('CH4 (ppb)')
    ax.set_xlabel('N2O (ppb)')
    ax.text(.05, .95, f'CH4 emissions: {ch4_flux: 2.4f} kg/hr\nR$^2$: {r2: 2.4f}', ha='left', va='top',
            transform=ax.transAxes, bbox=dict(boxstyle="round", ec='k', fc='lightcyan',))
    ax.set_title('Single-tracer regression (N2O)')

    if save:
        outtxt.write(f'Single-tracer regression flux (N2O): {ch4_flux: 2.4f} kg CH4/hr\n\n')
        outtxt.write(str(results.summary()))
        outtxt.write('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')

    if verbose:
        print('Single-tracer regression flux (N2O):')
        print(ch4_flux, 'kg CH4/hr\n')
        print(results.summary())
        print('\n\n')


    """
    Start single tracer regression (C2H2).
    """

    results = sm.OLS(ch4_1, c2h2_1).fit()

    beta1 = results.predict([1])[0]
    ch4_flux = beta1*c2h2_flux/22.4*16.043*60/1e3
    r2 = results.rsquared

    x = [c2h2_1.min(), c2h2_1.max()]
    y = results.predict(x)

    ax = fig.add_subplot(2, 2, 3)
    ax.scatter(c2h2_1, ch4_1)
    ax.plot(x, y, 'r--')
    ax.set_ylabel('CH4 (ppb)')
    ax.set_xlabel('C2H2 (ppb)')
    ax.text(.05, .95, f'CH4 emissions: {ch4_flux: 2.4f} kg/hr\nR$^2$: {r2: 2.4f}', ha='left', va='top',
            transform=ax.transAxes, bbox=dict(boxstyle="round", ec='k', fc='lightcyan',))
    ax.set_title('Single-tracer regression (C2H2)')
    if save:
        outtxt.write(f'Single-tracer regression flux (C2H2): {ch4_flux: 2.4f} kg CH4/hr\n\n')
        outtxt.write(str(results.summary()))
        outtxt.write('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')

    if verbose:
        print('Single-tracer regression flux (C2H2):')
        print(ch4_flux, 'kg CH4/hr\n')
        print(results.summary())
        print('\n\n')


    """
    Start dual tracer regression.
    """

    results = sm.OLS(n2o_1, c2h2_1).fit()
    nc_ratio = results.predict([1])[0]
    r2nc = results.rsquared_adj

    exo = pd.concat([c2h2_1, n2o_1], axis=1)
    results = sm.OLS(ch4_1, exo).fit()

    beta1 = results.predict([0, 1])[0]
    beta2 = results.predict([1, 0])[0]
    ch4_flux = (beta1*c2h2_flux + beta2*n2o_flux)/22.4*16.043*60/1e3
    r2 = results.rsquared_adj

    x = np.array([c2h2_1.min(), c2h2_1.max()])
    y = x*nc_ratio
    ax = fig.add_subplot(2, 2, 4)
    ax.scatter(c2h2_1, n2o_1)
    ax.plot(x, y, 'r--')
    ax.set_ylabel('N2O (ppb)')
    ax.set_xlabel('C2H2 (ppb)')
    ax.text(.05, .95, f'N2O/C2H2 ratio: {nc_ratio: 2.4f}\nN2O/C2H2 corr.: {r2nc: 2.4f}\nCH4 emissions: {ch4_flux: 2.4f} kg/hr\nR$^2$: {r2: 2.4f}',
            transform=ax.transAxes, bbox=dict(boxstyle="round", ec='k', fc='lightcyan',), ha='left', va='top')

    # x1 = np.linspace(n2o_1.min(), n2o_1.max(), 30)
    # x2 = np.linspace(c2h2_1.min(), c2h2_1.max(), 30)
    # xx, yy = np.meshgrid(x1, x2)
    # z = xx*beta2 + yy*beta1
    # ax = fig.add_subplot(2, 2, 4, projection='3d')
    # ax.scatter(n2o_1, c2h2_1, ch4_1)
    # ax.plot_surface(xx, yy, z, color='red', linewidth=0.1, shade=False)
    # ax.set_ylabel('C2H2 (ppb)')
    # ax.set_xlabel('N2O (ppb)')
    # ax.set_zlabel('CH4 (ppb)')
    # ax.text(0, 0, 1, f'CH4 emissions: {ch4_flux: 2.4f} kg/hr\n R$^2$: {r2: 2.4f}',
    #         transform=ax.transAxes, bbox=dict(boxstyle="round", ec='k', fc='lightcyan',))
    ax.set_title('Dual-tracer regression')

    plt.tight_layout(h_pad=1.5, w_pad=1.5)
    plt.subplots_adjust(top=0.93)
    fig.suptitle(f'{date} {sttime}-{endtime} transect')
    if save:
        plt.savefig(f'./figures/{date}/{date}_{sttime.replace(":", "")}T{endtime.replace(":", "")}_transect.png', bbox_inches='tight')
        outtxt.write(f'Dual-tracer regression flux: {ch4_flux: 2.4f} kg CH4/hr\n\n')
        outtxt.write(str(results.summary()))
        outtxt.close()
        plt.close()

    if verbose:
        print('Dual-tracer regression flux:')
        print(ch4_flux, 'kg CH4/hr\n')
        print(results.summary())
        print('\n\n')
        plt.show()
