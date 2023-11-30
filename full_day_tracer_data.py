#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Jeff Nivitanont, U. Wyoming 2022

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from glob import glob
from scipy.integrate import simps
import sys
import os
from data_io import mergePicarroFiles

makeplot = False
date = sys.argv[1]

if __name__ == "__main__":
    # read in TILDAS files
    year, month, day = date[:4], date[4:6], date[6:]
    til_files = glob(f"./data/{date}/TILDAS/*TILDAS*")+glob(f"./data/{date}/*TILDAS*")
    print('TILDAS files:')
    print(til_files)

    if not os.path.exists(f'./figures/{date}'):
        os.makedirs(f'./figures/{date}')

    df_list = []

    for til in til_files:
        try:
            dftemp = pd.read_csv(til, skiprows=3, usecols=['PC', 'CH4 Concentration (ppb)','C2H2 Concentration (ppb)', 'N2O Concentration (ppb)']) #, 'Lat', 'Lon' ])
    #         dfp = dftemp.dropna()
            df_list.append(dftemp)
        except:
            print(f'Error reading {til}')
    tildf = pd.concat(df_list, ignore_index=True)
    tildf['time'] = tildf['PC'].apply(dt.datetime.strptime, args=('%H%M%S*%Y%m%d',))
    tildf.set_index('time', inplace=True)
    tildf.sort_index(inplace=True)
    tildf.rename(columns={'CH4 Concentration (ppb)': 'CH4','C2H2 Concentration (ppb)': 'C2H2', 'N2O Concentration (ppb)': 'N2O'}, inplace=True)
    t0, t1 = tildf.index[0], tildf.index[-1]
    time = tildf.index

    # read in Airmar WX files
    airmar_headers = ['PC', 'UTC hhmmss', 'UTC Year', 'UTC Month', 'UTC Day',
        'Latitude (DD.ddd +N)', 'Longitude (DDD.ddd -W)', 'GPS Quality',
        'Altitude (m)', 'Air Temperature (C)', 'RH(%)', 'Dew Point (C)',
        'Wind Direction (Deg True)', 'Wind Direction (Deg Mag)',
        'Wind Speed (m/s)', 'Pressure (bar)', 'PCB1 Temperature (C)',
        'PCB2 Temperature (C)', 'Supply Voltage (VDC)', 'Heading(deg)',
        'GPSCorWindDirTrue (deg)', 'GPSCorWindDirMag (deg)',
        ' GPSCorWindSpeed (kts)', 'GPSCorWindSpeed (m/s)',
        'GPSGroundSpeed (m/s)']
    use_headers = ['PC', 'Latitude (DD.ddd +N)', 'Longitude (DDD.ddd -W)',
        'Altitude (m)', 'Air Temperature (C)', 'RH(%)', 'Dew Point (C)',
        'Pressure (bar)','Heading(deg)', 'GPSCorWindDirTrue (deg)',
        'GPSCorWindDirMag (deg)', 'GPSCorWindSpeed (m/s)', 'GPSGroundSpeed (m/s)']
    head_dict = {'Latitude (DD.ddd +N)': 'Lat', 'Longitude (DDD.ddd -W)': 'Lon'}
    wx_files = glob(f"./data/{date}/*WX*")+glob(f"./data/{date}/AIRMAR/*WX*")
    print('AIRMAR files:')
    print(wx_files)
    
    # read in Picarro files
    pic_files = glob(f"./data/{date}/PICARRO/BFADS*.dat")+glob(f"./data/{date}/BFADS*.dat")
    if not pic_files:
        print('No Picarro files found in directory.')
    else:
        print('Picarro files:')
        print(pic_files)
        picdf = mergePicarroFiles(pic_files)
    
    if not os.path.exists(f'./figures/{date}'):
        os.makedirs(f'./figures/{date}')

    df_list = []

    for wxf in wx_files:
        try:
            dftemp = pd.read_csv(wxf, skiprows=4, index_col=False, names=airmar_headers, usecols=use_headers)
    #         dfp = dftemp.dropna()
            df_list.append(dftemp)
        except:
            print(f'Error reading {wxf}')
    wxdf = pd.concat(df_list, ignore_index=True)
    wxdf['time'] = wxdf['PC'].apply(dt.datetime.strptime, args=('%H%M%S*%Y%m%d',))
    wxdf.set_index('time', inplace=True)
    wxdf.sort_index(inplace=True)
    wxdf.rename(columns=head_dict, inplace=True)

    if makeplot:
        fig, ax = plt.subplots(1,1, figsize=(16,4), dpi=100)
        fig.subplots_adjust(right=0.75)

        twin1 = ax.twinx()
        twin2 = ax.twinx()

        # Offset the right spine of twin2.  The ticks and label have already been
        # placed on the right by twinx above.
        twin2.spines['right'].set_position(("axes", 1.1))

        p1, = ax.plot(time, tildf['N2O'], 'b', label="N2O")
        p2, = twin1.plot(time, tildf['C2H2'], 'r',  label="C2H2")
        p3, = twin2.plot(time, tildf['CH4'], 'g', label="i12CH4")

        ax.set_xlabel("index")
        ax.set_ylabel("N2O (ppbv)")
        ax.set_ylim([330,360])
        ax.set_xlim([t0, t1])

        ax.grid(axis='y')
        twin1.set_ylabel("C2H2 (ppbv)")
        twin1.set_ylim([0, 150])
        twin2.set_ylabel("CH4 (ppbv)")
        twin2.set_ylim([1800,3000])

        ax.yaxis.label.set_color(p1.get_color())
        twin1.yaxis.label.set_color(p2.get_color())
        twin2.yaxis.label.set_color(p3.get_color())

        # plt.savefig(f'./figures/{date}/{date}_TILDAS_tracer.png', bbox_inches='tight')
        plt.show()
