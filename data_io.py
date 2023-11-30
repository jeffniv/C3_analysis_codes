#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Jeff Nivitanont, U. Wyoming 2022

import numpy as np
import pandas as pd
import os
from glob import glob
import datetime as dt

def read_picarro(fid):
    dfp = pd.read_table(fid, delim_whitespace=True)[['EPOCH_TIME', 'CH4_dry']]
    timestructs = list(map(dt.datetime.utcfromtimestamp, dfp['EPOCH_TIME']))
    dfp['time'] = pd.to_datetime(timestructs)
    dfp.set_index('time', inplace=True)
    dfp.rename(columns={'CH4_dry': 'CH4'}, inplace=True)
    dfp['CH4'] = dfp['CH4']*1e3
    return dfp

def mergePicarroFiles(pic_files):
    df_list = []
    for pic in pic_files:
        try:
            df_list.append(read_picarro(pic))
        except:
            print(f'Error reading {pic}')
    dfp = pd.concat(df_list)
    dfp.sort_index(inplace=True)
    return dfp