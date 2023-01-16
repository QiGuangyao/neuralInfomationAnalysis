#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:25:03 2022

@author: qiguangyao
"""

#%% import lib
import numpy as np
import os # os stands for "operating system" and includes read/write routines etc. 
from scipy import io # this is for importing Matlab data files
from scipy import stats # here we import a whole sub-library of stats functions
from matplotlib import pyplot as plt # all of our plotting is done with plt
import pickle as pkl
import seaborn as sns
import copy
import mat73
import scipy.io as sio
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import pandas as pd
import random
from scipy.stats import pearsonr
#%% import data
path1 = '/Volumes/10.10.49.10/gyqi/Light Sheet Data/'
cellRespLand20221029_3_1_5Bins = sio.loadmat(path1+'/cellRespLand20221029_3_1_5Bins.mat')
cell_resp_bins = cellRespLand20221029_3_1_5Bins['cell_resp_bins']
#%% functions
def corrfunc(x, y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, _ = pearsonr(x, y)
    ax = ax or plt.gca()
    ax.annotate(f'ρ = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)
def eta_squared(aov):
    aov['eta_sq'] = 'NaN'
    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    return aov
def omega_squared(aov):
    mse = aov['sum_sq'][-1]/aov['df'][-1]
    aov['omega_sq'] = 'NaN'
    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*mse))/(sum(aov['sum_sq'])+mse)
    return aov
def getAovANOVA(signTFR,oneWay):
    """
    signTFR: single-trial, single-unit acti
    """
    columns = ('context','location','acti')
    #condition
    FR = copy.deepcopy(signTFR)
    # FRCondition[:,-3] = fishLocaToClassfication(FRCondition[:,-3])
    dfFR = pd.DataFrame(FR)
    dfFR = dfFR.replace([np.inf, -np.inf], np.nan)
    dfFR = dfFR.dropna(axis=0,how='any')
    dfFR.columns = columns
    if oneWay == True:
        formular = 'acti' +' ~ C(location)'
    else:
        formular = 'acti' +' ~ C(context) +C(location)+C(context)*C(location)'
    modelANOVA = ols(formular, dfFR,missing='drop').fit()
    if len(list(modelANOVA.conf_int().index))>2:
        aov_table = anova_lm(modelANOVA)
        eta_squared(aov_table)
        omega_squared(aov_table)
    else:
        aov_table =  pd.DataFrame()
    return aov_table
# the fractional variable variance (FVV)
def fvv(PEV1,PEV2):
    """
    input: percentage of explained varience
            PEV1 (context), PEV2 (location)
    output: fractional variable variance (FVV)
            FVV = (PEV1-PEV2)/(PEV1+PEV2)
    """
    return (PEV1-PEV2)/(PEV1+PEV2)
#%% calculate overall PEV
loca1 = [i for i in range(5)]*3*(int(cell_resp_bins.shape[1]/15)-2)
cont1 = [0 for i in range(5)]+[1 for i in range(5)]+[2 for i in range(5)]
cont1 = cont1*(int(cell_resp_bins.shape[1]/15)-2)
contLocaPEV = np.full((cell_resp_bins.shape[0],4),np.nan)#P_value, eta C(context), C(location)
threContPEV = np.full((cell_resp_bins.shape[0],3,2),np.nan)#P_value, eta
contLocaFVV = np.full((cell_resp_bins.shape[0],1),np.nan)
for n in range(cell_resp_bins.shape[0]):
    print(n)
    activ = cell_resp_bins[n,range(15,cell_resp_bins.shape[1]-15)]
    signTFR = np.full((activ.shape[0],3),np.nan)
    signTFR[:,2] = activ
    signTFR[:,0] = cont1
    signTFR[:,1] = loca1
    ANOVAResu = getAovANOVA(signTFR)
    contLocaPEV[n,0] = ANOVAResu['PR(>F)'][0]
    contLocaPEV[n,1] = ANOVAResu['eta_sq'][0]
    contLocaPEV[n,2] = ANOVAResu['PR(>F)'][1]
    contLocaPEV[n,3] = ANOVAResu['eta_sq'][1]
    contLocaFVV[n,0] = fvv(contLocaPEV[n,1],contLocaPEV[n,3])
#%% calculate overall PEV for each context
loca1 = [i for i in range(5)]*3*(int(cell_resp_bins.shape[1]/15)-2)
cont1 = [0 for i in range(5)]+[1 for i in range(5)]+[2 for i in range(5)]
cont1 = cont1*(int(cell_resp_bins.shape[1]/15)-2)
threContPEV = np.full((cell_resp_bins.shape[0],3,2),np.nan)#P_value, eta
for n in range(cell_resp_bins.shape[0]):
    print(n)
    activ = cell_resp_bins[n,range(15,cell_resp_bins.shape[1]-15)]
    signTFR = np.full((activ.shape[0],3),np.nan)
    signTFR[:,2] = activ
    signTFR[:,0] = cont1
    signTFR[:,1] = loca1
    for c in range(3):
        seleSignTFR = signTFR[np.array(cont1)==c,:]
        ANOVAResu = getAovANOVA(seleSignTFR,oneWay = True)
        threContPEV[n,c,0] = ANOVAResu['PR(>F)'][0]
        threContPEV[n,c,1] = ANOVAResu['eta_sq'][0] 
#%% save data PEV for location in each context
pkl.dump({
    'threContPEV':threContPEV
    },
         open(path1+'threContPEV_20221029_3_1_5Bins.pickle','wb'))
#%% load data
threContPEV_20221029_3_1_5Bins = pkl.load(open(path1+'/threContPEV_20221029_3_1_5Bins.pickle','rb'))
threContPEV = threContPEV_20221029_3_1_5Bins['threContPEV']
#%% pair-plot
locaPEV = np.full((threContPEV.shape[0],3),np.nan)
for i in range(3):
    locaPEV[:,i] = threContPEV[:,i,1]
seleIndes = [i for i in range(threContPEV.shape[0]) if threContPEV[i,0,0]<0.05 or threContPEV[i,1,0]<0.05 or threContPEV[i,2,0]<0.05 ]
locaPEV = locaPEV[seleIndes,:]
locaPEVDf = pd.DataFrame(locaPEV)
columnsCont = ['Context_A','Context_B','Context_C']
locaPEVDf.columns = columnsCont
#%% shuffle data
loca1 = [i for i in range(5)]*3*(int(cell_resp_bins.shape[1]/15)-2)
cont1 = [0 for i in range(5)]+[1 for i in range(5)]+[2 for i in range(5)]
cont1 = cont1*(int(cell_resp_bins.shape[1]/15)-2)
loca2 = copy.deepcopy(loca1)
cont2 = copy.deepcopy(cont1)

contLocaPEV_Shuf = np.full((cell_resp_bins.shape[0],4),np.nan)#P_value, eta
threContPEV_Shuf = np.full((cell_resp_bins.shape[0],3,2),np.nan)#P_value, eta
contLocaFVV_Shuf = np.full((cell_resp_bins.shape[0],1),np.nan)
for n in range(cell_resp_bins.shape[0]):
    print(n)
    activ = cell_resp_bins[n,range(15,cell_resp_bins.shape[1]-15)]
    signTFR = np.full((activ.shape[0],3),np.nan)
    signTFR[:,2] = activ
    random.shuffle(cont2)
    random.shuffle(loca2)
    signTFR[:,0] = cont2
    signTFR[:,1] = loca2
    
    ANOVAResu = getAovANOVA(signTFR)
    contLocaPEV_Shuf[n,0] = ANOVAResu['PR(>F)'][0]
    contLocaPEV_Shuf[n,1] = ANOVAResu['eta_sq'][0]
    contLocaPEV_Shuf[n,2] = ANOVAResu['PR(>F)'][1]
    contLocaPEV_Shuf[n,3] = ANOVAResu['eta_sq'][1]
    contLocaFVV_Shuf[n,0] = fvv(contLocaPEV_Shuf[n,1],contLocaPEV_Shuf[n,3])
#%%
for n in range(cell_resp_bins.shape[0]):
    print(n)
    contLocaFVV_Shuf[n,0] = fvv(contLocaPEV_Shuf[n,1],contLocaPEV_Shuf[n,3])
#%% save PEV for location and context
pkl.dump({
    'contLocaPEV':contLocaPEV,
    'contLocaFVV':contLocaFVV,
    'contLocaPEV_Shuf':contLocaPEV_Shuf,
    'contLocaFVV_Shuf':contLocaFVV_Shuf,
    },
    open(path1+'/PEV_Context_Location_20221029_3_1_5Bins.pickle','wb'))




























