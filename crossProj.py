#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 11:00:11 2022

@author: qiguangyao
"""


#%% import lib
import numpy as np
from scipy import io # this is for importing Matlab data files
from sklearn.decomposition import PCA # check out all the other dimensionality reduction methods in "decomposition"
import pickle as pkl
import copy
import mat73
import scipy.io as sio
#%% import data
path1 = '/Volumes/10.10.49.10/gyqi/Light Sheet Data'
cellRespLand20221029_3_1_5Bins = sio.loadmat(path1+'/cellRespLand20221029_3_1_5Bins.mat')
cell_resp_bins = cellRespLand20221029_3_1_5Bins['cell_resp_bins']
#%% get averaged peristimulus time histograms (PSTHs) seperately
averTriaActi = np.full((cell_resp_bins.shape[0],15),np.nan)
for n in range(cell_resp_bins.shape[0]):
    print(n)
    for t in range(15):
        averTriaActi[n,t] = np.nanmean(cell_resp_bins[n,[i*15+t for i in range(1,54)]])
#%% get averaged peristimulus time histograms (PSTHs) seperately first/second half
averTriaActiFirs = np.full((cell_resp_bins.shape[0],15),np.nan)
averTriaActiSeco = np.full((cell_resp_bins.shape[0],15),np.nan)
badTria = [1,20,55]
for n in range(cell_resp_bins.shape[0]):
    print(n)
    for t in range(15):
        averTriaActiFirs[n,t] = np.nanmean(cell_resp_bins[n,[i*15+t for i in [k for k in range(28) if k !=0 and k != 19]]])
        averTriaActiSeco[n,t] = np.nanmean(cell_resp_bins[n,[i*15+t for i in [k for k in range(28,55) if k !=54]]])
#%% PCA
covMat = [[],[],[]]
featValue = [[],[],[]]
n_featVec = [[],[],[]]
#landmark 1
pcaLand1 = PCA(svd_solver = 'full')
averTriaActiPrep1 = averTriaActi[:,range(5)].T# - np.nanmean(averTriaActi[:,range(5)].T,0))/np.nanstd(averTriaActi[:,range(5)].T,0)
pcaLand1.fit(averTriaActiPrep1)
land1_tran = pcaLand1.transform(averTriaActiPrep1)
covMat[0] = np.cov(averTriaActiPrep1,rowvar=0)#计算协方差矩阵，rowvar=0表示数据的每一列代表一个feature
featValue[0] = pcaLand1.explained_variance_ratio_
n_featVec[0] = pcaLand1.components_
print('landmark 1')
#landmark 2
pcaLand2 = PCA(svd_solver = 'full')
averTriaActiPrep2 = averTriaActi[:,range(5,10)].T# - np.nanmean(averTriaActi[:,range(5,10)].T,0))/np.nanstd(averTriaActi[:,range(5,10)].T,0)
pcaLand2.fit(averTriaActiPrep2)
land2_tran = pcaLand2.transform(averTriaActiPrep2)
covMat[1] = np.cov(averTriaActiPrep2,rowvar=0)#计算协方差矩阵，rowvar=0表示数据的每一列代表一个feature
featValue[1] = pcaLand2.explained_variance_ratio_
n_featVec[1] = pcaLand2.components_
print('landmark 2')
#landmark 3
pcaLand3 = PCA(svd_solver = 'full')
averTriaActiPrep3 = averTriaActi[:,range(10,15)].T# - np.nanmean(averTriaActi[:,range(10,15)].T,0))/np.nanstd(averTriaActi[:,range(10,15)].T,0)
pcaLand3.fit(averTriaActiPrep3)
land3_tran = pcaLand3.transform(averTriaActiPrep3)
covMat[2] = np.cov(averTriaActiPrep3,rowvar=0)#计算协方差矩阵，rowvar=0表示数据的每一列代表一个feature
featValue[2] = pcaLand3.explained_variance_ratio_
n_featVec[2] = pcaLand3.components_
print('landmark 3')
# index = np.argsort(featValue[0])[::-1]
# n_featVec[0] = featVec[0][:, index][:,:PCs]
#%% save data
pkl.dump({
    'covMat':covMat,
    'featValue':featValue,
    'n_featVec':n_featVec,
    'pcaLand1':pcaLand1,
    'pcaLand2':pcaLand2,
    'pcaLand3':pcaLand3    
    },
         open(path1+'/PCACrossProjNeurCova.pickle','wb'))
#%% cross-projected variance
crosProjVari1 = [[],[],[]]
crosProjVari2 = [[],[],[]]
crosProjVari3 = [[],[],[]]

crosProjVari1[0] = featValue[0]
landB_on_landA = np.dot(np.dot(n_featVec[0], covMat[1]), n_featVec[0].T)
landB_on_landA = landB_on_landA.diagonal()/landB_on_landA.diagonal().sum()
crosProjVari1[1] = landB_on_landA[range(5)]
landC_on_landA = np.dot(np.dot(n_featVec[0], covMat[2]), n_featVec[0].T)
landC_on_landA = landC_on_landA.diagonal()/landC_on_landA.diagonal().sum()
crosProjVari1[2] = landC_on_landA[range(5)]

landA_on_landB = np.dot(np.dot(n_featVec[1], covMat[0]), n_featVec[1].T)
landA_on_landB = landA_on_landB.diagonal()/landA_on_landB.diagonal().sum()
crosProjVari2[0] = landA_on_landB[range(5)]
crosProjVari2[1] = featValue[1]
landC_on_landB = np.dot(np.dot(n_featVec[1], covMat[2]), n_featVec[1].T)
landC_on_landB = landC_on_landB.diagonal()/landC_on_landB.diagonal().sum()
crosProjVari2[2] = landC_on_landB[range(5)]

landA_on_landC = np.dot(np.dot(n_featVec[2], covMat[0]), n_featVec[2].T)
landA_on_landC = landA_on_landC.diagonal()/landA_on_landC.diagonal().sum()
crosProjVari3[0] = landA_on_landC[range(5)]
landB_on_landC = np.dot(np.dot(n_featVec[2], covMat[1]), n_featVec[2].T)
landB_on_landC = landB_on_landC.diagonal()/landB_on_landC.diagonal().sum()
crosProjVari3[1] = landB_on_landC[range(5)]
crosProjVari3[2] = featValue[2]
#%% save data
pkl.dump({
    'averTriaActi':averTriaActi,
    'pcaLand1':pcaLand1,
    'pcaLand2':pcaLand2,
    'pcaLand3':pcaLand3,
    'crosProjVari1':crosProjVari1,
    'crosProjVari2':crosProjVari2,
    'crosProjVari3':crosProjVari3,
    },
    open(path1+'/corssProjVariPCA_20221029_3_1_5Bins_Raw.pickle','wb'))
















