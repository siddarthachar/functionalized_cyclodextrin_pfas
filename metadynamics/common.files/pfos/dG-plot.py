import os
import os.path
import sys
import importlib
import yaml
from copy import deepcopy
import subprocess

from time import time
from datetime import datetime

import numpy as np
import pandas as pd
import math

import scipy
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib import rcParams
import matplotlib as mpl
from matplotlib  import cm
from PIL import Image
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from matplotlib.cbook import get_sample_data

import multiprocessing

from sklearn.metrics import mean_squared_error
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap
import copy

import importlib
#import calculateFES
#importlib.reload(calculateFES)
from numpy import ones,copy,cos,tan,pi,linspace

from scipy.interpolate import RegularGridInterpolator,LinearNDInterpolator, RectBivariateSpline,interpn, griddata
import warnings
warnings.filterwarnings("ignore")

############## PRODUCTION RUNNNNNN

C0 = 1/(1661*1e-3)
kb = 8.31446261815324*0.001 #kJ/molK
temperature = 300
kbt = kb*temperature  #kJ/mol
vmin = -70 # use this for min in contour
# probe= ['acd','bcd','gcd']
# analyte = ['pfos','sds','TCAA']

probe=[sys.argv[1]]
analyte=[sys.argv[2]]
#path_here=os.getcwd().split('/')[-1].split('-')
#probe = [path_here[0]]
#analyte = [path_here[-1]]
#print(probe)
print(analyte)

# probe = ['gcd']
# analyte = ['tcaa']

#choose number of blocks
preFactor = 2*kbt #kJ/mol for entropy correction
blocks = 1 #3
cutoff1 = np.sqrt(3.5)
cutoff2 = 2
sfes_mast = []
smast = []

for dimr in [200]:
    for cutoff in [1]:
        dimz = dimr
        pmf_probe_analyte_zero_ref_boundary = {} # collect pmfs
        pmf_probe_analyte_zero_ref_global_minima = {} # collect pmfs

        #for ii in probe:
        for kk in analyte:
            for ii in probe:
                sfes = []
                sval = []
                pmf_probe_analyte_zero_ref_global_minima[ii+'_'+kk] = {}
                pmf_probe_analyte_zero_ref_boundary[ii+'_'+kk]  = {}

                #create holders
                fep = []
                avgFE,sdFE = [],[]
                avgkb,sdkb = [],[]
                #load data 
                colvarFile = 'reweight/COLVAR_reweight' #simDir + systemDir + colvarFileName
                
                print(colvarFile)
                try:
                    colvarData = pd.read_csv(colvarFile, comment="#", 
                                             delim_whitespace=True,names=["time", "angle", "r", 
                                                                          "baxis","mag", "pr", "pz", "bias"],
                                             skip_blank_lines=True)
                    print(colvarData)
                except:
                    continue

                # discard initial few frames [frames where free energy is still being mapped]
                # subsample and randomize
                if kk=='pfos': #probe=='acd' or probe=='gcd' and analyte=='pfos':
                    s = 0
                    d = -1 #500001
                    colvarD = colvarData[s:d].sample(frac = 1)
                else:
                    s = 0
                    d = -1 #1000000
                    colvarD = colvarData[s:d:2].sample(frac = 1)        
                jmp = len(colvarD)//blocks

                print(ii,kk, len(colvarD), blocks, jmp)
                #\print(colvarD)
                for jj in range(blocks):
                    colvarData_init = colvarD[jmp*jj:jmp*(jj+1)] 
                    print(f'Before filter: {colvarData_init.shape}')


                    #colvarData = colvarData_init[(colvarData_init['r']**2 + colvarData_init['baxis']**2 < 3.5)]

                    #colvarData_2 = colvarData_init[(colvarData_init['r'] < cutoff1)]
                    #colvarData_3 = colvarData_2[(colvarData_2['baxis'] < cutoff1)]
                    #colvarData = colvarData_3[(colvarData_3['baxis'] >= -cutoff1)]
                    
                    #colvarData = colvarData_init

                    print(f'After filter (<3. spherical cut): {colvarData.shape}')

                    ############################################################################################################
                    #make a 2D heat map
                    bias = colvarData_init['bias'] 
                    r = colvarData_init['r']
                    z = colvarData_init['baxis']
                    print(len(bias), ' is bias length')
                    # calculate weights of each frame in colvarData - CHECK THIS APPROACH (SUMMING THE BIASES)
                    weights = np.exp((bias) /kbt)        
                    xedges = np.linspace(0,cutoff2,dimr)
                    yedges = np.linspace(-cutoff2,cutoff2,dimz)
                    dr,dz = abs(xedges[0]-xedges[1]),abs(yedges[0]-yedges[1])
                    # get unbiased probabilities
                    fe, xedges,yedges = np.histogram2d(r,z,range=[[0,cutoff2],[-cutoff2,cutoff2]],
                                                       bins=(xedges,yedges),
                                                       weights=weights,density=True)


                    # get unbiased probabilities -- 1 D along r, z
                    fe_r, xedges_r = np.histogram(r, range=[0,cutoff1], bins=xedges, weights=weights,density=True)
                    fe_z, xedges_z = np.histogram(z, range=[-cutoff1,cutoff1], bins=yedges, weights=weights,density=True)
                    fe_r = fe_r.T
                    fe_z = fe_z.T
                    fe_r = -kbt*np.log(fe_r)
                    fe_z = -kbt*np.log(fe_z)
                    
                    fe = fe.T
                    feDf = -kbt*np.log(fe) # Compute free energy
######## correction
                    #correction
                    for ki in range(xedges.shape[0]-1):
                        for kj in range(yedges.shape[0]-1):
        #                     if xedges[ki]**2 + np.abs(yedges[kj])**2 > 1.25: # 2.0:
                   
                            if np.isnan(feDf[kj,ki]) or np.isinf(feDf[kj,ki]):
                                continue
    
                            s = np.sqrt(xedges[ki]**2+yedges[kj]**2)
                            feDf[kj,ki] += preFactor*np.log(s)
####################################



                   # print(feDf)

                    feDf_zero_ref_global_minima = feDf - feDf.min() # UNCOMMENT THIS FOR SETTING GLOBAL FREE ENERGY MINIMA AS ZERO

                    _block_data_global_minima = {'fe': feDf_zero_ref_global_minima, 
                                                 'xedges': xedges, 'yedges': yedges}

                    pmf_probe_analyte_zero_ref_global_minima[ii+'_'+kk]['block'+str(jj)] = _block_data_global_minima

                    # Find values of free energy on the boundary of the surface
                    avg_fe_beyond_cutoff=[]
                    for ki in range(xedges.shape[0]-1):
                        for kj in range(yedges.shape[0]-1):
        #                     if xedges[ki]**2 + np.abs(yedges[kj])**2 > 1.25: # 2.0:
                   
                            if np.isnan(feDf[kj,ki]) or np.isinf(feDf[kj,ki]):
                                continue
                            
                            #cutoffs
                            s = np.sqrt(xedges[ki]**2+yedges[kj]**2)
                            sfes.append(feDf[kj,ki])
                            sval.append(s)
                            if s > 1.2 and s < np.sqrt(3.5): #xedges[ki] > cutoff or yedges[kj] < -cutoff or yedges[kj] > cutoff and xedges[ki] < cutoff1 and yedges[kj] < cutoff1 and yedges[kj] > cutoff1:
                                avg_fe_beyond_cutoff.append(feDf[kj,ki])

                           

                    avg_fe_beyond_cutoff = np.array(avg_fe_beyond_cutoff)

                    print(f'avg:{avg_fe_beyond_cutoff.mean()}, sd:{np.std(avg_fe_beyond_cutoff)},max:{avg_fe_beyond_cutoff.max()}, min:{avg_fe_beyond_cutoff.min()}')

                    #feDf_zero_ref_boundary = feDf - feDf.max() #collect_fe_k_max # UNCOMMENT THIS FOR SETTING BOUNDARY FREE ENERGY TO ZERO
                    feDf_zero_ref_boundary = feDf - avg_fe_beyond_cutoff.mean() #collect_fe_k_max # UNCOMMENT THIS FOR SETTING BOUNDARY FREE ENERGY TO ZERO
                    fe_r = fe_r - avg_fe_beyond_cutoff.mean()
                    fe_z = fe_z - avg_fe_beyond_cutoff.mean()

                    print(xedges.min(), xedges.max(), yedges.min(), yedges.max())

                    _block_data_boundary = {'fe': feDf_zero_ref_boundary,
                                            'xedges': xedges,
                                            'yedges': yedges,
                                            'fe_r': fe_r,
                                             'xedges_r':xedges_r,
                                            'fe_z': fe_z,
                                             'xedges_z':xedges_z}
                                            #'fe_points_on_boundary': fe_points_on_boundary}

                    pmf_probe_analyte_zero_ref_boundary[ii+'_'+kk]['block'+str(jj)] = _block_data_boundary
                    sfes_mast.append(sfes)
                    smast.append(sval)
                    print(ii, kk, jj)
                    print(dimr,dimz,cutoff)

#colors=['blue','green','red']
mt = []
iii = 0
for kk in analyte:
    mtot = []
    for ii in probe:
        bins = np.linspace(0, 2 + 1e-12, 51) # 10 bins, so 11 bin boundaries
        c = np.digitize(smast[iii], bins)
        c = list(c)
        indx = []
        lst =  list(enumerate(c))
        for ijk in range(len(bins)):
            indices = []
            for val in range(len(lst)):
                if lst[val][1]==ijk:
                    indices.append(lst[val][0])
            indx.append(indices)
        mval,stval = [],[]
        sfesp = np.array(sfes_mast[iii])
        for ijk in range(len(indx)):
            mval.append(np.mean(sfesp[indx[ijk]]))
           # if ijk>49 and ijk<95:
           # corrval.append())
        #    stval.append(np.std(sfesp[indx[ijk]]))
        mtot.append(mval)
        iii+=1
    mt.append(mtot)

        #plt.plot(bins[:-1], fe_r, color=colors[i], label=kk.upper())


    #plt.show()
   # plt.savefig('s_fedf_%s.png' %(ii.upper()),dpi=300)
   # plt.close()

colors=['blue','green','red']
iii = 0
for kk in analyte:
    plts = mt[iii]
    plt.figure(figsize=(16, 12))
    plt.xlabel('s [nm]',fontsize=24)
    plt.ylabel('PMF [kJ/mol]',fontsize=24)
    jjj = 0
    plt.xticks(np.linspace(0,50,5),np.linspace(0,2,5))
    for ii in ['$\\beta$-CD']:
        pltp = np.subtract(plts[jjj],np.mean(plts[jjj][30:47]))
        plt.plot(pltp, color=colors[jjj], label=ii)
        jjj+=1
    jjj = 0
    #plt.legend()
    plt.axvline(x=30,color='k',label='Unbounded')
    plt.axvline(x=47,color='k')    
    plt.axhline(y=0,color='k')
    plt.legend()
    # plt.show()
    # plt.savefig('s_fedf_%s.png' %(ii.upper()),dpi=300)
    plt.savefig('s_fed_PMF.png')
    plt.close()
    iii +=1

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# Set the figure size
fig = plt.figure(figsize=(12, 8))

# Adjust spacing
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.35)

# Create a 3x3 grid
grid = GridSpec(2, 2)

# probe= ['bcd'] #,'bcd','gcd']
# analyte = ['pfos']

selected_block = 0 # entire data is used.


colors=['blue','green','red']
for k, kk in enumerate(analyte):
    
    ax = plt.subplot(grid[0,k])
    ax.set_title(kk.upper())
    
    for i, ii in enumerate(probe):
        # Create subplots
        
        
        print(kk,ii)
        #selected_fe = pmf_probe_analyte_zero_ref_global_minima
        selected_fe = pmf_probe_analyte_zero_ref_boundary
        
        
        fe_r = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['fe_r']
        xedges_r = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['xedges_r']
        
        ax.plot(xedges_r[:-1], fe_r, color=colors[i], label=ii.upper())
        
        ax.set_xlabel('r [nm]')
        ax.set_ylabel('PMF [kJ/mol]')
        
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.tick_params(which='major', length=8)
        ax.tick_params(which='minor', length=2)
        
        ax.axvline(x=1.25)
        ax.axhline(y=0)
        
        
ax.legend(frameon=False)
        
colors=['blue','green','red']
for k, kk in enumerate(analyte):
    
    ax = plt.subplot(grid[1,k])
    
    for i, ii in enumerate(probe):
        # Create subplots
        
        
        print(kk,ii)
        #selected_fe = pmf_probe_analyte_zero_ref_global_minima
        selected_fe = pmf_probe_analyte_zero_ref_boundary
        
        
        fe_z = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['fe_z']
        xedges_z = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['xedges_z']
        
        ax.plot(xedges_z[:-1], fe_z, color=colors[i])
        
        ax.set_xlabel('z [nm]')
        ax.set_ylabel('PMF [kJ/mol]')
                  

#         ax.set_xlabel('r [nm]')
#         ax.set_ylabel('z [nm]')
#         ax.set_title(kk.upper()+':'+ii.upper())
        
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.tick_params(which='major', length=8)
        ax.tick_params(which='minor', length=2)
        
        ax.axvline(x=1.25)
        ax.axvline(x=-1.25)
        
        ax.axhline(y=0)
plt.savefig('pmf_r_z.png', dpi=400)

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# Set the figure size
fig = plt.figure(figsize=(18, 12))

# Adjust spacing
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.3)

# Create a 3x3 grid
grid = GridSpec(2, 3)

# probe= ['acd','bcd','gcd']
# probe = ['bcd']

# analyte = ['pfos']#,'TCAA']
# lbl = ['$\\beta$-CD']#,'$\\beta$-CD','$\\gamma$-CD']
lbl = probe
selected_block = 0 # entire data is used.

for k, kk in enumerate(analyte):
    for i, ii in enumerate(probe):
        
        # Create subplots
        ax = plt.subplot(grid[k, i])
        
        #selected_fe = pmf_probe_analyte_zero_ref_global_minima
        selected_fe = pmf_probe_analyte_zero_ref_boundary
       
        fe = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['fe']
        #fe = np.flip(fe,axis=0)
        xedges = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['xedges']
        yedges = selected_fe[ii+'_'+kk]['block'+str(selected_block)]['yedges']
        
        x, y = np.meshgrid(xedges, yedges)
        #contour_levels = np.linspace(fe.min(), 0, 15) 
        contour_levels = np.linspace(vmin, 0, 15) 
        
    
        ax.contour(x[:-1, :-1], y[:-1, :-1], fe, cmap='rainbow', levels=contour_levels,vmin=vmin,vmax=0)

        contour_plot = ax.contourf(x[:-1, :-1], y[:-1, :-1], fe, cmap='rainbow', levels=contour_levels, extend='max',vmin=vmin,vmax=0)
        cbar = fig.colorbar(contour_plot, ax=ax)
        cbar.set_label('PMF [kJ/mol]')  # You can change the label as needed

        # plt.show()
        #print(fe.max(), np.isinf(fe).any(), np.isnan(fe).any())
        
        # if i==2: #ldns
        # cbar_ax = fig.add_axes([0.25, -0.025, 0.5, 0.025]) #x-location, y-location, width, height of bar
        # # tick_bounds = [-70,-60,-50,-40,-30,-20,-10,0]
        # tick_bounds = list(np.linspace(-30,0,11))
        # # bounds = [-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0]
        # bounds = list(np.linspace(-30,0,11))
        # ncolors = get_cmap('rainbow')(np.linspace(0, 1, len(bounds)))
        # cmap = ListedColormap(ncolors)
        # norm = mpl.colors.BoundaryNorm(bounds, len(bounds))
        # cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap,
        #                                 boundaries=bounds, 
        #                                 norm=mpl.colors.Normalize(vmin=bounds[0], vmax=bounds[-1]),
        #                                 extend='min',
        #                                 spacing='uniform', orientation='horizontal')
        # #cb.ax.set_xticklabels(bounds)
        # cb.ax.set_xticks(tick_bounds)
        # cb.set_label('PMF [kJ/mol]')
        

        # fe_beyond_cutoff = []
        # print(xedges.min(), xedges.max(), yedges.min(), yedges.max())
        # for ki in range(xedges.shape[0]-1):
        #     for kj in range(yedges.shape[0]-1):
                
        #         if np.isnan(fe[kj,ki]) or np.isinf(fe[kj,ki]):
        #             continue

        #         if .99 <= xedges[ki] <= 1.01 and -0.94  <= yedges[kj] <= 0.95:
        #             fe_beyond_cutoff.append(fe[kj,ki])

        ax.set_xlabel('r [nm]')
        ax.set_ylabel('z [nm]')
        ax.set_title(kk.upper()+':'+lbl[i])
       
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.set_xlim(0.0,2.0001)
        ax.set_ylim(-2,2.0001)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.tick_params(which='major', length=8)
        ax.tick_params(which='minor', length=2)
        plt.savefig('pmf_2D.png',dpi=300)
        
# get dG values
V0 = 1661*1e-3
Vf = 32*np.pi/3
CP = V0/Vf
kb = 8.31446261815324*0.001 #kJ/molK
temperature = 300
kbt = kb*temperature  #kJ/mol
# probe= ['acd','bcd','gcd']
# analyte = ['pfos','sds','TCAA']
#choose number of blocks
blocks = 3
#define dimensions in r and z directions
#dimz = 100
#dimr = 100

pmf_probe_analyte_zero_ref_boundary = {} # collect pmfs
pmf_probe_analyte_zero_ref_global_minima = {} # collect pmfs
dimvec = [200]
cutvec = [0.8,1,1.2,1.4,1.6]

deltag_kb_probe_analyte_block = {}
cutoff1 = np.sqrt(3.5)
dgs = []
kbs = []
kbe = []
dge = []
iii = 0
for ii in probe:
    for kk in analyte:
        kbvec = []
        dgvec = []
        kberrs = []
        dgerrs = []
    
        for ll in range(len(dimvec)):
            for mm in range(len(cutvec)):
                dimz = dimvec[ll]
                cutoff = cutvec[mm]
                dimr = dimz
                pmf_probe_analyte_zero_ref_global_minima[ii+'_'+kk] = {}
                pmf_probe_analyte_zero_ref_boundary[ii+'_'+kk]  = {}

                deltag_kb_probe_analyte_block[ii+'_'+kk] = {}

                #create holders
                fep = []
                avgFE,sdFE = [],[]
                avgkb,sdkb = [],[]
                #load data 
                colvarFile = 'reweight/COLVAR_reweight' #simDir + systemDir + colvarFileName

                colvarData = pd.read_csv(colvarFile, comment="#", 
                                         delim_whitespace=True,names=["time", "angle", "r", 
                                                                      "baxis","mag", "pr", "pz", "bias"],
                                         skip_blank_lines=True)

                # discard initial few frames [frames where free energy is still being mapped]
                # subsample and randomize
                if kk=='pfos': 
                    # s = 5000
                    s = 0
                    # d = 50001 # changed it from 500000
                    d = -1 
                    colvarD = colvarData[s:d].sample(frac = 1)
                else:
                    # s = 10000
                    s = 0
                    # d = 100000 # changed it from 1000000
                    d = -1 
                    colvarD = colvarData[s:d:2].sample(frac = 1)        
                jmp = len(colvarD)//blocks

                lockb = []
                locdg = []
                for jj in range(blocks): 
                    colvarData_init = colvarD[jmp*jj:jmp*(jj+1)] 
                    print(jmp*jj,jmp*(jj+1))


                    colvarData_2 = colvarData_init[(colvarData_init['r'] < cutoff1)]
                    colvarData_3 = colvarData_2[(colvarData_2['baxis'] < cutoff1)]
                    colvarData = colvarData_3[(colvarData_3['baxis'] >= -cutoff1)]
                    ############################################################################################################
                    #make a 2D heat map
                    bias = colvarData['bias'] 
                    r = colvarData['r']
                    z = colvarData['baxis']
                    
                    # calculate weights of each frame in colvarData - CHECK THIS APPROACH (SUMMING THE BIASES)
                    weights = np.exp((bias) /kbt)        
                    xedges = np.linspace(0,cutoff1,dimr)
                    yedges = np.linspace(-cutoff1,cutoff1,dimz)
                    #print(xedges,yedges)
                    dr,dz = abs(xedges[0]-xedges[1]),abs(yedges[0]-yedges[1])
                    # get unbiased probabilities
                    fe, xedges,yedges = np.histogram2d(r,z,range=[[0,cutoff1],[-cutoff1,cutoff1]],
                                                       bins=(xedges,yedges),
                                                       weights=weights,density=True)

                    fe = fe.T
                    feDf = -kbt*np.log(fe) # Compute free energy

                    ######## correction
                    #correction
                    for ki in range(xedges.shape[0]-1):
                        for kj in range(yedges.shape[0]-1):
        #                     if xedges[ki]**2 + np.abs(yedges[kj])**2 > 1.25: # 2.0:
                   
                            if np.isnan(feDf[kj,ki]) or np.isinf(feDf[kj,ki]):
                                continue
    
                            s = np.sqrt(xedges[ki]**2+yedges[kj]**2)
                            feDf[ki,kj] += preFactor*np.log(s)
####################################

            
                    # Find values of free energy on the boundary of the surface
                    avg_fe_beyond_cutoff=[]
                    for ki in range(xedges.shape[0]-1):
                        for kj in range(yedges.shape[0]-1):

                            if np.isnan(feDf[kj,ki]) or np.isinf(feDf[kj,ki]):
                                continue
                           # print(bool(xedges[ki]>0.8),ki, xedges[ki])
                            #if 1.48 <= xedges[ki] < 1.5 and -1.48**2 <= yedges[kj]**2 < 1.5**2:
                            if xedges[ki] > cutoff or yedges[kj] < -cutoff or yedges[kj] > cutoff:
                                avg_fe_beyond_cutoff.append(feDf[kj,ki])
                        



                    avg_fe_beyond_cutoff = np.array(avg_fe_beyond_cutoff)
                 #   print(avg_fe_beyond_cutoff,len(avg_fe_beyond_cutoff))
                    cutmean = avg_fe_beyond_cutoff.mean()
                    cutstd = avg_fe_beyond_cutoff.std()
                    feDf_zero_ref_boundary = feDf - avg_fe_beyond_cutoff.mean() #collect_fe_k_max # UNCOMMENT THIS FOR SETTING BOUNDARY FREE ENERGY TO ZERO
                    _block_data_boundary = {'fe': feDf_zero_ref_boundary,
                                            'xedges': xedges,
                                            'yedges': yedges}

                    pmf_probe_analyte_zero_ref_boundary[ii+'_'+kk]['block'+str(jj)] = _block_data_boundary

                    selected_fe = pmf_probe_analyte_zero_ref_boundary[ii+'_'+kk]['block'+str(jj)]['fe']


                    ##Compute deltaG and kb
                    ##Also visualize the cutoff area for bounded region
                    print('Computing free energies ...')
                    dgh=0
                    dgh_unbound = 0
                    #simulation_volume = 0 # term 1/C in Kb. 1/(simulation concentration)

                    for ki in range(xedges.shape[0]-1):
                        for kj in range(yedges.shape[0]-1):

                           
                            # integrate below the upper cutoff of the above region.
                            if xedges[ki] < cutoff and (yedges[kj] > -cutoff and yedges[kj] < cutoff):
                
                                # convert selected_fe to probability (redundant step; done previously for visualization)
                                #dgh += np.exp(-selected_fe/kbt)[kj,ki]*dr*dz*xedges[ki]
                                dgh += np.exp(-selected_fe[kj,ki]/kbt)*dr*dz*xedges[ki]

                    dg = -kbt*np.log(2*np.pi*dgh*(1/V0))
                    kb = (1/V0)*2*np.pi*dgh

                    print(f'{ii.upper()}:{kk.upper()} Block:{jj} Kb:{kb} dG:{dg} kJ/mol cutmean:{cutmean} cutstd:{cutstd}')
                    print(f'{ii.upper()}:{kk.upper()}', f' dimz: {dimz} cutoff: {cutoff}')
                    deltag_kb_probe_analyte_block[ii+'_'+kk]['block'+str(jj)] = {'deltag': dg, 'kb': kb}
                    lockb.append(kb)
                    locdg.append(dg)

                kbvec.append(np.mean(lockb))
                kberrs.append(np.std(lockb))
                dgvec.append(np.mean(locdg))
                dgerrs.append(np.std(locdg))
            kbs.append(kbvec)
            kbe.append(kberrs)
            dgs.append(dgvec)
            dge.append(dgerrs)

outdata = [kbs, kbe, dgs, dge]
print(dgs, dge)
np.save('kb_dg.npy',outdata)
