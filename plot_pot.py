#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pyProbeParticle.GridUtils      as GU


# =============== arguments definition

dz=.1
df = np.load('FFLJ_z.npy')
lvec = np.load('FFLJ_vec.npy')

print lvec

namez = 'z_el_pot'

for k in range(df.shape[0]):
    dff = np.array(df[k,:,:]).copy()
    kk = k*lvec[3,2]/df.shape[0]
    print kk
    name_file=namez+'_'+str(kk)+'.dat'
    # ploting part here:
    plt.figure( figsize=( lvec[1,0]/10 , lvec[2,1]/10 ) )
    plt.imshow( dff, origin='image')#, extent=extent , cmap='gray')
    plt.xlabel(r' Tip_x $\AA$')
    plt.ylabel(r' Tip_y $\AA$')
    #plt.title(name_plot_df)
    plt.savefig( namez+'_s_'+str(kk)+'.png' , bbox_inches='tight' )
    #plt.show()



