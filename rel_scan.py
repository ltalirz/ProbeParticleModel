#!/usr/bin/python -u

import os
import numpy as np
#import matplotlib.pyplot as plt
import sys

'''
import basUtils
import elements
import GridUtils as GU
import ProbeParticle      as PP;    PPU = PP.PPU;
'''

import pyProbeParticle                as PPU     
import pyProbeParticle.GridUtils      as GU
import pyProbeParticle.core           as PPC
import pyProbeParticle.HighLevel      as PPH

#import PPPlot 		# we do not want to make it dempendent on matplotlib

print "Amplitude ", PPU.params['Amplitude']

# =============== arguments definition

from optparse import OptionParser
parser = OptionParser()
parser.add_option( "-k",       action="store", type="float", help="tip stiffenss [N/m]" )
parser.add_option( "--krange", action="store", type="float", help="tip stiffenss range (min,max,n) [N/m]", nargs=3)
parser.add_option( "-q",       action="store", type="float", help="tip charge [e]" )
parser.add_option( "--qrange", action="store", type="float", help="tip charge range (min,max,n) [e]", nargs=3)

parser.add_option( "--pos",    action="store_true", default=False, help="save probe particle positions" )
parser.add_option( "--disp",    action="store_true", default=False, help="save probe particle displacements")

(options, args) = parser.parse_args()
opt_dict = vars(options)

# =============== Setup

PPU.loadParams( 'params.ini' )

print opt_dict
# Ks
if opt_dict['krange'] is not None:
	Ks = np.linspace( opt_dict['krange'][0], opt_dict['krange'][1], opt_dict['krange'][2] )
elif opt_dict['k'] is not None:
	Ks = [ opt_dict['k'] ]
else:
	Ks = [ PPU.params['stiffness'][0] ]
# Qs

charged_system=False
if opt_dict['qrange'] is not None:
	Qs = np.linspace( opt_dict['qrange'][0], opt_dict['qrange'][1], opt_dict['qrange'][2] )
elif opt_dict['q'] is not None:
	Qs = [ opt_dict['q'] ]
else:
	Qs = [ PPU.params['charge'] ]

for iq,Q in enumerate(Qs):
        if ( abs(Q) > 1e-7):
                charged_system=True

print "Ks   =", Ks 
print "Qs   =", Qs 

print " ============= RUN  "

if ( charged_system == True):
        print " load Electrostatic Force-field "
        FFel, lvec = GU.loadVecFieldNpy( "FFel" )

print " load Lenard-Jones Force-field "
FFLJ, lvec = GU.loadVecFieldNpy( "FFLJ" )
PPU.lvec2params( lvec )
PPC.setFF( FFLJ )

xTips,yTips,zTips,lvecScan = PPU.prepareScanGrids( )

for iq,Q in enumerate( Qs ):
        if ( charged_system == True):
	        FF = FFLJ + FFel * Q
        else:
                FF = FFLJ
	PPC.setFF_Pointer( FF )
	for ik,K in enumerate( Ks ):
		dirname = "Q%1.2fK%1.2f" %(Q,K)
		print " relaxed_scan for ", dirname
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
		PPC.setTip( kSpring = np.array((K,K,0.0))/-PPU.eVA_Nm )
		fzs,PPpos = PPH.relaxedScan3D( xTips, yTips, zTips )
		GU.saveNpy( dirname+'/OutFz', fzs, lvecScan)

#                print "SHAPE", PPpos.shape, xTips.shape, yTips.shape, zTips.shape
                if opt_dict['disp']:
                    PPdisp=PPpos.copy()
                    nx=PPdisp.shape[2]
                    ny=PPdisp.shape[1]
                    nz=PPdisp.shape[0]
                    test=np.meshgrid(xTips,yTips,zTips)
#                    print "TEST SHAPE", np.array(test).shape
#                    print nx,ny,nz
                    i=0
                    while i<nx:
                        j=0
                        while j<ny:
                            k=0
                            while k<nz:
                                PPdisp[k][j][i]-=np.array([xTips[i],xTips[j],zTips[k]])+ np.array([PPU.params['r0Probe'][0],PPU.params['r0Probe'][1],-PPU.params['r0Probe'][2]])
                                k+=1
                            j+=1
                        i+=1
		    
                    GU.saveVecFieldNpy( dirname+'/PPdisp', PPdisp, lvec )
                    
		if opt_dict['pos']:
			GU.saveVecFieldNpy( dirname+'/PPpos', PPpos, lvec )
			GU.saveVecFieldNpy( dirname+'/PPpos', PPpos, lvecScan ) ### !!!! Changed, due to the need to have the right position to see non-relaxed scan !!!! ###
		# the rest is done in plot_results.py; For df, go to plot_results.py

print " ***** ALL DONE ***** "

