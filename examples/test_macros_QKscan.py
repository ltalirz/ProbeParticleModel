#!/usr/bin/python

#import matplotlib
#matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import os
import numpy as np
import matplotlib.pyplot as plt
import elements
#import XSFutils
import basUtils


from memory_profiler import profile




print " # ========== make & load  ProbeParticle C++ library " 

def makeclean( ):
	import os
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".so") ]
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".o") ]
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".pyc") ]

makeclean( )  # force to recompile 

import ProbeParticle as PP

print " >> WARNING!!! OVEWRITING SETTINGS by params.ini  "

PP.loadParams( 'params_carbox.ini' )

#Fx,lvec,nDim,head=XSFutils.loadXSF('Fx.xsf')
#Fy,lvec,nDim,head=XSFutils.loadXSF('Fy.xsf')
#Fz,lvec,nDim,head=XSFutils.loadXSF('Fz.xsf')

Fx,lvec,nDim,head=PP.loadXSF('Fx.xsf')
Fy,lvec,nDim,head=PP.loadXSF('Fy.xsf')
Fz,lvec,nDim,head=PP.loadXSF('Fz.xsf')

PP.params['gridA'] = lvec[ 1,:  ].copy()
PP.params['gridB'] = lvec[ 2,:  ].copy()
PP.params['gridC'] = lvec[ 3,:  ].copy()
PP.params['gridN'] = nDim.copy()

FFLJ  = np.zeros( (nDim[0],nDim[1],nDim[2],3) )
FFel  = np.zeros( np.shape( FFLJ ) )

FFel[:,:,:,0]=Fx
FFel[:,:,:,1]=Fy
FFel[:,:,:,2]=Fz

del Fx; del Fy; del Fz


cell =np.array([
PP.params['gridA'],
PP.params['gridB'],
PP.params['gridC'],
]).copy() 
gridN = PP.params['gridN']

print "cell", cell

PP.setFF( FFLJ, cell  )

print " # ============ define atoms "

atoms    = basUtils.loadAtoms('carboxylics.bas', elements.ELEMENT_DICT )
Rs       = np.array([atoms[1],atoms[2],atoms[3]]);  
iZs      = np.array( atoms[0])

if not PP.params['PBC' ]:
	print " NO PBC => autoGeom "
	PP.autoGeom( Rs, shiftXY=True,  fitCell=True,  border=3.0 )
	print " NO PBC => params[ 'gridA'   ] ", PP.params[ 'gridA' ] 
	print " NO PBC => params[ 'gridB'   ] ", PP.params[ 'gridB'   ]
	print " NO PBC => params[ 'gridC'   ] ", PP.params[ 'gridC'   ]
	print " NO PBC => params[ 'scanMin' ] ", PP.params[ 'scanMin' ]
	print " NO PBC => params[ 'scanMax' ] ", PP.params[ 'scanMax' ]

Rs = np.transpose( Rs, (1,0) ).copy() 
Qs = np.array( atoms[4] )

if PP.params['PBC' ]:
	iZs,Rs,Qs = PP.PBCAtoms( iZs, Rs, Qs, avec=PP.params['gridA'], bvec=PP.params['gridB'] )

print "shape( Rs )", np.shape( Rs ); 

print " # ============ define Scan and allocate arrays   - do this before simulation, in case it will crash "

dz    = PP.params['scanStep'][2]
zTips = np.arange( PP.params['scanMin'][2], PP.params['scanMax'][2]+0.00001, dz )[::-1];


PP.setTip()

xTips  = np.arange( PP.params['scanMin'][0], PP.params['scanMax'][0]+0.00001, 0.1 )
yTips  = np.arange( PP.params['scanMin'][1], PP.params['scanMax'][1]+0.00001, 0.1 )
extent=( xTips[0], xTips[-1], yTips[0], yTips[-1] )





lvecScan =np.array([
PP.params['scanMin'],
[        PP.params['scanMax'][0],0.0,0.0],
[0.0,    PP.params['scanMax'][1],0.0    ],
[0.0,0.0,PP.params['scanMax'][2]        ]
]).copy() 

headScan='''
ATOMS
 1   0.0   0.0   0.0

BEGIN_BLOCK_DATAGRID_3D                        
   some_datagrid      
   BEGIN_DATAGRID_3D 
'''


nslice = 10;

#quit()

# ==============================================
#   The costly part of simulation starts here
# ==============================================

print " # =========== Sample LenardJones "

computeLJ( xsfLJ = False )

#Ks = [ 0.125, 0.25, 0.5, 1.0 ]
#Qs = [ -0.4, -0.3, -0.2, -0.1, 0.0, +0.1, +0.2, +0.3, +0.4 ]


#@profile
def main():
	for iq,Q in enumerate( Qs ):
		FF = FFLJ + FFel * Q
		PP.setFF_Pointer( FF )
		for ik,K in enumerate( Ks ):
			dirname = "Q%1.2fK%1.2f" %(Q,K)
			os.makedirs( dirname )
			PP.setTip( kSpring = np.array((K,K,0.0))/-PP.eVA_Nm )
			fzs = relaxedScan3D( xTips, yTips, zTips )
			PP.saveXSF( dirname+'/OutFz.xsf', headScan, lvecScan, fzs )
			dfs = PP.Fz2df( fzs, dz = dz, k0 = PP.params['kCantilever'], f0=PP.params['f0Cantilever'], n=int(PP.params['Amplitude']/dz) )
			plotImages( dirname+"/df", dfs, slices = range( 0, len(dfs) ) )

main()

print " ***** ALL DONE ***** "

#plt.show()




