#!/usr/bin/python

#import matplotlib
#matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import os
import numpy as np
import matplotlib.pyplot as plt
import elements
#import XSFutils
import basUtils

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

FF   = np.zeros( (nDim[0],nDim[1],nDim[2],3) )
FFel = np.zeros( np.shape( FF ) )

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

PP.setFF( FF, cell  )

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

Rs     = np.transpose( Rs, (1,0) ).copy() 

Qs = np.array( atoms[4] )

if PP.params['PBC' ]:
	iZs,Rs,Qs = PP.PBCAtoms( iZs, Rs, Qs, avec=PP.params['gridA'], bvec=PP.params['gridB'] )

print "shape( Rs )", np.shape( Rs ); 
#print "Rs : ",Rs

print " # ============ define Scan and allocate arrays   - do this before simulation, in case it will crash "

dz    = PP.params['scanStep'][2]
zTips = np.arange( PP.params['scanMin'][2], PP.params['scanMax'][2]+0.00001, dz )[::-1];
ntips = len(zTips); 
print " zTips : ",zTips
rTips = np.zeros((ntips,3))
rs    = np.zeros((ntips,3))
fs    = np.zeros((ntips,3))

rTips[:,0] = 1.0
rTips[:,1] = 1.0
rTips[:,2] = zTips 

PP.setTip()

xTips  = np.arange( PP.params['scanMin'][0], PP.params['scanMax'][0]+0.00001, 0.1 )
yTips  = np.arange( PP.params['scanMin'][1], PP.params['scanMax'][1]+0.00001, 0.1 )
extent=( xTips[0], xTips[-1], yTips[0], yTips[-1] )
fzs    = np.zeros(( len(zTips), len(yTips ), len(xTips ) ));

nslice = 10;

#quit()

# ==============================================
#   The costly part of simulation starts here
# ==============================================

print " # =========== Sample LenardJones "

#xsfLJ       = True
xsfLJ       = False
recomputeLJ = True

if (     xsfLJ  and os.path.isfile('FFLJ_x.xsf')):
	recomputeLJ = False
if ((not xsfLJ) and os.path.isfile('FFLJ_y.npy')):
	recomputeLJ = False

if recomputeLJ:
	FFparams = PP.loadSpecies        ( 'atomtypes.ini'  )
	C6,C12   = PP.getAtomsLJ( PP.params['probeType'], iZs, FFparams )
	PP.setFF( FF, cell  )
	PP.setFF_Pointer( FF )
	PP.getLenardJonesFF( Rs, C6, C12 )
	if xsfLJ:
		PP.saveXSF('FFLJ_x.xsf', head, lvec, FF[:,:,:,0] )
		PP.saveXSF('FFLJ_y.xsf', head, lvec, FF[:,:,:,1] )
		PP.saveXSF('FFLJ_z.xsf', head, lvec, FF[:,:,:,2] )
	else:
		np.save('FFLJ_x.npy', FF[:,:,:,0] )
		np.save('FFLJ_y.npy', FF[:,:,:,1] )
		np.save('FFLJ_z.npy', FF[:,:,:,2] )
else:
	if xsfLJ:
		FF[:,:,:,0],lvec,nDim,head=PP.loadXSF('FFLJ_x.xsf')
		FF[:,:,:,1],lvec,nDim,head=PP.loadXSF('FFLJ_y.xsf')
		FF[:,:,:,2],lvec,nDim,head=PP.loadXSF('FFLJ_z.xsf')
	else:
		FF[:,:,:,0] = np.load('FFLJ_x.npy' )
		FF[:,:,:,1] = np.load('FFLJ_y.npy' )
		FF[:,:,:,2] = np.load('FFLJ_z.npy' )

# ======= plot 
#plt.figure(figsize=( 5*nslice,5 )); plt.title( ' FF LJ ' )
#for i in range(nslice):
#	plt.subplot( 1, nslice, i+1 )
#	plt.imshow( FF[i,:,:,2], origin='image', interpolation='nearest' )


withElectrostatics = ( abs( PP.params['charge'] )>0.001 )
if withElectrostatics: 
	print " # =========== Sample Coulomb "
	FF += FFel*PP.params['charge']
	PP.setFF_Pointer( FF )

del FFel

#plt.figure(figsize=( 5*nslice,5 )); plt.title( ' FF total ' )
#for i in range(nslice):
#	plt.subplot( 1, nslice, i+1 )
#	plt.imshow( FF[i,:,:,2], origin='image', interpolation='nearest' )

print " # ============  Relaxed Scan 3D "

for ix,x in enumerate( xTips  ):
	print "relax ix:", ix
	rTips[:,0] = x
	for iy,y in enumerate( yTips  ):
		rTips[:,1] = y
		itrav = PP.relaxTipStroke( rTips, rs, fs ) / float( len(zTips) )
		fzs[:,iy,ix] = fs[:,2].copy()
		#print itrav
		#if itrav > 100:
		#	print " bad convergence > %i iterations per pixel " % itrav
		#	print " exiting "
		#	break
		

print " # ============  convert Fz -> df "

dfs = PP.Fz2df( fzs, dz = dz, k0 = PP.params['kCantilever'], f0=PP.params['f0Cantilever'], n=int(PP.params['Amplitude']/dz) )

print " # ============  Plot Relaxed Scan 3D "

#slices = range( PP.params['plotSliceFrom'], PP.params['plotSliceTo'], PP.params['plotSliceBy'] )
#print "plotSliceFrom, plotSliceTo, plotSliceBy : ", PP.params['plotSliceFrom'], PP.params['plotSliceTo'], PP.params['plotSliceBy']
#print slices 
#nslice = len( slices )

slices = range( 0, len(dfs) )

for ii,i in enumerate(slices):
	print " plotting ", i
	plt.figure( figsize=( 10,10 ) )
	plt.imshow( dfs[i], origin='image', interpolation=PP.params['imageInterpolation'], cmap=PP.params['colorscale'], extent=extent )
	z = zTips[i] - PP.params['moleculeShift' ][2]
	plt.colorbar();
	plt.xlabel(r' Tip_x $\AA$')
	plt.ylabel(r' Tip_y $\AA$')
	plt.title( r"df Tip_z = %2.2f $\AA$" %z  )
	plt.savefig( 'df_%3.3i.png' %i, bbox_inches='tight' )



print " ***** ALL DONE ***** "


plt.show()




