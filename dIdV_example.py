#!/usr/bin/python

import os
import numpy as np

import pyProbeParticle.GridUtils      as GU
import pyProbeParticle.ProbeSTM       as PS

import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt

# --- specification of paths

path=''
path_pos='Q0.00K0.50/'
path_df = path_pos+'Amp0.40/'

# --- specification of PBC and electronics

pbc=(0,0)
lvs = None
#lvs = np.array([[15., 0.],[0.,15.]]
WorkFunction = 5.0 #more or less standart.
fermi=None	# the Fermi from GPAW ; means: -5.04612664712 eV
orbs= 'sp'	# only 'sp' works now
cut_min=-1.0	# HOMO -0.88 bellow the Fermi Level
cut_max=+1.0	# LUMO -0.88 above the Fermi Level
cut_at=-1	# All atoms of the molecule
eta = 0.01	# very low, to pronounce the single orbitals only
# WF_decay=1.0  # for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)

#eigEn, coefs, Ratin  = PS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=None, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_at=[0,1.]);
Ratin = PS.read_fire_atoms(path+'crazy_mol.xyz',pbc=pbc,cut_at=cut_at)
eigEn, coefs = PS.read_fire_coeffs(name = path+'phik_example_', fermi=-5.04612664712, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,cut_at=cut_at);

tip_r1, lvec = GU.loadVecFieldNpy( path_pos+"PPpos" )

dz=0.1
dx=dy =0.1

#ppheight = -4.0  # not needed any more, position of PP in LvecScan
xl = lvec[1,0]
yl = lvec[2,1]
zl = lvec[3,2]
extent = (lvec[0,0],lvec[0,0]+xl,lvec[0,1],lvec[0,1]+yl)

tip_r2 = PS.mkSpaceGrid(lvec[0,0],lvec[0,0]+xl,dx,lvec[0,1],lvec[0,1]+yl,dy,lvec[0,2],lvec[0,2]+zl,dz)

Voltages=[-0.88,+0.88]
namez=['HOMO','LUMO']

#Voltages=np.arange(-1.0,+1.0+0.01,0.1) # this part is important for scans over slabs at different voltages
#namez = []
#for V in Voltages:
#    namez.append(str(round(V,1)))

df, lvec2 = GU.loadNpy(path_df+'df')

for WorkFunction in [WorkFunction]:
    i=0;
    for V in Voltages:
	for eta in [eta]:
	    current0 = PS.dIdV( V, WorkFunction, eta, eigEn,  tip_r2 , Ratin, coefs, orbs=orbs , s=1.0, px=0.0, py=0.0, pz = 0.0)
	    np.save('didv_s-fixed_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'.npy',current0)
	    current1 = PS.dIdV( V, WorkFunction, eta, eigEn,  tip_r1 , Ratin, coefs, orbs=orbs , s=1.0, px=0.0, py=0.0, pz = 0.0)
	    np.save('didv_s_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'.npy',current1)
	    current2 = PS.dIdV( V, WorkFunction, eta, eigEn,  tip_r1 , Ratin, coefs, orbs=orbs , s=0.0, px=1.0, py=1.0, pz = 0.0)
	    np.save('didv_pxy_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'.npy',current2)
	    #current3 = PS.dIdV( V, WorkFunction, eta, eigEn,  tip_r1 , Ratin, coefs, orbs=orbs , s=0.0, px=0.0, py=0.0, pz = 1.0)
	    #np.save('didv_pz_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'.npy',current3)
	    #current3 = PS.STM( V, nV, WorkFunction, eta, eigEn,  tip_r1 , Ratin, coefs, orbs=orbs , pz=0.0 ,s=0.0, px=0.5, py=0.5, WF_decay=WF_decay)
	    print " plotting "
	
	    for k in range(df.shape[0]):
		dff = np.array(df[k,:,:]).copy()
		curr0 = np.array(current0[k,:,:]).copy()
		curr1 = np.array(current1[k,:,:]).copy()
		curr2 = np.array(current2[k,:,:]).copy()
		#curr3 = np.array(current3[k,:,:]).copy()
		
		name_file='didV-'+namez[i]+'_%03d.dat' %k
		name_plot_df='height:%03dA; df [Hz]' %k
		name_plot0=namez[i]+';height:%03dA; dIdV [G0] s-fixed-tip' %k
		name_plot1=namez[i]+';height:%03dA; dIdV [G0] s-tip' %k
		name_plot2=namez[i]+';height:%03dA; dIdV [G0] pxy-tip' %k
		#name_plot3=namez[i]+';height:%03dA; dIdV [G0] pz-tip' %k

		# ploting part here:
		plt.figure( figsize=(1.5* xl , 1.5*yl/4 ) )
		plt.subplot(1,4,1)
		plt.imshow( dff, origin='image', extent=extent , cmap='gray')
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot_df)

		# ploting part here:
		plt.subplot(1,4,2)
		plt.imshow( curr0, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot0)

		# ploting part here:
		plt.subplot(1,4,3)
		plt.imshow( curr1, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot1)

		# ploting part here:
		plt.subplot(1,4,4)
		plt.imshow( curr2, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot2)

		plt.savefig( 'didv_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
		#plt.show()
		#
		#tmp_curr=curr.flatten()
		#out_curr=np.zeros((len(tmp_curr),3))
		#out_curr[:,0]=tip_r[k,:,:,0].flatten()
		#out_curr[:,1]=tip_r[k,:,:,1].flatten()
		#out_curr[:,2]=tmp_curr.copy()
		#f=open(name_file,'w')
		#print >> f, "WSxM file copyright Nanotec Electronica"
		#print >> f, "WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al."
		#print >> f, "X[A]  Y[A]  dIdV[G0]"
		#print >> f, ""
 		#np.savetxt(f, out_curr)
		#f.close()
		#
	
	plt.show()
	i = i+1
	

print 
print
print "Done"