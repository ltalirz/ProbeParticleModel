#!/usr/bin/python

import os
import numpy as np
from   ctypes import c_int, c_double, c_char_p
import ctypes
import basUtils as bU
import elements

import cpp_utils

# this library has functions for reading STM coefficients and make a grid for non-relaxed 3D scan

# ==============================
# ============================== Pure python functions
# ==============================

LIB_PATH = os.path.dirname( os.path.realpath(__file__) )
print " ProbeParticle Library DIR = ", LIB_PATH

def mkSpaceGrid(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz):
	'''
	mkSpaceGridsxmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz):
	Give rectangular grid along the main cartesian axes for non-relaxed dI/dV or STM - 4D grid of xyz coordinates.
		'''
	h = np.mgrid[xmin:xmax+0.0001:dx,ymin:ymax+0.0001:dy,zmin:zmax+0.0001:dz]
	f = np.transpose(h)
	sh = f.shape
	print "Grid has dimensios: ", sh
	return f;	#, sh;

def	read_AIMS_all(name = 'KS_eigenvectors.band_1.kpt_1.out', geom='geometry.in', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[],header=False):
	'''
	read_AIMS_all(name = 'KS_eigenvectors.band_1.kpt_1.out', geom='geometry.in', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
	read eigen energies, coffecients, Fermi Level and geometry  from the GPAW  *.gpw file.
	orbs - only 'sp' works 	can read only sp structure of valence orbitals (hydrogens_has to be at the end !!!!)
	Fermi - set to zero by AIMS itself
	pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
	imaginary = False (other options for future k-points dependency
	cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
	cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
	lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first four atoms by 0.5
	'''
	assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
	assert (imaginary == False), "sorry imaginary version is under development" 	
	print "reading FHI-AIMS LCAO coefficients for basis: ",orbs	
	from ase import Atoms; from ase.io import read;
	slab=read(geom)
	num_at = slab.get_number_of_atoms()
	at_num = slab.get_atomic_numbers()
	if not ((cut_at == -1)or(cut_at == num_at)):
		print "cutting attoms"
		Ratin = slab.get_positions()[:cut_at,:]
		num_at = cut_at
		at_num = at_num[:cut_at]
	else:
		#print "NOT! cutting attoms"
		Ratin = slab.get_positions()
	if (pbc != ((0,0)or(0.,0.))):
		print "Applying PBC"
		atoms = Atoms('H%d' % len(Ratin), positions=Ratin)
		cell = slab.get_cell()
		atoms.set_cell(cell)
		if (pbc == (0.5,0.5)):
			atoms *= (2,2,1)
		else:
			atoms *= ( (int(2*pbc[0])+1),(int(2*pbc[1])+1),1 )
			atoms.translate( [( -int(pbc[0])*cell[0,0]-int(pbc[1])*cell[1,0] ) , (-int(pbc[0])*cell[0,1]-int(pbc[1])*cell[1,1]),0] )
		#from ase.visualize import view
		#view(atoms)
		Ratin = atoms.get_positions()
		print " Number of atoms after PBC: ", len(Ratin)
	print "geometry read"

	filein = open(name )
	for i in range(4):
		tmp=filein.readline()
	if header :
		tmp=filein.readline()
	del tmp;
	pre_eig = filein.readline().split()
	filein.close()
	pre_eig=np.delete(pre_eig,[0,1,2],0)
	n_bands = len(pre_eig)
	eig = np.zeros(n_bands)
	for i in range(n_bands):
		eig[i] = float(pre_eig[i])
	if (fermi!=None):
		eig += -fermi
	del pre_eig;
	print "Fermi Level set to:", fermi, "AIMS set Fermi to zero by itself"
	n_min = -1
	n_max = -1
	j = 0
	for i in eig:
		if (i < cut_min):
			n_min = j
		if (i < cut_max):
			n_max = j
		j += 1
	if (n_min and n_max != -1):
		assert (n_min < n_max), "no orbitals left for dI/dV"
	print "eigenenergies read"
	
	if header :
		tmp = np.genfromtxt(name,skip_header=6, usecols=(1,2,3,4,5),dtype=None)
	else:
		tmp = np.genfromtxt(name,skip_header=5, usecols=(1,2,3,4,5),dtype=None)
	if (orbs == 'sp'):
		orb_pos=np.zeros((num_at,4), dtype=np.int)
		orb_sign=np.zeros((num_at,4), dtype=np.int)
	orb_pos += -1
	el = elements.ELEMENTS
	for j in range(num_at):
		Z = at_num[j];
		per = el[Z][2]
		temp=int((np.mod(2,2)-0.5)*2)	# phase of radial function in long distance for l=0: if n even - +1, if odd - -1
		orb_sign[j]=[temp,-1*temp,-1*temp,temp]		# {1, 1, 1, -1};(*Dont change, means - +s, +py +pz -px*) but l=1 has opposite phase than l=0 ==>  sign[s]*{1, -1, -1, 1};
	for i in range(len(tmp)):
		for j in range(num_at):
			Z = at_num[j];
			per = el[Z][2]
			if ((tmp[i][0]==j+1)and(tmp[i][1]=='atomic')):
				if (tmp[i][2]==per):
					if 	(tmp[i][3]=='s'):
						orb_pos[j,0]=i
					elif (tmp[i][3]=='p'):
						if (tmp[i][4]==-1):
							orb_pos[j,1]=i
						elif (tmp[i][4]==0):
							orb_pos[j,2]=i
						elif (tmp[i][4]==1):
							orb_pos[j,3]=i
	del tmp; del temp;
	if header :
		tmp = np.genfromtxt(name,skip_header=6, usecols=tuple(xrange(6, n_bands*2+6, 2))) #tmp = np.genfromtxt(name,skip_header=5)#, usecols=(6,))
	else:
		tmp = np.genfromtxt(name,skip_header=5, usecols=tuple(xrange(6, n_bands*2+6, 2))) #tmp = np.genfromtxt(name,skip_header=5)#, usecols=(6,))
	if (orbs == 'sp'):
		coef = np.zeros((n_bands,num_at,4))
		for j in range(num_at):
			for l in range(4):
				if (orb_pos[j,l]!=-1):
					coef[:,j,l] = tmp[orb_pos[j,l]]
					coef[:,j,l] *= orb_sign[j,l]
	del tmp;
	if (lower_atoms != []):
		print 'lowering atoms hoppings for atoms:', lower_atoms
		i_coef = 0;
		for j in lower_atoms:
			coef[:,j,:] *= lower_coefs[i_coef]
			i_coef +=1
	coeff = coef.flatten()
	coeffs = coeff.reshape((len(eig),num_at*4)).copy()
	del coef; del coeff;
	#now removing non-effective orbitals:
	if (n_max != -1):
		print "cutting orbitals with too high eigen energies"
		coeffs = np.delete(coeffs,range(n_max+1,len(eig)),0)
		eig = np.delete(eig,range(n_max+1,len(eig)),0)
	if (n_min != -1):
		print "cutting orbitals with too low eigen energies"
		coeffs = np.delete(coeffs,range(n_min+1),0)
		eig = np.delete(eig,range(n_min+1),0)
	# applying PBC
	if ((pbc != (0,0))or(pbc != (0.0,0.0))) :
		print "applying pbc"
		coeff =np.repeat(coeffs,int(pbc[0]*2+1)*int(pbc[1]*2+1),0).flatten()
		num_at *=int(pbc[0]*2+1)*int(pbc[1]*2+1)
		coeffs = coeff.reshape((len(eig),num_at*4)) if (orbs == 'sp') else coeff.reshape((len(eig),num_at*9));
	print "coefficients read; now we are getting geometry, cut_at:", cut_at
	print "coefficients read"
	return eig.copy(), coeffs.copy(), Ratin.copy();

def	read_GPAW_all(name = 'OUTPUT.gpw', fermi = None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[] ):
	'''
	read_GPAW_all(name = 'OUTPUT.gpw', fermi = None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
	read eigen energies, coffecients, Fermi Level and geometry  from the GPAW  *.gpw file.
	If fermi = None then Fermi comes from the GPAW calculation
	orbs - only 'sp' works 	can read only sp structure of valence orbitals (hydrogens_has to be at the end !!!!)
	pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
	imaginary = False (other options for future k-points dependency
	cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
	cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
	lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first four atoms by 0.5
	'''
	assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
	assert (imaginary == False), "sorry imaginary version is under development" 	
	print "reading GPAW LCAO coefficients for basis: ",orbs	
	from ase import Atoms
	from gpaw import GPAW
	calc = GPAW(name)
	slab = calc.get_atoms()
	num_at = slab.get_number_of_atoms()
	n_bands = calc.get_number_of_bands()
	eig = calc.get_eigenvalues(kpt=0, spin=0, broadcast=True)
	at_num = slab.get_atomic_numbers()
	efermi = calc.get_fermi_level()
	if (fermi == None):
		fermi = efermi
	else:
		fermi += efermi
	print "Fermi Level: ", fermi, " eV"
	eig -=fermi
	n_min = -1
	n_max = -1
	j = 0
	for i in eig:
		if (i < cut_min):
			n_min = j
		if (i < cut_max):
			n_max = j
		j += 1
	if (n_min and n_max != -1):
		assert (n_min < n_max), "no orbitals left for dI/dV"
	coef = np.zeros((n_bands,num_at,4))
	for i in range(n_bands):
		h=0
		for j in range(num_at):
			coef[i,j,0] = calc.wfs.kpt_u[0].C_nM[i,h]
			coef[i,j,1] = calc.wfs.kpt_u[0].C_nM[i,h+1]
			coef[i,j,2] = calc.wfs.kpt_u[0].C_nM[i,h+2]
			coef[i,j,3] = calc.wfs.kpt_u[0].C_nM[i,h+3]
			#if (at_num[j] == lower_at[0]):
			    #coef[i,j,:] *= lower_at[1]
			    #print "atomic hoppings from at:", lower_at[0], " l0wered by coef.:", lower_at[1]
			h += calc.wfs.setups[j].nao
	if (lower_atoms != []):
		print 'lowering atoms hoppings for atoms:', lower_atoms
		i_coef = 0;
		for j in lower_atoms:
			coef[:,j,:] *= lower_coefs[i_coef]
			i_coef +=1
	coeff = coef.flatten()
	coeffs = coeff.reshape((len(eig),num_at*4))
	if (cut_at != -1):
		coeffs=np.delete(coeffs,range(cut_at*4,num_at*4),1)
		num_at = cut_at
	#now removing non-effective orbitals:
	if (n_max != -1):
		print "cutting orbitals with too high eigen energies"
		coeffs = np.delete(coeffs,range(n_max+1,len(eig)),0)
		eig = np.delete(eig,range(n_max+1,len(eig)),0)
	if (n_min != -1):
		print "cutting orbitals with too low eigen energies"
		coeffs = np.delete(coeffs,range(n_min+1),0)
		eig = np.delete(eig,range(n_min+1),0)
	# applying PBC
	if ((pbc != (0,0))or(pbc != (0.0,0.0))) :
		print "applying pbc"
		coeff =np.repeat(coeffs,int(pbc[0]*2+1)*int(pbc[1]*2+1),0).flatten()
		num_at *=int(pbc[0]*2+1)*int(pbc[1]*2+1)
		coeffs = coeff.reshape((len(eig),num_at*4)) if (orbs == 'sp') else coeff.reshape((len(eig),num_at*9));
	print "coefficients read; now we are getting geometry, cut_at:", cut_at
	num_at = slab.get_number_of_atoms()
	if not ((cut_at == -1)or(cut_at == num_at)):
		print "cutting attoms"
		Ratin = slab.get_positions()[:cut_at,:]
	else:
		#print "NOT! cutting attoms"
		Ratin = slab.get_positions()
	print " Number of atoms: ", len(Ratin)
	if (pbc != ((0,0)or(0.,0.))):
		print "Applying PBC"
		atoms = Atoms('H%d' % len(Ratin), positions=Ratin)
		cell = slab.get_cell()
		atoms.set_cell(cell)
		if (pbc == (0.5,0.5)):
			atoms *= (2,2,1)
		else:
			atoms *= ( (int(2*pbc[0])+1),(int(2*pbc[1])+1),1 )
			atoms.translate( [( -int(pbc[0])*cell[0,0]-int(pbc[1])*cell[1,0] ) , (-int(pbc[0])*cell[0,1]-int(pbc[1])*cell[1,1]),0] )
		#from ase.visualize import view
		#view(atoms)
		Ratin = atoms.get_positions()
		print " Number of atoms after PBC: ", len(Ratin)
	return eig.copy(), coeffs.copy(), Ratin.copy();

#def	read_AIMS_all(name = 'KS_eigenvectors.band_1.kpt_1.out', geom='geometry.in', fermi=0.0, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[],header=False):

def	read_FIREBALL_all(name = 'phi_' , geom='answer.bas', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
	'''
	read_FIREBALL_all(name = 'phi_' , geom='answer.bas', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
	read coffecients and eigen numbers from Fireball made (iwrtcoefs = -2) files phik_0001_s.dat, phik_0001_py.dat ....
	fermi - If None the Fermi Level from the Fireball calculations (in case of molecule and visualising some molecular orbitals it can be move to their energy by putting there real value)
	orbs = 'sp' read only sp structure of valence orbitals (spd works, but calculator is done for sp only for now)
	pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
	imaginary = False (other options for future k-points dependency
	cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
	cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
	lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first for atoms by 0.5
	note: sometimes hydrogens have to have hoppings lowered by 0.5 this is under investigation
	'''
	assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 	
	assert (imaginary == False), "sorry imaginary version is under development" 	
	print " # ============ define atoms "
	atoms    = bU.loadAtoms(geom, elements.ELEMENT_DICT )
	assert (cut_at <= len(atoms[1])), "wrong cut for atoms"
	if not ((cut_at == -1)or(cut_at == len(atoms[1]))):
		atoms2 = [atoms[0][:cut_at],atoms[1][:cut_at],atoms[2][:cut_at],atoms[3][:cut_at]]
	else:
		atoms2 = atoms
	n_atoms= len(atoms2[1])
	print " Number of atoms: ", n_atoms
	if (pbc != ((0,0)or(0.,0.))):
		print "Applying PBC"
		if (pbc == (0.5,0.5)):
			atoms2 = bU.multCell( atoms2, [[lvs[0,0],lvs[0,1],0.],[lvs[1,0],lvs[1,1],0.],[0.,0.,100.]], m=(2,2,1) )
			Rs = np.array([atoms2[1],atoms2[2],atoms2[3]])
		else:
			atoms2 = bU.multCell( atoms2, ((lvs[0,0],lvs[0,1],0.),(lvs[1,0],lvs[1,1],0.),(0.,0.,100.)), m=( (int(2*pbc[0])+1),(int(2*pbc[1])+1),1 ) )
			Rs = np.array([atoms2[1],atoms2[2],atoms2[3]]); 
			Rs[0] -= int(pbc[0])*lvs[0,0]+int(pbc[1])*lvs[1,0]
			Rs[1] -= int(pbc[0])*lvs[0,1]+int(pbc[1])*lvs[1,1]
		print " Number of atoms after PBC: ", len(Rs[0])
	else:
		Rs = np.array([atoms2[1],atoms2[2],atoms2[3]]) 
	Ratin    = np.transpose(Rs).copy()
	del Rs;
	print "atomic geometry read"

	print "reading fireball LCAO coefficients for basis: ",orbs	
	filein = open(name+'s.dat' )
	pre_eig = filein.readline().split()
	filein.close()
	num_at=int(pre_eig[0]); n_bands= int(pre_eig[1]);
	if (fermi==None):
		fermi = float(pre_eig[2]);
	else:
		fermi += float(pre_eig[2]);
	del pre_eig;
	assert (n_atoms==num_at), "different number of atoms in geometry file"
	del n_atoms;
	eig = np.loadtxt(name+'s.dat',skiprows=1, usecols=(0,))
	assert (len(eig)==n_bands), "number of bands wrongly specified"
	eig -= fermi
	n_min = -1
	n_max = -1
	j = 0
	for i in eig:
		if (i < cut_min):
			n_min = j
		if (i < cut_max):
			n_max = j
		j += 1
	if (n_min and n_max != -1):
		assert (n_min < n_max), "no orbitals left for dI/dV"
	if (orbs == 'sp'):
		coef = np.zeros((n_bands,num_at,4))
		if (num_at > 1):
			coef[:,:,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
		else:
			coef[:,0,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
		if (lower_atoms != []):
			print 'lowering atoms hoppings for atoms:', lower_atoms
			i_coef = 0;
			for j in lower_atoms:
				coef[:,j,:] *= lower_coefs[i_coef]
				i_coef +=1
		coeff = coef.flatten()
		coeffs = coeff.reshape((n_bands,num_at*4))
		if (cut_at != -1):
			coeffs=np.delete(coeffs,range(cut_at*4,num_at*4),1)
			num_at = cut_at
	else:
		coef = np.zeros((n_bands,num_at,9))
		if (num_at > 1):
			coef[:,:,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,:,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
		else:
			coef[:,0,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
			coef[:,0,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(xrange(1, num_at*2+1, 2)) )
		coeff = coef.flatten() #xy yz z2 xz x2-y2
		coeffs = coeff.reshape((n_bands,num_at*9))
		if (cut_at != -1):
			coeffs=np.delete(coeffs,range(cut_at*9,num_at*9),1)
			num_at = cut_at
	#now removing non-effective orbitals:
	if (n_max != -1):
		print "cutting orbitals with too high eigen energies"
		coeffs = np.delete(coeffs,range(n_max+1,len(eig)),0)
		eig = np.delete(eig,range(n_max+1,len(eig)),0)
	if (n_min != -1):
		print "cutting orbitals with too low eigen energies"
		coeffs = np.delete(coeffs,range(n_min+1),0)
		eig = np.delete(eig,range(n_min+1),0)
	#now pbc applied
	if ((pbc != (0,0))or(pbc != (0.0,0.0))) :
		print "applying pbc"
		coeff =np.repeat(coeffs,int(pbc[0]*2+1)*int(pbc[1]*2+1),0).flatten()
		num_at *=int(pbc[0]*2+1)*int(pbc[1]*2+1)
		coeffs = coeff.reshape((len(eig),num_at*4)) if (orbs == 'sp') else coeff.reshape((len(eig),num_at*9));
	print "coefficients read"
	return eig.copy(), coeffs.copy(), Ratin.copy();

############## END OF LIBRARY ##################################
