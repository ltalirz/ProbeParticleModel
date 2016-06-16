#!/usr/bin/python

import os
import numpy as np
from   ctypes import c_int, c_double, c_char_p
import ctypes
import basUtils as bU
import elements

import cpp_utils


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

def	read_GPAW_all(name = 'OUTPUT.gpw', fermi = None,orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_at=[0,1.] ):
	'''
	read_GPAW_coeffs(name = 'OUTPUT.gpw' ,orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_at=[0,1.] ):
	read eigen energies, coffecients, Fermi Level and geometry  from the GPAW  *.gpw file.
	If fermi = None then Fermi comes from the GPAW calculation
	orbs - only 'sp' works 	can read only sp structure of valence orbitals (hydrogens_has to be at the end !!!!)
	pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
	imaginary = False (other options for future k-points dependency
	cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
	cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
	lower_at=[0,1.] ... do nothing; lower_at=[8,0.5] lower coefficients (=hoppings) for all oxygens (=8) by factor 0.5 -- necessarry for some systems
	'''
	assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
	assert (imaginary == False), "sorry imaginary version is under development" 	
	print "reading GPAW LCAO coefficients for basis: ",orbs	
	from ase import *
	from gpaw import GPAW
	calc = GPAW(name)
	slab = calc.get_atoms()
	num_at = slab.get_number_of_atoms()
	n_bands = calc.get_number_of_bands()
	eig = calc.get_eigenvalues(kpt=0, spin=0, broadcast=True)
	at_num = slab.get_atomic_numbers()
	if (fermi == None):
		fermi = calc.get_fermi_level()
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
			if (at_num[j] == lower_at[0]):
			    coef[i,j,:] *= lower_at[1]
			    #print "atomic hoppings from at:", lower_at[0], " l0wered by coef.:", lower_at[1]
			h += calc.wfs.setups[j].nao
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



def	read_fire_coeffs(fermi, name = 'phik_' ,orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
	'''
	read_fire_coeffs(fermi, name = 'phik_' ,orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[] ):
	read coffecients and eigen numbers from Fireball made (iwrtcoefs = -2) files phik_0001_s.dat, phik_0001_py.dat ....
	fermi - the Fermi Level from the Fireball calculations (in case of molecule and visualising some molecular orbitals it can be move to their energy
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
	print "reading fireball LCAO coefficients for basis: ",orbs	
	eig = np.loadtxt(name+'s.dat',skiprows=0, usecols=(0,))
	num_at = int(eig[0])
	eig = np.delete(eig,0)
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
		coef = np.zeros((len(eig),num_at,4))
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
		coeffs = coeff.reshape((len(eig),num_at*4))
		if (cut_at != -1):
			coeffs=np.delete(coeffs,range(cut_at*4,num_at*4),1)
			num_at = cut_at
	else:
		coef = np.zeros((len(eig),num_at,9))
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
		coeffs = coeff.reshape((len(eig),num_at*9))
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
	return eig.copy(), coeffs.copy();

def read_fire_atoms(name, cut_at=-1, pbc=(1,1), lvs=np.array([[0.0,0.0],[0.0,0.0]]) ):
	'''
	read_fire_atoms(name, cut=-1, pbc=(1,1), lvs=np.array([[0.0,0.0],[0.0,0.0]]):
	Give you the only array with atomic coordinates needed for the STM procedures
	reads .xyz or .bas files.(name)
	pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
	lvs = 2x2 array of x and y lattice vectors. !!! BE AWARE no checking for non-zero lattice !!!!
	cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
	'''
	#if(pbc != (0,0)or(0.,0.)):
	#	assert (lvs != np.array([[0.0,0.0],[0.0,0.0]])), "wrongly specified PBC and lattice vectors"
	print " # ============ define atoms "
	atoms    = bU.loadAtoms(name, elements.ELEMENT_DICT )
	assert (cut_at <= len(atoms[1])), "wrong cut for atoms"
	if not ((cut_at == -1)or(cut_at == len(atoms[1]))):
		atoms2 = [atoms[0][:cut_at],atoms[1][:cut_at],atoms[2][:cut_at],atoms[3][:cut_at]]
	else:
		atoms2 = atoms
	print " Number of atoms: ", len(atoms2[1])
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
	return Ratin ;#, Natin;

def dIdV( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0):
	'''
	dIdV( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0):
	V - voltage = (energy vs. the Fermi Level in eV);
	WF - workfunction (normally ~5 eV gives reasonable results),
	eta - energy smearing (0.5-0.30 eV) deppending on system (for single orbital very low number
	eig - eigenenergies of sample states (=molecular orbitals)
	R input of points in whish you calculate dI/dV (relaxed via PP afm, or nonrelaxed via mkSpaceGrid)
	coes -- LCAO coefficients from read_fire_coes (Fireball, maybe FHI-AIMS & mathematica) or read_GPAW_all
	orbs = 'sp' orbitals of the sample (spd don't work at the moment
	s and/or px and/or py and/or pz orbitals at the PP
	unification of all the predefined dI/dV procedures from C++, you can choose, whatever PP orbital you want
	'''
	#assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 
	assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
	assert ((s > 0.0)or(px > 0.0)or(py > 0.0)or(pz > 0.0)), "all tip orbitals are zero"
	assert ((s >= 0.0)and(px >= 0.0)and(py >= 0.0)and(pz >= 0.0)), "you cannot have negative current"
	bo=False
	if (orbs == 'sp'):
		if ((px>0.0)and(px==py)):
			cur = px*dIdV_pxy_sp( V, WF, eta ,eig, R, Rat, coes)
			bo=True
			if (s>0.0):
				cur += s*dIdV_s_sp( V, WF, eta ,eig, R, Rat, coes)
			if (pz>0.0):
				cur += pz*dIdV_pz_sp( V, WF, eta ,eig, R, Rat, coes)
		else:
			if (s>0.0):
				cur = s*dIdV_s_sp( V, WF, eta ,eig, R, Rat, coes)
				bo=True
			if (px>0.0):
				curpx = px*dIdV_px_sp( V, WF, eta ,eig, R, Rat, coes)
				if bo:
					cur += curpx
				else:
					cur = curpx
					bo=True
			if (py>0.0):
				curpy = py*dIdV_py_sp( V, WF, eta ,eig, R, Rat, coes)
				if bo:
					cur += curpy
				else:
					cur = curpy
					bo=True
			if (pz>0.0):
				curpz = pz*dIdV_pz_sp( V, WF, eta ,eig, R, Rat, coes)
				if bo:
					cur += curpz
				else:
					cur = curpz
					bo=True
	'''
	else:
		if (s>0.0):
			cur = s*dIdV_s_spd( V, WF, eta ,eig, R, Rat, coes)
			bo=True
		if (px>0.0):
			curpx = px*dIdV_px_spd( V, WF, eta ,eig, R, Rat, coes)
			if bo:
				cur += curpx
			else:
				cur = curpx
				bo=True
		if (py>0.0):
			curpy = py*dIdV_py_spd( V, WF, eta ,eig, R, Rat, coes)
			if bo:
				cur += curpy
			else:
				cur = curpy
				bo=True
		if (pz>0.0):
			curpz = pz*dIdV_pz_spd( V, WF, eta ,eig, R, Rat, coes)
			if bo:
				cur += curpz
			else:
				cur = curpz
				bo=True
	'''
	return cur;


def dIdV_para( V, WF, eta ,eig, R, Rat, coes, orbs, s, px, py, pz, name):
	'''
	dIdV( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, cur_in):
	V - voltage = (energy vs. the Fermi Level in eV);
	WF - workfunction (normally 5eV),
	eta - energy smearing (0.5-0.15), eig - eigenenergies of sample states (=molecular orbitals)
	R input of points in whish you calculate dI/dV,...
	unification of all the predefined dI/dV procedures from C++, you can choose, sheather you want
	sp or spd orbitals of the sample, and s and/or px and/or py and/or pz orbitals of the tip apex
	'''
	#assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 
	assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
	assert ((s > 0.0)or(px > 0.0)or(py > 0.0)or(pz > 0.0)), "all tip orbitals are zero"
	assert ((s >= 0.0)and(px >= 0.0)and(py >= 0.0)and(pz >= 0.0)), "you cannot have negative current"
	bo=False
	if (orbs == 'sp'):
		if ((px>0.0)and(px==py)):
			cur = px*dIdV_pxy_sp( V, WF, eta ,eig, R, Rat, coes)
			bo=True
			if (s>0.0):
				cur += s*dIdV_s_sp( V, WF, eta ,eig, R, Rat, coes)
			if (pz>0.0):
				cur += pz*dIdV_pz_sp( V, WF, eta ,eig, R, Rat, coes)
		else:
			if (s>0.0):
				cur = s*dIdV_s_sp( V, WF, eta ,eig, R, Rat, coes)
				bo=True
			if (px>0.0):
				curpx = px*dIdV_px_sp( V, WF, eta ,eig, R, Rat, coes)
				if bo:
					cur += curpx
				else:
					cur = curpx
					bo=True
			if (py>0.0):
				curpy = py*dIdV_py_sp( V, WF, eta ,eig, R, Rat, coes)
				if bo:
					cur += curpy
				else:
					cur = curpy
					bo=True
			if (pz>0.0):
				curpz = pz*dIdV_pz_sp( V, WF, eta ,eig, R, Rat, coes)
				if bo:
					cur += curpz
				else:
					cur = curpz
					bo=True
	#print "cur:", cur
	np.save(name,cur)
	# The end of parallel procedure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def STM( V, nV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, WF_decay=1.0):
	'''
	STM( V, nV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, WF_decay=1.0):
	summing more dI/dV via rectangle integration, be aware Work Function is changing with Voltage!
	'''
	assert (float(V) != 0.0),"you cannot have zero Voltage"
	print "STM simulation via more dI/dV calculations"
	print "Be aware Work Function is changing with Voltage by a factor: ",  WF_decay
	ii=1;
	for v_ in np.linspace(0.,V,nV):
		print "Start to calculate voltage step %d of %d in total." %(ii, nV)
		ii +=1
		assert (WF-v_*WF_decay> 0.1), "Non-physical Work Function or Voltage, together WF <= 0.1 eV	"
		print "WF for this step is: " , WF-v_*WF_decay, " eV"
		i_ = dIdV( v_, WF-v_*WF_decay, eta ,eig, R, Rat, coes, orbs=orbs, s=s, px =px, py=py, pz=pz)
		#print "maximal dI/dV: " , max(i_)
		if (v_ == 0):
			cur = i_
		else:
			cur += i_
	cur *= abs(V)*7.7480917346E-05 # rescaling into Amper
	print "All dI/dV steps done, current rescalled into Ampers"
	return cur;

# ==============================
# ============================== interface to C++ core 
# ==============================

# ============================== interface to C++ core 

cpp_name='ProbeSTM_0.0.10sp'
#cpp_utils.compile_lib( cpp_name  )
cpp_utils.make("STM")
lib    = ctypes.CDLL(  cpp_utils.CPP_PATH + "/" + cpp_name + cpp_utils.lib_ext )     # load dynamic librady object using ctypes 

# define used numpy array types for interfacing with C++

array1i = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array1d = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
array4d = np.ctypeslib.ndpointer(dtype=np.double, ndim=4, flags='CONTIGUOUS')

# ========
# ======== Python warper function for C++ functions
# ========

#************* s *************
# void proc_dIdVssp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur)
lib.proc_dIdVssp.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d ]
lib.proc_dIdVssp.restype  = None
def dIdV_s_sp( V, WF, eta ,eig, R, Rat, coes):
	print "Entering the dI/dV (s-sp) procedure"
	NoAt = len(Rat)
	NoOrb = len(eig)
	sh = R.shape
	cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
	Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
	assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
	if (len(coes) != 0):
		assert (NoOrb == len(coes)*len(coes[0])/(4*NoAt)), "Different eigennumbers, than basis"	
	print "We're going to C++"
	lib.proc_dIdVssp( NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

#************* px *************
# void proc_dIdVpxsp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur)
lib.proc_dIdVpxsp.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d ]
lib.proc_dIdVpxsp.restype  = None
def dIdV_px_sp( V, WF, eta ,eig, R, Rat, coes):
	print "Entering the dI/dV (px-sp) procedure"
	NoAt = len(Rat)
	NoOrb = len(eig)
	sh = R.shape
	cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
	Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
	assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
	if (len(coes) != 0):
		assert (NoOrb == len(coes)*len(coes[0])/(4*NoAt)), "Different eigennumbers, than basis"	
	print "We're going to C++"
	lib.proc_dIdVpxsp( NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

#************* py *************
# void proc_dIdVpysp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur)
lib.proc_dIdVpysp.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d ]
lib.proc_dIdVpysp.restype  = None
def dIdV_py_sp( V, WF, eta ,eig, R, Rat, coes):
	print "Entering the dI/dV (py-sp) procedure"
	NoAt = len(Rat)
	NoOrb = len(eig)
	sh = R.shape
	cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
	Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
	assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
	if (len(coes) != 0):
		assert (NoOrb == len(coes)*len(coes[0])/(4*NoAt)), "Different eigennumbers, than basis"	
	print "We're going to C++"
	lib.proc_dIdVpysp( NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

#************* pxy - fast pxy procedure *************
# void proc_dIdVpxysp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur)
lib.proc_dIdVpxysp.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d ]
lib.proc_dIdVpxysp.restype  = None
def dIdV_pxy_sp( V, WF, eta ,eig, R, Rat, coes):
	print "Entering the dI/dV (pxy-sp) procedure"
	NoAt = len(Rat)
	NoOrb = len(eig)
	sh = R.shape
	cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
	Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
	assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
	if (len(coes) != 0):
		assert (NoOrb == len(coes)*len(coes[0])/(4*NoAt)), "Different eigennumbers, than basis"	
	print "We're going to C++"
	lib.proc_dIdVpxysp( NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

#************* pz *************
# void proc_dIdVpzsp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur)
lib.proc_dIdVpzsp.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d ]
lib.proc_dIdVpzsp.restype  = None
def dIdV_pz_sp( V, WF, eta ,eig, R, Rat, coes):
	print "Entering the dI/dV (pz-sp) procedure"
	NoAt = len(Rat)
	NoOrb = len(eig)
	sh = R.shape
	cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
	Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
	assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
	if (len(coes) != 0):
		assert (NoOrb == len(coes)*len(coes[0])/(4*NoAt)), "Different eigennumbers, than basis"	
	print "We're going to C++"
	lib.proc_dIdVpzsp( NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();


############## END OF LIBRARY ##################################
