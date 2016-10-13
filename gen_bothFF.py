#!/usr/bin/python
import sys
import numpy as np
import os
import __main__ as main


import pyProbeParticle                as PPU     
from   pyProbeParticle            import basUtils
from   pyProbeParticle            import elements   
import pyProbeParticle.GridUtils      as GU
#import pyProbeParticle.core          as PPC
import pyProbeParticle.HighLevel      as PPH
import pyProbeParticle.fieldFFT       as fFFT

HELP_MSG="""Use this program in the following way:
"""+os.path.basename(main.__file__) +""" -i <filename> [ --sigma <value> ]

Supported file fromats are:
   * cube
   * xsf """

from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-i", "--input", action="store", type="string", help="format of input file")
parser.add_option( "-g", "--geometry", action="store", type="string", help="format of input file")
parser.add_option( "-t", "--tip", action="store", type="string", help="tip model (multipole)", default='s')
parser.add_option( "-w", "--sigma", action="store", type="float",  help="gaussian width for convolution in Electrostatics [Angstroem]", default=0.7)
(options, args) = parser.parse_args()

print options


if options.input==None:
    sys.exit("ERROR!!! Please, specify the input file with the '-i' option \n\n"+HELP_MSG)

is_xyz  = options.input.lower().endswith(".xyz")
is_cube = options.input.lower().endswith(".cube")
is_xsf  = options.input.lower().endswith(".xsf")

gs_xyz  = options.input.lower().endswith(".xyz")
gs_cube = options.input.lower().endswith(".cube")
gs_xsf  = options.input.lower().endswith(".xsf")


print " >> OVEWRITING SETTINGS by params.ini  "
PPU.loadParams( 'params.ini' )

print " ========= get electrostatic forcefiled from hartree "

# TODO with time implement reading a hartree potential generated by different software
print " loading Hartree potential from disk "

if( is_xsf ):
    print "Use loadXSF"
    V, lvec, nDim, head = GU.loadXSF(options.input)
elif( is_cube ):
    print "Use loadCUBE"
    V, lvec, nDim, head = GU.loadCUBE(options.input)
    V*=27.211396132
else:
    sys.exit("ERROR!!! Unknown format of the input file\n\n"+HELP_MSG)
rho = None
multipole = None
if options.tip in {'s','px','py','pz','dx2','dy2','dz2','dxy','dxz','dyz'}:
    rho = None
    multipole={options.tip:1.0}
elif options.tip.endswith(".xsf"):
    rho, lvec_tip, nDim_tip, tiphead = GU.loadXSF(options.tip)
    if any(nDim_tip != nDim):
        sys.exit("Error: Input file for tip charge density has been specified, but the dimensions are incompatible with the Hartree potential file!")    

print " computing convolution with tip by FFT "
Fel_x,Fel_y,Fel_z = fFFT.potential2forces(V, lvec, nDim, rho=rho, sigma = options.sigma, multipole = multipole)
Fel = GU.packVecGrid(Fel_x,Fel_y,Fel_z)

print " saving electrostatic forcefiled "

GU.saveVecFieldNpy( 'FFel', Fel, lvec)

del Fel_x,Fel_y,Fel_z,V, Fel

if (options.noLJ):
	print "Computing LJ potential with parameters from electrostatic grid"
	#lvec from Electrostatic !!!!
	PPU.params['gridA'] =    lvec[ 1,:  ].copy() 
	PPU.params['gridB'] =    lvec[ 2,:  ].copy()
	PPU.params['gridC'] =    lvec[ 3,:  ].copy()
	PPU.params['gridN'] = np.array([nDim[2],nDim[1],nDim[0]])

	if options.geom==None:
		if(is_cube):
			atoms = basUtils.loadAtomsCUBE(options.input,elements.ELEMENT_DICT)
		elif(is_xsf):
			atoms, nDim, lvec = basUtils.loadXSFGeom( options.input )
		else:
			sys.exit("ERROR!!! Unknown format of geometry system. Supported formats: .xyz, .cube \n\n")
	else:
		if not (gs_xyz or gs_cube or gs_xsf):
			sys.exit("ERROR!!! Unknown format of the geometry input file\n\n"+HELP_MSG)

		print "--- Compute Lennard-Jones Force-field ---"
		if(gs_xyz):
				atoms = basUtils.loadAtoms(options.geom, elements.ELEMENT_DICT )
		elif(gs_cube):
				atoms = basUtils.loadAtomsCUBE(options.geom,elements.ELEMENT_DICT)
		elif(gs_xsf):
			atoms, nDim, lvec = basUtils.loadXSFGeom( options.input )
		else:
			sys.exit("ERROR!!! Unknown format of geometry system. Supported formats: .xyz, .cube \n\n")

		FFparams=None
		if os.path.isfile( 'atomtypes.ini' ):
				print ">> LOADING LOCAL atomtypes.ini"  
				FFparams=PPU.loadSpecies( 'atomtypes.ini' ) 
		iZs,Rs,Qs = PPH.parseAtoms( atoms, autogeom = False, PBC = options.noPBC )
		FFLJ      = PPH.computeLJ( Rs, iZs, FFLJ=None, FFparams=FFparams )

		GU.limit_vec_field( FFLJ, Fmax=10.0 ) # remove too large valuesl; keeps the same direction; good for visualization 

	print "--- Save  ---"
	GU.saveVecFieldNpy( 'FFLJ', FFLJ, lvec)

print "--- Force-field(s) saved "
