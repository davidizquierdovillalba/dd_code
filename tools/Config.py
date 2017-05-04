import os
import sys

##############################
#                            #
#           FOLDERS          #
#                            #
##############################

home = os.path.expanduser('~')+'/'              # current home folder (automatic)
here = home + 'Documentos/thesis/dd_code/'      # folder where the code is
LG_dir = home + 'Documentos/thesis/LGalaxies/'  # LGalaxies code root folder
LGout_dir = LG_dir + 'output/'                  # LGalaxies outputs folder
AuxCode_dir = LG_dir + 'AuxCode/Python/'        # folder in which useful python code is
plots_dir = here + 'plots/'                     # (this code) output plots folder
if not os.path.exists(plots_dir):
	os.makedirs(plots_dir)


#################################
#                               #
#          PARAMETERS           #
#                               #
#################################

zlist = [0, 1, 2, 3]
GALTREE = False
MII = False
Plotlog = True
Factor10 = [True,True]
removeh =  False
if Plotlog:
	sizebin = 0.1
	if Factor10[0]:
		MoL_func_xlim = [4.0,9.9]
	else:
		MoL_func_xlim = [-6.0,-1.9]
else:
	sizebin = 1/20
        if Factor10[0]:
                MoL_func_xlim = [10**4.0,10**9.9]
        else:
                MoL_func_xlim = [10**-6.0,10**-1.9]
MoL_func_ylim  = [10**-5, 10**-1.0]
firstfile = 0
lastfile = 11
DMap_RES = 500 # Resolution (number of bins in each axis) for density maps
ns = plots_dir + 'check.png'

wtp = [0, 3]  # properties to plot  0:BH Mass, 1:StellarMass, 2:Mvir, 3:Sfr, 4:BulgeMass, 5:DiskMass


sys.path.append(AuxCode_dir + 'Misc/')
from read_lgal import * 
if GALTREE:
	fileprefx = 'SA_galtree_'
else:
	fileprefx = 'SA_z'
from dictionaries import *
P1 = Prop(wtp[0]) 
P2 = Prop(wtp[1])

# COSMOLOGIES & DARK MATTER SIMS #

WMAP1=0
PLANCK=1

if WMAP1: 
	FullRedshiftList=[0.00,0.41,0.99,2.07,3.06,3.87] 
	FullSnapshotList=[63,50,41,32,27,24]  
	if MII:
		BoxSize  = 100. #full MRII      
	else:
		BoxSize  = 500. #full MR 
	Hubble_h      = 0.73
	Omega_M       = 0.25 
	Omega_Lambda  = 0.75
	MaxTreeFiles  = 512

if PLANCK: 
	FullRedshiftList=[0.00,0.11,0.40,1.04,2.07,3.11,3.95] 
	FullSnapshotList=[58,53, 47,38,30,25,22]  
	if MII:
                BoxSize  = 100. * 0.960558 #full MRII      
        else:
                BoxSize  = 500. * 0.960558 #full MR
	Hubble_h      = 0.673
	Omega_M       = 0.315 
	Omega_Lambda  = 0.683
	MaxTreeFiles  = 512

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = Hubble_h*100, Om0 = Omega_M)


#################################
#                               #
#     LG OUTPUT STRUCTURES      #
# (decomment only what is used) #
#                               #
#################################

sys.path.append(AuxCode_dir)
if(GALTREE == False):
	from LGalaxies_Henriques2015a_struct import LGalaxiesStruct 
	from LGalaxies_Henriques2015a_struct import PropertiesToRead
else:
	from LGalaxies_tree_Henriques2015a_struct import LGalaxiesStruct 
	from LGalaxies_tree_Henriques2015a_struct import PropertiesToRead_tree

