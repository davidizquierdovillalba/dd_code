import os
import sys

from dictionaries import *
from Useful_func import *


#################################
#                               #
#          PARAMETERS           #
#                               #
#################################


GALTREE = False          # if 'True' reads galtree mode outputs and structures
MII = True              # if 'True' MII boxsize is used (see below)
wtp = [6, 1]             # properties to plot on x and y axis respectively. Available quatities are:  0:BH Mass,  1:StellarMass,  2:Mvir,  3:Sfr,  4:BulgeMass,  5:DiskMass, 6:seedTime
ztoplot = [0,1,2.5,6]     # list of 4 redshift to be plotted (ONLY 4!!!!!)

params_to_read = ['FileNameGalaxies', 'FirstFile', 'LastFile', 'BlackHoleSeedMass', 'AccretionModel', 'OutputDir']    # parameters/names to be read from LGalaxies input.par

plot_last_run = False     # if 'True' output_dir name gets read from param_file variable (see below, section "RUN LGalaxies")
fastcheck = False        # if 'True' plots only up to a fixed LastFile (see LGparams definition below)
Plotlog = False           # if 'True' log-scale plots are produced
Factor10 = [True,True]   # if 'True' uses 1.**10 units in plots (x and y axis, respectively)
removeh =  True         # if 'True' plotted quantities are in "h-free" units (for data comparison)
DMap_RES = 200           # Resolution (number of bins in each axis) for density maps (optimal: 100)

sizebin = 0.1		 # Size bin in MoL_func when Plotlog = True; otherwise other sizebin is automatically selected





###################################
#                                 #
#             FOLDERS             #
#                                 #
###################################
home = os.path.expanduser('~')+'/'
user = home[12:-1]

if(user == 'dspinoso'):
	here = home + 'works/phd/dd_code/'              # folder where the code is
	LG_dir = home + 'LGalaxies/'                    # LGalaxies code root folder
        if MII == True:
                LG_inParFile = LG_dir + 'input/input_dani_MRII_W1_PLANCK.par'
                final = 'MRII/'
        else:
                LG_inParFile = LG_dir+'/input/input_dani_MR_W1_PLANCK.par'
                final = 'MR/'
        LG_output_z = LG_dir + 'input/desired_output_redshifts_dani.txt'

	#################################
	#                               #
	#        RUN LGalaxies          #
	#                               #
	#################################
	
	ptc = np.array(['BlackHoleSeedMass', 'AgnEfficiency','OutputDir'])  # run-dependent parameters to be written in the output folder name
	ptc_alis = np.array(['seedM','eff'])              #parameter-alias to be written in the output folder name
	par_changed = [1,0,1]   #The last one MUST be ALWAYS = 1


elif(user == 'dizquierdo'):
	here = home + '/Documentos/thesis/dd_code/'     # folder where the code is
	LG_dir = home + 'Documentos/thesis/LGalaxies/'  # LGalaxies code root folder
	if MII == True:
                LG_inParFile = LG_dir + 'input/input_david.par'
                final = 'MRII/'
        else:
                LG_inParFile = LG_dir + 'input/input_david.par'
                final = 'MR/'
        LG_output_z = LG_dir + 'input/david_redshift.txt'

	#################################
	#                               #
	#        RUN LGalaxies          #
	#                               #
	#################################
	
	ptc = np.array(['BlackHoleSeedMass', 'AgnEfficiency', 'AccretionModel', 'OutputDir'])  # run-dependent parameters to be written in the output folder name
	ptc_alis = np.array(['seedM','eff', 'accmod'])              #parameter-alias to be written in the output folder name
	par_changed = [1,0,1,1]   #The last one MUST be ALWAYS = 1
        
else:
	print 'Error: User not found in the Config.py'

LGout_dir = LG_dir + 'output/' + final                  # LGalaxies outputs folder
AuxCode_dir = LG_dir + 'AuxCode/Python/'        # folder in which Bruno's useful python code is

plots_dir = here + 'plots/'                     # directory to store plots made by this code
if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        
sys.path.append(AuxCode_dir + 'Misc/')          # import read_snap and read_tree functions
from read_lgal import *


LGparams = read_LG_inParamFile(LG_inParFile, LG_output_z, params_to_read)
if fastcheck == True:
        LGparams['LastFile'] = 1 # Comment this



###################################
#                                 #
#      FILES and PLOTS NAMES      #
#                                 #
###################################
seedMass = '1.e3'  # should be = BlackHoleSeedMass in LGal input files. This goes in the file prefix
fileprefx = get_prefix(GALTREE, MII, LGparams)  # automatic file prefix. Always call with all keywords

png = True              # How to save the plots: In png (png = True) or in pdf (png = False)
if(png == True):
	dens_ns = plots_dir + Prop(wtp[0])+'_vs_'+Prop(wtp[1])+'density_plot.png'  # density plot name
	func_ns = plots_dir + Prop(wtp[0])+'_function.png'        # luminosity/mass function plot name
        if MII == True:
                dens_ns = plots_dir + Prop(wtp[0])+'_vs_'+Prop(wtp[1])+'density_plot_MII_seed'+fileprefx[-7:-2]+'.png'  # density plot name
                func_ns = plots_dir + Prop(wtp[0])+'_function_MII_seed'+fileprefx[-7:-2]+'.png'        # luminosity/mass function plot name
                
else:
	dens_ns = plots_dir + Prop(wtp[0])+'_vs_'+Prop(wtp[1])+'density_plot.pdf'  # density plot name
	func_ns = plots_dir + Prop(wtp[0])+'_function.pdf'        # luminosity/mass function plot name
        if MII == True:
                dens_ns = plots_dir + Prop(wtp[0])+'_vs_'+Prop(wtp[1])+'density_plot_MII_seed'+fileprefx[-7:-2]+'.pdf'  # density plot name
                func_ns = plots_dir + Prop(wtp[0])+'_function_MII_seed'+fileprefx[-7:-2]+'.pdf'        # luminosity/mass function plot name


                


########################################
#                                      #
#    COSMOLOGIES & DARK MATTER SIMS    #
#                                      #
########################################

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





        
