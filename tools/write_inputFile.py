import sys
import numpy as np
from Config import *
from dictionaries import *
from Useful_func import *


#THIS CODE CREATE A TAILORED INPUT FILE
#BASED ON CONFIG.PY OPTIONS

def do_and_move_parFile():
	if MII == True:
	    temp_fold = ['./output/MRII/']
	else:
	    temp_fold = ['./output/MR/']

	#extract info from input parameter file
	f = open(LG_inParFile, 'r')
	lines = f.readlines()
	f.close()
	g = open(LG_inParFile, 'w')

	par_changed2 = np.array(par_changed, dtype=bool)
	ptu = ptc[par_changed2]
	ptu_alis = ptc_alis[par_changed2[:-1]]
	for s,kk in zip(ptu, range(len(ptu))): # The last one will always be OutputDir
	    for i in range(len(lines)):
		if s in lines[i] and s == ptu[kk] and s != 'OutputDir':
		    aa = str(lines[i].split()[1])
                    if(s == 'BHGrowthInDiskInstabilityModel'): aa_val = ptu_alis[kk]+str(int(lines[i].split()[1])) 
                    if(s == 'BlackHoleCutoffVelocity'): aa_val = ptu_alis[kk]+str(float(lines[i].split()[1]))
                    if(s == 'BlackHoleGrowthRate'): aa_val = ptu_alis[kk]+str(float(lines[i].split()[1]))
                    if(s == 'BlackHoleSeedMass'): aa_val = ptu_alis[kk]+str( np.round(np.log10( float(lines[i].split()[1])),2) )
                    if(s == 'AccretionModel'): aa_val = ptu_alis[kk]+str(int(lines[i].split()[1]))
                    if(s == 'DIduringMerger'): aa_val = ptu_alis[kk]+str(int(lines[i].split()[1])) 
                    if(s == 'BHGrRaDI'): aa_val = ptu_alis[kk]+str(float(lines[i].split()[1])) 
		    else: aa_val = ptu_alis[kk]+str(float(lines[i].split()[1]))
		    #else: aa_val = ptu_alis[kk]+str( np.round(np.log10( float(lines[i].split()[1])),2) )
		    temp_fold.append(aa_val)
		if s in lines[i] and s == 'OutputDir':
		    dummy1 = str(lines[i].split()[0])
		    dummy2 = temp_fold[0]
		    for j in temp_fold[1:]:
			dummy2 += j+'_'
		    dummy2 = dummy2[:-1]+'\n'
		    lines[i] = dummy1+'               '+dummy2
	g.writelines(lines)
	g.close()
	if MII == True:
	    jjjj = dummy2[14:-1]
	else:
	    jjjj = dummy2[12:-1]

	if not os.path.exists(LGout_dir+jjjj):
		os.makedirs(LGout_dir+jjjj)

	os.chdir(LGout_dir+jjjj)
	copy_to_outFolder = 'cp '+LG_inParFile+' .'
	os.system(copy_to_outFolder)
