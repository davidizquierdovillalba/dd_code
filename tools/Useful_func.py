

def update_print(phrase, appendix=' ', percent=False):
    import sys
    import time
    from time import sleep
    sys.stdout.write('\r')
    if percent==True:
        sys.stdout.write(phrase + '%d%%' % appendix )
        if appendix == 100:
            print '\n'
    if percent==False:
        sys.stdout.write(phrase + appendix )
        if appendix == 'done':
            print '\n'
    sys.stdout.flush()
    sleep(1.e-6)



    
def hist(values,sizebins):
	import numpy as np
        mm = min(values)
        mx = max(values)
        bins=np.arange(mm,mx,sizebins)
        hist, b_edges = np.histogram(values,bins)
        centers = (b_edges[:-1] + b_edges[1:]) / 2.0
        hist = np.array(hist, dtype=np.float64)
        return centers, hist




def get_prefix(gltree, mill2, seed):
    if gltree == True:
	pfx = 'SA_galtree_'
    elif mill2 == True:
	pfx = 'SA_MII_BHseed1e3_z'#'SA_MII__z'
    else:
        pfx = 'SA_z'
    return pfx



def get_outFolder_names(plot_last_run, LGparams, LGout_dir, here):
    import os
    import numpy as np
    if plot_last_run == True:
        out_LGfiles_dir = LGparams['OutputDir']
    else:
        os.chdir(LGout_dir)
        os.system('ls -1d */ > out_folders_list.txt')
        os.chdir(here)
        fold_names = np.loadtxt(LGout_dir+'out_folders_list.txt', dtype=np.str)
        print 'Output folders available for plots:\n', fold_names
        fold_mask = []
        for i in fold_names:
            a = raw_input('Plot files from  '+i+'  folder? (1=selection, 0=no):\n')
            fold_mask.append(int(a))
        fold_mask = np.array(fold_mask, dtype=bool)
        out_LGfiles_dir = []
        for i in fold_names[fold_mask]:
            out_LGfiles_dir.append(LGout_dir + i)
        out_LGfiles_dir = np.array(out_LGfiles_dir, dtype=np.str)

    return out_LGfiles_dir




def read_LG_inParamFile(inputs, z_list, ptr):
    params = {}
    
    #extract info from input parameter file
    f = open(inputs, 'r')
    lines = f.readlines()
    f.close()
    for s in ptr:
        for i in range(len(lines)):
            if s in lines[i] and s != 'FileNameGalaxies' and s != 'OutputDir':
                if s != 'BlackHoleSeedMass':
                    params[s] = int(lines[i].split()[1])
                else:
                    params[s] = float(lines[i].split()[1])
            if s in lines[i] and s == 'FileNameGalaxies':
                params[s] = lines[i].split()[1]+'_z'
            if s in lines[i] and s == 'OutputDir':
                params[s] = lines[i].split()[1]
    #extract info from redshift list (without importing numpy)
    f = open(z_list, 'r')
    params['zlist'] = [float(i) for i in f.read().split()]
    f.close()

    return params



    
def get_sb(macs, minim):
    szbin = 1./20
    szbin = (macs - minim)*szbin
        
    return szbin

