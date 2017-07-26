import numpy as np

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




def get_prefix(gltree,mill2, seed):
    if gltree == True:
	pfx = 'SA_galtree_'
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
                if s != 'BlackHoleSeedMass' and s!= 'BlackHoleGrowthRate' and s!='BlackHoleCutoffVelocity' and s!='AgnEfficiency' and s!='BHGrRaDI':
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



def Hopkins_fit(phi_star, L_star, gamma1, gamma2, array_Lum_erg):
	#phi = phi_star * ((array_Lum_erg/L_star)**(-gamma1)) * np.exp(-((array_Lum_erg/L_star)**(gamma2)))
	#phi = phi_star/((array_Lum_erg/L_star)**gamma1 + (array_Lum_erg/L_star)**gamma2)
	phi = phi_star/((array_Lum_erg/L_star)**gamma1 + (array_Lum_erg/L_star)**gamma2)
	
	return phi
	
def Hopkins_fit_err(phi_star, errphi_star, L_star, err_L_star, gamma1, err_gamma1,gamma2, err_gamma2, array_Lum_erg):


	der_phi_star = errphi_star/((array_Lum_erg/L_star)**gamma1 + (array_Lum_erg/L_star)**gamma2)
	der_phi_star = (der_phi_star)**2

	der_L_star = (phi_star * err_L_star * L_star**-2) * ((array_Lum_erg/L_star)**(gamma1-1) + (array_Lum_erg/L_star)**(gamma2-1)) / ((array_Lum_erg/L_star)**gamma1 + (array_Lum_erg/L_star)**gamma2)**2
	der_L_star = (der_L_star)**2

	der_gamma1 = -phi_star* err_gamma1 * (((array_Lum_erg/L_star)**gamma1 + (array_Lum_erg/L_star)**gamma2)**(-2)) * (gamma1*(array_Lum_erg/L_star)**(gamma1)) * np.log((array_Lum_erg/L_star))
	der_gamma1 = (der_gamma1)**2

	der_gamma2 = -phi_star* err_gamma2 * (((array_Lum_erg/L_star)**gamma1 + (array_Lum_erg/L_star)**gamma2)**(-2)) * (gamma2*(array_Lum_erg/L_star)**(gamma2)) * np.log((array_Lum_erg/L_star))
	der_gamma2 = (der_gamma2)**2

	err_tot = (der_phi_star+der_L_star+der_gamma1+der_gamma2)**0.5
	phi = Hopkins_fit(phi_star, L_star, gamma1, gamma2, array_Lum_erg)

	phi_max= phi + err_tot
	phi_min = phi - err_tot
	return phi_max,phi_min



def Add_Quasar_LF_data(redshift, data=False, fit=False):
        from Config import *
        if(data == True):
                zz, Lum, phi_Lum, phi_Lum_err, band, author = np.loadtxt(data_Dir + "/Quasar_LF", unpack = True)
                Lum = Lum + 33.6 # erg/s
                _sel = np.where(zz == redshift)
                if(len(_sel[0])==0):
                        return 0, 0 , 0
                else:
                        return Lum[_sel], phi_Lum[_sel], phi_Lum_err[_sel]
        if(fit == True):
                zz,phi,err_phi,L_s, err_L_s, gam1,err_gam1,gam2,err_gam2= np.loadtxt(data_Dir + "/Quasar_LF_Hopkins_fit", unpack = True)
                err_phi = err_phi * np.log(10) * 10**phi
                L_s = L_s + 33.6
                err_L_s = err_L_s * np.log(10) * 10**L_s
		L_s = 10**L_s
                phi = 10**phi
                _sel = np.where(zz == redshift)
                if(len(_sel[0])==0):
                        return 0, 0
                else:
                        luminosities = np.arange(40,50,0.05)
                        luminosities = 10 ** luminosities
                        phi_luminosities = Hopkins_fit(phi[_sel][0],L_s[_sel][0],gam1[_sel][0],gam2[_sel][0],luminosities)
                        phi_max, phi_min = Hopkins_fit_err(phi[_sel][0],err_phi[_sel][0],L_s[_sel][0],err_L_s[_sel][0],gam1[_sel][0],err_gam1[_sel][0],gam2[_sel][0],err_gam2[_sel][0],luminosities)
                        return luminosities, phi_luminosities,  phi_max, phi_min




def Plot_Fit_Haring_Rix_2004():
	print 'Haring & Rix (2004) is a fit for 30 galaxies in the local universe'
	alpha = 8.20
	err_alpha = 0.10
	beta = 1.12
	err_beta = 0.06
	Bulge_mass = np.arange(0.006,40,0.001)
	Bulge_mass = np.log10(Bulge_mass)
	MBH_log = alpha + beta * Bulge_mass
	#MBH_log_lower_lim = (alpha-err_alpha) + (beta-err_beta) * Bulge_mass
	#MBH_log_up_lim = (alpha+err_alpha) + (beta+err_beta) * Bulge_mass
	MBH_lim = ((MBH_log * err_beta)**2 + (err_alpha)**2)**0.5
	MBH_log_up_lim = MBH_log - MBH_lim
	MBH_log_lower_lim = MBH_log + MBH_lim
	Bulge_mass = Bulge_mass + 11
	return Bulge_mass, MBH_log, MBH_log_lower_lim, MBH_log_up_lim


def Plot_Fit_McConell():
	print 'Nicholas J. McConnell and Chung-Pei Ma fit --> arxive: http://iopscience.iop.org/article/10.1088/0004-637X/764/2/184/pdf'
	alpha = 8.46
	err_alpha = 0.08
	beta = 1.05
	err_beta = 0.11
	Bulge_mass = np.arange(0.006,40,0.001)	# /1e11 M_sum
	Bulge_mass = np.log10(Bulge_mass)
	MBH_log =  alpha + beta*Bulge_mass
        MBH_lim = ((MBH_log * err_beta)**2 + (err_alpha)**2)**0.5
        MBH_log_up_lim = MBH_log - MBH_lim
        MBH_log_lower_lim = MBH_log + MBH_lim
        Bulge_mass = Bulge_mass + 11
        return Bulge_mass, MBH_log, MBH_log_lower_lim, MBH_log_up_lim
	
	
def Median_cloud_points(x,y):
	nbins = 30.
	sb = (11.5 - 9.0)/nbins
	x_points = []
	medians = []
	p84 = []
	p16 = []
	for i in np.arange(0,nbins,1):
		v_min = 9 + i*sb
		v_max = 9 + (i+1)*sb
		dummy = np.where((x>=v_min) & (x<v_max))
		x_points.append(v_max - (sb/2.0))
		medians.append(np.percentile(y[dummy],50))
		p16.append(np.percentile(y[dummy],16))
		p84.append(np.percentile(y[dummy],84))
	x_points = np.array(x_points)
	medians = np.array(medians)
	p84 = np.array(p84)
	p16 = np.array(p16)
	return x_points, medians, p16, p84



def plot_BH_MassFunction_data(redshift_desired,Ma2004=False, Me2008=False, Sha2004 = False):
	from Config import *
        if(Ma2004==True):
                BHmassMarconi2004, phiMarconi2004, phi_maxMarconi2004, phi_minMarconi2004 = np.loadtxt( data_Dir + '/Marconi2004', unpack =True)
                BHmassUeda_et_al, phiUeda_et_al, phi_maxUeda_et_al, phi_minUeda_et_al = np.loadtxt( data_Dir + '/Ueda_et_al', unpack =True)

                return BHmassMarconi2004, phiMarconi2004, phi_maxMarconi2004, phi_minMarconi2004, BHmassUeda_et_al, phiUeda_et_al, phi_maxUeda_et_al, phi_minUeda_et_al
        elif(Me2008==True):
                #zz = np.array([0.10,0.30,0.60,1.0,1.5,2.0,3.0,4.0,5.0])
                zz = np.array([0.0,0.30,0.60,1.0,1.5,2.0,3.0,4.0,5.0]) # Be carefullly is not z = 0, is z = 0.1 but for plots it's right
                name_data = '/mf_mh08.dat'
                BHmass, phi_max, BHmass_dummy, phi_min = np.loadtxt(data_Dir + name_data, unpack =True)
                #import pdb as pdb
                #pdb.set_trace()
                BHmass = np.split(BHmass,9)
                phi_max = np.split(phi_max,9)
                phi_min = np.split(phi_min,9)
                pos_redshift_desired = np.where(zz == redshift_desired)
                if(len(pos_redshift_desired[0])==0):
                        pos_redshift_desired_min = np.where(zz<redshift_desired)
                        pos_redshift_desired_max = np.where(zz>redshift_desired)
                        zz_min = zz[pos_redshift_desired_min[0][-1]]
                        zz_max = zz[pos_redshift_desired_max[0][0]]
                        diff_max = abs(zz_max-redshift_desired)
                        diff_min = abs(zz_min-redshift_desired)
                        if(diff_max<diff_min):
                                return BHmass[pos_redshift_desired_max[0][0]], phi_max[pos_redshift_desired_max[0][0]], phi_min[pos_redshift_desired_max[0][0]]
                        else:
                                return BHmass[pos_redshift_desired_min[0][-1]], phi_max[pos_redshift_desired_min[0][-1]], phi_min[pos_redshift_desired_min[0][-1]]

                else:
                        print 'pos_redshift_desired', pos_redshift_desired[0], 'redshift_desired', redshift_desired, 'zz[pos_redshift_desired]', zz[pos_redshift_desired]
                        return BHmass[pos_redshift_desired[0]], phi_max[pos_redshift_desired[0]], phi_min[pos_redshift_desired[0]]

        elif(Sha2004==True):
                BHmassShankar2004, phi_minShankar2004, phi_Shankar2004, phi_maxShankar2004 = np.loadtxt( data_Dir + '/Shankar_2004', unpack =True)
                return BHmassShankar2004, phi_minShankar2004, phi_Shankar2004, phi_maxShankar2004
        else:
                return -1,-1,-1,-1


def HistFixedMhalo(mstellar, mhalo):
        print mstellar
        masahalobin=[]
        inter=[]
        nBines=np.arange(9,14.75,0.25)
        MedianasHalo = np.empty(len(nBines)-1)
        HaloMass = np.empty(len(nBines)-1)
        MedianasHalo.fill(-19)
        for i in np.arange(1,len(nBines),1):
                mini = nBines[i-1]
                maxi = nBines[i]
                _sel = np.where((mhalo>=mini) & (mhalo<maxi))
                MedianasHalo[i-1] = np.median(mstellar[_sel])
                HaloMass[i-1] = nBines[i-1] + (0.25/2.)
        return HaloMass[np.where(MedianasHalo>-19)], MedianasHalo[np.where(MedianasHalo>-19)]

