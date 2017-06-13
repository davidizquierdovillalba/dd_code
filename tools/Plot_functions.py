import matplotlib.colors as colors
import pylab as plt
import matplotlib.lines as mlines

from Config import *
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
from dictionaries import *
from Useful_func import *

def PvsP_densityplot():

        P1 = Prop(wtp[0])    # property to be plotted on x-axis
        P2 = Prop(wtp[1])    # property to be plotted on y-axis
        
        print '\n###########################################################################\n'
        print '    PLOT:  ', P1, 'vs', P2
        print '    Redshifts used in plots:', ztoplot
        print '    Output prefix:', LGparams['FileNameGalaxies']
        print '\n############################################################################\n'

        Folders_to_do = get_outFolder_names(plot_last_run, LGparams, LGout_dir, here) # Files we have to plot
        Plotlog = True    # to be removed?
        nx = Labels(wtp[0],Plotlog,removeh)
	ny = Labels(wtp[1],Plotlog,removeh)
	row = [0,0,1,1]
        col = [0,1,0,1]
        f, axarr = plt.subplots(2, 2, figsize = (12,8))
	if(GALTREE == True):
		a = read_tree(out_LGfiles_dir,LGparams['FileNameGalaxies'],LGparams['FirstFile'],LGparams['LastFile'],PropertiesToRead_tree,LGalaxiesStruct)
		gg = a[1]
		del a
	maximX = np.empty(len(ztoplot))
	maximY = np.empty(len(ztoplot))
	minimX = np.empty(len(ztoplot))
	minimY = np.empty(len(ztoplot))

	for i,j in zip(ztoplot,np.arange(0,len(ztoplot),1)):

		if(GALTREE == False):

			filepref = LGparams['FileNameGalaxies'] + np.str(zdict(i, MII))
			a = read_snap(out_LGfiles_dir,filepref,LGparams['FirstFile'],LGparams['LastFile'],PropertiesToRead,LGalaxiesStruct)
			gg = a[3]
			p1 = gg[P1]
			p2 = gg[P2]
		else:
			_select = np.where(gg['SnapNum'] == snapdict(i))
			p1 = np.log10(gg[_select][P1])
			p2 = np.log10(gg[_select][P2])

		if(Plotlog == True):
			_select = np.where((p1>0) & (p2>0))
			if(Factor10[0]==True):
				p1 = 10 + np.log10(p1[_select])
			else:
				p1 = np.log10(p1[_select])
				
			if(Factor10[1]==True):
				p2 = 10 + np.log10(p2[_select])
			else:
				p2 = np.log10(p2[_select])

			if(removeh == True):
				if(P1 != 'Sfr'):
					p1 = p1 - np.log10(cosmo.h)
				if(P2 != 'Sfr'):
					p2 = p2 - np.log10(cosmo.h)

		else:
			_select = np.where((p1>0) & (p2>0))
			if(Factor10[0]==True):
				p1 = p1[_select] * 1.**10
			else:   
				p1 = p1[_select]
				
			if(Factor10[1]==True):
				p2 = p2[_select] * 1.**10
			else:
				p2 = p2[_select]
			
			if(removeh == True):
				if(P1 != 'Sfr'):
					p1 = p1/cosmo.h
				if(P2 != 'Sfr'):	
					p2 = p2/cosmo.h


		H , y, x = np.histogram2d(p2,p1,bins = DMap_RES, normed = True)
		d = axarr[row[j],col[j]].pcolor(x, y, H, norm = colors.LogNorm(),cmap='jet', label = '$z$ = ' + zdict(i, MII))
		d.set_clim(vmin=0.01, vmax=2)
		if(row[j] == 1 and col[j]==1):
			position=f.add_axes([0.90,0.12,0.02,0.75])
			f.colorbar(d,cax=position, label = 'N/$\mathrm{N_{max}}$')
		leg = axarr[row[j],col[j]].legend(loc = 'lower right',fontsize = 18, handlelength=0, handletextpad=0, fancybox=True)
		for item in leg.legendHandles:
			item.set_visible(False)
	
		maximX[j] = max(p1)
		minimX[j] = min(p1)
		maximY[j] = max(p2)
		minimY[j] = min(p2)
	
		if(row[j] == 0):
			axarr[row[j],col[j]].set_xticklabels([''])
		if(col[j]>0):
			axarr[row[j],col[j]].set_yticklabels([''])
	
	for i in np.arange(0,len(ztoplot),1):
		axarr[row[i],col[i]].set_ylim(min(minimY),max(maximY))
		axarr[row[i],col[i]].set_xlim(min(minimX),max(maximX))
	
	seedMasses = np.array(seedMasses)
	acc_Model = np.array(acc_Model, int)
	for kk in range(len(seedMasses)):
		c = mlines.Line2D([], [], color=sfcol[kk], linestyle = '-', linewidth = 2, label = 'BH seed = %2.1e'%float(seedMasses[kk]) + r'$\rm [M_{\odot}]$' ' Model ' + np.str(acc_Model[kk]))
	handles.append(c)
	labels = [h.get_label() for h in handles]

	f.canvas.draw()
	labels = [item.get_text() for item in axarr[1,0].get_xticklabels()]
	labels[-1] = ' '
	axarr[1,0].set_xticklabels(labels)
	f.subplots_adjust(wspace=0.0)
	f.subplots_adjust(hspace=0.0)
	f.subplots_adjust(left=0.11)
	f.subplots_adjust(right=0.88)
	f.subplots_adjust(bottom=0.09)
	f.text(0.0325, 0.5, ny, va='center', rotation='vertical', fontsize = 22)
	f.text(0.45, 0.025, nx, va='center', rotation='horizontal', fontsize = 22)
	del p1,p2,gg
	plt.savefig(dens_ns, aspect='equal')
	plt.show()


def MoL_func():
        P1 = Prop(wtp[0])    # property to be used for luminosity/mass function
        
        print '\n###########################################################################\n'
        print '    PLOT:  ', P1, ' function'
        print '    Redshifts used in plots:', ztoplot
        print '    Output prefix:', LGparams['FileNameGalaxies']
        print '\n############################################################################\n'

	Folders_to_do = get_outFolder_names(plot_last_run, LGparams, LGout_dir, here) # Files we have to plot
        
	nx, ny = MoL_labels(wtp[0],Plotlog,removeh)
	
	maximX = np.empty((len(Folders_to_do),len(ztoplot)))
	minimX = np.empty((len(Folders_to_do),len(ztoplot)))
	maximY = np.empty((len(Folders_to_do),len(ztoplot)))
	minimY = np.empty((len(Folders_to_do),len(ztoplot)))
	f, axarr = plt.subplots(2, 2, figsize = (10,6))
	row = [0,0,1,1]
	col = [0,1,0,1]
	sfcol = plt.cm.spectral(np.linspace(0.2,0.9,len(Folders_to_do)))
        seedMasses = []
        acc_Model = []
        handles = []


	for nams,kk in zip(Folders_to_do, range(len(Folders_to_do))):
                loc_parFile = nams + LG_inParFile[LG_inParFile.find('input/input')+6:]
                LGparams_loc = read_LG_inParamFile(loc_parFile, LG_output_z, params_to_read)
                seedMasses.append(str(LGparams_loc['BlackHoleSeedMass']))
                acc_Model.append(str(LGparams_loc['AccretionModel']))
		if(GALTREE == True):
			a = read_tree(nams,LGparams_loc['FileNameGalaxies'],LGparams_loc['FirstFile'],LGparams_loc['LastFile'],PropertiesToRead_tree,LGalaxiesStruct)
			gg = a[1]

		if(removeh == True):
			Volume = ((BoxSize/cosmo.h)**3.0) * (LGparams_loc['LastFile'] - LGparams_loc['FirstFile'] + 1) / MaxTreeFiles # Mpc^3
		else:
			Volume = (BoxSize**3.0) * (LGparams_loc['LastFile'] - LGparams_loc['FirstFile'] + 1) / MaxTreeFiles # Mpc^3 h^-3


		for i,j in zip(ztoplot,np.arange(0,len(ztoplot),1)):
			if(GALTREE == False):
				filepref = LGparams_loc['FileNameGalaxies'] + np.str(zdict(i, MII))
				a = read_snap(nams,filepref,LGparams_loc['FirstFile'],LGparams_loc['LastFile'],PropertiesToRead,LGalaxiesStruct)
				gg = a[3]

				if(Factor10[0] == True):
					if(Plotlog == True):
						pp = (10+ np.log10(gg[np.where(gg[P1]>0)][P1]))
					else:
						pp = gg[np.where(gg[P1]>0)][P1] * 1.**10
				else:
					if(Plotlog == True):
						pp = np.log10(gg[np.where(gg[P1]>0)][P1])
					else:
						pp = gg[np.where(gg[P1]>0)][P1]

				if(removeh == True):
					if(Plotlog == True):
						pp = pp - np.log10(cosmo.h)
					else:
						pp = pp/cosmo.h

			else:
				_select = np.where((gg['SnapNum'] == snapdict(i)) & (gg[P1]>0))

				if(Factor10[0] == True):
					if(Plotlog == True):
						pp = (10+ np.log10(gg[_select][P1]))
					else:
						pp = gg[_select][P1] * 1.**10
				else:
					if(Plotlog == True):
						pp = np.log10(gg[_select][P1])
					else:
						pp = gg[_select][P1] 

				if(removeh == True):
					if(Plotlog == True):
						pp = pp - np.log10(cosmo.h)
					else:
						pp = pp/cosmo.h

			if Plotlog == False:
				sb = get_sb(max(pp), min(pp))
			else:
				sb = sizebin

			pp_c, phi = hist(pp,sb)
			phi = phi/(Volume * sb)
			if(plot_obs_data==True and kk == 0):
				if(wtp[0] == 0):
					BHmass, phi_max, phi_min = plot_BH_MassFunction_data(ztoplot[j])
					phi_data = (phi_max - phi_min)/2.
					#axarr[row[j],col[j]].errorbar(BHmass, phi_data, yerr = [phi_min,phi_max], color = 'k', fmt = 'o', markersize=3.5, label = '$z$ = ' + zdict(i, MII))
					axarr[row[j],col[j]].fill_between(BHmass, phi_min,phi_max, olor = 'k', alpha = 0.5, label = '$z$ = ' + zdict(i, MII))
			maximX[kk,j] = max(pp_c[pp_c > 0])
			minimX[kk,j] = min(pp_c[pp_c > 0])
			maximY[kk,j] = max(phi[phi != 0])
			minimY[kk,j] = min(phi[phi != 0])
			
			axarr[row[j],col[j]].plot(pp_c,phi,color = sfcol[kk], linewidth=2, label = '$z$ = ' + zdict(i, MII))
			axarr[row[j],col[j]].set_yscale('log')

			if(kk == 0):
                                leg = axarr[row[j],col[j]].legend(loc = 'upper left',fontsize = 18, handlelength=0, handletextpad=0, fancybox=True)
                                for item in leg.legendHandles:
                                        item.set_visible(False)
                                if (row[j] == 0) and (col[j] == 0):
                                        axarr[row[j],col[j]].add_artist(leg)
                                        
                                        
                        if(row[j] == 0):
                                axarr[row[j],col[j]].set_xticklabels([''])
                        if(col[j]>0):
                                axarr[row[j],col[j]].set_yticklabels([''])


	for i in np.arange(0,len(ztoplot),1):
                axarr[row[i],col[i]].set_xlim( np.min(np.min(minimX, axis=1), axis=0), np.max(np.max(maximX, axis=1), axis=0) )
                axarr[row[i],col[i]].set_ylim(2*np.min(np.min(minimY, axis=1), axis=0), 1.15*np.max(np.max(maximY, axis=1), axis=0) ) # Check the factor 5

        seedMasses = np.array(seedMasses)
	acc_Model = np.array(acc_Model, int)
        for kk in range(len(seedMasses)):
                c = mlines.Line2D([], [], color=sfcol[kk], linestyle = '-', linewidth = 2, label = 'BH seed = %2.1e'%float(seedMasses[kk]) + r'$\rm [M_{\odot}]$' ' Model ' + np.str(acc_Model[kk]))
                handles.append(c)
        labels = [h.get_label() for h in handles]
        #axarr[row[0],col[0]].legend(handles, labels, bbox_to_anchor=(1.55, 0.85),loc = "upper right", fontsize = 11)
        axarr[row[0],col[0]].legend(handles, labels, loc = "upper right", fontsize = 11)
                
	fig = plt.tight_layout()
	f.subplots_adjust(wspace=0)
	f.subplots_adjust(hspace=0.0)
	f.subplots_adjust(left=0.11)
	f.subplots_adjust(bottom=0.09)
	f.text(0.0175, 0.5, ny, va='center', rotation='vertical', fontsize = 20)
	f.text(0.45, 0.025, nx, va='center', rotation='horizontal', fontsize = 20)
        #f.canvas.draw()
        #labels = [item.get_text() for item in axarr[1,0].get_xticklabels()]
        #labels[-1] = ' '
        #axarr[1,0].set_xticklabels(labels)
        #f.canvas.draw()
        #labels = [item.get_text() for item in axarr[0,0].get_yticklabels()]
        #labels[-1] = ' '
        #axarr[0,0].set_yticklabels(labels)

	plt.savefig(func_ns)
	plt.show()


def generic_histo(): #SHOULD BE GENERALIZED TO PLOT MULTIPLE HISTOGRAMS
        P1 = Prop(wtp[0])    # property to be used for luminosity/mass function
        
        print '\n###########################################################################\n'
        print '    PLOT:  ', P1, ' histogram'
        print '    Redshifts used in plots:', ztoplot
        print '    Output prefix:', LGparams['FileNameGalaxies']
        print '\n############################################################################\n'

        # get_outFolder_names routine asks the user to select the folders
        # containing input files needed for the desired plots (see Useful_func.py)
	Folders_to_do = get_outFolder_names(plot_last_run, LGparams, LGout_dir, here)

        # we set, Plotlog == False if, by chance, it is not (only for seedTime plots)
        Plotlog = False
        removeh = False
	nx = Labels(wtp[0], Plotlog, removeh)
	
	maximX = np.empty((len(Folders_to_do),len(ztoplot)))
	minimX = np.empty((len(Folders_to_do),len(ztoplot)))
	maximY = np.empty((len(Folders_to_do),len(ztoplot)))
	minimY = np.empty((len(Folders_to_do),len(ztoplot)))
	f, axarr = plt.subplots(2, 2, figsize = (12,8))
	row = [0,0,1,1]
	col = [0,1,0,1]
	sfcol = plt.cm.jet(np.linspace(0.2,0.9,len(Folders_to_do)))
        seedMasses = []
        handles = []

        # loop on all the desired folders
	for nams, kk in zip(Folders_to_do, range(len(Folders_to_do))):
                loc_parFile = nams + LG_inParFile[LG_inParFile.find('input/input')+6:]       # parameter file name, saved in the corresponding folder
                LGparams_loc = read_LG_inParamFile(loc_parFile, LG_output_z, params_to_read) # reads the parameter file corresponding to the right folder
                seedMasses.append(str(LGparams_loc['BlackHoleSeedMass']))                    # each folder should have LG ouputs obtained with different seed mass
		if(GALTREE == True):
			a = read_tree(nams,LGparams_loc['FileNameGalaxies'],LGparams_loc['FirstFile'],LGparams_loc['LastFile'],PropertiesToRead_tree,LGalaxiesStruct)
			gg = a[1]

		for i, j in zip( ztoplot, np.arange(0,len(ztoplot),1) ):  # loops over the 4 selected redshift values to be used for the plots
                        if(GALTREE == False):
				filepref = LGparams_loc['FileNameGalaxies'] + np.str(zdict(i, MII))
				a = read_snap(nams, filepref, LGparams_loc['FirstFile'], LGparams_loc['LastFile'], PropertiesToRead, LGalaxiesStruct)
				gg = a[3]

                                gg['BlackHoleMass'] = gg['BlackHoleMass'] / (Hubble_h*1.e-10)
                                _select = ((gg[P1] > 0) & (gg['BlackHoleMass'] == LGparams_loc['BlackHoleSeedMass']) )
                                _NOTselect = ((gg[P1] > 0) & (gg['BlackHoleMass'] > LGparams_loc['BlackHoleSeedMass']) )
				if(Factor10[0] == True):
                                        #pp = np.log10( gg[_select][P1] * 1.**10 )       ###########
                                        #pp1 = np.log10( gg[_NOTselect][P1] * 1.**10 )   ###########
                                        pp = gg[_select][P1] * 1.**10       ###########
                                        pp1 = gg[_NOTselect][P1] * 1.**10   ###########
				else:
                                        #pp = np.log10( gg[_select][P1] )       ###########
                                        #pp1 = np.log10( gg[_NOTselect][P1] )   ###########
                                        pp = gg[_select][P1]       ###########
                                        pp1 = gg[_NOTselect][P1]   ###########
			else:
				_select = np.where( (gg['SnapNum'] == snapdict(i)) & (gg[P1] > 0))

				if(Factor10[0] == True):
                                        pp = gg[_select][P1] * 1.**10
				else:
                                        pp = gg[_select][P1]
                                        
                        sb = get_sb(np.log10(max(pp)), np.log10(min(pp)))   #automatic size bins (in Useful_func.py)                                        
                        sb1 = get_sb(np.log10(max(pp1)), np.log10(min(pp1)))   #automatic size bins (in Useful_func.py)

			pp_c, phi = hist(pp,sb)
			pp_c1, phi1 = hist(pp1,sb)
			maximX[kk,j] = max(pp_c[pp_c > 0])
			minimX[kk,j] = min(pp_c[pp_c > 0])
			maximY[kk,j] = max(phi[phi != 0])
			minimY[kk,j] = min(phi[phi != 0])
			
			axarr[row[j],col[j]].plot(pp_c, phi, color = sfcol[kk], alpha=0.8, linewidth=2, label = '$z$ = ' + zdict(i, MII))
			axarr[row[j],col[j]].plot(pp_c1, phi1, color = 'r', alpha=0.8, linewidth=2)
			axarr[row[j],col[j]].set_yscale('log')

                        # The plot is done. Now it's just a matter of limits, labels and legends 
			if(kk == 0):
                                leg = axarr[row[j],col[j]].legend(loc = 'upper left',fontsize = 18, handlelength=0, handletextpad=0, fancybox=True)
                                for item in leg.legendHandles:
                                        item.set_visible(False)
                                if (row[j] == 0) and (col[j] == 0):
                                        axarr[row[j],col[j]].add_artist(leg)
                        if(row[j] == 0):
                                axarr[row[j],col[j]].set_xticklabels([''])
                        if(col[j]>0):
                                axarr[row[j],col[j]].set_yticklabels([''])

        for i in np.arange(0,len(ztoplot),1):
                axarr[row[i],col[i]].set_xlim( [np.min(np.min(minimX, axis=1), axis=0), np.max(np.max(maximX, axis=1), axis=0)] )
                
        # plot a line in the legend for each seedMass
        seedMasses = np.array(seedMasses)
        for kk in range(len(seedMasses)):
                c = mlines.Line2D([], [], color=sfcol[kk], linestyle = '-', linewidth = 2, label = 'BHM = seed M = %2.1e'%float(seedMasses[kk]))
                d = mlines.Line2D([], [], color='r', linestyle = '-', linewidth = 2, label = 'BHM != seed M = %2.1e'%float(seedMasses[kk]))
                handles.append(c)
                handles.append(d)
        labels = [h.get_label() for h in handles]
        axarr[row[0],col[0]].legend(handles, labels, loc = "upper right", fontsize = 8)

	fig = plt.tight_layout()
	f.subplots_adjust(wspace=0)
	f.subplots_adjust(hspace=0.0)
	f.subplots_adjust(left=0.11)
	f.subplots_adjust(bottom=0.09)
	f.text(0.0175, 0.5, '# counts', va='center', rotation='vertical', fontsize = 22)
	f.text(0.45, 0.025, nx, va='center', rotation='horizontal', fontsize = 22)
        
	f.canvas.draw()
        labels = axarr[row[i],col[i]].get_xticklabels()
	labels[-1] = ' '
	axarr[1,0].set_xticklabels(labels)
	labels = [item.get_text() for item in axarr[0,0].get_yticklabels()]
	labels[-1] = ' '
	axarr[0,0].set_yticklabels(labels)
        
	plt.savefig(func_ns)
	plt.show()


def Quasar_LumFunction(): #Generate the quasars luminosity functions 

        print '\n###########################################################################\n'
        print '    PLOT:  Quasar_LumFunction function'
        print '    Redshifts used in plots:', ztoplot
        print '    Output prefix:', LGparams['FileNameGalaxies']
        print '\n############################################################################\n'

        Folders_to_do = get_outFolder_names(plot_last_run, LGparams, LGout_dir, here) # Files we have to plot


        f, axarr = plt.subplots(2, 2, figsize = (10,7))
        row = [0,0,1,1]
        col = [0,1,0,1]
        #sfcol = plt.cm.spectral(np.linspace(0.2,0.9,len(Folders_to_do)))
	col_names = ['magenta','blue', 'red', 'green']
        seedMasses = []
        acc_Model = []
        handles = []
	legend_names = []
        legData = []

        for nams,kk in zip(Folders_to_do, range(len(Folders_to_do))):
                loc_parFile = nams + LG_inParFile[LG_inParFile.find('input/input')+6:]
                LGparams_loc = read_LG_inParamFile(loc_parFile, LG_output_z, params_to_read)
                seedMasses.append(str(LGparams_loc['BlackHoleSeedMass']))
                acc_Model.append(str(LGparams_loc['AccretionModel']))
		Volume = ((BoxSize/cosmo.h)**3.0) * (LGparams_loc['LastFile'] - LGparams_loc['FirstFile'] + 1) / MaxTreeFiles # Mpc^3

                for i,j in zip(ztoplot,np.arange(0,len(ztoplot),1)):
			filepref = LGparams_loc['FileNameGalaxies'] + np.str(zdict(i, MII))
			a = read_snap(nams,filepref,LGparams_loc['FirstFile'],LGparams_loc['LastFile'],PropertiesToRead,LGalaxiesStruct)
			gg = a[3]
			gg = gg[np.where(gg['QuasarLum']>0.0)]
			gg['QuasarLum'] = np.log10(gg['QuasarLum']) + 40 # The code give us L_bol / 10^40 [erg/s]
			bins=np.arange(min(gg['QuasarLum']),max(gg['QuasarLum']),0.25)	
			hist, b_edges = np.histogram(gg['QuasarLum'],bins)
			log_Lum = (b_edges[:-1] + b_edges[1:]) / 2.0
			hist = np.array(hist, dtype=np.float64)
                        phi = hist/(Volume * abs(bins[1]- bins[0])) # dex^-1 Mpc^-3
			if(Quasar_LumFunctions == True):
				if(Quasar_data == True):
					Lum, phi_L, phi_L_err = Add_Quasar_LF_data(i,data=True)
					#Lum = np.log10(Lum)
					phi_L_err_min = phi_L - phi_L_err
					phi_L_err_max = phi_L + phi_L_err
					legHop = axarr[row[j],col[j]].errorbar(Lum, 10**phi_L, yerr =[(10**phi_L - 10**phi_L_err_min), (10** phi_L_err_max - 10**phi_L)],color = 'k', fmt = 'o', markersize=4.5)
				if(Quasar_data_fit == True):
					Lum, phi_L,phi_mx, phi_mi = Add_Quasar_LF_data(i,fit=True)
					Lum = np.log10(Lum)
					legHop = axarr[row[j],col[j]].plot(Lum, phi_L, color = 'k', linestyle = '-')
					axarr[row[j],col[j]].fill_between(Lum,phi_mi,phi_mx, color = 'k', alpha = 0.15)
				if(kk==0 and j==0):
					legData.append(legHop)
					legend_names.append('From Hopkins et al. 2007')
			if(int(str(LGparams_loc['AccretionModel'])) == 0):
				cc = col_names[0] 
			elif(int(str(LGparams_loc['AccretionModel'])) == 1):
				cc = col_names[1]
			elif(int(str(LGparams_loc['AccretionModel'])) == 2):
				cc = col_names[2]
			else:
				cc = col_names[3]
                        axarr[row[j],col[j]].plot(log_Lum,phi,color = cc, linewidth=2.25, label = '$z$ = ' + zdict(i, MII))
                        axarr[row[j],col[j]].set_yscale('log')

                        if(kk == 0):
                                leg = axarr[row[j],col[j]].legend(loc = 'upper right',fontsize = 18, handlelength=0, handletextpad=0, fancybox=True)
                                for item in leg.legendHandles:
                                        item.set_visible(False)
                                if (row[j] == 0) and (col[j] == 0):
                                        axarr[row[j],col[j]].add_artist(leg)

                        if(row[j] == 0):
                                axarr[row[j],col[j]].set_xticklabels([''])
                        if(col[j]>0):
                                axarr[row[j],col[j]].set_yticklabels([''])


        for i in np.arange(0,len(ztoplot),1):
                axarr[row[i],col[i]].set_xlim(42,47.9)
                axarr[row[i],col[i]].set_ylim(1.1e-10,1e-2) # Check the factor 5


	leg = axarr[row[0],col[0]].legend(legData,legend_names, loc = 'center left', fontsize = 13.5)
        axarr[row[0],col[0]].add_artist(leg)

        seedMasses = np.array(seedMasses)
        acc_Model = np.array(acc_Model, int)
        for kk in range(len(acc_Model)):
		if(acc_Model[kk] == 0):
			cc = col_names[0]
		elif(acc_Model[kk] == 1):
			cc = col_names[1]
		elif(acc_Model[kk] == 2):
			cc = col_names[2]
		else:
			cc = col_names[3]
                c = mlines.Line2D([], [], color=cc, linestyle = '-', linewidth = 2, label = 'BH seed = %2.1e'%float(seedMasses[kk]) + r'$\rm [M_{\odot}]$' ' Model ' + np.str(acc_Model[kk]))
                handles.append(c)
        labels = [h.get_label() for h in handles]
        axarr[row[0],col[0]].legend(handles, labels, loc = "lower left", fontsize = 11)

        fig = plt.tight_layout()
        f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace=0.0)
        f.subplots_adjust(left=0.11)
        f.subplots_adjust(bottom=0.09)
        f.text(0.0175, 0.5, r'$\phi \rm[dex^{-1} Mpc^{-3}]$', va='center', rotation='vertical', fontsize = 20)
        f.text(0.45, 0.025, r'$\mathrm{log_{10}(L_{bol}[erg/s])}$', va='center', rotation='horizontal', fontsize = 20)

        plt.savefig(plots_dir + 'Quasar_LF.pdf')
        plt.show()


def BH_mass_Function():

        print '\n###########################################################################\n'
        print '    PLOT:  Black Hole mass function'
        print '    Redshifts used in plots:', ztoplot
        print '    Output prefix:', LGparams['FileNameGalaxies']
        print '\n############################################################################\n'

        Folders_to_do = get_outFolder_names(plot_last_run, LGparams, LGout_dir, here) # Files we have to plot


        f, axarr = plt.subplots(2, 2, figsize = (10,7))
        row = [0,0,1,1]
        col = [0,1,0,1]
        sfcol = plt.cm.spectral(np.linspace(0.2,0.9,len(Folders_to_do)))
        seedMasses = []
        acc_Model = []
        handles = []
	legend_names = []
	legData = []
	BHg = []

        for nams,kk in zip(Folders_to_do, range(len(Folders_to_do))):
                loc_parFile = nams + LG_inParFile[LG_inParFile.find('input/input')+6:]
                LGparams_loc = read_LG_inParamFile(loc_parFile, LG_output_z, params_to_read)
                seedMasses.append(str(LGparams_loc['BlackHoleSeedMass']))
                acc_Model.append(str(LGparams_loc['AccretionModel']))
		BHg.append(str(LGparams_loc['BlackHoleGrowthRate']))
		Volume = ((BoxSize/cosmo.h)**3.0) * (LGparams_loc['LastFile'] - LGparams_loc['FirstFile'] + 1) / MaxTreeFiles # [Mpc^3]

                for i,j in zip(ztoplot,np.arange(0,len(ztoplot),1)):
			filepref = LGparams_loc['FileNameGalaxies'] + np.str(zdict(i, MII))
			a = read_snap(nams,filepref,LGparams_loc['FirstFile'],LGparams_loc['LastFile'],PropertiesToRead,LGalaxiesStruct)
			gg = a[3]
			pp = gg['BlackHoleMass'] * 1e10 / cosmo.h
			pp = pp[np.where(pp>0.)]
			pp = np.log10(pp) # BHMass [M_sun]
			bins=np.arange(min(pp),max(pp),0.20)
                        hist, b_edges = np.histogram(pp,bins)
                        pp_c= (b_edges[:-1] + b_edges[1:]) / 2.0
                        phi = hist/(Volume * abs(bins[1]- bins[0])) # [dex^-1 Mpc^3]

			################################## DATA (config file details)#########################################################
                        if(plot_obs_data==True and kk == 0):
				if(Marconi2004==True):
					if(ztoplot[j]==0 and kk ==0): # I only add this plot in the redshift z = 0. (local universe)
						BHmassM, phi_M, phi_maxM, phi_minM, BHmassU, phi_U, phi_maxU, phi_minU = plot_BH_MassFunction_data(ztoplot[j], Ma2004=True)
						#legMarconi = axarr[row[j],col[j]].errorbar(BHmassM, 10**phi_M, yerr =[(10**phi_M - 10**phi_minM), (10** phi_maxM - 10**phi_M)],color = 'k', fmt = 'o', markersize=4.5)
						legMarconi = axarr[row[j],col[j]].plot(BHmassM,10**phi_M, color = 'black', linewidth = 2.25, linestyle = '--')
						axarr[row[j],col[j]].fill_between(BHmassM,10**phi_minM, 10**phi_maxM, color = 'k', alpha = 0.15)
						legData.append(legMarconi[0])
						#legUeda = axarr[row[j],col[j]].errorbar(BHmassU, 10**phi_U, yerr =[(10**phi_U - 10**phi_minU), (10** phi_maxU - 10**phi_U)],color = 'g', fmt = 's', markersize=4.5)
						#legData.append(legUeda)
						legend_names.append('Marconi et al. 2004')
						#legend_names.append('Ueda et al. 2003')
				if(Shankar2004 == True):
                                        if(ztoplot[j]==0 and kk ==0): # I only add this plot in the redshift z = 0. (local universe)
                                                BHmassS, phi_minS, phi_S, phi_maxS = plot_BH_MassFunction_data(ztoplot[j],Sha2004 = True)
                                                legSa2004 = axarr[row[j],col[j]].errorbar(BHmassS, 10**phi_S, yerr =[(10**phi_S - 10**phi_minS), (10** phi_maxS - 10**phi_S)],color = 'brown', fmt = 'p', markersize=4.5)
                                                legData.append(legSa2004)
                                                legend_names.append('Shankar et al. 2004')
				if(Merloni2008==True):
					BHmass, phi_max, phi_min = plot_BH_MassFunction_data(ztoplot[j],Me2008=True)
					phi_data = (phi_max - phi_min)/2.
					legMerloni = axarr[row[j],col[j]].fill_between(BHmass, phi_min,phi_max, color = 'k', alpha = 0.35) 
					if(kk ==0):
						legData.append(legMerloni)
						legend_names.append('Merloni & Heinz et al. 2008')
			################################## End DATA #########################################################################

                        axarr[row[j],col[j]].plot(pp_c,phi,color = sfcol[kk], linewidth=2, label = '$z$ = ' + zdict(i, MII))
                        axarr[row[j],col[j]].set_yscale('log')
                        if(kk == 0):
                                leg = axarr[row[j],col[j]].legend(loc = 'upper right',fontsize = 14, handlelength=0, handletextpad=0, fancybox=True)
                                for item in leg.legendHandles:
                                        item.set_visible(False)
                                if (row[j] == 0) and (col[j] == 0):
                                        axarr[row[j],col[j]].add_artist(leg)
                        if(row[j] == 0):
                                axarr[row[j],col[j]].set_xticklabels([''])
                        if(col[j]>0):
                                axarr[row[j],col[j]].set_yticklabels([''])


        for i in np.arange(0,len(ztoplot),1):
                axarr[row[i],col[i]].set_xlim(6,10)
                axarr[row[i],col[i]].set_ylim(1e-8,1e-2) 

	############# legend names ###########
	print legData,legend_names
        leg = axarr[row[0],col[0]].legend(legData,legend_names, loc = 'center left', fontsize = 13.5)
        axarr[row[0],col[0]].add_artist(leg)

        seedMasses = np.array(seedMasses)
        acc_Model = np.array(acc_Model, int)
	BHg = np.array(BHg)
        for kk in range(len(seedMasses)):
                c = mlines.Line2D([], [], color=sfcol[kk], linestyle = '-', linewidth = 2, label = 'BH seed = %2.1e'%float(seedMasses[kk]) + r'$\rm [M_{\odot}]$' ' Model ' + np.str(acc_Model[kk]) + ' Growth rate: '+ np.str(BHg[kk]))
                handles.append(c)
        labels = [h.get_label() for h in handles]
        axarr[row[0],col[0]].legend(handles, labels, loc = "lower left", fontsize = 11)
	############# end legend names #######

        fig = plt.tight_layout()
        f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace=0.0)
        f.subplots_adjust(left=0.11)
        #f.subplots_adjust(right=1.)
        f.subplots_adjust(bottom=0.09)
        f.text(0.0175, 0.5, r'$\phi(\rm{M_{BH}}) \rm[dex^{-1} Mpc^{-3}]$', va='center', rotation='vertical', fontsize = 20)
        f.text(0.45, 0.025, r'$\rm log_{10}(M_{BH}[M_{\odot}])$', va='center', rotation='horizontal', fontsize = 20)
	f.canvas.draw()
        labels = [item.get_text() for item in axarr[0,0].get_yticklabels()]
        labels[1] = ' '
        axarr[0,0].set_yticklabels(labels)
	labels = [item.get_text() for item in axarr[1,0].get_xticklabels()]
        labels[-1] = ' '
        axarr[1,0].set_xticklabels(labels)

        plt.savefig(plots_dir + 'BHole_Mass_function.pdf')
        plt.show()





def Bulge_MBH_local_universe():
        print '\n###########################################################################\n'
        print '    PLOT:  Bulge vs BH mass local universe'
        print '    Redshifts used in plots: $z$ = 0'
        print '    Output prefix:', LGparams['FileNameGalaxies']
        print '\n############################################################################\n'

        Folders_to_do = get_outFolder_names(plot_last_run, LGparams, LGout_dir, here) # Files we have to plot

        sfcol = plt.cm.spectral(np.linspace(0.2,0.9,len(Folders_to_do)))
        seedMasses = []
        acc_Model = []
        handles = []
	BHg = []
        legend_names = []
        legData = []
	redshift = 0.0
	f, ax = plt.subplots()
        for nams,kk in zip(Folders_to_do, range(len(Folders_to_do))):
                loc_parFile = nams + LG_inParFile[LG_inParFile.find('input/input')+6:]
                LGparams_loc = read_LG_inParamFile(loc_parFile, LG_output_z, params_to_read)
                seedMasses.append(str(LGparams_loc['BlackHoleSeedMass']))
                acc_Model.append(str(LGparams_loc['AccretionModel']))
                BHg.append(str(LGparams_loc['BlackHoleGrowthRate']))

		filepref = LGparams_loc['FileNameGalaxies'] + np.str(zdict(redshift, MII))
		a = read_snap(nams,filepref,LGparams_loc['FirstFile'],LGparams_loc['LastFile'],PropertiesToRead,LGalaxiesStruct)
		gg = a[3]
		BH = gg['BlackHoleMass'] * 1e10 / cosmo.h # M_sun
		Bulge = gg['BulgeMass'] * 1e10 / cosmo.h # M_sun 
		BH = np.log10(BH) # log(BHMass) [M_sun]
		print 'Bulge', Bulge
		print 'BH', BH
		Bulge = np.log10(Bulge) # log10(Bulge) [M_sun]
		dummy = np.where((BH>0) & (Bulge>0))
		print BH[dummy]
		print Bulge[dummy]
		#ax.plot(Bulge,BH, color = sfcol[kk], marker = ',', linestyle = ' ')#, label= 'BH seed = %2.1e'%float(seedMasses[kk]) + r'$\rm [M_{\odot}]$' ' Model ' + np.str(acc_Model[kk]))
		pos_Bulge, med_BH, p16_BH, p84_BH = Median_cloud_points(Bulge[dummy],BH[dummy])
		ax.plot(pos_Bulge,med_BH,linestyle = '--', linewidth = 2.25, color = sfcol[kk])
		#ax.fill_between(pos_Bulge,p16_BH, p84_BH,color = sfcol[kk], alpha = 0.2)
		if(kk ==0): 
			Bulge_Fit, BH_fit, low, up = Plot_Fit_Haring_Rix_2004()
			ax.fill_between(Bulge_Fit, low, up,color = 'grey', alpha = 0.5)
			ax.plot(Bulge_Fit,BH_fit, color = 'k', linewidth = 2.25, linestyle = '-', label  = 'Haring & Rix (2004)')
	
        seedMasses = np.array(seedMasses)
        acc_Model = np.array(acc_Model, int)
        BHg = np.array(BHg)
        for kk in range(len(seedMasses)):
                c = mlines.Line2D([], [], color=sfcol[kk], marker = 'o', linestyle = ' ', linewidth = 2, label = 'BH seed = %2.1e'%float(seedMasses[kk]) + r'$\rm [M_{\odot}]$' ' Model ' + np.str(acc_Model[kk]) + ' Growth rate: '+ np.str(BHg[kk]))
                handles.append(c)
        labels = [h.get_label() for h in handles]
        sim_data_leg = ax.legend(handles, labels, loc = "lower right", fontsize = 9.5)
	ax.add_artist(sim_data_leg)
	ax.set_xlim(9,12.5)
	ax.set_ylim(6,10)
	ax.legend(loc = 'upper left',fontsize = 12.5)
	ax.text(11.5,8.0, '$z =$'+ np.str(redshift),fontsize = 22)
        fig = plt.tight_layout()
        f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace=0.0)
        f.subplots_adjust(left=0.11)
	f.subplots_adjust(bottom=0.09)
        #f.subplots_adjust(right=1.)
	f.text(0.0175, 0.5, r"$\rm log(M_{BH} [M_{\odot}])$", va='center', rotation='vertical', fontsize = 20)
        f.text(0.45, 0.025, r"$\rm log(M_{bulge} [M_{\odot}])$", va='center', rotation='horizontal', fontsize = 20)
        plt.savefig(plots_dir + 'Bulge_BH.pdf')
        plt.show()







