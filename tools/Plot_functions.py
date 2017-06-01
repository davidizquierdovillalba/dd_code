import matplotlib.colors as colors
import pylab as plt
import matplotlib.lines as mlines

from Config import *
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)
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

	Folders_to_do = tools.get_outFolder_names() # Files we have to plot
        
	nx, ny = MoL_labels(wtp[0],Plotlog,removeh)
	
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


	for nams,kk in zip(Folders_to_do, range(len(Folders_to_do))):
                loc_parFile = nams + LG_inParFile[LG_inParFile.find('input/input')+6:]
                LGparams_loc = read_LG_inParamFile(loc_parFile, LG_output_z, params_to_read)
                seedMasses.append(str(LGparams_loc['BlackHoleSeedMass']))
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
        for kk in range(len(seedMasses)):
                c = mlines.Line2D([], [], color=sfcol[kk], linestyle = '-', linewidth = 2, label = 'seed M = %2.1e'%float(seedMasses[kk]))
                handles.append(c)
        labels = [h.get_label() for h in handles]
        #axarr[row[0],col[0]].legend(handles, labels, bbox_to_anchor=(1.55, 0.85),loc = "upper right", fontsize = 11)
        axarr[row[0],col[0]].legend(handles, labels, loc = "upper right", fontsize = 11)
                
	f.canvas.draw()
	labels = [item.get_text() for item in axarr[1,0].get_xticklabels()]
	labels[-1] = ' '
	axarr[1,0].set_xticklabels(labels)
	labels = [item.get_text() for item in axarr[0,0].get_yticklabels()]
	labels[-1] = ' '
	axarr[0,0].set_yticklabels(labels)

	fig = plt.tight_layout()
	f.subplots_adjust(wspace=0)
	f.subplots_adjust(hspace=0.0)
	f.subplots_adjust(left=0.11)
	f.subplots_adjust(bottom=0.09)
	f.text(0.0175, 0.5, ny, va='center', rotation='vertical', fontsize = 22)
	f.text(0.45, 0.025, nx, va='center', rotation='horizontal', fontsize = 22)
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


