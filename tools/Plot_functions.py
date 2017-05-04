import matplotlib.colors as colors
import pylab as plt
from Config import *
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)

from dictionaries import *
from Useful_func import *

def PvsP_densityplot(LGout_dir,fileprefx,firstfile,lastfile,P1,P2,ns, Plotlog = True, Factor10 = [True,True], removeh = False, GALTREE=False):

        print '###########################################################################'
        print '    The redshifts that are going to be used are:', zlist
        print '    Properties:', P1, 'vs', P2
        print '############################################################################'

        nx = Labels(wtp[0],Plotlog,removeh)
	ny = Labels(wtp[1],Plotlog,removeh)
	row = [0,0,1,1]
        col = [0,1,0,1]
        f, axarr = plt.subplots(2, 2, figsize = (12,8))

        if(GALTREE == True):
                a = read_tree(LGout_dir,fileprefx,firstfile,lastfile,PropertiesToRead_tree,LGalaxiesStruct)
                gg = a[1]
                del a
        maximX = np.empty(len(zlist))
        maximY = np.empty(len(zlist))
        minimX = np.empty(len(zlist))
        minimY = np.empty(len(zlist))

        for i in zlist:
                if(GALTREE == False):

			filepref = fileprefx + np.str(zdict(i))
                        a = read_snap(LGout_dir,filepref,firstfile,lastfile,PropertiesToRead,LGalaxiesStruct)
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
                                p1 = p1 - np.log10(cosmo.h)

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
		d = axarr[row[i],col[i]].pcolor(x, y, H, norm = colors.LogNorm(),cmap='jet', label = '$z$ = ' + zdict(i))#cmap='RdGy')
		d.set_clim(vmin=0.01, vmax=2)
		if(row[i] == 1 and col[i]==1):
			position=f.add_axes([0.90,0.12,0.02,0.75])
			f.colorbar(d,cax=position, label = 'N/$\mathrm{N_{max}}$')
			#cbar = f.colorbar(d, cmap='jet', ax = axarr[row[i],col[i]], label = 'N/$\mathrm{N_{max}}$')

                #axarr[row[i],col[i]].text(min(p1),max(p2)*0.90,'$z$ = ' + zdict(i),fontsize = 20, weight =  'bold')
		leg = axarr[row[i],col[i]].legend(loc = 'lower right',fontsize = 18, handlelength=0, handletextpad=0, fancybox=True)
		for item in leg.legendHandles:
    			item.set_visible(False)

                maximX[i] = max(p1)
                minimX[i] = min(p1)
                maximY[i] = max(p2)
                minimY[i] = min(p2)

                if(row[i] == 0):
                        axarr[row[i],col[i]].set_xticklabels([''])
                if(col[i]>0):
                        axarr[row[i],col[i]].set_yticklabels([''])

	for i in zlist:
		axarr[row[i],col[i]].set_ylim(min(minimY),max(maximY))
                axarr[row[i],col[i]].set_xlim(min(minimX),max(maximX))
        f.canvas.draw()
        labels = [item.get_text() for item in axarr[1,0].get_xticklabels()]
        labels[-1] = ' '
        axarr[1,0].set_xticklabels(labels)
        #fig = plt.tight_layout()
        f.subplots_adjust(wspace=0.0)
        f.subplots_adjust(hspace=0.0)
        f.subplots_adjust(left=0.11)
        f.subplots_adjust(right=0.88)
        f.subplots_adjust(bottom=0.09)
        f.text(0.0325, 0.5, ny, va='center', rotation='vertical', fontsize = 22)
        f.text(0.45, 0.025, nx, va='center', rotation='horizontal', fontsize = 22)
        del p1,p2,gg
        plt.savefig(ns,aspect='equal')
        plt.show()






def MoL_func(LGout_dir,fileprefx,firstfile,lastfile,MaxTreeFiles,BoxSize,P1,ns,Plotlog = True, Factor10 = [True,True], removeh = False, GALTREE=False):

	nx, ny = MoL_labels(wtp[0],Plotlog,removeh)



        if(GALTREE == True):
                a = read_tree(LGout_dir,fileprefx,firstfile,lastfile,PropertiesToRead_tree,LGalaxiesStruct)
                gg = a[1]

        if(removeh == True):
                Volume = ((BoxSize/cosmo.h)**3.0) * (lastfile - firstfile + 1) / MaxTreeFiles # Mpc^3
        else:
                Volume = (BoxSize**3.0) * (lastfile - firstfile + 1) / MaxTreeFiles # Mpc^3 h^-3

        row = [0,0,1,1]
        col = [0,1,0,1]
        f, axarr = plt.subplots(2, 2, figsize = (12,8))

	for i in zlist:
                if(GALTREE == False):
			filepref = fileprefx + np.str(zdict(i))
                        a = read_snap(LGout_dir,filepref,firstfile,lastfile,PropertiesToRead,LGalaxiesStruct)
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

                if(Plotlog == True):
                       sb = sizebin
                else:
                        sb = (max(pp) - min(pp))*sizebin
                pp_c, phi = hist(pp,sb)
                phi = phi/(Volume * sb)
                axarr[row[i],col[i]].plot(pp_c,phi,color='red', linewidth=2, label = '$z$ = ' + zdict(i))
                axarr[row[i],col[i]].set_yscale('log')
                axarr[row[i],col[i]].set_ylim(MoL_func_ylim)
		leg = axarr[row[i],col[i]].legend(loc = 'upper right',fontsize = 18, handlelength=0, handletextpad=0, fancybox=True)
                for item in leg.legendHandles:
                        item.set_visible(False)		

		axarr[row[i],col[i]].set_xlim(MoL_func_xlim)

                if(row[i] == 0):
                        axarr[row[i],col[i]].set_xticklabels([''])
                if(col[i]>0):
                        axarr[row[i],col[i]].set_yticklabels([''])

	f.canvas.draw()
        labels = [item.get_text() for item in axarr[0,0].get_yticklabels()]
        labels[1] = ' '
        axarr[0,0].set_yticklabels(labels)

        fig = plt.tight_layout()
        f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace=0.0)
        f.subplots_adjust(left=0.11)
        f.subplots_adjust(bottom=0.09)
        f.text(0.0325, 0.5, ny, va='center', rotation='vertical', fontsize = 22)
        f.text(0.45, 0.025, nx, va='center', rotation='horizontal', fontsize = 22)
        plt.savefig(ns)
        plt.show()



