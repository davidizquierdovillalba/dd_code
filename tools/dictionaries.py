import numpy as np


def zdict(x, mill2):
    if mill2 == True:
        vals = {
            0.0 : '0.00',
            0.5 : '0.51',
            1.0 : '1.04',
            1.5 : '1.48',
            2.0 : '2.07',
            2.5 : '2.44',
            3.0 : '3.11',
            3.5 : '3.37',
            4.0 : '3.95',
            4.5 : '4.64',
            5.0 : '5.03',
            5.5 : '5.46',
            6.0 : '5.92',
            6.5 : '6.42',
            7.0 : '6.97',
            7.5 : '7.57',
            8.0 : '8.22',
            9.0 : '8.93',
            10.0 : '9.72'
        }.get(x, ' ')
    else:
        vals = {
            0.0 : '0.00',
            0.25: '0.26',
            0.5 : '0.51',
            1.0 : '1.04',
            1.25: '1.25',
            1.5 : '1.48',
            2.0 : '2.07',
            2.25: '2.25',
            2.5 : '2.44',
            3.0 : '3.11',
            3.25: '3.37',
            3.5 : '3.65',
            4.0 : '3.95',
            4.25: '4.28',
            4.5 : '4.64',
            5.0 : '5.03',
            5.25: '5.46',
            5.5 : '5.45',
            6.0 : '5.92',
            6.5 : '6.42',
            7.0 : '6.97',
            7.5 : '7.57',
            8.0 : '8.22',
            8.5 : '8.22',
            9.0 : '8.93',
            9.5 : '9.72',
        }.get(x, ' ')
    if vals == ' ' :
        raise ValueError("Unknown redshift value: %s. You may want to add it in the dictionary"%x)
    return vals

    

def snapdict(x):
    vals = {
        0 : '63',
        1 : '40',
        2 : '27',
        3 : '23'
    }.get(x, ' ')
    if vals == ' ' :
        raise ValueError("Invalid redshift value: %s. You may want to add it in the dictionary"%x)
    return vals




def Prop(x):
    val = {
        0:'BlackHoleMass', 
        1:'StellarMass',
        2:'Mvir',
        3:'Sfr',
        4:'BulgeMass', 
        5:'DiskMass'
    }.get(x, ' ')
    if val == ' ':
        raise ValueError('Property %d not found'%x)
    return val




def Labels(x, log, noH):
    if log == True:
        if noH == False:
            lab = {
                0:r'$\mathrm{log_{10}\left(M_{BH}[M_{\odot}/h]\right)}$',
                1:r'$\mathrm{log_{10}\left(M_{stellar}[M_{\odot}/h]\right)}$',
                2:r'$\mathrm{log_{10}\left(M_{vir}[M_{\odot}/h]\right)}$',
                3:r'$\mathrm{log_{10}\left(Sfr[M_{\odot}/yr]\right)}$',
                4:r'$\mathrm{log_{10}\left(M_{Bulge}[M_{\odot}/h]\right)}$',
                5:r'$\mathrm{log_{10}\left(M_{Disk}[M_{\odot}/h]\right)}$'
            }.get(x, ' ')
        else:
            lab = {
                0:r'$\mathrm{log_{10}\left(M_{BH}[M_{\odot}]\right)}$',
                1:r'$\mathrm{log_{10}\left(M_{stellar}[M_{\odot}]\right)}$',
                2:r'$\mathrm{log_{10}\left(M_{vir}[M_{\odot}]\right)}$',
                3:r'$\mathrm{log_{10}\left(Sfr[M_{\odot}/yr]\right)}$',
                4:r'$\mathrm{log_{10}\left(M_{Bulge}[M_{\odot}]\right)}$',
                5:r'$\mathrm{log_{10}\left(M_{Disk}[M_{\odot}]\right)}$'
            }.get(x, ' ')
        
    else:
        if noH == False:
            lab = {
                0:r'$\mathrm{M_{BH}[M_{\odot}/h]}$',
                1:r'$\mathrm{M_{stellar}[M_{\odot}/h]}$',
                2:r'$\mathrm{M_{vir}[M_{\odot}/h]}$',
                3:r'$\mathrm{Sfr[M_{\odot}/yr]}$',
                4:r'$\mathrm{M_{Bulge}[M_{\odot}/h]}$',
                5:r'$\mathrm{M_{Disk}[M_{\odot}/h]}$'
            }.get(x, ' ')
        else:
            lab = {
                0:r'$\mathrm{M_{BH}[M_{\odot}]}$',
                1:r'$\mathrm{M_{stellar}[M_{\odot}]}$',
                2:r'$\mathrm{M_{vir}[M_{\odot}]}$',
                3:r'$\mathrm{Sfr[M_{\odot}/yr]}$',
                4:r'$\mathrm{M_{Bulge}[M_{\odot}]}$',
                5:r'$\mathrm{M_{Disk}[M_{\odot}]}$'
            }.get(x, ' ')

    if lab == ' ':
        raise ValueError('Property %d not found'%x)
    return lab


def MoL_labels(x, log, noH):
    if log == True:
        if noH == False:
            lab = {
                0:(r'$\mathrm{log_{10}\left(M_{BH}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3\,dex^{-1}]}$'),
                1:(r'$\mathrm{log_{10}\left(M_{stellar}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3\,dex^{-1}]}$'),
                2:(r'$\mathrm{log_{10}\left(M_{vir}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3\,dex^{-1}]}$'),
                3:(r'$\mathrm{log_{10}\left(Sfr[M_{\odot}/yr]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3\,dex^{-1}]}$'),
                4:(r'$\mathrm{log_{10}\left(M_{Bulge}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3\,dex^{-1}]}$'),
                5:(r'$\mathrm{log_{10}\left(M_{Disk}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3\,dex^{-1}]}$')
            }.get(x, ' ')
        else:
            lab = {
                0:(r'$\mathrm{log_{10}\left(M_{BH}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,dex^{-1}]}$'),
                1:(r'$\mathrm{log_{10}\left(M_{stellar}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,dex^{-1}]}$'),
                2:(r'$\mathrm{log_{10}\left(M_{vir}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,dex^{-1}]}$'),
                3:(r'$\mathrm{log_{10}\left(Sfr[M_{\odot}/yr]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,dex^{-1}]}$'),
                4:(r'$\mathrm{log_{10}\left(M_{Bulge}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,dex^{-1}]}$'),
                5:(r'$\mathrm{log_{10}\left(M_{Disk}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,dex^{-1}]}$')
            }.get(x, ' ')
    else:
        if noH == False:
            lab = {
                0:(r'$\mathrm{\left(M_{BH}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3]}$'),
                1:(r'$\mathrm{\left(M_{stellar}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3]}$'),
                2:(r'$\mathrm{\left(M_{vir}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3]}$'),
                3:(r'$\mathrm{\left(Sfr[M_{\odot}/yr]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3]}$'),
                4:(r'$\mathrm{\left(M_{Bulge}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3]}$'),
                5:(r'$\mathrm{\left(M_{Disk}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}\,h^3]}$')
            }.get(x, ' ')
        else:
            lab = {
                0:(r'$\mathrm{\left(M_{BH}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}]}$'),
                1:(r'$\mathrm{\left(M_{stellar}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}]}$'),
                2:(r'$\mathrm{\left(M_{vir}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}]}$'),
                3:(r'$\mathrm{\left(Sfr[M_{\odot}/yr]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}]}$'),
                4:(r'$\mathrm{\left(M_{Bulge}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}]}$'),
                5:(r'$\mathrm{\left(M_{Disk}[M_{\odot}/h]\right)}$',r'$\rm{\Phi\ [Mpc^{-3}]}$')
            }.get(x, ' ')
	
    if lab == ' ':
        raise ValueError('Property %d not found'%x)

    return lab 
