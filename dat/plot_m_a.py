#This is the plotting module going with newdiff.f90
#Called in the directory of dat/ as:
#
# python plot_m_a.py Hist*.dat or
# python plot_m_a.py Tran*.dat
#
#Note: 'xxxx*.dat' uses wildcard '*'; so make sure there are no
#irrelevant data files in the data directory beginning with
#'k','n','a','l' (check the fucntion 'do_plot' code for the
# reason)

import matplotlib.pyplot as plt
import numpy as np
import sys
plt.hold(True)

hwidth = 1.7
acx = 0.5
acy = 1.05
uc = 'upper center'
lr = 'lower right'
ul = 'upper left'

#smeta, sma, smdelta, smddelta = np.loadtxt('malist_smcos.dat', unpack = True, 
#        usecols = [0, 1, 2, 3], skiprows = 1)
ksm, smeta, lcdma, lcdmdelta, lcdmddelta = np.loadtxt('malist_lcdm.dat', unpack = True, 
        usecols = [0, 1, 2, 3, 4], skiprows = 1)
knc, nceta, nca, ncdelta, ncddelta = np.loadtxt('malist_nocouple.dat', unpack = True, 
        usecols = [0, 1, 2, 3, 4], skiprows = 1)

def do_histplot():
    plt.xlim(xmax = 1.0)
    plt.ylim(ymax = 1.2)
#    plt.plot(sma, smdelta, linestyle = '--', label = r'$\Lambda$'+'CDM')
    plt.plot(nca, ncdelta, linestyle = '-.', label = 'Non-Coupling',lw=hwidth)
    plt.plot(lcdma, lcdmdelta, linestyle = '--', label = r'$\Lambda$'+'CDM',lw=hwidth)
    plt.legend(loc=ul, ncol = 1, fancybox=True, shadow=True, frameon = False)
    plt.xlabel('a', fontsize='x-large')
    plt.ylabel(r'$\delta_m$', fontsize='x-large')

def do_tranplot():
    plt.xlim((0.002, 0.5))
#    plt.ylim(ymin = 0.99)
#    plt.plot(sma, smdelta, linestyle = '--', label = r'$\Lambda$'+'CDM')
    plt.legend(loc=ul, ncol = 1, fancybox=True, shadow=True, frameon = False)
    plt.xlabel(r'$k [h\ Mpc^{-1}]$', fontsize='x-large')
    plt.ylabel(r'$T(k)$ (normalized)', fontsize='x-large')

haven = False
havel = False
havea = False

isTran = False
isHist = False

for datafile in sys.argv[1:]:
    if(datafile[0] == 'H'):
	isHist = True
    	k,eta,a,delta,ddelta = np.loadtxt(datafile, unpack=True,
        	usecols=[0,1,2,3,4], skiprows=1)
    if(datafile[0] == 'T'):
	isTran = True
	k, a, delta =  np.loadtxt(datafile, unpack = True,
		usecols=[0,1,2], skiprows=1)

    if(datafile[5] == 'n'):
        haven = True
        paralabel = 'n='+datafile[7:12] #str(n[1])
        plt.figure(2)
    if(datafile[5] == 'l'):
        havel = True
        paralabel = r'$\lambda$='+datafile[7:13]+'E-6'
        plt.figure(3)
    if(datafile[5] == 'a'):
        havea = True
        paralabel = r'$\alpha$='+datafile[7:13]+r'$H_0^2$'
        plt.figure(4)

    if (isHist):
	plt.plot(a, delta, label=paralabel, lw=hwidth)
    if (isTran):
	plt.semilogx(k, delta/delta[0], label=paralabel, lw=hwidth)

if(haven):
    plt.figure(2)
    if(isHist):
	do_histplot()
    	plt.savefig('hist-n.eps', bbox_inches='tight')
    if(isTran):
	do_tranplot()
    	plt.savefig('tran-n.eps', bbox_inches='tight')

if(havel):
    plt.figure(3)
    if(isHist):
	do_histplot()
    	plt.savefig('hist-l.eps', bbox_inches='tight')
    if(isTran):
	do_tranplot()
    	plt.savefig('tran-l.eps', bbox_inches='tight')

if(havea):
    plt.figure(4)
    if(isHist):
	do_histplot()
    	plt.savefig('hist-a.eps', bbox_inches='tight')
    if(isTran):
	do_tranplot()
    	plt.savefig('tran-a.eps', bbox_inches='tight')
#plt.show()
