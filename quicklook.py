#!/usr/bin/env python

from archivehandler import ArchiveHandler
import numpy as np
from matplotlib.pyplot import *
import matplotlib.gridspec as gridspec
import time
import warnings 
import argparse



def fitfuncDM(p,f):
    return 4.149e9*p[0]/f**2 + p[1] #DMshift in mus, f in MHz, p[1] = t_inf
    return 4.149e3*p[0]/f**2 + p[1]#DMshift in s, f in MHz
def errfuncDM(p,f,y,err):
    return (y-fitfuncDM(p,f))/err
pinitDM = [0.0,0.0]


'''
This class takes a ArchiveHandler and does the very basic plotting
This is the main part of the code that calls everything else, where all of the graphics can be written
'''
class Quicklook:
    def __init__(self,FILE, templatefilename=None, NCHAN=16, save=False, disp = False, minutes = False, ext = None):
        warnings.filterwarnings("ignore") #probably should remove later
        start = time.time()
        self.templatefilename = templatefilename
        self.save = save
        self.disp = disp
        self.minutes = minutes
        self.ext = ext
        self.handler = ArchiveHandler(FILE,templatefilename = templatefilename, NCHAN=NCHAN)
        
        self.fig = figure(figsize= (15,9.5))
        self.fig.subplots_adjust(top=0.93)
        gs1 = gridspec.GridSpec(13,2)
        gs1.update(left = .05, right = .33, wspace = .05, hspace = .2)
        
        
        ax1 = subplot(gs1[2:,0])
        self.handler.plotProfiles(ax = ax1, labels = True)
        
        #ax2 = subplot(gs1[2:,1], sharey = ax1)
        ax2 = subplot(gs1[2:,1])
        self.handler.plotProfiles(ax = ax2, difference = True, labels = False)
        ax2.set_ylim(ax1.get_ylim())
        
        ax3 = subplot(gs1[0:2,0])
        self.handler.plotAverageProfile(ax = ax3)
        
        ax4 = subplot(gs1[0:2,1])
        self.handler.plotTemplate(ax = ax4)
        
        ax3.axes.get_xaxis().set_visible(False)
        ax4.axes.get_xaxis().set_visible(False)
        
        gs2 = gridspec.GridSpec(13,2)
        #gs2.update(left = .37, right = .61)
        gs2.update(left = .34, right = .58)
        
        #ax5 = subplot(gs2[2:,:], sharey = ax1)
        ax5 = subplot(gs2[2:,:])
        self.handler.plotDynamicSpectrum(ax = ax5, minutes = minutes, labels = True)
        ax5.set_ylim(ax1.get_ylim())
        
        
        gs3 = gridspec.GridSpec(26,12)
        #gs3.update(left = .65, right = .99, wspace = .02, hspace = .02)
        gs3.update(left = .62, right = .99, wspace = .02, hspace = .02)
        
        ax7 = subplot(gs3[2:9,1:11])
        self.handler.plotAcf2d(ax = ax7, minutes = minutes, labels = True)
        
        ax8 = subplot(gs3[11:20,1:11])
        self.handler.plotSecondarySpectrum(ax = ax8, minutes = minutes, labels = True)
        
        ax9 = subplot(gs3[22:26,0:5])
        if not args.nodm:
            print 'Now calculating SNR vs DM. This is slow; to skip this step use quicklook.py -nodm'
            self.handler.DMCalculator(ax9, initial=.99, end = 1.01, iters = np.abs(args.iters), max_depth = np.abs(args.depth))
        
        ax10 = subplot(gs3[22:26,7:12])
        self.handler.plotBasicHistogram(ax = ax10, bins = 48, labels = True)
    
        self.plotHeader()
        
        if self.handler.templatefilename is None:
            textstr = 'NONE'
        else:
            textstr = 'TMPL'
        
        props = dict(boxstyle='round', facecolor= 'none', edgecolor = 'black', alpha=0.5)
        self.fig.text(0.3, 0.92, textstr, fontsize=8, verticalalignment='top', bbox=props)
        
        #hide stuff
        ax3.axes.get_xaxis().set_visible(False)
        ax4.axes.get_xaxis().set_visible(False)
        ax4.axes.get_yaxis().set_visible(False)
        ax2.axes.get_xaxis().set_visible(False)
        ax2.axes.get_yaxis().set_visible(False)
        ax5.axes.get_yaxis().set_visible(False)
        
        
        end = time.time()
        print 'Run time: %0.2f s'%(end - start) 
            
        if ext is not None:
            self.save_to_pdf(ext=ext)
        
        if save:
            self.handler.save_to_npz(name = None)

        if disp:
            show()

            
    def plotHeader(self):
        """Plot the header for Quicklook"""
        calcDM = self.handler.calculated_DM
        if calcDM is None:
            calcDM = "Not calculated"
        else:
            calcDM = "%0.6f"%calcDM
        duration = self.handler.getDuration()
        if self.minutes:
            duration = "%0.2f min"%(duration/60)
        else:
            duration = "%0.2f s"%(duration)
        header = ('Name:  %s    Telescope:  %s    MJD:  %i\n'%(self.handler.getName(full=True),self.handler.getTelescope(),self.handler.getMJD())+
                  'PSR:  %s    Period:  %0.3f ms    Nbins:  %i    DM:  %0.6f    Calc DM:  %s\n'%(self.handler.ar.getName(),self.handler.getPeriod()*1e3,self.handler.getNbin(),self.handler.getDM(),calcDM)
                  'BW:  %0.2f MHz    Nchan:  %i    Duration  %s\n'%(abs(self.handler.getBandwidth()),np.size(self.handler.F),duration)+
                  'Time/Bin:  %0.2f us'%(self.handler.getPeriod()/self.handler.getNbin()*1e6))
                
        props = dict(boxstyle='round', facecolor= 'none', edgecolor = 'black', alpha=0.5)
        self.fig.text(0.35, 0.95, header, fontsize=13, verticalalignment='top')
            
    
    def save_to_pdf(self, name=None, ext = None):
        """Print the Quicklook plot to a file with extension ext"""
        if name is None:
            name = self.handler.name
        name = name + '.' + ext
        self.fig.savefig(name)
        
        
            
'''
Test in the terminal with: python quicklook.py
OR
in ipython: %run quicklook.py

'''

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Run quicklook with commands')
    parser.add_argument('-nchan', help = 'Number of frequency channels. Deafult is 16',default = 16, type = int)
    parser.add_argument('-save', '--save',help = 'Save to .npz file?', action = 'store_true')
    parser.add_argument('inputfile',help = 'File (string) to be opened using Quicklook', type = str)
    parser.add_argument('-noshow', '--noshow',help = 'Display the defualt Quicklook', action= 'store_true')
    parser.add_argument('-min', '--min',help = 'Plot with minutes as defualt', action= 'store_true', default = False)
    parser.add_argument('-ext','--ext', help = 'Print Quicklook to file with extension',default = None,choices = ['emf','png', 'pdf', 'ps','eps','svg','raw','rgba','svgz'])
    parser.add_argument('-iters','--iters', help = 'Number of DM iteration to go through when calculating optimal DM',default = 11,type = int)
    parser.add_argument('-depth','--depth', help = 'The max depth to reach when calculating optimal DM',default = 3,type = int)
    parser.add_argument('-nodm','--nodm', help = 'Do not calculate the optimal DM', action = 'store_true',default = False)
    parser.add_argument('-template', '--template',help = 'Give location of template file for pulsar',default = None, type = str)
    args = parser.parse_args()



    if args.inputfile is None:
        raise UserWarning('You have not supplied a file!')
    else:
        filename = args.inputfile


    
    disp = (not args.noshow)



    ql = Quicklook(FILE=filename,templatefilename=args.template, NCHAN=args.nchan , save=args.save, disp=disp, minutes = args.min, ext=args.ext)
