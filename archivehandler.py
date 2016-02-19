import os
import numpy as np
from matplotlib.pyplot import *
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import scipy.optimize as optimize
import sys
import glob 
#import barycenter as bary



import pypulse.utils as u
from pypulse.archive import *
from pypulse.singlepulse import *
from pypulse.dynamicspectrum import *





def fitfuncDM(p,f):
    return 4.149e9*p[0]/f**2 + p[1] #DMshift in mus, f in MHz, p[1] = t_inf
    return 4.149e3*p[0]/f**2 + p[1]#DMshift in s, f in MHz
def errfuncDM(p,f,y,err):
    return (y-fitfuncDM(p,f))/err
pinitDM = [0.0,0.0]


NCHAN = 16

'''
This class builds off of the Archive class, see archive.py

'''

class ArchiveHandler:
    def __init__(self,filename,NCHAN = NCHAN,templatefilename=None,windowsize=256,prepare = True,**kwargs):
        name , file_extension = os.path.splitext(filename)
        path , self.name = os.path.split(name)
        self.fullname = "%s%s" % (self.name,file_extension)
        if file_extension == '.npz':
            self.npz = np.load(filename)
        #if file_extension == '.fits' or file_extension == '.zap':
        else:
            self.npz = None
        #else:
        #    raise RuntimeError('File has improper extension for Quicklook, and this code will not work')
        
        name, file_extension = os.path.splitext(filename)
        self.filename = filename
        self.templatefilename = templatefilename
        self.windowsize = windowsize
        
        self.ds = None #is all this necessary?
        self.ss = None
        self.acf2d = None
        self.profiles = None
        self.difference_profiles = None
        self.SS_xaxis = None
        self.SS_yaxis = None
        self.ACF2D_xaxis = None
        self.ACF2D_yaxis = None
        self.peak_DM = None
        self.DM_arr = None
        self.calculated_DM = None
        
        
        if self.templatefilename is not None:
                artemp = Archive(self.templatefilename, prepare=False)
                artemp.pscrunch()
                temp = u.normalize(u.center_max(artemp.getData()),simple=True)
                self.sptemp = SinglePulse(temp,windowsize=windowsize)
        
        if self.npz is None:
            self.ar = Archive(self.filename, prepare=prepare)
            self.F = self.ar.getAxis('F')
            # temporary holdover from archive.py frequency issues
            #if self.F[0] > self.F[-1]:
            #    self.F = self.F[::-1]
            self.T = self.ar.getAxis('T')
            self.data = self.ar.getData()
    
            #self.calculateAverageProfile()
            
            #calculateAverageProfile
            self.ar.tscrunch()
            self.ar.fscrunch()
            avgprof = self.ar.getData()
            imax = np.argmax(avgprof)
            self.average_profile = u.center_max(avgprof) #not normalized!
            
            if self.templatefilename is None:
                self.spavgprof = SinglePulse(self.average_profile,windowsize=windowsize)
                self.sptemp = SinglePulse(u.normalize(self.average_profile),mpw=self.spavgprof.mpw,opw=self.spavgprof.opw)
            else:
                self.spavgprof = SinglePulse(self.average_profile,mpw=self.sptemp.mpw,opw=self.sptemp.opw)
        
            self.spavgprof.remove_baseline()
            
            #Reset archive
            self.ar.reset()
    
            self.nbins = len(self.average_profile)
            alignval = self.nbins/2 - imax
            self.data = np.roll(self.data,alignval)
            self.DM0 = self.ar.getDM()
        
        else:
            self.F = self.npz['frequency']
            self.T = self.npz['time']
            self.average_profile = self.npz['average_profile'][0,:] #is it bad that we need to do this?
            if self.templatefilename is None:
                self.spavgprof = SinglePulse(self.average_profile,windowsize=windowsize)
                self.sptemp = SinglePulse(u.normalize(self.average_profile),mpw=self.spavgprof.mpw,opw=self.spavgprof.opw)
            else:
                self.spavgprof = SinglePulse(self.average_profile,mpw=self.sptemp.mpw,opw=self.sptemp.opw)
            self.spavgprof.remove_baseline()
            self.nbins = len(self.average_profile)
            self.DM0 = self.npz['DM']
            self.peak_DM = self.npz['peak_DM']
            self.DM_arr = self.npz['DM_arr']
        
        if len(self.F) % NCHAN != 0:
            str(len(self.F)) + '%' + str(NCHAN) + '=' + str( len(self.F) % NCHAN )
            raise UserWarning('Number of frequency channels in file is not multiple of your provided NCHAN')
    
    
    def getFlux(self,mode='max'): 
        """
        Return various statistics of the flux of the average profile
        """
        if self.templatefilename is None:
            if mode == 'max':
                return np.max(self.average_profile)
            if mode == 'mean':
                return np.mean(self.average_profile)
        else:
            if mode == 'max':
                return self.spavgprof.fitPulse(self.sptemp.data)[2]
            if mode == 'mean':
                return self.spavgprof.fitPulse(self.sptemp.data)[2]/self.nbins


    def getSN(self): 
        """
        Return the S/N of the average profile
        """
        if self.templatefilename is None:
            return np.max(self.average_profile)/self.spavgprof.getOffpulseNoise()
        else:
            return self.spavgprof.fitPulse(self.sptemp.data)[-2]


    def getDeltaDM(self,nchan=NCHAN,barycenter=True): 
        """
        Calculate the DM offset 
        """
        data = np.mean(self.data,axis=0)
        newdata = np.zeros((nchan,self.nbins))
        window = len(data)/nchan
        for i in range(nchan):
            newdata[i,:] = np.mean(data[i*window:(i+1)*window],axis=0)
        F = u.decimate(self.F,window)



        freqs = []
        tauhats = []
        sigma_taus = []
        for i,prof in enumerate(newdata):
            spprof = SinglePulse(prof,mpw=self.spavgprof.mpw,opw=self.spavgprof.opw)                
            #plot(u.normalize(prof))
            #plot(self.sptemp.data)
            #show()
            retval = spprof.fitPulse(self.sptemp.data) #WHAT IF NO TEMPLATE?

            if retval is None:
                continue
            tauhat,sigma_tau = retval[1],retval[3]
            freqs.append(F[i])
            tauhats.append(tauhat)
            sigma_taus.append(sigma_tau)

        freqs = np.array(freqs)
        tauhats = np.array(tauhats)
        sigma_taus = np.array(sigma_taus)

        out = optimize.leastsq(errfuncDM,pinitDM,args=(freqs,tauhats,sigma_taus),full_output=1)
        residuals = tauhats-fitfuncDM(out[0],freqs)
        s_sq = np.sum(residuals**2)/(len(tauhats)-len(pinitDM))

        #print "DM",out[0][0],np.sqrt(out[1][0,0]*s_sq)
        #errorbar(freqs,tauhats,yerr=sigma_taus,fmt='k.')
        #plot(freqs,fitfuncDM(out[0],freqs))
        #show()
        #errorbar(freqs,residuals,yerr=sigma_taus,fmt='k.')
        #show()

        DM, DMerr = out[0][0],np.sqrt(out[1][0,0]*s_sq) #check the inversion of DM?

        if barycenter:
            #print "DM",DM
            return DM, DMerr
            #DM = bary.convertDMtopo(DM,self.ar.getTelescopeCoords(),self.ar.getPulsarCoords(parse=False),self.ar.getMJD(full=True))
            #print "DMnew",DM
            #err?


        return DM, DMerr


    def getCorrectedDM(self): #UY
        """
        Return the header DM and calculated DM offset
        """
        DeltaDM,DeltaDMerr = self.getDeltaDM()
        return self.DM0 + DeltaDM,DeltaDMerr


    def getDynamicSpectrum(self,fast=True, reset = False):
        """
        Calculate the dynamic spectrum via the DynamicSpectrum class
        """
        if self.ds is None:
            reset = True
        if reset:
            if self.npz is None:
                if self.templatefilename is None or fast:
                    if np.ndim(self.data) <= 2: #quick patch
                        return
                    ds = np.transpose(np.mean(self.data,axis=2)) #no need to reset archive
                    #this is the average over all 2048 phase bins, not the peak!
                else:
                    self.ar.tscrunch()
                    self.ar.fscrunch()
                    alignval = self.nbins/2 - np.argmax(self.ar.getData())
                    self.ar.reset(prepare=True)
                                                                 
                    ds,dsoff,dserr = self.ar.getDynamicSpectrum(template=self.sptemp.data,mpw=self.sptemp.mpw,align=alignval,verbose=True)
            else:
                ds = self.npz['DynamicSpectrum']
            self.ds = ds
        else:
            ds = self.ds
            
        return ds


    def getProfiles(self,nchan = NCHAN,difference=False, reset = False): #may run several times due to there deing difference option
        """
        Return profiles of a given frequency channelization
        """
        if self.profiles is None:
            reset = True
        if self.difference_profiles is None:
            reset = True
        if reset:
            if self.npz is None:
                if np.ndim(self.data) <= 2: #patch
                    data = self.data
                else:
                    data = np.mean(self.data,axis=0)
                newdata = np.zeros((nchan,self.nbins))
                window = len(data)/nchan
                for i in range(nchan):
                    newdata[i,:] = np.mean(data[i*window:(i+1)*window],axis=0)
                F = u.decimate(self.F,window)
                self.profiles = newdata
        
                if difference:
                    if self.templatefilename is None:
                        sptemp = self.spavgprof                
                    else:
                        sptemp = self.sptemp
                    for i,prof in enumerate(newdata):
                        spprof = SinglePulse(prof,mpw=sptemp.mpw,opw=sptemp.opw)
                        retval = spprof.fitPulse(sptemp.data)
                        if retval is None:
                            continue
                        tauhat,bhat = retval[1],retval[2]
        
                        newdata[i] = prof - bhat*sptemp.shiftit(tauhat)
                    self.diff_profiles = newdata 
                return newdata
            
            else:
                if difference:
                    return self.npz['differenceprofiles']
                else:
                    return self.npz['profiles']
        else:
            if difference:
                return self.difference_profiles
            else:
                return self.profiles 


    def getSecondarySpectrum(self,fast=True, reset = False):
        """
        Calculate the secondary spectrum
        """
        if self.ss is None:
            reset = True
        if reset:
            if self.npz is None:
                if self.ds is None:
                    self.ds = self.getDynamicSpectrum(fast=fast)
                ss = np.log10(np.abs(np.fft.fftshift(np.fft.fft2(self.ds)))**2)
            else:
                ss = self.npz['SecondarySpectrum']
            
            self.ss = ss
        else:
            ss = self.ss
        return ss


    def plotTemplate(self,ax=None):
        """
        Plot the template shape, useful if on a specific axis
        """
        if self.templatefilename is None:
            return
        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
        #ax.plot(self.sptemp.data,'k')
        ax.plot(self.sptemp.mpw,self.sptemp.data[self.sptemp.mpw],'g')
        ax.plot(self.sptemp.opw,self.sptemp.data[self.sptemp.opw],'r')
        ax.set_xlim(0,np.max(self.sptemp.bins))
        ax.set_ylim(-0.05,1.05)
        ax.set_xlabel("Phase bins")
        ax.set_ylabel("Normalized Intensity")
        if doshow:
            show()
        return ax


    def plotAverageProfile(self,ax=None,labels = False,spavgprof_mpw=None,spavgprof_mpw_data=None,spavgprof_opw=None,spavgprof_opw_data=None,spavgprof=None):
        """
        Plot the average profile, useful if on a specific axis
        """
        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
        #ax.plot(self.average_profile,'k')
        if spavgprof_mpw is None:
            spavgprof_mpw = self.spavgprof.mpw
            spavgprof_mpw_data = self.spavgprof.data[self.spavgprof.mpw]
        if spavgprof_opw is None:
            spavgprof_opw = self.spavgprof.opw
            spavgprof_opw_data = self.spavgprof.data[self.spavgprof.opw]
        if spavgprof is None:
            spavgprof = self.average_profile 
        ax.plot(spavgprof_mpw,spavgprof_mpw_data,'g')
        ax.plot(spavgprof_opw,spavgprof_opw_data,'r')
        dy = np.ptp(spavgprof)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_xlim(0,len(spavgprof))
        ax.set_ylim(np.min(spavgprof)-0.05*dy,np.max(spavgprof)+0.05*dy)
        if labels:
            ax.set_xlabel("Phase bins")
            ax.set_ylabel("mJy")
        if doshow:
            show()
        return ax
    
             
    def plotDynamicSpectrum(self,ax=None,fast=True, minutes = False,labels = False, ds = None, T = None, F = None):
        """
        Plot the dynamic spectrum, useful if on a specific axis
        """

        if np.ndim(self.data) <= 2:
            return
        
        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
        if ds is None:
            ds = self.getDynamicSpectrum(fast=fast)
        
        if T is None:
            T = self.T
        if F is None:
            F = self.F
            
        if minutes:
            T = T / 60.0
        
        D = DynamicSpectrum(ds,F = F,T = T)
        D.remove_baseline()
        ax = D.imshow(ax=ax,cmap=cm.jet,alpha=False,show=False)#,cbar=True)
            
        if labels:
            ax.axes.set_ylabel('MHz')
            if minutes:
                ax.axes.set_xlabel('Minutes')
            else:
                ax.axes.set_xlabel('Seconds')
        if doshow:
            show()
        return ax

        
    def plotSecondarySpectrum(self,ax=None,fast=True,minutes=False,labels=False, ss= None, SS_xaxis = None, SS_yaxis = None):
        """
        Plot the secondary spectrum, useful if on a specific axis
        """

        if np.ndim(self.data) <= 2:
            return

        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
            
        if ss is None:
            ss = self.getSecondarySpectrum(fast=fast)
        
        if None in (SS_xaxis,SS_yaxis): #WHAT IF ONE AXIS GIVEN?
            SS_xaxis,SS_yaxis = self.getSecondarySpectrumAxes()
        
        if minutes: 
            SS_xaxis = SS_xaxis*60
        
        extent = [SS_xaxis[0],SS_xaxis[-1],SS_yaxis[0],SS_yaxis[-1]]
        im=u.imshow(ss, ax=ax, extent = extent)
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        if labels:
            ax_right = ax.twinx()   
            ax_right.set_yticks([])
            ax.axes.set_ylabel('Conjugate Frequency (1/MHz)')
            ax_right.axes.set_ylabel("Secondary Spectrum")
            if minutes:
                ax.axes.set_xlabel('Conjugate Time (1/Minutes)')
            else:
                ax.axes.set_xlabel('Conjugate Time (1/Seconds)')
        
        if doshow:
            show()
        return ax
    
    
    def getSecondarySpectrumAxes(self, reset = False): 
        """
        Get the Secondary Spectrum Axes
        """
        if None in (self.SS_xaxis, self.SS_yaxis):
            reset = True
        if reset:
            if self.npz is None:
                taxis = self.T
                Ttransform = np.fft.fft(taxis)
                size = Ttransform.size
                d = taxis[1]-taxis[0] #get time step better
                SS_xaxis = np.fft.fftfreq(size,d)
                SS_xaxis = np.fft.fftshift(SS_xaxis)
                
                faxis = self.F
                ftransform = np.fft.fft(faxis)
                size = ftransform.size
                d =  faxis[1]-faxis[0] #get frequency step better
                SS_yaxis = np.fft.fftfreq(size,d)
                SS_yaxis = np.fft.fftshift(SS_yaxis)
            else:
                SS_xaxis = self.npz['SS_xaxis']
                SS_yaxis = self.npz['SS_yaxis']
            
            self.SS_xaxis = SS_xaxis
            self.SS_yaxis = SS_yaxis
        else:
            SS_xaxis = self.SS_xaxis
            SS_yaxis = self.SS_yaxis

        # Holdover from negative bandwidth issues in pypulse.archive
        if SS_yaxis[0] > SS_yaxis[-1]:
            SS_yaxis = SS_yaxis[::-1]

        return SS_xaxis, SS_yaxis

         
    def plotProfiles(self,ax=None,nchan=NCHAN,difference=False, labels = False, profiles = None, reset = False, BW = None, F = None):
        """
        Plot the data profiles, useful if on a specific axis
        """
    
        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
        if profiles is None:
            profiles = self.getProfiles(nchan=nchan,difference=difference, reset = reset)
        #for profile in profiles:
        #    plt.clf()
        #    plt.plot(profile)
        #    plt.show()
        #raise SystemExit
        
        if BW is None:
            BW = self.getBandwidth()#np.abs(self.getBandwidth())
        if F is None:
            F = self.F #nchan is number of averaged channels
        
        step_size = BW/np.size(F)
        fsize = np.size(F)
        
        scale_diff = np.amax(profiles) - np.amin(profiles)
        for i,prof in enumerate(profiles):
            #gets middle for frequency range 
            index = (i*fsize/nchan + fsize/nchan/2)*step_size + F[0]
            if difference:
                scale_diff = np.amax(prof) - np.amin(prof)
                ax.plot( (prof/scale_diff) * fsize/nchan*abs(step_size) + index, color = 'k')
            else:
                ax.plot(u.normalize(prof,simple=True)*fsize/nchan*abs(step_size) + index, color = 'k')
        ax.set_xlim(0,len(prof))
        if labels:
            ax.axes.set_ylabel('Frequency (MHz)')
            ax.axes.set_xlabel('Phase Bins')
        if doshow:
            show()
        return ax


    def getAcf2d(self,fast=True, reset = False): 
        """
        Calculate the acf2D 
        """
        if self.acf2d is None:
            reset = True
            
        if reset:
            if self.npz is None:
                if self.ds is None:
                    self.ds = self.getDynamicSpectrum(fast=fast)
                D = DynamicSpectrum(self.ds,F=self.F,T=self.T)
                acf2d = D.acf2d()
            else:
                acf2d = self.npz['ACF2D']
        else:
            acf2d = self.acf2d
            
        self.acf2d = acf2d
        return acf2d
    
 
    def plotAcf2d(self,ax=None,fast=True,minutes=False,labels=False,acf2d=None,tlags=None,flags=None):
        """
        Plot the acf2d, useful if on a specific axis
        """

        if np.ndim(self.data) <= 2:
            return
        
        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
            
        if None in (acf2d, tlags, flags):
            tlags, flags = self.getAcf2dAxes()
            acf2d = self.getAcf2d(fast=fast) 
        if minutes:
            extent = [tlags[0]/60.0,tlags[-1]/60.0,flags[0],flags[-1]]
        else:
            extent = [tlags[0],tlags[-1],flags[0],flags[-1]]
            
        if labels:
            ax_right = ax.twinx()  
            ax_right.set_yticks([])
            ax_right.axes.set_ylabel("ACF2D")
            ax.axes.set_ylabel('Lag (MHz)') #input dynamic units here
            if minutes:
                ax.axes.set_xlabel('Lag (Minutes)')
            else:
                ax.axes.set_xlabel('Lag (Seconds)')
        im=u.imshow(acf2d,ax=ax,extent=extent)

        if doshow:
            show()
        return ax


    def getAcf2dAxes(self, reset=False): 
        """
        Get the secondary spectrum axes
        """
        if None in (self.ACF2D_xaxis, self.ACF2D_yaxis):
            reset = True
        if reset:
            if self.npz is None:
                tlags = u.lagaxis(self.T)
                flags = u.lagaxis(self.F)
            else:
                tlags = self.npz['ACF2D_xaxis']
                flags = self.npz['ACF2D_yaxis']
            
            self.ACF2D_xaxis = tlags
            self.ACF2D_yaxis = flags
        else:
            tlags = self.ACF2D_xaxis
            flags = self.ACF2D_yaxis
            
        # Holdover from negative bandwidth issues in pypulse.archive
        if flags[0] > flags[-1]:
            flags = flags[::-1]
            
        return tlags, flags
    
    
    def plotBasicHistogram(self, ax = None, bins = 20, labels = False):
        """
        Calculates the maximum pulse amplitude  at each frequency, and plots a histogram
        of occureences vs. pulse amplitudes 
        """

        if np.ndim(self.data) <= 2:
            return
        
        doshow = False
        if ax is None:
            doshow = True
            fig = figure()
            ax = fig.add_subplot(111)
        
        arr = self.getDynamicSpectrum()
        abs_max = np.amax(arr)
        abs_min = np.amin(arr)

        step_size = (abs_max - abs_min) / bins
        #remove dead channels
        arr = arr[arr != 0]
        histogram, bin_edges = np.histogram(arr, bins)
        left_edges = bin_edges[0:len(bin_edges)-1]
        ax.bar(left_edges,histogram, width= step_size)
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
        if labels:
            ax.axes.set_ylabel('Occurrences')
            ax.axes.set_xlabel('Amplitude')
            
        if doshow:
            show()
        return ax
 
   
    def FitSN(self,iters = 10,initial = .99,end = 1.01, ax = None,labels = False, plot = True, showfit = False, recalculate = True, DM = None, DM_arr = None):
        """
        Calculate the best fit DM of the data set by finding the peak S/N, and plot it NEEDS WORK
        """
        if DM is None:
            DM = self.getDM()
          
        initial = initial * DM
        end = end * DM
        dms = np.linspace(initial, end, iters)

        sep = dms[1]-dms[0]
    
        if DM_arr is not None:
            DM_arr = DM_arr
            
        elif (self.DM_arr is not None) and (not recalculate):
            DM_arr = self.DM_arr
            
        
        else:
            DM_arr = np.zeros( (iters,4 ) )
            self.ar.dedisperse(DM = (DM-initial), reverse= True)
            for i,dm in enumerate(dms):
                if i >0:
                    self.ar.dedisperse(DM = sep)
                sn = self.ar.getSN()
                DM_arr[i] = [dm,sn,np.max(self.ar.spavg.data),self.ar.spavg.getOffpulseNoise()]
            self.ar.dedisperse(DM = (dm-DM), reverse=True) 
        
        if plot:
            doshow = False
            if ax is None:
                doshow = True
                fig = figure()
                ax = fig.add_subplot(111)
                
            ax.scatter(DM_arr[:,0],DM_arr[:,1])
        
            if showfit:    
                ax.plot(newarr[:,0],newarr[:,1], color = 'g')
                ax.plot(x_new, y_new, color ='r')
        
            min_tick = np.amin(DM_arr[:,0]) - .005*np.amin(DM_arr[:,0])
            max_tick = np.amax(DM_arr[:,0]) + .005*np.amin(DM_arr[:,0])
            ax.xaxis.set_ticks(np.arange(min_tick,max_tick,(max_tick-min_tick)/5.0))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.02f'))
            if labels:
                ax.axes.set_ylabel('Peak SNR')
                ax.axes.set_xlabel('DM')
        
            if doshow:
                show()
         
        
        return DM_arr
   
    def DMCalculator(self, ax=None,plot=True, labels = True, DM = None, iters = 10, initial = .999, end = 1.001, max_depth = 3, arr=None):
        """
        Plot the calculated optimal for each measurement for the timeseries
        """
        if max_depth == 0:
            if arr is None:
                raise UserWarning('Need max depth > 0 if no arr is given ')
            max_value_index = np.argmax(arr[:,1])
            newarr = arr[max_value_index-1:max_value_index + 2]
            output = np.polyfit(newarr[:,0],newarr[:,1],deg = 2) #change back to sn once that's working 
            f = np.poly1d(output)
            x_new = np.linspace(newarr[0][0], newarr[-1][0], num = 10000)
            y_new = f(x_new)
            new_dm = x_new[np.argmax(y_new)]
            peak_SN = np.amax(y_new)
            self.DM_arr = arr
            
            if plot is True:
                doshow = False 
                if ax is None:
                    doshow = True
                    fig = figure()
                    ax = fig.add_subplot(111)
                ax.scatter(arr[:,0],arr[:,1])
                ax.axvline(self.getDM(), color = 'r')
                min_tick = np.amin(arr[:,0]) - .005*np.amin(arr[:,0])
                max_tick = np.amax(arr[:,0]) + .005*np.amin(arr[:,0])
                ax.xaxis.set_ticks(np.arange(min_tick,max_tick,(max_tick-min_tick)/5.0))
                ax.xaxis.set_major_formatter(FormatStrFormatter('%.02f'))
                if labels:
                    ax.axes.set_ylabel('Peak SNR')
                    ax.axes.set_xlabel('DM')
                if doshow:
                    show()
            
            #try:
            #    fwhm, L, R = u.FWHM(DM_arr[:,1],notcentered=True) #check not centered tag!!!!
            #except:
            #    fwhm = 0
            #fwhm = sep * fwhm #does this give the correct width???
            #fwhm, L, R = u.FWHM(dm_arr[:,1], simple = True, notcentered = True)
            #lindex = list(dm_arr[:,1]).index(L)
            #rindex = list(dm_arr[:,1]).index(R)
            #
            #fwhm = dm_arr[rindex][0] - dm_arr[lindex][0]
            fwhm = 0
            
            self.calculated_DM = new_dm
            return ax, arr, new_dm, peak_SN, fwhm
        
        
        if DM is None:
            DM = self.getDM()
        dm_arr = self.FitSN(iters = iters, initial= initial, end = end, DM = DM, plot = False)
        if arr is not None:
            dm_arr = np.concatenate((arr,dm_arr), axis = 0)
            dm_arr = dm_arr[dm_arr[:,0].argsort()]
        opt_dm_pos = np.argmax(dm_arr[:,1])
        new_dm = dm_arr[opt_dm_pos,0]
        sub_arr =  dm_arr[opt_dm_pos-1:opt_dm_pos+2,0] 
        initial = sub_arr[0]/new_dm
        end = sub_arr[-1]/new_dm
        next_depth = max_depth - 1
        self.DMCalculator(ax=ax, plot=plot, labels = labels, DM = new_dm, iters = iters, initial=initial, end = end, arr = dm_arr, max_depth = next_depth)
        
   
    def save_to_npz(self, name = None):
        if name is None:
            name = self.name
        if self.npz is None:
            SS_xaxis, SS_yaxis = self.getSecondarySpectrumAxes()
            ACF2D_xaxis, ACF2D_yaxis = self.getAcf2dAxes()
            average_profile = self.average_profile,
            spavgprof_mpw = self.spavgprof.mpw,
            spavgprof_mpw_data = self.spavgprof.data[self.spavgprof.mpw]
            spavgprof_opw = self.spavgprof.opw, 
            spavgprof_opw_data = self.spavgprof.data[self.spavgprof.opw]
            np.savez(name,
                     telescope = self.getTelescope(),
                     MJD = self.getMJD(),
                     BW = self.getBandwidth(),
                     DM = self.getDM(),
                     DUR = self.getDuration(),
                     PER = self.getPeriod(),
                     NBIN = self.getNbin(),
                     SUB = self.getNsubint(),
                     NCHAN = self.getNchan(),
                     NPOL = self.getNpol(),
                     time = self.T,
                     frequency = self.F,
                     profiles = self.getProfiles(difference = False), 
                     differenceprofiles = self.getProfiles(difference = True),
                     DynamicSpectrum = self.getDynamicSpectrum(),
                     ACF2D = self.getAcf2d(),
                     ACF2D_xaxis = ACF2D_xaxis, ACF2D_yaxis = ACF2D_yaxis,
                     SecondarySpectrum = self.getSecondarySpectrum(),
                     SS_xaxis = SS_xaxis, SS_yaxis = SS_yaxis,
                     average_profile = average_profile,
                     spavgprof_mpw = spavgprof_mpw, spavgprof_mpw_data = spavgprof_mpw_data,
                     spavgprof_opw = spavgprof_opw, spavgprof_opw_data = spavgprof_opw_data,
                     peak_DM = self.calculated_DM,
                     DM_arr = self.DM_arr)
        else:
            raise UserWarning('Cannot save from .npz to .npz, this can cuase errors')
        
     
    ###GETTERS####  
    def getName(self,full=False):
        if full:
            return self.fullname
        return self.name
            
    def getTelescope(self):
        if self.npz is None:
            return self.ar.getTelescope()
        else:
            return self.npz['telescope']
            
    def getMJD(self):
        if self.npz is None:
            return self.ar.getMJD()
        else:
            return self.npz['MJD']
    
    def getBandwidth(self):
        if self.npz is None:
            return self.ar.getBandwidth()
        else:
            return self.npz['BW']
    
    def getDM(self):
        if self.npz is None:
            return self.ar.getDM()
        else:
            return self.npz['DM']
    
    def getDuration(self):
        if self.npz is None:
            return self.ar.getDuration()
        else:
            return self.npz['DUR']
       
    def getPeriod(self):
        if self.npz is None:
            return self.ar.getPeriod()
        else:
            return self.npz['PER']
         
    def getNbin(self):
        if self.npz is None:
            return self.ar.getNbin()
        else:
            return self.npz['NBIN']

    def getNsubint(self):
        if self.npz is None:
            return self.ar.getNsubint()
        else:
            return self.npz['NSUB']

    def getNchan(self):
        if self.npz is None:
            return self.ar.getNchan()
        else:
            return self.npz['NCHAN']

    def getNpol(self):
        if self.npz is None:
            return self.ar.getNpol()
        else:
            return self.npz['NPOL']

