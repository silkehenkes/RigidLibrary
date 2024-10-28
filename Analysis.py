#!/usr/bin/python

# Silke Henkes, 29.08.18: 
# Created Analysis class to actually analyse results and plot them
# Detangled code parts, contains: 
# - Stress statistics
# - Cluster statistics
# - Cluster and stress cross-correlations
# - plotStresses, and plotPebbles detailed plots
# - Rigidity analysis for modes, and plot
# - Correlations between modes and pebble rigidity, and plot
# - Correlations between displacements and pebble rigidity, and plot
# - Analysis and plotting of D2_min



import random
import sys, os, glob

# Each analysis has a Configuration, a Pebbles and a Hessian 
# (of which Pebbles and Hessian are allowed to be empty shells if you're not using these functions)
from Configuration import *
from Pebbles import *
from Hessian import *
from Tiling import *

# Note: changed to simple pickle on move to python 3. Compatibility problems with loading old pickle files expected!
#import cPickle as pickle
import pickle as pickle
import copy as cp
import numpy as np

# This does need all the plotting libraries
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.lines as lne
from matplotlib.colors import LinearSegmentedColormap

#import matplotlib
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 4
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=16.0
matplotlib.rcParams['legend.fontsize']=14.0


class Analysis:

    ### =========== Analysis constructor: Give it a configuration (mandatory), and Pebbles and Hessian, which are allowed to be empty shells
    def __init__(self,conf0,pebbles0,hessian0,tiling0='skip',fgiven0=0.001,verbose0=False):
        self.conf=conf0
        self.pebbles=pebbles0
        self.hessian=hessian0
        self.tiling=tiling0
        self.verbose=verbose0
        self.fgiven=fgiven0
        # basic dimensions ...
        self.N=self.conf.N
        self.ncon=self.conf.ncon
        self.Lx=self.conf.Lx
        self.Ly=self.conf.Ly
        # Data for the Maxwell-Cremona tiling plotting
        self.I = self.conf.I
        self.J = self.conf.J
        self.fn = self.conf.fnor
        self.ft = self.conf.ftan
        self.fullmobi = self.conf.fullmobi
        self.ftot = np.sqrt(np.square(self.fn)+np.square(self.ft)) #Use this for color scheme forces?
        if self.tiling != 'skip':
            self.tiles = self.tiling.tiles
        
        # this is the distance beween lines in the double contact plots
        if self.conf.datatype == 'experiment_square':
                self.small = 8.0
        elif self.conf.datatype == 'experiment_annulus':
                self.small = 12.0
        elif self.conf.datatype == 'simulation':
                self.small = 0.2
        elif self.conf.datatype =='lattice':
            self.small = 0.1
        else:
            self.small = 0.2

    
    # ====== Color helper function, essentially give jet color map in a convenient form ====
    def color_init(self,Fmax=3.):
        # Maximum set by default to 3 times the average
        # else use input as max
        # approximation to color map jet
        cl=( [ 0   , (0.,0. ,0.5) ],
                [ 0.5, (0.,0. ,1.) ],
                [ 1. , (0.,0.5,1.) ],
                [ 1.25 , (0.5,0.75,0.5) ],
                [ 1.7 , (1.,0.5 ,0.) ],
                [ 2.0 , (1.,0.,0.) ],
                [ Fmax , (0.5,0.,0.) ])
        n=len(cl)
        red=interp1d([l[0] for l in cl],[l[1][0] for l in cl])
        green=interp1d([l[0] for l in cl],[l[1][1] for l in cl])
        blue=interp1d([l[0] for l in cl],[l[1][2] for l in cl])
        def Fcolor(p):
            return float(red(p)),float(green(p)),float(blue(p))

        cm= { 'red'  : [(l[0]/Fmax,l[1][0],l[1][0]) for l in cl] ,
                'green': [(l[0]/Fmax,l[1][1],l[1][1]) for l in cl] ,
                'blue' : [(l[0]/Fmax,l[1][2],l[1][2]) for l in cl] }

        Fcmap=LinearSegmentedColormap('F_cm',cm,N=256)
        
        return Fcolor, Fcmap
    
    #================================= PLOTTING ======================================
    
    # ===== Plotting force and Contact related pieces
    def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
        if 'figure' in kwargs:
            fig=plt.figure(figsize=(20,20))
            axval= fig.add_subplot(1, 1, 1)
        elif 'axis' in kwargs:
            axval=kwargs['axis']
        else:
            fig=plt.figure()
            axval= fig.add_subplot(1, 1, 1)
        if 'timestamp' in kwargs:
            axval.text(0.0,0.4*self.Ly,'T= ' +str(kwargs['timestamp']),fontsize=18)
        for k in range(self.N):
            # Particle circles
            if (plotCir):
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.5, 0.5, 0.5),fc='none', linewidth=2)
                axval.add_patch(cir1)
                #angline=lne.Line2D([self.x[k],self.x[k]+self.rad[k]*np.cos(self.alpha[k])],[self.y[k],self.y[k]+self.rad[k]*np.sin(self.alpha[k])],color=(0.5, 0.5, 0.5))
                #plt.gca().add_line(angline)
                if self.verbose:
                    axval.text(self.conf.x[k],self.conf.y[k],str(k),fontsize=8)
            # Velocity vectors
            if (plotVel):
                scale=2*np.mean(self.conf.rad)/np.sqrt(np.mean(self.conf.dx*self.conf.dx)+np.mean(self.conf.dy*self.conf.dy))
                velvect=ptch.Arrow(self.conf.x[k],self.conf.y[k],scale*self.conf.dx[k],scale*self.conf.dy[k],width=0.3,color='b')
                axval.add_patch(velvect)
                
        # Plotting contacts
        if (plotCon):
            for k in range(len(self.conf.I)):
                x0,x1,y0,y1=self.conf.getConPos(k)
                if (self.fullmobi[k]==1):
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.7,0,0),lw=3.0)
                    axval.add_line(conlink)
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0,0,0),lw=3.0)
                    axval.add_line(conlink)
                if self.verbose:
                    if (k>=self.ncon):
                        axval.text(x0+0.5*(x1-x0)-self.small*np.sin(ang),y0+0.5*(y1-y0)+self.small*np.cos(ang),str(k),fontsize=8)
                    else:
                        axval.text(x0+0.5*(x1-x0),y0+0.5*(y1-y0),str(k),fontsize=8)
                        
        if (plotF):
            #fscale=np.mean(self.fnor)
            fscale=self.fgiven
            Fcolor,Fmap=self.color_init(np.sqrt(np.amax(self.conf.fnor))/fscale)
            for k in range(len(self.conf.I)):
                fval=np.sqrt(self.conf.fnor[k]/fscale)
                x0,x1,y0,y1=self.conf.getConPos(k)
                conlink=lne.Line2D([x0,x1],[y0,y1],color=Fcolor(fval),lw=2*fval)
                axval.add_line(conlink)
#				if self.verbose:
#					if (k>=self.ncon):
#						axval.text(x0+0.5*(x1-x0)-self.small*np.sin(ang),y0+0.5*(y1-y0)+self.small*np.cos(ang),str(k),fontsize=8)
#					else:
#						axval.text(x0+0.5*(x1-x0),y0+0.5*(y1-y0),str(k),fontsize=8)
                        
        if (plotStress):
            print("Plot Stress is still under construction.")
        plt.title('Forces and contacts')				
        axval.axis('scaled')
        if self.conf.datatype=='simulation':
            axval.set_xlim(-self.Lx/2,self.Lx/2)
            axval.set_ylim(-self.Ly/2,self.Ly/2)
        elif self.conf.datatype=='experiment':
            xmin=np.amin(self.conf.x)-20
            ymin=np.amin(self.conf.y)-20
            axval.set_xlim(xmin,xmin+self.Lx+120)
            axval.set_ylim(ymin,ymin+self.Ly+120)
        # returning figure pointer for use outside (like plotting or saving)
        if 'figure' in kwargs:
            return fig
        elif 'axis' in kwargs:
            return axval
        else:
            return fig
        
    def cluster_sorter(self):
        return labels
    
    # =============== pebble plotting ========================================
    # boolean arguments are plotting options (what to add or not)
    #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,axval,running):
    def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):
        if 'figure' in kwargs:
            fig=plt.figure(figsize=(20,20))
            #fig=plt.figure()
            axval= fig.add_subplot(1, 1, 1)
        elif 'axis' in kwargs:
            axval=kwargs['axis']
        else:
            fig=plt.figure()
            axval= fig.add_subplot(1, 1, 1)
        if 'timestamp' in kwargs:
            axval.text(0.0,0.4*self.Ly,'T= ' +str(kwargs['timestamp']),fontsize=18)
        for k in range(self.N):
            t=0
            # Particle circles
            if (plotCir):
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.5, 0.5, 0.5),fc='none', linewidth=2)
                axval.add_patch(cir1)
                #angline=lne.Line2D([self.x[k],self.x[k]+self.rad[k]*np.cos(self.alpha[k])],[self.y[k],self.y[k]+self.rad[k]*np.sin(self.alpha[k])],color=(0.5, 0.5, 0.5))
                #axval.add_line(angline)
                if self.verbose:
                    axval.text(self.conf.x[k],self.conf.y[k],str(k),fontsize=8)
            # The pebbles: one red, two green, three red
            if (plotPeb):
                for j in range(3):
                    if(self.pebbles.pebbles[k,j]==-1):
                        t=t+1;
                        if(t==1):
                            color='r'
                        if(t==2):
                            color='g'
                        if(t==3):
                            color='b'
                        cir2=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=0.75*self.conf.rad[k],ec=(0.5, 0.5, 0.5),fc=color,linewidth=2)
                        axval.add_patch(cir2)

            
        # Plotting clusters, using the external (global!) random_read color pattern
        if (plotClus):
            self.rgbclus=np.empty((self.pebbles.ncon2,3))
            if self.pebbles.cidx>0:
                rgbpattern=np.random.rand(self.pebbles.cidx+1,3)
            else:
                # no clusters found, allocate dummy matrix for smooth run of doing nothing
                rgbpattern=np.random.rand(1,3)
            print(self.pebbles.cidx)
            rgbpattern[0,:]=np.array([0,0,1])
            for k in range(self.pebbles.ncon2):
                if (self.pebbles.cluster[k]>=0):
                    # check for global colormap; else use random pattern
                    try:
                        self.rgbclus[k,:]=self.random_read[self.pebbles.cluster[k].astype(int),:]
                    except:
                        self.rgbclus[k,:]=rgbpattern[self.pebbles.cluster[k].astype(int),:]
                        #print(rgbclus[k,:])
                else:
                    self.rgbclus[k,:]=np.array([0,0,0])
            for i in range(self.N):
                if len(self.pebbles.pcluster[i])>0:
                    try:
                        color=self.random_read[self.pebbles.pluster[i][0],:]
                    except:
                        color=rgbpattern[self.pebbles.pcluster[i][0],:]
                        #print(rgbclus[k,:])
                    cir2=ptch.Circle((self.conf.x[i],self.conf.y[i]),radius=self.conf.rad[i],ec=(0.5, 0.5, 0.5),fc=color,linewidth=2)
                    axval.add_patch(cir2)

            for k in range(len(self.pebbles.Ifull)):
                # this version depends on pebble numbering, so use 2nd version
                x0,x1,y0,y1=self.conf.getConPos2(self.pebbles.Ifull[k],self.pebbles.Jfull[k])
                if (self.pebbles.cluster[k]!=-1):
                    if (k>=self.ncon):
                        ang=np.arctan((y1-y0)/(x1-x0))
                        conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=3,color=(self.rgbclus[k,0],self.rgbclus[k,1],self.rgbclus[k,2]))
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=3,color=(self.rgbclus[k,0],self.rgbclus[k,1],self.rgbclus[k,2]))
                    axval.add_line(conlink)
                else:
                    if (k>=self.ncon):
                        ang=np.arctan((y1-y0)/(x1-x0))
                        conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=3,color=(0.65,0.65,0.65))
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=3,color=(0.65,0.65,0.65))
                    axval.add_line(conlink)
                # contact label for debugging
                if self.verbose:
                    if (k>=self.ncon):
                        axval.text(x0+0.5*(x1-x0)-self.small*np.sin(ang),y0+0.5*(y1-y0)+self.small*np.cos(ang),str(k),fontsize=8)
                    else:
                        axval.text(x0+0.5*(x1-x0),y0+0.5*(y1-y0),str(k),fontsize=8)
                        
        # Plotting the overconstrained contacts
        if (plotOver):
            for k in range(len(self.pebbles.Ifull)):
                x0,x1,y0,y1=self.conf.getConPos2(self.pebbles.Ifull[k],self.pebbles.Jfull[k])
                if (self.pebbles.Overcon[k]==True):
                    if (k>=self.ncon):
                        ang=np.arctan((y1-y0)/(x1-x0))
                        conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=3,color=(0.8,0.2,0.2))
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=3,color=(0.8,0.2,0.2))
                    axval.add_line(conlink)
                else:
                    if (k>=self.ncon):
                        ang=np.arctan((y1-y0)/(x1-x0))
                        conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=3,color=(0.65,0.65,0.65))
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=3,color=(0.65,0.65,0.65))
                    axval.add_line(conlink)
                # contact label for debugging
                if self.verbose:
                    if (k>=self.ncon):
                        axval.text(x0+0.5*(x1-x0)-self.small*np.sin(ang),y0+0.5*(y1-y0)+self.small*np.cos(ang),str(k),fontsize=8)
                    else:
                        axval.text(x0+0.5*(x1-x0),y0+0.5*(y1-y0),str(k),fontsize=8)
        # Plotting the contacts that were added during the pebble game
        if (plotPebCon):
            for k in range(len(self.pebbles.Ifull)):
                x0,x1,y0,y1=self.conf.getConPos2(self.pebbles.Ifull[k],self.pebbles.Jfull[k])
                if (self.pebbles.ptr[k]!=-1):
                    if (k>=self.ncon):
                        ang=np.arctan((y1-y0)/(x1-x0))
                        conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=3,color=(0.6,0,0.0))
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=3,color=(0.6,0,0))
                    axval.add_line(conlink)
                else:
                    if (k>=self.ncon):
                        ang=np.arctan((y1-y0)/(x1-x0))
                        conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=3,color=(0.5,0.5,0.5))
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=2,color=(0.5,0.5,0.5))
                    axval.add_line(conlink)
                if self.verbose:
                    if (k>=self.ncon):
                        axval.text(x0+0.5*(x1-x0)-self.small*np.sin(ang),y0+0.5*(y1-y0)+self.small*np.cos(ang),str(k),fontsize=8)
                    else:
                        axval.text(x0+0.5*(x1-x0),y0+0.5*(y1-y0),str(k),fontsize=8)  
        axval.axis('scaled')
        if self.conf.datatype=='simulation':
            axval.set_xlim(-self.Lx/2,self.Lx/2)
            axval.set_ylim(-self.Ly/2,self.Ly/2)
        elif self.conf.datatype=='experiment_square':
            xmin=np.amin(self.conf.x)-20
            ymin=np.amin(self.conf.y)-20
            axval.set_xlim(xmin,xmin+self.Lx+120)
            axval.set_ylim(ymin,ymin+self.Ly+120)
        axval.set_title('System with ' + str(self.pebbles.freepeb) + ' free pebbles and ' + str(self.pebbles.fail) + ' failed contacts')
        # returning figure pointer for use outside (like plotting or saving)
        if 'figure' in kwargs:
            return fig
        elif 'axis' in kwargs:
            return axval
        else:
            return fig
                    

        # ========================== ANALYSIS =========================================
        
    #========================== Conventional force, stress and contact statistics ===========
    # This one lives in conf since it depends wholly on conf variables, i.e. positions and forces
    def getStressStat(self):
        zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx=self.conf.getStressStat()
        print ("Mean coordination number: " + str(zav))
        print ("Mean mobilisation: " + str(nm))
        print ("Mean pressure: " + str(pres))
        print ("Force balance: " +str(fxbal)+ "x, " + str(fybal) + "y")
        print ("Torque balance: " + str(torbal))
        return zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx
    
    # ============== Cluster statistics ===================
    def clusterStatistics(self):
        # Cluster statistics, contact and particle based
        # Total contacts, fraction part of cluster
        incluster=[val for val in self.pebbles.cluster if val>-1]
        frac=len(incluster)/(1.0*self.pebbles.ncon2)
        print('Cluster fraction ' + str(frac))
        # Largest cluster, fraction of contacts
        islargest=np.nonzero(self.pebbles.cluster==self.pebbles.maxidx)[0]
        fracmax=len(islargest)/(1.0*self.pebbles.ncon2)
        print('Largest Cluster fraction ' + str(fracmax))
        # Largest cluster, size in x and y direction. Sort and concatenate all the x values. Look for largest gap. lenx is L-gap
        # Rinse and repeat for y. Disregard gaps the size of particle distance
        if (len(islargest)>0):
            gapmax=2*np.amax(self.conf.rad)
            Iarr=np.array(self.pebbles.Ifull)
            Jarr=np.array(self.pebbles.Jfull)
            xval=np.sort(np.union1d(self.conf.x[Iarr[islargest]],self.conf.x[Jarr[islargest]]))
            # Differences between consecutive elements of an array
            xdiff=np.ediff1d(xval)
            addl=xval[0]-xval[-1]+self.conf.Lx
            #xdiff[0]+=self.L
            gap=max(np.amax(xdiff),addl)
            if gap>gapmax:
                lenx=self.Lx-gap
            else:
                lenx=self.Lx
            print('cluster length x ' +str(lenx))
            yval=np.sort(np.union1d(self.conf.y[Iarr[islargest]],self.conf.y[Jarr[islargest]]))
            ydiff=np.ediff1d(yval)
            addl=yval[0]-yval[-1]+self.Ly
            gap=max(np.amax(ydiff),addl)
            if gap>gapmax:
                leny=self.Ly-gap
            else:
                leny=self.Ly
            print('cluster length y ' +str(leny))
        else:
            lenx=0.0
            leny=0.0
        return frac,fracmax,lenx,leny
        
    ###### ======= Cross-correlation between between clusters and local sliding (see also displacement decomposition below)
    def crossCorrelate(self):
        # Correlation between contact force and being rigid
        # Note: going only to ncon. The further (double) bonds up to ncon2 are not represented in the force list, and by construction
        # they belong to the same cluster as the original bond.
        incluster = np.nonzero(self.pebbles.cluster[0:self.ncon]>-1)[0]
        outcluster = np.nonzero(self.pebbles.cluster[0:self.ncon]==-1)[0]
        favn_in=np.mean(self.conf.fnor[incluster])
        favn_out=np.mean(self.conf.fnor[outcluster])
        mobav_in=np.mean(self.conf.fullmobi[incluster])
        mobav_out=np.mean(self.conf.fullmobi[outcluster])
        # sliding at the contact between frames (normalize at some point!)
        # convert into some non-affine measure some time ...
        drijx=self.conf.dx[self.conf.J]-self.conf.dx[self.conf.I]
        drijy=self.conf.dy[self.conf.J]-self.conf.dy[self.conf.I]
        drijnor=self.conf.nx*drijx+self.conf.ny*drijy
        drijtan=-self.conf.ny*drijx+self.conf.nx*drijy
        dnor_in=np.mean(np.abs(drijnor[incluster]))
        dnor_out=np.mean(np.abs(drijnor[outcluster]))
        dtan_in=np.mean(np.abs(drijtan[incluster]))
        dtan_out=np.mean(np.abs(drijtan[outcluster]))
        
        #if (self.verbose):
        print('normal force in cluster ' + str(favn_in))
        print('normal force outside cluster ' + str(favn_out))
        print('mobilization in cluster ' + str(mobav_in))
        print('mobilization outside cluster ' + str(mobav_out))
        print('normal sliding in cluster ' + str(dnor_in))
        print('normal sliding outside cluster ' + str(dnor_out))
        print('tangential sliding in cluster ' + str(dtan_in))
        print('tangential sliding outside cluster ' + str(dtan_out))
        
        # Contact sliding and orientation correlations (setting to rest the DOS issue if fully mobilized contacts slide,
        # and if yes, in which direction)
        ismobi=np.nonzero(self.conf.fullmobi==1)[0]
        notmobi=np.nonzero(self.conf.fullmobi==0)[0]
        # This only works if angles are defined, i.e. not in the experiment
        if self.conf.hasAngles:
            # See equation 2 of the Henkes et al. 2010 EPL paper:
            # For our packings with Hertz-Mindlin forces the effective potential is given by [V]
            # where deltat = (delta r.t) - (Ri delta alpha_i + Rj delta alpha_j) is the tangential 
            # displacement of the contact point of the two particles,
            # as illustrated in fig. 2.
            dijtang = drijtan - (self.conf.rad[self.conf.I]*self.conf.dalpha[self.conf.I]+self.conf.rad[self.conf.J]*self.conf.dalpha[self.conf.J])
            
            dtang_in=np.mean(np.abs(dijtang[incluster]))
            dtang_out=np.mean(np.abs(dijtang[outcluster]))
            dtang_mob=np.mean(np.abs(dijtang[ismobi]))
            dtang_not=np.mean(np.abs(dijtang[notmobi]))
            dtang_mob_signed=np.mean(dijtang[ismobi]*np.sign(self.conf.ftan[ismobi]))
            dtang_not_signed=np.mean(dijtang[notmobi]*np.sign(self.conf.ftan[notmobi]))
            
            #if (self.verbose):
            print('tang sliding in cluster ' + str(dtang_in))
            print('tang sliding outside cluster ' + str(dtang_out))
            print('tang sliding mobilized ' + str(dtang_mob))
            print('signed version ' + str(dtang_mob_signed))
            print('tang sliding not mobilized ' + str(dtang_not))
            print('signed version ' + str(dtang_not_signed))
        else:
            print('No angles in data set, tang variables are set to 0')
            dtang_in=0
            dtang_out=0
            dtang_mob=0
            dtang_mob_signed=0
            dtang_not=0
            dtang_not_signed=0
        
        return favn_in, favn_out, mobav_in, mobav_out, dnor_in, dnor_out, dtan_in, dtan_out, \
            dtang_in, dtang_out, dtang_mob, dtang_mob_signed, dtang_not, dtang_not_signed
            

                
    
    # ================ Compute which bits of the clusters are rigid, as per normal modes ==========
    def ModeClusters(self,mode,thresh=2e-5):
        # Attempt to identify rigid clusters by several methods
        fig=plt.figure()
        axval=fig.add_subplot(1, 1, 1)
        
        # Level 0: threshold on square displacements on particle (particle clusters)
        # Get the zero modes according to the hessian first
        #ModeThresh(self,thresh,plotRigid=False):	
        self.hessian.ModeThresh(thresh)
        # have now determined which of these are rigid and which not
        for k in self.hessian.isrigid:
            cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.0, 0.0, 0.0),fc='none', linewidth=2)
            axval.add_patch(cir1)
        for k in self.hessian.notrigid:
            cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.65, 0.65, 0.65),fc='none', linewidth=2)
            axval.add_patch(cir1)
        axval.plot(self.conf.x[self.hessian.isrigid],self.conf.y[self.hessian.isrigid],'.k')

        # Level 1: Work with relative displacements. That still leaves global rotations possible
        # Start to work on contacts
        # normal, tangential and rotational (eff sliding for a little, look for gearing in a bit)
        # This is now done within the Hessian class through ModeContacts
        eign, eigt, eigr, eiggear = self.hessian.ModeContacts()
        # Replace by internal threshold instead!
        thresh=0.25*np.mean(eign+eigt)

        # Make proper threshold plots
        for k in range(self.ncon):
            
            x0,x1,y0,y1=self.conf.getConPos(k)
            if mode == 'translations':
                if (eign[k]+eigt[k])<thresh:
                        if eigr[k]>thresh:
                            conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.7,0,0),lw=3.0)
                        else:
                            conlink=lne.Line2D([x0,x1],[y0,y1],color=(0,0.0,1.0),lw=3.0)
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.65,0.65,0.65),lw=3.0)
                axval.add_line(conlink)
                axval.set_title('translations')
            elif mode == 'rotations':
                if (eigr[k])<thresh:
                    if (eign[k]+eigt[k])>thresh:
                        conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.7,0,0),lw=3.0)
                    else:
                        conlink=lne.Line2D([x0,x1],[y0,y1],color=(0,0.0,1.0),lw=3.0)
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.65,0.65,0.65),lw=3.0)
                axval.add_line(conlink)
                axval.set_title('rotations')
            elif mode == 'strict':
                if ((eigr[k])<thresh) and ((eign[k]+eigt[k])<thresh):
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0,0.0,1.0),lw=3.0)
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.65,0.65,0.65),lw=3.0)
                axval.add_line(conlink)
                axval.set_title('strict')
            elif mode == 'sliding':
                if (eiggear[k]<thresh):
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(1.0,0.0,0.0),lw=3.0)
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],color=(0.65,0.65,0.65),lw=3.0)
                axval.add_line(conlink)
                axval.set_title('not sliding')
            else:
                print ("Unknown mode, doing nothing!")
        axval.axis('scaled')
        if self.conf.datatype=='simulation':
            axval.set_xlim(-self.Lx/2,self.Lx/2)
            axval.set_ylim(-self.Ly/2,self.Ly/2)
        elif self.conf.datatype=='experiment':
            xmin=np.amin(self.conf.x)-20
            ymin=np.amin(self.conf.y)-20
            axval.set_xlim(xmin,xmin+self.Lx+120)
            axval.set_ylim(ymin,ymin+self.Ly+120)
        return fig
    
    # ========== Correlation between pebble rigidity and modes, with optional plot ========= 
    def RigidModesCorrelate(self,thresh,plotDisp=True):
    # Correlations
        fig=plt.figure()
        Corr_PGDM=np.zeros(4)
        P_eig_if_pebble=0
        P_pebble_if_eig=0
        # Compute correlations between being part of a rigid cluster and contact displacements
        eign, eigt, eigr, eiggear = self.hessian.ModeContacts()
        for k in range(self.ncon):
            if self.pebbles.cluster[k]>-1:
                Corr_PGDM[0]+=1.0
                if (eigr[k])<thresh and (eign[k]+eigt[k])<thresh:
                    Corr_PGDM[1]+=1.0
            if (eigr[k])<thresh and (eign[k]+eigt[k])<thresh:
                    Corr_PGDM[2]+=1.0
                    if self.pebbles.cluster[k]>-1:
                        Corr_PGDM[3]+=1.0
                            
        if plotDisp:
            axval= fig.add_subplot(1, 1, 1)
            for k in range(self.N):
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.7, 0.7, 0.7),fc='none', linewidth=2)
                axval.add_patch(cir1)
            Fcolor,Fmap=self.color_init(5.0)
            intThresh=np.mean(eign+eigt)
            for k in range(self.ncon):
                fval=np.sqrt((eign[k]+eigt[k])/intThresh)
                if fval<0.2:
                    fval=0.2
                x0,x1,y0,y1=self.conf.getConPos(k)
                conlink=lne.Line2D([x0,x1],[y0,y1],lw=1.0/fval,color=Fcolor(1.0/fval))
                axval.add_line(conlink)
            axval.axis('scaled')
            if self.conf.datatype=='simulation':
                    axval.set_xlim(-self.Lx/2,self.Lx/2)
                    axval.set_ylim(-self.Ly/2,self.Ly/2)
            elif self.conf.datatype=='experiment':
                    xmin=np.amin(self.conf.x)-20
                    ymin=np.amin(self.conf.y)-20
                    axval.set_xlim(xmin,xmin+self.Lx+120)
                    axval.set_ylim(ymin,ymin+self.Ly+120)
            plt.title('Relative lack of motion from zero modes')
        # Conditional probabilities: Eigenvalue rigid when pebble rigid
        if Corr_PGDM[0]>0:
            P_eig_if_pebble=Corr_PGDM[1]/Corr_PGDM[0]
        else:
            P_eig_if_pebble=0.0
        if Corr_PGDM[2]>0:
            P_pebble_if_eig=Corr_PGDM[3]/Corr_PGDM[2]
        else:
            P_pebble_if_eig=0.0
        if plotDisp:
            return P_eig_if_pebble,P_pebble_if_eig, fig
        else:
            return P_eig_if_pebble,P_pebble_if_eig, fig
 
    # ========= Corrlation between displacements and pebble rigidity, with optinal plot ======
    def RigidDisplacementsCorrelate(self,minThresh,plotDisp=True):
        # Correlations
        Corr_PGDP=np.zeros(4)
        P_disp_if_pebble=0
        P_pebble_if_disp=0
        Corr_DPDM=np.zeros(4)
        P_disp_if_eig=0
        P_eig_if_disp=0
        eign, eigt, eigr, eiggear = self.hessian.ModeContacts()
        # Compute correlations between being part of a rigid cluster and contact displacements
        if self.conf.hasAngles:
            disp2n, disp2t, disp2r, disp2gear,intThresh = self.conf.Disp2Contacts(minThresh,True)
        else:
            disp2n, disp2t, intThresh = self.conf.Disp2Contacts(minThresh,True)
        for k in range(self.ncon):
            if self.pebbles.cluster[k]>-1:
                Corr_PGDP[0]+=1.0
                if self.conf.hasAngles:
                    if (disp2r[k])<intThresh and (disp2n[k]+disp2t[k])<intThresh:
                        Corr_PGDP[1]+=1.0
                else:
                    if (disp2n[k]+disp2t[k])<intThresh:
                        Corr_PGDP[1]+=1.0
            if (eigr[k])<minThresh and (eign[k]+eigt[k])<minThresh:
                Corr_DPDM[0]+=1.0
                if self.conf.hasAngles:
                    if (disp2r[k])<intThresh and (disp2n[k]+disp2t[k])<intThresh:
                        Corr_DPDM[1]+=1.0
                else:
                    if (disp2n[k]+disp2t[k])<intThresh:
                        Corr_DPDM[1]+=1.0
            if self.conf.hasAngles:
                if (disp2r[k])<intThresh and (disp2n[k]+disp2t[k])<intThresh:
                    Corr_PGDP[2]+=1.0
                    Corr_DPDM[2]+=1.0
                    if self.pebbles.cluster[k]>-1:
                        Corr_PGDP[3]+=1.0
                    if (eigr[k])<minThresh and (eign[k]+eigt[k])<minThresh:
                        Corr_DPDM[3]+=1.0

            else:
                if (disp2n[k]+disp2t[k])<intThresh:
                    Corr_PGDP[2]+=1.0
                    Corr_DPDM[2]+=1.0
                    if self.pebbles.cluster[k]>-1:
                        Corr_PGDP[3]+=1.0
                    if (eigr[k])<minThresh and (eign[k]+eigt[k])<minThresh:
                        Corr_DPDM[3]+=1.0
                                
        if plotDisp:
            fig=plt.figure()
            axval= fig.add_subplot(1, 1, 1)
            plt.quiver(self.conf.x,self.conf.y,self.conf.dx,self.conf.dy)
            for k in range(self.N):
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.65, 0.65, 0.65),fc='none', linewidth=2)
                axval.add_patch(cir1)
            Fcolor,Fmap=self.color_init(5.0)
            for k in range(self.ncon):
                fval=np.log(disp2n[k]+disp2t[k])+3
                if fval<0.000001:
                    fval=0.000001
                if fval>4.999999:
                    fval=4.999999
                x0,x1,y0,y1=self.conf.getConPos(k)
                conlink=lne.Line2D([x0,x1],[y0,y1],lw=5-fval,color=Fcolor(5-fval))
                axval.add_line(conlink)
            axval.axis('scaled')
            if self.conf.datatype=='simulation':
                axval.set_xlim(-self.Lx/2,self.Lx/2)
                axval.set_ylim(-self.Ly/2,self.Ly/2)
            elif self.conf.datatype=='experiment':
                xmin=np.amin(self.conf.x)-20
                ymin=np.amin(self.conf.y)-20
                axval.set_xlim(xmin,xmin+self.Lx+120)
                axval.set_ylim(ymin,ymin+self.Ly+120)
            plt.title('Actual relative lack of motion')
        # Conditional probabilities: Eigenvalue rigid when pebble rigid
        if Corr_PGDP[0]>0:
            P_disp_if_pebble=Corr_PGDP[1]/Corr_PGDP[0]
        else:
            P_disp_if_pebble=0.0
        if Corr_PGDP[2]>0:
            P_pebble_if_disp=Corr_PGDP[3]/Corr_PGDP[2]
        else:
            P_pebble_if_disp=0.0

        if Corr_DPDM[0]>0:
            P_disp_if_eig=Corr_DPDM[1]/Corr_DPDM[0]
        else:
            P_disp_if_eig=0.0
        if Corr_DPDM[2]>0:
            P_eig_if_disp=Corr_DPDM[3]/Corr_DPDM[2]
        else:
            P_eig_if_disp=0.0
            
        if plotDisp:
            return P_disp_if_pebble,P_pebble_if_disp,P_disp_if_eig,P_eig_if_disp, fig
        else:
            return P_disp_if_pebble,P_pebble_if_disp,P_disp_if_eig,P_eig_if_disp
        
    # ==== Computes D2min, and plots it together with the rigid cluster ===
    # At this point unclear if it is a useful measure ...
    def DisplacementCorrelateD2min(self,plotdis=True,threshold_rad=10.0):
        fig=plt.figure()
        axval=fig.add_subplot(1, 1, 1)
        for k in range(self.N):
            if len(self.pebbles.pcluster[k])>0:
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0,0,0),fc='none',linewidth=2)
            else:
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.7,0.7,0.7),fc='none',linewidth=2)
            axval.add_patch(cir1)

        if (plotdis):
            
            D2_min=self.conf.getD2min(threshold_rad)
            #plt.figure()
            #n, bins, patches = plt.hist(np.log(D2_min), 50, facecolor='g')
            # use now a color scheme not unlike the ones above.
            # This has a rather odd distribution ...
            D2_scale=np.mean(np.log(D2_min))-np.amin(np.log(D2_min))
            radscale=np.mean(self.conf.rad)
            Fcolor,Fmap=self.color_init(3.0)
            for k in range(self.N):
                duse=np.log(D2_min[k])-np.amin(np.log(D2_min))
                if (duse>3.0*D2_scale):
                    duse=3.0*D2_scale
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=0.5*radscale*duse/D2_scale,ec=Fcolor(duse/D2_scale),fc=Fcolor(duse/D2_scale), linewidth=2)
                axval.add_patch(cir1)
        
        for k in range(len(self.pebbles.Ifull)):
            # this version depends on pebble numbering, so use 2nd version
            x0,x1,y0,y1=self.conf.getConPos2(self.pebbles.Ifull[k],self.pebbles.Jfull[k])
            if (self.pebbles.cluster[k]!=-1):
                if (k>=self.ncon):
                    ang=np.arctan((y1-y0)/(x1-x0))
                    conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=2,color=(0,0,1))
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=2,color=(0,0,1))
                axval.add_line(conlink)
            else:
                if (k>=self.ncon):
                    ang=np.arctan((y1-y0)/(x1-x0))
                    conlink=lne.Line2D([x0-self.small*np.sin(ang),x1-self.small*np.sin(ang)],[y0+self.small*np.cos(ang),y1+self.small*np.cos(ang)],linewidth=2,color=(0.65,0.65,0.65))
                else:
                    conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=2,color=(0.65,0.65,0.65))
                axval.add_line(conlink)

        axval.axis('scaled')
        #axval.axis('off')
        if self.conf.datatype=='simulation':
            axval.set_xlim(-self.Lx/2,self.Lx/2)
            axval.set_ylim(-self.Ly/2,self.Ly/2)
        elif self.conf.datatype=='experiment':
            xmin=np.amin(self.conf.x)-20
            ymin=np.amin(self.conf.y)-20
            axval.set_xlim(xmin,xmin+self.Lx+120)
            axval.set_ylim(ymin,ymin+self.Ly+120)
        #plt.colorbar(cpick,label="deformation field")
        plt.title(' D2 (log-scale)')
        return fig
    
    #Function for plotting the Maxwell-Cremona tiling
    def tileplotter(self, colorscheme, filled):
        #Create figure
        fig = plt.figure()
        axval = fig.add_subplot(1, 1, 1)    
        
        #Generate these colorschemes before the loop
        if colorscheme == 'random':
            colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(self.tiles))]
        
        elif colorscheme == 'colorblind':
            colors = ['#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB'] #Tol's color blind friendly color scheme
        
        #Only allow filled tiling for random and colorblind colorscheme
        if (colorscheme =='force' or colorscheme =='cluster') and filled:
            sys.exit('Only random and colorblind colorschemes are possible for filled tiling, exiting program.')
            
        #Loop over all tiles
        for n in range(len(self.tiles)):
            #Retrieve data
            data = self.tiles[n]
            i = data[0]
            vertices = data[1]
            
            #Plot filled tiles
            if filled:
                #Remove contact point and return only vertex coordinates
                for vertex in vertices:
                    del vertex[0]
                axval.add_patch(ptch.Polygon(vertices, facecolor=colors[n%len(colors)], edgecolor=None, alpha=0.7))
            
            #Else plot edges
            if not filled:
                for k in range(len(vertices)-1):
                    #Extract vertices
                    vertex1 = vertices[k]
                    vertex2 = vertices[k+1]
                    
                    #Extract contact and find corresponding index
                    j = vertex2[0]
                    arg = np.argwhere((self.I == i) & (self.J == j)).flatten()
                    if len(arg) == 0:
                        arg = np.argwhere((self.I == j) & (self.J == i)).flatten()
                    try:
                        idx = arg[0]
                    except:
                        idx = 0 # give color of cluster 0 if all fails
                    
                    #Define colors and plot
                    if colorscheme == 'cluster':
                        #colors = self.cluster_colors
                        label = int(self.pebbles.cluster[idx])
                        #color=(self.rgbclus[k,0],self.rgbclus[k,1],self.rgbclus[k,2])
                        if label != -1:
                            conlink = lne.Line2D([vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color=(self.rgbclus[label,0],self.rgbclus[label,1],self.rgbclus[label,2]), marker='o')
                            axval.add_line(conlink)
                        else:
                            conlink = lne.Line2D([vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color=(0.65, 0.65, 0.65), marker='o')
                            axval.add_line(conlink)
                    
                    elif colorscheme == 'force':
                        fscale = self.fgiven
                        Fcolor, Fmap = self.color_init(np.sqrt(np.amax(self.fn))/fscale) #Use fn or ftot?
                        fval = np.sqrt(self.fn[idx]/fscale)
                        conlink = lne.Line2D([vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color=Fcolor(fval), lw=2*fval, marker='.')
                        axval.add_line(conlink)
                        
                    elif colorscheme == 'random':
                        conlink = lne.Line2D([vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color=colors[k], marker='o')
                        axval.add_line(conlink)
                        
                    elif colorscheme == 'colorblind':
                        conlink = lne.Line2D([vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color=colors[k%len(colors)], marker='o')
                        axval.add_line(conlink)
                    
                    else: 
                        sys.exit('Not a valid color scheme, exiting program.')
                    
                
        plt.title('Maxwell-Cremona Tiling')
        plt.xlabel(r'$F_x\ [N]$')
        plt.ylabel(r'$F_y\ [N]$')
        plt.axis('equal')
        return fig

    #Function that return some statistics of the tiling
    def tiling_statistics(self):
        err = []
        f = []
        fvertices = []
        for n in range(len(self.tiles)):
            #Retrieve data
            data = self.tiles[n]
            vertices = data[1]
  
            #First and last vertex of a tile
            vertex_start = vertices[0]
            vertex_end = vertices[-1]
            
            #Coordintates of first and last vertex
            vertex_start_x = vertex_start[1]
            vertex_start_y = vertex_start[2]
            vertex_end_x = vertex_end[1]
            vertex_end_y = vertex_end[2]

            #Compute force imbalance
            errx = np.abs(vertex_start_x - vertex_end_x)
            erry = np.abs(vertex_start_y - vertex_end_y)
            err.append(np.sqrt(np.square(errx) + np.square(erry)))
            
            #Set force on a vertex equal to zero
            fver = 0
            
            #Loop over all vertices of a tile
            for k in range(len(vertices)-1):
                #Extract vertices
                vertex1 = vertices[k]
                vertex2 = vertices[k+1]
                fx = vertex1[1] - vertex2[1]
                fy = vertex1[2] - vertex2[2]
                #Compute sum of all magnitudes of forces
                fmag = np.sqrt(np.square(fx) + np.square(fy))
                f.append(fmag)
                fver += fmag
            
            fvertices.append(fver)
                
        #Average over all tiles
        err_mean = np.mean(err)
        #Average over all contacts
        f_mean = np.mean(f)
        #Compute off-force-balance ratio
        err_force_ratio = err_mean/f_mean       
        
        #Print and return
        print("mean force imbalance = ", err_mean)
        print("mean force magnitude = ", f_mean)
        print("mean off-force balance ratio = ", err_force_ratio)
        return err, fvertices, err_mean, f_mean, err_force_ratio

        




