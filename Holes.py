#!/usr/bin/python

# Silke Henkes, 11.10.19: 
# Created Holes clase to produce half-edge tesselation of the rigid clusters
# Detangled code parts, contains: 
# - Contacts only in rigid clusters
# - MakeLoop to create half-edge tesselation
# - Cluster and stress cross-correlations
# - plotHoles to visualise tesselation

import random
import sys, os, glob

# Each analysis has a Configuration, a Pebbles and a Hessian 
# (of which Pebbles and Hessian are allowed to be empty shells if you're not using these functions)
from Configuration import *
from Pebbles import *

# Note: changed to simple pickle on move to python 3. Compatibility problems with loading old pickle files expected!
#import cPickle as pickle
import pickle as pickle
import copy as cp
import numpy as np

# This does need all the plotting libraries
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.lines as lne
#from matplotlib.colors import LinearSegmentedColormap

#import matplotlib
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 4
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=16.0
matplotlib.rcParams['legend.fontsize']=14.0


class Holes:
    
        ### =========== Analysis constructor: Give it a configuration (mandatory), and Pebbles and Hessian, which are allowed to be empty shells
        def __init__(self,conf0,pebbles0,verbose0=False):
                self.conf=conf0
                self.pebbles=pebbles0
                self.verbose=verbose0
                # basic dimensions ...
                self.N=self.conf.N
                self.Lx=self.conf.Lx
                self.Ly=self.conf.Ly
                self.conf.y[-4]+=500
                self.conf.y[-3]-=500
                self.conf.x[-2]-=500
                self.conf.x[-1]+=500
        
        # Make a half-edge structure on a given connectivity set
        # The options here are either all of the connectivity, the rigid clusters, or the rigid regions
        def makeHalfEdge(self):
            
            # We need to loop over all the bonds and identify if they are part of a rigid cluster
            # if yes, keep them, if not, throw them out
            # we need to also throw out the second bonds in case the first one has already been registered
            # recall they are in the region ncon < k < ncon2
            Rigbonds = []
            # There can only be one bond between two particles
            # and the stupid
            ifull = np.array(self.pebbles.Ifull)
            jfull = np.array(self.pebbles.Jfull)
            idone=[]
            jdone=[]
            for k in range(self.pebbles.ncon2):
                if self.pebbles.cluster[k]>=0:
                    i = self.pebbles.Ifull[k]
                    j = self.pebbles.Jfull[k]
                    #print "working on bond " + str(k) + " between particles " + str(i) + " and " + str(j)
                    # Unfortunately, a bit messy to look for the pair
                    isi = [index for index,value in enumerate(idone) if value==i]
                    isj = [index for index,value in enumerate(jdone) if value==i]
                    seen=False
                    if len(isi) > 0:
                        hmm=np.array(jdone)
                        if j in hmm[isi]:
                            if self.verbose:
                                print "This contact " + str(k) + " has been seen already, not adding it!"
                            seen = True
                    if len(isj)>0 :
                        hmm=np.array(idone)
                        if j in hmm[isj]:
                            if self.verbose:
                                print "This contact " + str(k) + " has been seen already, not adding it!"
                            seen = True
                    #print seen
                    if not seen:
                        Rigbonds.append(k)
                        idone.append(i)
                        jdone.append(j)
                        
            # Now that we have them, extract their connectivity
            # and because stupid
            self.Ibonds = ifull[Rigbonds]
            self.Jbonds = jfull[Rigbonds]
            
            # this is juct still bonds that are single-counted!
            self.doneR = [ False for k in range(len(Rigbonds))]
            self.doneL = [ False for k in range(len(Rigbonds))]
            
            # Here is our actual halfedge structure. Start the counting fresh ...
            # keep the vertex (particle) labels, but nothing else
            # holes will be our faces
            self.I = []
            self.J = []
            self.Next = []
            self.Prev = []
            self.zFace = []
            self.nbond = 0
            self.nface = 0
            self.facecon = []
            for k in range(len(Rigbonds)):
                if self.verbose:
                    print "Working on bond " + str(k) + " linking particles " + str(self.Ibonds[k]) + " and " + str(self.Jbonds[k])
                if (not self.doneR[k]):
                    # starting a new face
                    # stick the current bond into the structure, labeled by the current bond label
                    self.I.append(self.Ibonds[k])
                    self.J.append(self.Jbonds[k])
                    self.facecon.append(self.nbond)
                    self.doneR[k]=True # this one is cooked
                    # can't know previous yet, keep the label 
                    self.Prev.append(-1)
                    bini = self.nbond
                    self.MakeLoop(k,1,bini)
                    
                if (not self.doneL[k]):
                    # stick the current bond into the structure, labeled by the current bond label
                    self.I.append(self.Jbonds[k])
                    self.J.append(self.Ibonds[k])
                    self.facecon.append(self.nbond)
                    self.doneL[k]=True # this one is cooked
                    # can't know previous yet, keep the label 
                    self.Prev.append(-1)
                    bini = self.nbond
                    # next counterclockwise bond
                    self.MakeLoop(k,-1,bini)
                    
            # and that should be it!
            print self.zFace
                
        def MakeLoop(self,k,orient,bini):
            # next counterclockwise bond
            z=1
            ccn,lr=self.ccnext(k,orient)
            # going round the loop now
            while (not ccn==k):
                # this new bond's previous is the current counter
                self.Prev.append(self.nbond)
                self.nbond+=1
                # and the previous bond's next is the incremented counter
                self.Next.append(self.nbond)
                # stick into the structure in the correct orientation
                if lr==1:
                    self.I.append(self.Ibonds[ccn])
                    self.J.append(self.Jbonds[ccn])
                    self.doneR[ccn]=True
                else:
                    self.I.append(self.Jbonds[ccn])
                    self.J.append(self.Ibonds[ccn])
                    self.doneL[ccn]=True
                # we're still on the same face
                #self.Face.append(self.nface)
                z+=1
                # now keep going
                ccn,lr=self.ccnext(ccn,lr)
                #debug
                #if z>2000:
                    #ccn=k
            self.Prev[bini]=self.nbond
            self.Next.append(bini)
            self.zFace.append(z)
            self.nbond+=1
            if self.verbose:
                print "Done with face " + str(self.nface) + " which has " + str(z) + " edges."
            self.nface+=1
        
        # computes the counter-clockwise next bond of an oriented bond, and return its own orientation   
        # orientation is either 1 or -1, which makes me go look for the next bond via either I or J
        def ccnext(self,k,orient ):
            if orient == 1:
                particle = self.Jbonds[k]
                # angle of the original bond (signed), and between -pi to pi
                x0,x1,y0,y1=self.conf.getConPos3(self.Ibonds[k],particle)
            else:
                particle = self.Ibonds[k]
                # angle of the original bond (signed), and between -pi to pi
                x0,x1,y0,y1=self.conf.getConPos3(self.Jbonds[k],particle)
            ang0 = np.arctan2(y1-y0,x1-x0)
            if self.verbose:
                print "Looking for turning angle"
                print "Locating neighbours of particle " + str(particle)
                print "Initial angle: " + str(ang0)
            # labels of the next bonds
            # This lacks elegance.
            nnei1 = [index for index,value in enumerate(self.Ibonds) if value==particle]
            nnei2 = [index for index,value in enumerate(self.Jbonds) if value==particle]
            # list addition = concatenation
            nnei = nnei1+nnei2
            if self.verbose:
                print "My neighbour bonds are: "
                print nnei
            # throw out bond itself
            # But careful: Only when there is more than one! Otherwise allow it to double back ...
            if len(nnei)>1:
                nnei.remove(k)
            # get the angle of the oriented bond
            angmax=-np.pi-0.0001
            ccn = -1
            lr = 1
            for l in nnei:
                lr0=1
                if self.Ibonds[l]==particle:
                    x0,x1,y0,y1=self.conf.getConPos3(particle,self.Jbonds[l])
                    lr0 = 1
                elif self.Jbonds[l]==particle:
                    x0,x1,y0,y1=self.conf.getConPos3(particle,self.Ibonds[l])
                    lr0 = -1
                else:
                    print "Problem!!"
                ang = np.arctan2(y1-y0,x1-x0)
                dang = ang-ang0
                if (dang < -np.pi):
                    dang += 2*np.pi 
                if (dang > np.pi):
                    dang -= 2*np.pi
                if (dang > angmax): # there appear to be stray double bonds ...
                    ccn = l
                    angmax = dang
                    lr = lr0
            if self.verbose:
                print "Chose contact " + str(ccn) + " originally pointing from " + str(self.Ibonds[ccn]) + " to " + str(self.Jbonds[ccn])
                print "Turning angle " + str(angmax) + " and lr = " + str(lr)
            return ccn, lr
        
        # Analyse the holes: Which ones are outer face? Which ones contain particles? What are their areas / perimeters?
        def sortHoles(self):
            ha0 = []
            hp0 = []
            zp0 = []
            self.conf.y[-4]-=500
            self.conf.y[-3]+=500
            self.conf.x[-2]+=500
            self.conf.x[-1]-=500
            # labels of outer faces 
            self.outer = []
            areaouter = 0.0
            perimouter = 0.0
            zouter = 0
            self.isHole = [ False for k in range(self.nface)]
            areascale = self.Lx*self.Ly / self.N
            for f in range(self.nface):
                # starting bond
                bstart = self.facecon[f]
                bend=bstart
                done=False
                c=0
                # get coordinates of polygon
                polyx=[]
                polyy=[]
                while not done:
                    x0,x1,y0,y1=self.conf.getConPos3(self.I[bend],self.J[bend])
                    polyx.append(x0)
                    polyy.append(y0)
                    bend = self.Next[bend]
                    if bend ==bstart:
                        done = True
                    if c>2000:
                        done = True
                    c+=1
                parea = self.polyArea(np.array(polyx),np.array(polyy))
                pperim = self.polyPerim(np.array(polyx),np.array(polyy))
                if parea <0:
                    self.outer.append(f)
                    areaouter -= parea
                    perimouter += pperim
                    zouter +=self.zFace[f]
                else:
                    # Empirically: 5 and below cannot enclose a particle
                    if self.zFace[f]>5 and parea/areascale>2:
                        ha0.append(parea)
                        hp0.append(pperim)
                        zp0.append(self.zFace[f])
                        self.isHole[f]=True
                    # elif (self.zFace[f]>6) and (pperim/parea**0.5 < 4.5):
                    #     ha0.append(parea)
                    #     hp0.append(pperim)
                    #     zp0.append(self.zFace[f])
                    #     self.isHole[f]=True
                    # elif (self.zFace[f]>5) and (pperim/parea**0.5 < 4.0):
                    #     ha0.append(parea)
                    #     hp0.append(pperim)
                    #     zp0.append(self.zFace[f])
                    #     self.isHole[f]=True
            # Let's compute the remaining, uncovered area
            ha0.append(self.Lx*self.Ly-areaouter)
            hp0.append(perimouter)
            zp0.append(zouter)
            # Compute our output data now
            # mean area per particle
            self.HoleZ=np.array(zp0)
            self.HoleArea=np.array(ha0)/areascale
            self.HolePerim=np.array(hp0)/np.sqrt(areascale)
            self.nHole=len(zp0)
            # Basic statistics
            # plt.figure()
            # plt.plot(self.HoleArea,self.HolePerim/np.sqrt(self.HoleArea),'o')
            # plt.xlabel('Hole area (in particle units)')
            # plt.ylabel('Shape factor (3.54 for a circle)')
            
            return self.nHole,self.HoleZ,self.HolePerim,self.HoleArea
            
        def polyArea(self,x,y):
            # will return positive number for a ccw non-intersecting polygon
            return -0.5*(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
        def polyPerim(self,x,y):
            return np.sum(np.sqrt((x-np.roll(x,1))**2+(y-np.roll(y,1))**2))
        
        def plotFaces(self,**kwargs):
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
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.5, 0.5, 0.5),fc='none', linewidth=2)
                axval.add_patch(cir1)
                if self.verbose:
                    axval.text(self.conf.x[k],self.conf.y[k],str(k),fontsize=8)
                
                            
            # now plotting the actual holes
            # as polygons
            rgbpattern=np.random.rand(self.nface+1,3)
            for f in range(self.nface):
                # starting bond
                bstart = self.facecon[f]
                bend=bstart
                done=False
                c=0
                polyx=[]
                polyy=[]
                while not done:
                    try:
                        x0,x1,y0,y1=self.conf.getConPos3(self.I[bend],self.J[bend])
                        polyx.append(x0)
                        polyy.append(y0)
                        conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=1,color=(rgbpattern[f,0],rgbpattern[f,1],rgbpattern[f,2]))
                        axval.add_line(conlink)
                        bend = self.Next[bend]
                        if bend ==bstart:
                            done = True
                        #if c>2000:
                            #done = True
                        c+=1
                    except:
                        done = True
                newface = ptch.Polygon(np.transpose(np.array([polyx, polyy])),facecolor=(rgbpattern[f,0],rgbpattern[f,1],rgbpattern[f,2]),alpha=0.5)
                axval.add_patch(newface)
                    
            axval.axis('scaled')
            if self.conf.datatype=='simulation':
                    axval.set_xlim(-self.Lx/2,self.Lx/2)
                    axval.set_ylim(-self.Ly/2,self.Ly/2)
            elif self.conf.datatype=='experiment':
                    xmin=np.amin(self.conf.x)-20
                    ymin=np.amin(self.conf.y)-20
                    axval.set_xlim(xmin,xmin+self.Lx+120)
                    axval.set_ylim(ymin,ymin+self.Ly+120)
            axval.set_title('System with ' + str(self.nface) + " faces")
            # returning figure pointer for use outside (like plotting or saving)
            if 'figure' in kwargs:
                return fig
            elif 'axis' in kwargs:
                return axval
            else:
                return fig
            
        def plotHoles(self,**kwargs):
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
            plt.xticks([], [])
            plt.yticks([], [])
            for k in range(self.N):
                t=0
                # Particle circles
                cir1=ptch.Circle((self.conf.x[k],self.conf.y[k]),radius=self.conf.rad[k],ec=(0.65, 0.65, 0.65),fc='none', linewidth=2)
                axval.add_patch(cir1)
                if self.verbose:
                    axval.text(self.conf.x[k],self.conf.y[k],str(k),fontsize=8)
                
            for k in range(len(self.Ibonds)):
                x0,x1,y0,y1=self.conf.getConPos2(self.Ibonds[k],self.Jbonds[k])
                conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=1.5,color=(0,0,0))
                axval.add_line(conlink)
                if self.verbose:
                    axval.text(x0+0.5*(x1-x0),y0+0.5*(y1-y0),str(k),fontsize=8)
                            
            # now plotting the actual holes
            # as polygons
            rgbpattern=np.random.rand(self.nface+1,3)
            for f in range(self.nface):
                if self.isHole[f]:
                    # starting bond
                    bstart = self.facecon[f]
                    bend=bstart
                    done=False
                    c=0
                    polyx=[]
                    polyy=[]
                    while not done:
                        x0,x1,y0,y1=self.conf.getConPos2(self.I[bend],self.J[bend])
                        polyx.append(x0)
                        polyy.append(y0)
                        #conlink=lne.Line2D([x0,x1],[y0,y1],linewidth=3,color=(rgbpattern[f,0],rgbpattern[f,1],rgbpattern[f,2]))
                        #axval.add_line(conlink)
                        bend = self.Next[bend]
                        if bend ==bstart:
                            done = True
                        #if c>2000:
                            #done = True
                        c+=1
                    # if(len(polyx)<40):
                    newface = ptch.Polygon(np.transpose(np.array([polyx, polyy])),edgecolor=(rgbpattern[f,0],rgbpattern[f,1],rgbpattern[f,2]), linewidth=1, facecolor=(rgbpattern[f,0],rgbpattern[f,1],rgbpattern[f,2]),alpha=0.5)
                    axval.add_patch(newface)
                    
            axval.axis('scaled')
            if self.conf.datatype=='simulation':
                    axval.set_xlim(-self.Lx/2,self.Lx/2)
                    axval.set_ylim(-self.Ly/2,self.Ly/2)
            elif self.conf.datatype=='experiment':
                    xmin=np.amin(self.conf.x)#-20
                    ymin=np.amin(self.conf.y)
                    axval.set_xlim(xmin,xmin+self.Lx+80)
                    axval.set_ylim(ymin,ymin+self.Ly+80)
            axval.spines['bottom'].set_color('b')
            axval.spines['top'].set_color('r')
            axval.spines['left'].set_color('r')
            axval.spines['right'].set_color('b')
            for axis in ['top','bottom','left','right']:
                axval.spines[axis].set_linewidth(6)
            # axval.set_title('System with ' + str(self.nHole) + " potential holes")
            # returning figure pointer for use outside (like plotting or saving)
            if 'figure' in kwargs:
                return fig
            elif 'axis' in kwargs:
                return axval
            else:
                return fig
                        
                    
            
            
            
        
        
        
