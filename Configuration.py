# Silke Henkes 11.07.2013: Accelerated version of pebble game and rigid cluster code
# 
# # Silke Henkes, 29.08.18: 
# Created Configuration class to read both experimental and simulation data
# Detangled code parts, contains: 
# - Read-in for all circumstances
# - Analysis functions which depend on particle positions only, called by Analysis and others

# Silke Henkes, 29.10.23
# Refactored Configuration class to deal with 4 (and potentially more) setups: simulation, square experiment, annulus and lattice
# These are now derived classes of a base Configuration class, which allows for a much less confusing / entagled setup
# There is now a unified interface for initialising a configuration: Constructor which takes an (arbitrary) param dictionary, and a ReadData which is pointed to the right files via
# posfile, confile, and step information where appropriate
# This unfortunately breaks backwards compatibility with existing runscripts, but it was getting unavoidable 

#!/usr/bin/python

import sys, os, glob
import numpy as np

# Global list of exisiting geometries
# Include in script, not here
#setups={'simulation':ConfigurationSimulation,'lattice':ConfigurationLattice,'experiment_square':ConfigurationExpSquare,'experiment_annulus':ConfigurationExpAnnulus}


class Configuration:
        
    def __init__(self,param0):
        self.params=param0
        self.folder = self.params["folder"]
        self.experiment = self.params['experiment']
        self.periodic = self.params['periodic']
        self.periodicx = self.params['periodicx']
        self.periodicy = self.params['periodicy']
        self.addBoundarySquare = False
        self.addBoundaryAnnulus = False
        
    def ReadData(self,posfile,confile,step,verbose=False):
        if self.experiment:
            self.ReadDataExperiment(posfile,confile,step,verbose)    
        else:
            if verbose:
                print("Somehow defaulted to Configuration base class. Doing nothing.")

    def AddBoundaryContacts(self,verbose=False):
        if verbose:
            print("Defaulted to Configuration base class. Doing nothing.")


    # Read in data from posfile and confile, which is full file names set by the running script
    def ReadDataExperiment(self, posfile,confile,step,verbose=False):    
        fname = self.folder+posfile
        print(fname)
        self.isPosdata=True
        coords=np.loadtxt(self.folder+posfile, delimiter=',')
        # fudging placeholder data, delete and revert
        #coords=coords[coords[:,0] == step]
        self.id=coords[:,0]
        self.x=coords[:,1]
        self.y=coords[:,2]
        self.rad=coords[:,3]
        self.boundary=coords[:,4]
        self.N=len(self.rad)
        print("Experiment with " +str(self.N)+ " particles")
        self.Lx=np.amax(self.x)-np.amin(self.x)
        self.Ly=np.amax(self.y)-np.amin(self.y)
        del coords
        # try:
        #     coords=np.loadtxt(self.folder+posfile, delimiter=',')
        #     coords=coords[coords[:,0] == step]
        #     self.id=coords[:,1]
        #     self.x=coords[:,2]
        #     self.y=coords[:,3]
        #     self.rad=coords[:,4]
        #     self.boundary=coords[:,5]
        #     self.N=len(self.rad)
        #     print("Experiment with " +str(self.N)+ " particles")
        #     self.Lx=np.amax(self.x)-np.amin(self.x)
        #     self.Ly=np.amax(self.y)-np.amin(self.y)
        #     del coords
        # except:
        #     self.isPosdata=False
        #     self.x=0
        #     self.y=0
        #     self.rad=1 # purely so that we can divide by it ...
        #     self.N=0
        #     print("Error: there is no position data here")
        #Load in contact data
        #confile = self.folder +'/Adjacency_list.txt'
        self.isCondata=True
        fname = self.folder+confile
        print(fname)
        try:
            condata=np.loadtxt(self.folder+confile, delimiter=',')
            condata=condata[condata[:,0] == step]
        except:
            self.isCondata=False
            print ("Error: there is no contact data here")
            return 1
        print(condata)
        
        #Create empty lists
        #Lists with particle ids
        self.I=[]
        self.J=[]
        
        #Lists with forces
        fn0=[]
        ft0=[]
        
        #Frictional or sliding list
        fm0=[]
        
        #Drop duplicates
        for k in range(len(condata[:,0])):
            if condata[:,1][k] > condata[:,2][k]:
                #Add contact id's
                i = condata[:,1][k].astype(int)
                j = condata[:,2][k].astype(int)
                argi = np.argwhere(self.id == i).flatten()[0]
                argj = np.argwhere(self.id == j).flatten()[0]
                                    
                #For now we identify particle id's with indices in the list. Not ideal for debugging, but it works.                   
                self.I.append(argi)
                self.J.append(argj)
                
                #Extract forces
                fn = condata[:,4][k]
                ft = condata[:,3][k]
                fn0.append(fn)
                ft0.append(ft)
                
                #Determine frictional or sliding
                if (abs(ft)/fn>self.mu):
                    fm0.append(1)
                else:
                    fm0.append(0)
        
        #Final arrays
        self.ncon=len(fm0)
        self.I = np.array(self.I)
        self.J = np.array(self.J)
        self.fnor=np.array(fn0)
        self.ftan=np.array(ft0)
        self.fullmobi=np.array(fm0)
        
        self.nx=np.zeros(self.ncon)
        self.ny=np.zeros(self.ncon)
        for k in range(self.ncon):
            x1=self.x[self.I[k]]
            y1=self.y[self.I[k]]
            x2=self.x[self.J[k]]
            y2=self.y[self.J[k]]
            rij=np.sqrt((x1-x2)**2+(y1-y2)**2)
            self.nx[k]=(x2-x1)/rij
            self.ny[k]=(y2-y1)/rij

        print("Config frame #" +str(step)+ " created")      
        return 0
    #### ======================== Analysis helper functions	 =================================================
    # computes basic contact, force, torque, and stress statistics
    def getStressStat(self):
        fsumx=np.empty((self.N,))
        fsumy=np.empty((self.N,))
        torsum=np.empty((self.N,))
        #------- those useful for plotting ----
        self.prepart=np.empty((self.N,))
        self.sxxpart=np.empty((self.N,))
        self.syypart=np.empty((self.N,))
        self.sxypart=np.empty((self.N,))
        self.syxpart=np.empty((self.N,))
        #-------------------------------------
        for u in range(self.N):
            # problem - self.I doesn't seem to behave as an integer ...
            # apparently this only works with arrays, not with lists??
            c1=np.nonzero(np.array(self.I)==u)
            c2=np.nonzero(np.array(self.J)==u)
    
            fsumx[u]=np.sum(-self.fnor[c1[0]]*self.nx[c1[0]])+np.sum(self.fnor[c2[0]]*self.nx[c2[0]])+np.sum(self.ftan[c1[0]]*(-self.ny[c1[0]]))+np.sum(-self.ftan[c2[0]]*(-self.ny[c2[0]]))
            fsumy[u]=np.sum(-self.fnor[c1[0]]*self.ny[c1[0]])+np.sum(self.fnor[c2[0]]*self.ny[c2[0]])+np.sum(self.ftan[c1[0]]*self.nx[c1[0]])+sum(-self.ftan[c2[0]]*self.nx[c2[0]])
            torsum[u]=self.rad[u]*(np.sum(self.ftan[c1[0]])+np.sum(self.ftan[c2[0]]))
            self.prepart[u]=np.sum(self.fnor[c1[0]]*self.rad[u])+np.sum(self.fnor[c2[0]]*self.rad[u])
            self.sxxpart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.nx[c1[0]]*self.nx[c1[0]])+np.sum(self.fnor[c2[0]]*self.nx[c2[0]]*(-self.nx[c2[0]]))+np.sum(self.ftan[c1[0]]*(-self.ny[c1[0]])*(-self.ny[c1[0]]))+np.sum(-self.ftan[c2[0]]*(-self.ny[c2[0]])*(self.ny[c2[0]])))
            self.syypart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.ny[c1[0]]*(self.ny[c1[0]]))+np.sum(self.fnor[c2[0]]*self.ny[c2[0]]*(-self.ny[c2[0]]))+np.sum(self.ftan[c1[0]]*(self.nx[c1[0]])*self.nx[c1[0]])+sum(-self.ftan[c2[0]]*self.nx[c2[0]]*(-self.nx[c2[0]])))
            self.sxypart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.nx[c1[0]]*(-self.ny[c1[0]]))+np.sum(self.fnor[c2[0]]*self.nx[c2[0]]*(self.ny[c2[0]]))+np.sum(self.ftan[c1[0]]*(-self.ny[c1[0]])*(self.nx[c1[0]]))+np.sum(-self.ftan[c2[0]]*(-self.ny[c2[0]])*(-self.nx[c2[0]])))
            self.syxpart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.ny[c1[0]]*(self.nx[c1[0]]))+np.sum(self.fnor[c2[0]]*self.ny[c2[0]]*(-self.nx[c2[0]]))+np.sum(self.ftan[c1[0]]*(self.nx[c1[0]])*(-self.ny[c1[0]]))+sum(-self.ftan[c2[0]]*self.nx[c2[0]]*(self.ny[c2[0]])))
            del c1 
            del c2
        
        # statistical output: don't do histograms here just yet; the problem of consistent and appropriate binning is too messy. 
        # return the whole thing.
        # except for the mobilization distribution which is unique
        # just get a few signature parameters here
        zav=2*len(self.I)/(1.0*self.N)
        nm=len(np.nonzero(self.fullmobi>0)[0])/(1.0*self.N)
        pres=np.sum(self.prepart)/(self.Lx*self.Ly)
        sxx=np.sum(self.sxxpart)/(self.Lx*self.Ly)
        syy=np.sum(self.syypart)/(self.Lx*self.Ly)
        sxy=np.sum(self.sxypart)/(self.Lx*self.Ly)
        syx=np.sum(self.syxpart)/(self.Lx*self.Ly)
        fxbal=np.mean(abs(fsumx))/np.mean(self.fnor)
        fybal=np.mean(abs(fsumy))/np.mean(self.fnor)
        torbal=np.mean(abs(torsum))/(np.mean(self.fnor)*np.mean(self.rad)) # correct units; do *not* normalize with ftan; consider a mu though
        mobin=np.linspace(-1,1,101)
        mohist,bin_edges=np.histogram(self.ftan/(self.mu*self.fnor),mobin)
        
        return zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx	
    
    ##================ Finally, pure helper functions: positions of both ends of a contact ===============
    # get the cleaned up, periodic boundary conditions sorted out positions corresponding to two ends of a contact. 
    # Basic plotting helper function
    def getConPos(self,k,nobound=False):
        x0=self.x[self.I[k]]
        x1=self.x[self.J[k]]
        y0=self.y[self.I[k]]
        y1=self.y[self.J[k]]
        if self.periodic:
                x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
                yover=np.round((y1-y0)/self.Ly)
                if (yover!=0):
                        y1=y1-self.Ly*yover
                        x1-=self.Lx*self.strain*yover
                        x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
        if self.periodicx:
            x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
        if self.periodicy:
            y1=y1-self.Ly*np.round((y1-y0)/self.Ly)
        if not nobound:
            if self.addBoundarySquare:
                ival=self.I[k]
                if ((ival==self.bindices[0]) or (ival==self.bindices[1])): #top or bottom
                    x0=x1
                if ((ival==self.bindices[2]) or (ival==self.bindices[3])): #left or right
                    y0=y1
            if self.addBoundaryAnnulus:
                ival=self.I[k]
                l = 100 #Length in pixels of the contacts with the boundary
                #If in contact with the inner boundary
                if (ival==self.bindices[0]):
                    #Compute position vector from midpoint
                    r = np.array([x1-self.mid_x, y1-self.mid_y])
                    #Compute length of r
                    r_len = np.linalg.norm(r)
                    #Compute angle between position vector from midpoint and (1,0) vector
                    ang = np.arctan2(r[1], r[0])
                    #Compute coordinates of outward inward coordinates
                    x0 = self.mid_x + (r_len - l)*np.cos(ang)
                    y0 = self.mid_y + (r_len - l)*np.sin(ang)
                #If in contact with the outer boundary 
                if (ival==self.bindices[1]):
                    #Compute position vector from midpoint
                    r = np.array([x1-self.mid_x, y1-self.mid_y])
                    #Compute length of r
                    r_len = np.linalg.norm(r)
                    #Compute angle between position vector from midpoint and (1,0) vector
                    ang = np.arctan2(r[1], r[0])
                    #Compute coordinates of outward pointing coordinates
                    x0 = self.mid_x + (r_len + l)*np.cos(ang)
                    y0 = self.mid_y + (r_len + l)*np.sin(ang)
        return x0,x1,y0,y1
    
    # same, but based on existing particle labels (in case those come from elsewhere)
    def getConPos2(self,k1,k2,nobound=False):
        x0=self.x[k1]
        x1=self.x[k2]
        y0=self.y[k1]
        y1=self.y[k2]
        if self.periodic:
                x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
                yover=np.round((y1-y0)/self.Ly)
                if (yover!=0):
                        y1=y1-self.Ly*yover
                        x1-=self.Lx*self.strain*yover
                        x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
        if self.periodicx:
            x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
        if self.periodicy:
            y1=y1-self.Ly*np.round((y1-y0)/self.Ly)
        if not nobound:
            if self.addBoundarySquare:
                if ((k1==self.bindices[0]) or (k1==self.bindices[1])): #top or bottom
                    x0=x1
                if ((k1==self.bindices[2]) or (k1==self.bindices[3])): #left or right
                    y0=y1
            if self.addBoundaryAnnulus:
                l = 100 #Length in pixels of the contacts with the boundary
                if (k1==self.bindices[0] and k2==self.bindices[1]) or (k1==self.bindices[1] and k2==self.bindices[0]):
                    return x0,x1,y0,y1
                #If in contact with the inner boundary
                if (k1==self.bindices[0]):
                    #Compute position vector from midpoint
                    r = np.array([x1-self.mid_x, y1-self.mid_y])
                    #Compute length of r
                    r_len = np.linalg.norm(r)
                    #Compute angle between position vector from midpoint and (1,0) vector
                    ang = np.arctan2(r[1], r[0])
                    #Compute coordinates of outward inward coordinates
                    x0 = self.mid_x + (r_len - l)*np.cos(ang)
                    y0 = self.mid_y + (r_len - l)*np.sin(ang)
                #If in contact with the outer boundary 
                if (k1==self.bindices[1]):
                    #Compute position vector from midpoint
                    r = np.array([x1-self.mid_x, y1-self.mid_y])
                    #Compute length of r
                    r_len = np.linalg.norm(r)
                    #Compute angle between position vector from midpoint and (1,0) vector
                    ang = np.arctan2(r[1], r[0])
                    #Compute coordinates of outward pointing coordinates
                    x0 = self.mid_x + (r_len + l)*np.cos(ang)
                    y0 = self.mid_y + (r_len + l)*np.sin(ang)
        return x0,x1,y0,y1

########## Experimental Square Configuration #################
class ConfigurationExpSquare(Configuration):

    def __init__(self,params):
        print("ConfigurationExpSquare: Created new Configuration Experimental Square")
        super(ConfigurationExpSquare,self).__init__(params)
  
        #self.setup = param["setup"] # redundant as needs to be same 'experiment_annulus' as in runscript
        self.setup = 'experiment_square'
        #print(self.param)
        #self.experiment = self.param.experiment
        #print(self.param.experiment)
        # These get set to true if add boudary functions are called
        # self.addBoundarySquare = True
        # debug - remove
        self.addBoundarySquare = False
        

        #Set friction coefficient
        self.mu = self.params["mu"]
    
        # Experimental data does not have periodic boundary conditions and there is no angle data either
        #self.periodic = False
        #self.hasAngles = False

        ### Geometric data
        # Get parameters for boundary distance from dictionary
        try:
            self.threshold = self.params["threshold"]
            self.Brad = self.params["Brad"]
        except:
            self.threshold=20,
            self.Brad=20.0

        #### Mechanical data
        # density of the material in kg/m^3
        self.density = self.params["density"]
        # Prefactor of the stiffness coefficient, k_n = \frac{\pi}{6} E h = 6490 kg /s^2
        self.stiffness = self.params["stiffness"]
        # radius conversion factor from pixel to m.
        self.rconversion = self.params["pixel_convert"]
        # height of the disks in m
        self.height = self.params["thickness"]
        # boundary width
        self.width = self.params["width"]

    #### ======================== Boundary integration =======================================================
    def AddBoundaryContactsSquare(self,):
        self.addBoundarySquare=True

        # Threshold to check if a particle is close enough to walls.
        upidx=np.argmax(self.y)
        downidx=np.argmin(self.y)
        leftidx=np.argmin(self.x)
        rightidx=np.argmax(self.x)
        
        # Boundary posiitons:
        # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
        up=self.y[upidx]
        yup = up+self.rad[upidx]
        down=self.y[downidx]
        ydown = down-self.rad[downidx]
        left=self.x[leftidx]
        xleft=left-self.rad[leftidx]
        right=self.x[rightidx]
        xright=right+self.rad[rightidx]
        
        # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
        Boundaries=np.zeros((4,3)) # Four boundary particles with their x,y and rad
        Boundaries[0,:]=[(left+right)*0.5,yup+self.Brad,self.Brad]
        Boundaries[1,:]=[(left+right)*0.5,ydown-self.Brad,self.Brad]
        Boundaries[2,:]=[xleft-self.Brad,(up+down)*0.5,self.Brad]
        Boundaries[3,:]=[xright+self.Brad,(up+down)*0.5,self.Brad]
        
        # Find the particles in contact with the boundary, and label correctly
        self.bindices=[self.N,self.N+1,self.N+2,self.N+3]
        padd=[]
        labels=[]
        pup =  np.nonzero(np.abs(self.y+self.rad-yup)<self.threshold)[0]
        padd.extend(pup)
        labels.extend([0 for k in range(len(pup))])
        pdown =  np.nonzero(np.abs(self.y-self.rad-ydown)<self.threshold)[0]
        padd.extend(pdown)
        labels.extend([1 for k in range(len(pdown))])
        pleft = np.nonzero(np.abs(self.x-self.rad-xleft)<self.threshold)[0]
        padd.extend(pleft)
        labels.extend([2 for k in range(len(pleft))])
        pright = np.nonzero(np.abs(self.x+self.rad-xright)<self.threshold)[0]
        padd.extend(pright)
        labels.extend([3 for k in range(len(pright))])
        
        fullmobi_add=[]
        fnor_add=[]
        ftan_add=[]
        nx_add=[]
        ny_add=[]
        for k in range(len(padd)):
            # does this guy have neighbours?
            neii=np.nonzero(self.I[:self.ncon]==padd[k])[0]
            neij=np.nonzero(self.J[:self.ncon]==padd[k])[0]
            # if yes add the boundary contacts
            if (len(neii)>0 or len(neij)>0):
                self.I.append(self.bindices[labels[k]])
                self.J.append(padd[k])
                if (labels[k])==0:
                    nx0=0
                    ny0=-1
                elif (labels[k]==1):
                    nx0=0
                    ny0=1
                elif (labels[k]==2):
                    nx0=1
                    ny0=0
                else:
                    nx0=-1
                    ny0=0
                # compute force on this contact by force balance
                # two minus signs on second part cancel out
                ftotx=np.sum(self.fnor[neii]*self.nx[neii]-self.ftan[neii]*self.ny[neii])-np.sum(self.fnor[neij]*self.nx[neij]-self.ftan[neij]*self.ny[neij])
                ftoty=np.sum(self.fnor[neii]*self.ny[neii]+self.ftan[neii]*self.nx[neii])-np.sum(self.fnor[neij]*self.ny[neij]+self.ftan[neij]*self.nx[neij])
                # (fx*nx+fy*ny)
                fnor0=ftotx*nx0+ftoty*ny0
                # (fx*(-ny)+fy*nx)
                ftan0=ftotx*(-ny0)+ftoty*nx0
                #print (ftan0)
                if (abs(ftan0)/fnor0>self.mu):
                    fullmobi_add.append(1)
                else:
                    fullmobi_add.append(0)
                fnor_add.append(fnor0)
                ftan_add.append(ftan0)
                nx_add.append(nx0)
                ny_add.append(ny0)
        # Finally stick it at the end of the existing data
        self.x=np.concatenate((self.x,Boundaries[:,0]))
        self.y=np.concatenate((self.y,Boundaries[:,1]))
        self.rad=np.concatenate((self.rad,Boundaries[:,2]))
        self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
        self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
        self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
        self.nx=np.concatenate((self.nx,np.array(nx_add)))
        self.ny=np.concatenate((self.ny,np.array(ny_add)))
        self.ncon=len(self.I)
        self.N+=4
        print ("Added boundaries!")         

########## Annulus Configuration #################
class ConfigurationExpAnnulus(Configuration):
    def __init__(self,params):
        print("ConfigurationExpAnnulus: Created new Configuration Experimental Annulus")
        super(ConfigurationExpAnnulus,self).__init__('experiment_annulus')
  
        #self.folder = param["folder"]
        #self.setup = param["setup"] # redundant as needs to be same 'experiment_annulus' as in runscript
        self.setup = 'experiment_annulus'
        # These get set to true if add boudary functions are called
        #self.addBoundarySquare = False
        #self.addBoundaryAnnulus = False
        #print ("Reading experimental data on annulus!")

        #Set friction coefficient
        self.mu = self.params["mu"]
    
        # Experimental data does not have periodic boundary conditions and there is no angle data either
        #self.periodic = False
        #self.hasAngles = False

        ### Geometric data
        # Radius R1 of the inner ring in experiment
        # Radius R2 of the outer ring in experiment
        try:
            self.R1 = self.params["R1"]
            self.R2 = self.params["R2"]
        except:
            self.R1 = 1548
            self.R2 = 2995

        #Coordinates of midpoint
        try:
            self.mid_x = self.params["mid_x"]
            self.mid_y = self.params["mid_y"]
        except:
            self.mid_x = 3134
            self.mid_y = 3028

        #Size of boundary texture
        try:
            self.texture = self.params["texture"]
        except:
            self.texture = 30
        #Strainstep of the experiment
        #self.step = strainstep

        #### Mechanical data
        try:
            # density of the material in kg/m^3
            self.density = self.params["density"]
            # Prefactor of the stiffness coefficient, k_n = \frac{\pi}{6} E h = 6490 kg /s^2
            self.stiffness = self.params["stiffness"]
            # radius conversion factor from pixel to m.
            self.rconversion = self.params["rconversion"]
            # height of the disks in m
            self.height = self.params["height"]
            # boundary width
            self.width = self.params["width"]
        
        except:
            self.density = 1060.0
            self.stiffness = 6490.0
            self.rconversion = 2.7e-4
            self.height = 3.1e-3
            self.width = 20e-3
    
     #### ======================== Boundary integration (Annulus) =======================================================
    def AddBoundaryContacts(self,verbose=False):
        # For getting positions
        self.addBoundaryAnnulus=True
        
        #Set radius of boundary particles
        #This needs some tweaking, since we are working with innies instead of outies. 
        Brad = self.texture
        Brad = 300
        
        # Boundary positions:
        # Coordinates of virtual boundary particles: at mid_x and y one Brad off from radi
        b1 = self.mid_y + self.R1 - Brad #Inner
        b2 = self.mid_y + self.R2 + Brad #Outer
                    
        # Coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
        Boundaries=np.zeros((2,3)) # Two boundary particles with their x,y and rad
        Boundaries[0,:]=[self.mid_x,b1,Brad] #Inner
        Boundaries[1,:]=[self.mid_x,b2,Brad] #Outer
        
        #Generate two new id's for boundary particles
        self.bindices=[self.N, self.N +1]
        padd=[]
        labels=[]
        
        #If in contact with inner boundart
        p1 =  np.nonzero(self.boundary == -1)[0]
        padd.extend(p1)
        labels.extend([0 for k in range(len(p1))])
        #If in contact with outer boundary
        p2 =  np.nonzero(self.boundary == 1)[0]
        padd.extend(p2)
        labels.extend([1 for k in range(len(p2))])
        
        #Some arrays for the information of the newly introduced contacts
        fullmobi_add=[]
        fnor_add=[]
        ftan_add=[]
        nx_add=[]
        ny_add=[]
        for k in range(len(padd)):
            self.I = np.append(self.I, self.bindices[labels[k]])
            self.J = np.append(self.J, padd[k])
            neii=np.nonzero(self.I[:self.ncon]==padd[k])[0]
            neij=np.nonzero(self.J[:self.ncon]==padd[k])[0]

            #Always add double bound
            fullmobi_add.append(0)
            
            #Compute forces if neighbours are present
            if (len(neii)>0 or len(neij)>0):
                #Compute angle of particle
                ang = np.arctan2(self.y[padd[k]], self.x[padd[k]])
                #If in contact with inner boundary
                if (labels[k])==0:
                    nx0=-np.cos(ang)
                    ny0=-np.sin(ang)
                
                #If in contact with outer boundary
                elif (labels[k]==1):
                    nx0=np.cos(ang)
                    ny0=np.sin(ang)

                # Compute force on this contact by force balance
                # Two minus signs on second part cancel out
                ftotx=np.sum(self.fnor[neii]*self.nx[neii]-self.ftan[neii]*self.ny[neii])-np.sum(self.fnor[neij]*self.nx[neij]-self.ftan[neij]*self.ny[neij])
                ftoty=np.sum(self.fnor[neii]*self.ny[neii]+self.ftan[neii]*self.nx[neii])-np.sum(self.fnor[neij]*self.ny[neij]+self.ftan[neij]*self.nx[neij])
                fnor0=ftotx*nx0+ftoty*ny0
                ftan0=ftotx*(-ny0)+ftoty*nx0

            else: fnor0 = ftan0 = nx0 = ny0 = 0 #if no neighbours, forces are zero
            
            #Add data to list
            fnor_add.append(fnor0)
            ftan_add.append(ftan0)
            nx_add.append(nx0)
            ny_add.append(ny0)
        
        #Add slipping contact between the two boundary particles
        self.I = np.append(self.I, self.bindices[0])
        self.J = np.append(self.J, self.bindices[1])
        fullmobi_add.append(1)

        #Forces between boundary particles are zero?
        fnor0 = ftan0 = nx0 = ny0 = 0

        #Add to lists
        fnor_add.append(fnor0)
        ftan_add.append(ftan0)
        nx_add.append(nx0)
        ny_add.append(ny0)

        # Finally stick it at the end of the existing data
        self.x=np.concatenate((self.x,Boundaries[:,0]))
        self.y=np.concatenate((self.y,Boundaries[:,1]))
        self.id=np.concatenate((self.id, self.bindices))
        self.rad=np.concatenate((self.rad,Boundaries[:,2]))
        self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
        self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
        self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
        self.nx=np.concatenate((self.nx,np.array(nx_add)))
        self.ny=np.concatenate((self.ny,np.array(ny_add)))
        self.ncon=len(self.I)
        self.N+=2 #Since we have two boundary particles
        print ("Added boundaries!")
    



########## Lattice Configuration #################
class ConfigurationLattice(Configuration):

    def __init__(self,param):
        print('ConfigurationLattice: Created new Configuration Lattice')
        super(ConfigurationLattice,self).__init__('lattice')

        #self.folder = param["folder"]
        #self.setup = param["setup"] # redundant as needs to be same 'experiment_annulus' as in runscript
        self.setup = 'lattice'
        #self.experiment = False

        print ("Reading lattic data!")
        # # Use mu to do the boundary condition ...
        # bc = param["bc"]
        # self.periodic=False
        # if bc =='open':
        #     self.periodic=False
        #     self.periodicx = False
        #     self.periodicy = False
        #     print('open lattice')
        # elif bc =='xy':
        #     self.periodic=True
        #     self.periodicx = True
        #     self.periodicy = True
        # elif bc =='x':
        #     self.periodicx = True
        #     self.periodicy = False
        #     print('x periodic lattice')
        # elif bc =='y':
        #     self.periodicy = True
        #     self.periodicx = False
        #     print('y periodic lattice')
        # else:
        #    self.periodic = False

        # and the system size
        nhex1 = self.params["nhex1"]
        nhex2 = self.params["nhex2"]
        self.Lx = nhex1*np.sqrt(3)/2
        self.Ly = nhex2*1.5

        # For ever hardcoded because maths
        self.hasAngles=False
        # basis for mass computations
        self.density = 1.0
        # Prefactor of the stiffness coefficients
        self.stiffness = 1.0
        # radius conversion factor
        self.rconversion = 1.0
        self.height=1.0
        self.width=1.0
        self.strain=0.0

    # Read in data from posfile and confile, which is full file names set by the running script
    def ReadData(self, posfile,confile,step=0,verbose=False):  
    #######=========Lattice data read in ==========
    #def readLatticedata(self,lattice,bc,verbose0):
        #self.verbose=verbose0
        #filename = self.folder +'/particle_positions' + lattice + bc + '.txt'
        if (verbose):
            print(posfile)
        self.isPosdata=False
        # Only if: 1) The file exists 2) is not empty 
        try:
            data=np.loadtxt(posfile, delimiter=',')
            if (len(data.shape)>1):
                self.isPosdata=True	
        except:
            pass
                    
        if (not self.isPosdata):
            print('error: no position data')
            self.x=0.0
            self.y=0.0
            self.N=0
        else:
            self.x = data[:,1]
            self.y = data[:,2]
            self.N = len(self.x)
            self.rad = 0.5*np.ones((self.N,))
            del data

        # contact data
        #filename = self.folder +'/Adjacency_list' + lattice + bc + '.txt'
        if (verbose):
            print(confile)
        self.isCondata=False
        # Only if: 1) The file exists 2) is not empty 3) is longer than 1 contact (garbage avoiding device .. hope for the best. Some of these got cut off mid-writing)
        try:
            data=np.loadtxt(confile, delimiter=',')
            if (len(data.shape)>1):
                self.isCondata=True	
        except:
            pass
        if ((not self.isCondata) or (not self.isPosdata)):
            print('error: no contact data')
            self.noConf=True
            self.I=-1
            self.J=-1
            self.fullmobi=0
            self.Ifull=-1
            self.Jfull=-1
        else:
            self.noConf=False
            self.I=list(data[:,0].astype(int)-1)
            self.J=list(data[:,1].astype(int)-1)
            self.fullmobi=data[:,4].astype(int)
            self.nx=data[:,2]
            self.ny=data[:,3]
            self.fnor = 0*self.nx
            self.ftan = 0*self.nx
            del data
        if self.periodicx:
            self.x-=self.Lx*np.round(self.x/self.Lx)
        if self.periodicy:
            self.y-=self.Ly*np.round(self.y/self.Ly)
        self.ncon=len(self.I)

        return 0

########## Simulation Configuration #################
class ConfigurationSimulation(Configuration):

    def __init__(self,param):
        print("ConfigurationSimulation: Created new Configuration Simulation")
        super(ConfigurationLattice,self).__init__('simulation')

        # Read in parameter file, sets N, mu k_rep, L, Lx, Ly, rad, xi, phi and gammadot
        self.getParameters(self.params["folder"])
        #self.folder = param["folder"]
        #self.setup = param["setup"] # redundant as needs to be same 'experiment_annulus' as in runscript
        self.setup = 'simulation'
        #self.experiment = False   

        try:
            self.distSnapshot=self.params['distSnapshot']
        except:
            self.distSnapshot=100.0
        
        # Simulation data has periodic boundary conditions and also angles as output
        self.periodic=True
        self.hasAngles=True
        # basis for mass computations
        self.density = 1.0
        # Prefactor of the stiffness coefficients
        self.stiffness = 1.0
        # radius conversion factor
        self.rconversion = 1.0
        self.height=1.0
        self.width=1.0


    #======= Simulation parameter read-in =====================
    # use Parameter read-in as function called by the constructor
    # eventually get rid of this in favour of a proper parameter file
    def getParameters(self,folder):
        parname=folder + 'Parameters_final.dat'
        if os.path.exists(parname):	column=1
        else:
            parname=folder +'Parameters_LE.par'
            if os.path.exists(parname):	
                column=2
            else:	
                print('Error: Couldn\'t find parameter file.')
        parfile = open(parname,'r')
        for line in parfile:
            s = line.split()
            if s == []: s = ['']
            if s[0] == 'np': self.N = eval(s[column])		# number of particles
            elif s[0] == 'k_rep': self.krep = eval(s[column])	# spring constant
            elif s[0] == 'boxVec1x': self.L = eval(s[column])	# system size
            elif s[0] == 'mu': self.mu = eval(s[column])		# friction coefficient
            elif s[0] == 'xivisc': self.xi = eval(s[column]) # viscous damping coefficient
            elif s[0] == 'phi': self.phi = eval(s[column]) # packing fraction
            elif s[0] == 'strain_rate': self.gammadot = eval(s[column]) # packing fraction
        parfile.close()
        self.Lx=self.L 
        self.Ly=self.L
        self.rad=np.loadtxt(folder +'posinit.dat',usecols=(3,))

    # Read in data from posfile and confile, which is full file names set by the running script
    def ReadData(self, posfile,confile,snap,verbose=False):  
    #######========== Simulation data read-in ================== 
    #def readSimdata(self,snap,verbose0,distSnapshot0=100.0):
        #self.verbose=verbose0
        # snap is the label, and distSnapshot tells me how many time units they are apart (in actual time, not steps)
        self.strain=self.gammadot*(1.0*snap)*self.distSnapshot
        if (verbose):
            print('L=' + str(self.L))
            print('Np=' + str(self.N))
            print('mu=' + str(self.mu))
            print('xi=' + str(self.xi))
            print('phi=' + str(self.phi))
            print('strain=' + str(self.strain))
        #filename=self.folder + '/PosVel' + str(snap) + '.dat'
        filename = posfile + str(snap) + '.dat'
        if (verbose):
            print(filename)
        self.isPosdata=False
        # Only if: 1) The file exists 2) is not empty 
        try:
            data=np.loadtxt(filename)
            if (len(data.shape)>1):
                self.isPosdata=True	
        except:
            pass
                
        if (not isPosdata):
            print('error: no position data at snapshot' + str(snap))
            self.x=0.0
            self.y=0.0
            self.dx=0.0
            self.dy=0.0
            self.alpha=0.0
            self.N=0
        else:
            self.x = data[:,0]
            self.y = data[:,1]
            self.alpha = data[:,2]
            self.dx = data[:,3]
            self.dy = data[:,4]
            self.dalpha = data[:,5]
            #print(self.x)
            del data
        
        #filename=self.folder + '/Contact' + str(snap) + '.dat'
        filename = confile + str(snap) + '.dat'
        if (verbose):
            print(filename)
        self.isCondata=False
        # Only if: 1) The file exists 2) is not empty 3) is longer than 1 contact (garbage avoiding device .. hope for the best. Some of these got cut off mid-writing)
        try:
            data=np.loadtxt(filename)
            if (len(data.shape)>1):
                self.isCondata=True	
        except:
            pass
        if ((not self.isCondata) or (not self.isPosdata)):
            print('error: no position data at snapshot' + str(snap))
            self.noConf=True
            self.I=-1
            self.J=-1
            self.fullmobi=0
            self.Ifull=-1
            self.Jfull=-1
        else:
            self.noConf=False
            self.I=list(data[:,0].astype(int))
            self.J=list(data[:,1].astype(int))
            self.fullmobi=data[:,4].astype(int)
            self.nx=data[:,2]
            self.ny=data[:,3]
            self.fnor=data[:,5]
            self.ftan=data[:,6]+data[:,7]
            del data
            self.x-=self.L*np.round(self.x/self.L)
            self.y-=self.L*np.round(self.y/self.L)
            self.ncon=len(self.I)
            



 
         
        
                    
          
       


       
        
        
            
            
        
        
