import sys, os, glob
import numpy as np

class Dynamics:

    __init__(self,configurations0):
        self.configurations=configurations0
    

           # same kind of procedure, but get only the next positions and don't redefine things
        # This is for comparison of the displacements compared to the the next dataset to pebbles and modes
        def ReadExpdataNextSquare(self,numlabel,scale=False):
            
                prefix = self.prefix1 + numlabel + self.prefix2
                print ("Reading experimental step " + prefix + " as next data set.")
                # Let's get the positions first
                self.isPosdataNext=True
                try:
                    coords=np.loadtxt(self.folder+prefix + 'ParticleData.dlm',delimiter=',')
                    self.xnext=coords[:,1]
                    self.ynext=coords[:,2]
                    self.radnext=coords[:,3]
                    self.Nnext=len(self.rad)
                    self.Lxnext=np.amax(self.x)-np.amin(self.x)
                    self.Lynext=np.amax(self.y)-np.amin(self.y)
                    del coords
                except:
                    self.isPosdataNext=False
                    self.xnext=0
                    self.ynext=0
                    self.radnext=1 # purely so that we can divide by it ...
                    self.Nnext=0
                # Slight overkill, but then this enforces unified syntax with simulation, making our job below easier
                if len(self.xnext)==len(self.x):
                    self.dx=self.xnext-self.x
                    self.dy=self.ynext-self.y


     def AddNextBoundaryContactsSquare(self,threshold=15,Brad=20.0):
            # Threshold to check if a particle is close enough to walls.
            upidx=np.argmax(self.ynext)
            downidx=np.argmin(self.ynext)
            leftidx=np.argmin(self.xnext)
            rightidx=np.argmax(self.xnext)
            
            # Boundary posiitons:
            # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            up=self.ynext[upidx]
            yup = up+Brad+self.radnext[upidx]
            down=self.ynext[downidx]
            ydown = down-Brad-self.radnext[downidx]
            left=self.xnext[leftidx]
            xleft=left-Brad-self.radnext[leftidx]
            right=self.xnext[rightidx]
            xright=right+Brad+self.radnext[rightidx]
            
            # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            Boundaries=np.zeros((4,3)) # Four boundary particles with their x,y and rad
            Boundaries[0,:]=[(left+right)*0.5,yup,Brad]
            Boundaries[1,:]=[(left+right)*0.5,ydown,Brad]
            Boundaries[2,:]=[xleft,(up+down)*0.5,Brad]
            Boundaries[3,:]=[xright,(up+down)*0.5,Brad]
            
            
            self.xnext=np.concatenate((self.xnext,Boundaries[:,0]))
            self.ynext=np.concatenate((self.ynext,Boundaries[:,1]))
            self.radnext=np.concatenate((self.radnext,Boundaries[:,2]))

            self.dx=self.xnext-self.x
            self.dy=self.ynext-self.y
            self.Nnext+=4



    #### =================== Computes the D2_min, amount of non-affine motion around a particle ==============
        def getD2min(self,threshold_rad):
                D2_min=np.zeros(len(self.x))
                for k in range(len(self.x)):
                    # Boundary is dubious
                    # the 10^5 is bad, but then they tend to be that value for experiment
                    # just a way to get the color scale to not f up
                    if (k in self.bindices):
                        D2_min[k]=1e5
                    else:
                        temp=[]
                        dist2=(self.x-self.x[k])**2+(self.y-self.y[k])**2
                        rad2=(self.rad[k]+np.mean(self.rad)+threshold_rad)**2
                        # yes, that includes myself, but that term drops out
                        temp=np.nonzero(dist2<rad2)
                        #for i in range(len(self.x)):
                                #if np.sqrt((self.x[k]-self.x[i])**2+(self.y[k]-self.y[i])**2)<(self.rad[k]+threshold_rad):
                                        #temp.append(i)
                        X=np.zeros((2,2))
                        Y=np.zeros((2,2))
                        epsilon=np.zeros((2,2))
                        if len(temp[0])>0:
                                
                                for neighbor in temp[0]:
                                        dx_next=self.xnext[neighbor]-self.xnext[k]
                                        dy_next=self.ynext[neighbor]-self.ynext[k]
                                        dx=self.x[neighbor]-self.x[k]
                                        dy=self.y[neighbor]-self.y[k]
                                        X[0,0]+=dx_next*dx
                                        X[0,1]+=dx_next*dy
                                        X[1,0]+=dy_next*dx
                                        X[1,1]+=dy_next*dy
                                        Y[0,0]+=dx*dx
                                        Y[0,1]+=dx*dy
                                        Y[1,0]+=dy*dx
                                        Y[1,1]+=dy*dy
                                epsilon[0,0]+=(X[0,0]/Y[0,0]+X[0,1]/Y[0,1])
                                epsilon[0,1]+=(X[0,0]/Y[1,0]+X[0,1]/Y[1,1])
                                epsilon[1,0]+=(X[1,0]/Y[0,0]+X[1,1]/Y[0,1])
                                epsilon[1,1]+=(X[1,0]/Y[1,0]+X[1,1]/Y[1,1])

                                for neighbor in temp[0]:
                                        dx_next=self.xnext[neighbor]-self.xnext[k]
                                        dy_next=self.ynext[neighbor]-self.ynext[k]
                                        dx=self.x[neighbor]-self.x[k]
                                        dy=self.y[neighbor]-self.y[k]
                                        D2_min[k]+=((dx_next- (epsilon[0,0]*dx+epsilon[0,1]*dy))**2+ (dy_next-(epsilon[1,0]*dx+epsilon[1,1]*dy))**2)
                        
                                                                    
                return D2_min
                
        
        # ====== Experimental or simulation displacements, decomposed into normal, tangential and potentiallt rotational displecements
        def Disp2Contacts(self,minThresh,debug=False):
                disp2n=np.zeros(self.ncon)
                disp2t=np.zeros(self.ncon)
                if self.hasAngles:
                    disp2r=np.zeros(self.ncon)
                    disp2gear=np.zeros(self.ncon)
                for k in range(self.ncon):
                        i=self.I[k]
                        j=self.J[k]
                        nx0=self.nx[k]
                        ny0=self.ny[k]
                        tx0=-self.ny[k]
                        ty0=self.nx[k]
                        disp2n[k] = ((self.dx[j]-self.dx[i])*nx0+(self.dy[j]-self.dy[i])*ny0)**2
                        disp2t[k] = ((self.dx[j]-self.dx[i])*tx0+(self.dy[j]-self.dy[i])*ty0)**2
                        if self.hasAngles:
                            disp2r[0]=(self.rad[j]*self.dalpha[j]-self.rad[i]*self.dalpha[i])**2
                            disp2gear[0]=((self.dx[j]-self.dx[i])*tx0+(self.dy[j]-self.dy[i])*ty0 - (self.rad[i]*self.dalpha[i]+self.rad[j]*self.dalpha[j]))**2
                thresh=np.mean(disp2n+disp2t)
                print ("Internally generated threshold of " + str(thresh))
                if thresh<minThresh:
                    thresh=minThresh
                if self.hasAngles:
                    return disp2n, disp2t, disp2r, disp2gear,thresh
                else:
                    return disp2n, disp2t, thresh