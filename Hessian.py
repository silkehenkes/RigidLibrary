#!/usr/bin/python

# Silke Henkes, 29.08.18: 
# Created Hessian class to compute and analyse the Dynamical matrix
# Detangled code parts, contains: 
# - Hessian matrix definition
# - Mode computation
# - Diagnostic plotting functions (see Analysis for publication-ready plotting)
# - Analysis function to decompose modes, called from Analysis class


from Configuration import*
 
# Here lives the diagonalisation algorithm
from numpy import linalg as LA

# Some amount of plotting library
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 4
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=16.0
matplotlib.rcParams['legend.fontsize']=14.0

class Hessian:
    
    #### === Hessian constructor, just giving it elematary information, does nothing with it ====
    def __init__(self,conf0):
        self.conf=conf0
        # do throw it a bone - basic dimensional information ...
        self.N=self.conf.N
        self.ncon=self.conf.ncon

    #=========================== Hessian calculation ===============================
    def makeHessian(self,frictional,recomputeFnor=False,stabilise=1e-8,verbose=False):
        # This matrix is 3N by 3N.
        print "Hessian: Info - allocating the " + str(3*self.N) + " by " + str(3*self.N) + " Hessian matrix."
        self.Hessian=np.zeros((3*self.N,3*self.N))
        # For friction, this needs to be a loop over all contacts. Not the Ifull one, the actual contact one!
        # compute average off equilibrium force
        favx=0
        favy=0
        frotav=0
        mrad=np.mean(self.conf.rad)
        mscale=self.conf.density*self.conf.height*np.pi*(self.conf.rconversion*mrad)**2
        self.eigscale=self.conf.stiffness/mscale
        print "Estimating eigenvalue scale at " + str(self.eigscale) + " kg/s^2"
        # Careful now: We are only single-counting contacts these days
        for k in range(self.ncon):
            i=self.conf.I[k]
            j=self.conf.J[k]
            # particle masses
            mi = self.conf.density*self.conf.height*np.pi*(self.conf.rconversion*self.conf.rad[i])**2
            mj = self.conf.density*self.conf.height*np.pi*(self.conf.rconversion*self.conf.rad[j])**2
            Ai = 1.0/2**0.5
            Aj = 1.0/2**0.5
            # The boundary particle is always i, per construction (see configuration)
            # see notes for detailed scaling 
            if self.conf.addBoundary:
                if (i in self.conf.bindices):
                    mi = self.conf.density*self.conf.height*self.conf.width*self.conf.rconversion*0.5*(self.conf.Lx+self.conf.Ly)
                    Ai = (12.0*mrad**2/((0.5*(self.conf.Lx+self.conf.Ly))**2+(2*mrad)**2))**0.5
                            
            nx0=self.conf.nx[k]
            ny0=self.conf.ny[k]
            tx0=-self.conf.ny[k]
            ty0=self.conf.nx[k]
            # x, y, R alpha components (note how R alpha has same dimensions as x and y)
            # see notes
            # Briefly: V = 1/2 \sum_<ij? [ k_n (dr.n)^2 -f_n/r (dt.t)^2 + kt Dt^2 ] , Dt = (dr.t) - (R_i da_i + R_j da_j)
            # We have here k_n = k_t = 1
            # See actual notes for derivation!!

                        # Simulation default repulsion constant
                        # Never make experiment recompute forces
            kn=self.conf.stiffness
            # get the r betwen particles
            dx=self.conf.x[j]-self.conf.x[i]
            if self.conf.periodic:
                dx-=self.conf.Lx*round(dx/self.conf.Lx)
            dy=self.conf.y[j]-self.conf.y[i]
            if self.conf.periodic:
                dy-=self.conf.Ly*round(dy/self.conf.Ly)
            rval=np.sqrt(dx**2+dy**2)
            if recomputeFnor:
                if self.conf.datatype=='experiment':
                    print "Warning: recomputing normal forces for the experiment. This is a bad idea!"
                # that's just overlap delta
                fn = kn*(self.conf.rad[i]+self.conf.rad[j]-rval)*self.conf.rconversion
                if fn<0:
                    print "Funny contact with negative overlap between " + str(i)+ " and " + str(j)
                    fn=0.0
            else:
                fn=self.conf.fnor[k]
            if frictional:
                if self.conf.fullmobi[k]==0:
                    kt=self.conf.stiffness
                else:
                    kt=0.0
            else:
                kt = 0.0
            # # This is our litte square in local coordinates (where nonzero)
            subsquare=np.zeros((3,3))
            subsquare[0,0]=-kn
            subsquare[1,1]=fn/(self.conf.rconversion*rval)-kt
            subsquare[1,2]=kt*Aj
            # note asymmetric cross-term
            subsquare[2,1]=-kt*Ai
            subsquare[2,2]=kt*Ai*Aj
            #collect our little forces
            favx+=fn*nx0+self.conf.ftan[k]*tx0
            favy+=fn*ny0+self.conf.ftan[k]*ty0
            frotav+=self.conf.rad[i]*self.conf.ftan[k]

            # Stick this into the appropriate places after rotating it away from the (n,t) frame
            Hij=np.zeros((3,3))
            Hij[0,0]=subsquare[0,0]*nx0**2+subsquare[1,1]*tx0**2
            Hij[0,1]=subsquare[0,0]*nx0*ny0+subsquare[1,1]*tx0*ty0
            Hij[1,0]=subsquare[0,0]*ny0*nx0+subsquare[1,1]*ty0*tx0
            Hij[1,1]=subsquare[0,0]*ny0**2+subsquare[1,1]*ty0**2
            Hij[0,2]=subsquare[1,2]*tx0
            Hij[1,2]=subsquare[1,2]*ty0
            Hij[2,0]=subsquare[2,1]*tx0
            Hij[2,1]=subsquare[2,1]*ty0
            Hij[2,2]=subsquare[2,2]

            # And put it into the Hessian, with correct elasticity prefactor
            # once for contact ij
            self.Hessian[3*i:(3*i+3),3*j:(3*j+3)]=Hij/(mi*mj)**0.5
            
            # see notes for the flip one corresponding to contact ji
            # both n and t flip signs. Put in here explicitly. Essentially, angle cross-terms flip sign
            # Yes, this is not fully efficient, but it's clearer. Diagonalisation is rate-limiting step, not this.
            # careful with the A's 
            subsquare[1,2]=subsquare[1,2]*Ai/Aj
            subsquare[2,1]=subsquare[2,1]*Aj/Ai
            Hji=np.zeros((3,3))
            Hji[0,0]=subsquare[0,0]*(-nx0)**2+subsquare[1,1]*(-tx0)**2
            Hji[0,1]=subsquare[0,0]*(-nx0)*(-ny0)+subsquare[1,1]*(-tx0)*(-ty0)
            Hji[1,0]=subsquare[0,0]*(-ny0)*(-nx0)+subsquare[1,1]*(-ty0)*(-tx0)
            Hji[1,1]=subsquare[0,0]*(-ny0)**2+subsquare[1,1]*(-ty0)**2
            Hji[0,2]=subsquare[1,2]*(-tx0)
            Hji[1,2]=subsquare[1,2]*(-ty0)
            Hji[2,0]=subsquare[2,1]*(-tx0)
            Hji[2,1]=subsquare[2,1]*(-ty0)
            Hji[2,2]=subsquare[2,2]

            # And put it into the Hessian
            # now for contact ji
            self.Hessian[3*j:(3*j+3),3*i:(3*i+3)]=Hji/(mi*mj)**0.5
            
            
            # Careful, the diagonal bits are not just minus because of the rotations
            diagsquare=np.zeros((3,3))
            diagsquare[0,0]=kn
            diagsquare[1,1]=-fn/(self.conf.rconversion*rval)+kt
            diagsquare[1,2]=kt*Ai
            diagsquare[2,1]=kt*Ai
            diagsquare[2,2]=kt*Ai**2
            
            # Stick this into the appropriate places:
            Hijdiag=np.zeros((3,3))
            Hijdiag[0,0]=diagsquare[0,0]*nx0**2+diagsquare[1,1]*tx0**2
            Hijdiag[0,1]=diagsquare[0,0]*nx0*ny0+diagsquare[1,1]*tx0*ty0
            Hijdiag[1,0]=diagsquare[0,0]*ny0*nx0+diagsquare[1,1]*ty0*tx0
            Hijdiag[1,1]=diagsquare[0,0]*ny0**2+diagsquare[1,1]*ty0**2
            Hijdiag[0,2]=diagsquare[1,2]*tx0
            Hijdiag[1,2]=diagsquare[1,2]*ty0
            Hijdiag[2,0]=diagsquare[2,1]*tx0
            Hijdiag[2,1]=diagsquare[2,1]*ty0
            Hijdiag[2,2]=diagsquare[2,2]
            
            # And then *add* it to the diagnual
            self.Hessian[3*i:(3*i+3),3*i:(3*i+3)]+=Hijdiag/mi
            
            #And once more for the jj contribution, which is the same whizz with the flipped sign of n and t 
            # and adjusted A's
            diagsquare=np.zeros((3,3))
            diagsquare[0,0]=kn
            diagsquare[1,1]=-fn/(self.conf.rconversion*rval)+kt
            diagsquare[1,2]=kt*Aj
            diagsquare[2,1]=kt*Aj
            diagsquare[2,2]=kt*Aj**2
            
            Hjidiag=np.zeros((3,3))
            Hjidiag[0,0]=diagsquare[0,0]*(-nx0)**2+diagsquare[1,1]*(-tx0)**2
            Hjidiag[0,1]=diagsquare[0,0]*(-nx0)*(-ny0)+diagsquare[1,1]*(-tx0)*(-ty0)
            Hjidiag[1,0]=diagsquare[0,0]*(-ny0)*(-nx0)+diagsquare[1,1]*(-ty0)*(-tx0)
            Hjidiag[1,1]=diagsquare[0,0]*(-ny0)**2+diagsquare[1,1]*(-ty0)**2
            Hjidiag[0,2]=diagsquare[1,2]*(-tx0)
            Hjidiag[1,2]=diagsquare[1,2]*(-ty0)
            Hjidiag[2,0]=diagsquare[2,1]*(-tx0)
            Hjidiag[2,1]=diagsquare[2,1]*(-ty0)
            Hjidiag[2,2]=diagsquare[2,2]
            
            # And then *add* it to the diagnual
            self.Hessian[3*j:(3*j+3),3*j:(3*j+3)]+=Hjidiag/mj
            


        # add a very small, diagonal bit to stabilise zero modes. Normalise by standard units
        self.Hessian[range(3*self.N),range(3*self.N)]+=stabilise*self.eigscale
        favx/=self.ncon
        favy/=self.ncon
        frotav/=self.ncon
        #print "Hessian: Estimating distance from mechanical equilibrium of initial configuration "
        
        if verbose:
            # This bit is commented since it eats up 10 gigs of memory for experimental-size Hessians
            #plt.figure()
            #plt.pcolor(self.Hessian)
            print "Scaled force sum is x:" + str(favx) + "  y:" + str(favy) + " rot:" + str(frotav)


        ## ==== Diagonalisation routine, with optional eigenvalue plotting ====
    def getModes(self,debug=False):
        # Let's have a look if what we get is in any way reasonable
        # Eigenvalues and eigenvectors
        # Only symmetrise to calculate - for clarity and debugging above
        HessianSym=0.5*(self.Hessian+np.transpose(self.Hessian))
        # Use routines for hermitian eigenvector decomposition
        # Default is ascending order, which suits us
        print "Starting Diagonalisation!"
        self.eigval, self.eigvec = LA.eigh(HessianSym)
        if debug:
            # start with some debugging output
            fig=plt.figure()
            eigrank=np.linspace(0,3*self.N,3*self.N)
            plt.plot(eigrank,self.eigval,'o')
            plt.xlabel('rank')
            plt.ylabel('eigenvalue')
            plt.title('Eigenvalues of the Hessian')
        return fig
        
        ##=============== Diagnostic output =============================================
        # Debugging output to check on usepts modes
    def plotModes(self,usepts):
        for u in usepts:
            fig=plt.figure()
            plt.quiver(self.conf.x,self.conf.y,self.eigvec[0:3*self.N:3,u],self.eigvec[1:3*self.N:3,u])
            # dimensional contributions
            wx=np.sum(self.eigvec[0:3*self.N:3,u]**2)
            wy=np.sum(self.eigvec[1:3*self.N:3,u]**2)
            wrot=np.sum(self.eigvec[2:3*self.N:3,u]**2)
            plt.title('mode #' + str(u) + ', eigval ' + str(self.eigval[u])+ ' ' + str(np.round(wx+wy,2)) + ' trans, ' + str(np.round(wrot,2)) + ' rot.' )
        return fig
        # first indication if any of the large displacement regions make sense
        # purely qualiative, use correlation plots in Analysis to do this properly
    def plotZeroModes(self,thresh=2e-8,simple=True):
        fig=plt.figure()
        if simple:
            nzero=0
            eigx=np.zeros((self.N,))
            eigy=np.zeros((self.N,))
            for u in range(3*self.N):
                if self.eigval[u]<thresh*self.eigscale:
                    nzero+=1
                    eigx+=self.eigvec[0:3*self.N:3,u]
                    eigy+=self.eigvec[1:3*self.N:3,u]
            eigx/=nzero
            eigy/=nzero
            plt.quiver(self.conf.x,self.conf.y,eigx,eigy)
            plt.title('Location of ' + str(nzero)+' zero eigenmodes')
            plt.gca().axis('scaled')
            if self.conf.datatype=='simulation':
                plt.xlim(-self.conf.Lx/2,self.conf.Lx/2)
                plt.ylim(-self.conf.Ly/2,self.conf.Ly/2)
            elif self.conf.datatype=='experiment':
                xmin=np.amin(self.conf.y)-20
                ymin=np.amin(self.conf.x)-20
                plt.xlim(xmin,xmin+self.conf.Lx+40)
                plt.ylim(ymin,ymin+self.conf.Ly+40)
        else:
            for u in range(3*self.N):
                if self.eigval[u]<thresh*self.eigscale:
                    plt.quiver(self.conf.x,self.conf.y,self.eigvec[0:3*self.N:3,u],self.eigvec[1:3*self.N:3,u])
                    # dimensional contributions
                    wx=np.sum(self.eigvec[0:3*self.N:3,u]**2)
                    wy=np.sum(self.eigvec[1:3*self.N:3,u]**2)
                    wrot=np.sum(self.eigvec[2:3*self.N:3,u]**2)
                    plt.title('mode #' + str(u) + ', eigval ' + str(self.eigval[u])+ ' ' + str(np.round(wx+wy,2)) + ' trans, ' + str(np.round(wrot,2)) + ' rot.' )
        return fig    
    
    ## =========================== Analysis helper functions ===================================
    # Computes total displacements of modes below threshold, based on particles
    def ModeThresh(self,thresh):    
                # Level 0: threshold on square displacements on particle (particle clusters)
        self.nzero=0
        eigx=np.zeros((self.N,))
        eigy=np.zeros((self.N,))
        eigr=np.zeros((self.N,))
        self.iszero=[]
        for u in range(3*self.N):
            if self.eigval[u]<thresh*self.eigscale:
                self.iszero.append(u)
                self.nzero+=1
                eigx+=self.eigvec[0:3*self.N:3,u]**2
                eigy+=self.eigvec[1:3*self.N:3,u]**2
                eigr+=self.eigvec[2:3*self.N:3,u]**2
        eigx/=self.nzero
        eigy/=self.nzero
        eigr/=self.nzero
        self.isrigid=[index for index,value in enumerate(eigx+eigy+eigr) if value<thresh*self.eigscale]
        self.notrigid=[x for x in range(self.N) if x not in self.isrigid]
        
         
        # Decompose zero modes onto normal translations, tangential translations, rotations and sliding
    def ModeContacts(self):
        # Level 1: Work with relative displacements. That still leaves global rotations possible
        # Start to work on contacts
        # normal, tangential and rotational (eff sliding for a little, look for gearing in a bit)

        eign=np.zeros(self.ncon)
        eigt=np.zeros(self.ncon)
        eigr=np.zeros(self.ncon)
        eiggear=np.zeros(self.ncon)
        for k in range(self.ncon):
            i=self.conf.I[k]
            j=self.conf.J[k]
            nx0=self.conf.nx[k]
            ny0=self.conf.ny[k]
            tx0=-self.conf.ny[k]
            ty0=self.conf.nx[k]
            eign[k] = np.sum(((self.eigvec[3*j,self.iszero]-self.eigvec[3*i,self.iszero])*nx0+(self.eigvec[3*j+1,self.iszero]-self.eigvec[3*i+1,self.iszero])*ny0)**2)/self.nzero
            eigt[k] = np.sum(((self.eigvec[3*j,self.iszero]-self.eigvec[3*i,self.iszero])*tx0+(self.eigvec[3*j+1,self.iszero]-self.eigvec[3*i+1,self.iszero])*ty0)**2)/self.nzero
            # note there is a plus here, other sign is a gearing motion
            eigr[k] = np.sum((self.eigvec[3*j+2,self.iszero]+self.eigvec[3*i+2,self.iszero])**2)/self.nzero
            # gearing
            eiggear[k]=np.sum(((self.eigvec[3*j,self.iszero]-self.eigvec[3*i,self.iszero])*tx0+(self.eigvec[3*j+1,self.iszero]-self.eigvec[3*i+1,self.iszero])*ty0-(self.eigvec[3*i+2,self.iszero]+self.eigvec[3*j+2,self.iszero]))**2)/self.nzero # note dimensions of rotational modes are lengths now, with precisely the radii as scale factors
        return eign, eigt, eigr, eiggear





