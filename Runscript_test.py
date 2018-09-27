import random
import sys, os, glob
import numpy as np
import matplotlib.pyplot as plt

# Location of the library files: insert into python path
pebblepath='/home/sh18581/Documents/Friction/RigidLibrary/'
sys.path.insert(1,pebblepath)

# Load pebble game and dynamical matrix libraries
import Configuration as CF
import Pebbles as PB
import Hessian as HS
import Analysis as AN


# Location of data to analyse
#topdir=''
topdir='/Users/kuangliu/Dropbox/Shear/'
# experimental friction coefficient
mu=0.3
# More information about data set
directions=['forward','reverse']
experiment_nums=['02','03','05','06','07','08','09','10','11','12','13','18','19','20','21','22','23','24','25','26','28','29','32','33','34','35']
#inilabel=5384
inilabel=3446
dlabel=2
#nsteps=12
nsteps=1000
AllCorrelations=[[],[],[],[],[],[],[],[],[],[],[]]
# explicit prefixes: These are the default values, so not necessary here unless this changes
prefix1='DSC'
prefix2='_solved_Tracke_'

# Loop over experiment
for experiment in experiment_nums:
        # and shear directions
        for direction in directions:
                foldername=topdir+experiment+'/'+direction+'/'
                print foldername
                # Create a new configuratin; the data in there will get replaced every time we read a new configuration
                #def __init__(self,folder,datatype,mu0=0.2,prefix10='DSC',prefix20='_solved_Tracke_'):
                ThisConf = CF.Configuration(foldername,'experiment',mu,prefix1,prefix2)
                # And finally loop over strain steps for a given experiment
                for u in range(nsteps-1):
                        ########## Reading in the configuration and the next step
                        # read in a specific step, careful with format of label
                        numlabel0 = inilabel+dlabel*u
                        numlabel = "%05d" %numlabel0
                        #def ReadExpdata(self,numlabel,scale=False):
                        ThisConf.ReadExpdata(numlabel)
                        ThisConf.AddBoundaryContacts()
                        #ThisConf.AddSomeContacts(1) # parameter is percentage of extra bonds.
                        # also read in the next data at this point: dlabel further along in the numbering 
                        #def ReadExpdataNext(self,numlabel,scale=False):
                        numlabel2 =  "%05d" %(numlabel0+dlabel)
                        ThisConf.ReadExpdataNext(numlabel2)
                        ThisConf.AddNextBoundaryContacts()
                        if ThisConf.isPosdata==True and ThisConf.isPosdataNext==True:
                        	if ThisConf.ncon>1:
		                        ########## Setting up and playing the pebble game
		                        # def __init__(self,conf, game10,game20,modifier='nothing',verbose0=False):
		                        ThisPebble = PB.Pebbles(ThisConf,3,3,'nothing',False)
		                        # play game
		                        ThisPebble.play_game()
		                        # compute rigid clusters
		                        cidx, clusterall, clusterallBonds, clusteridx, BigCluster=ThisPebble.rigid_cluster()
		                        

		                        ########### Setting up the dynamical matrix and getting eigenmodes
		                        # This itself does very little, just creates an empty Hessian class
		                        # __init__(self,conf0):
		                        ThisHessian = HS.Hessian(ThisConf)
		                        
		                        ########## Have a look at some analysis functions of the rigid clusters
		                        #def __init__(self,conf0,pebbles0,hessian0,verbose=False):
		                        ThisAnalysis=AN.Analysis(ThisConf,ThisPebble,ThisHessian,0.01,False)
		                        # stress statistics
		                        zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx=ThisAnalysis.getStressStat()
		                        # cluster statistics
		                        frac,fracmax,lenx,leny=ThisAnalysis.clusterStatistics()
		                        #def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
		                        fig1 = ThisAnalysis.plotStresses(True,False,False,True,False)
		                        #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):
		                        #ThisAnalysis.plotPebbles(True,True,True,False,False)
		                        #fig2 = ThisAnalysis.plotPebbles(True,True,False,True,False)
		                        
		                        
		                        ######### continuing with the Hessian now 
		                        # constructing the matrix
		                        #  makeHessian(self,frictional,recomputeFnor,stabilise,verbose=False):
		                        ThisHessian.makeHessian(True,False,0,True)
		                        # diagonalising the matrix
		                        # def getModes(self,debug=False):
		                        ThisHessian.getModes(True)
		                        
		                        ##### a couple of checks on the modes (optional)
		                        #plotModes(self,usepts):
		                        #usepts=[3*ThisConf.N-5,3*ThisConf.N-4,3*ThisConf.N-3,3*ThisConf.N-2,3*ThisConf.N-1]
		                        #ThisHessian.plotModes(usepts)
		                        #plotZeroModes(self,thresh=2e-8,simple=True):
		                        ThisHessian.plotZeroModes()
		                        
		                        ############ Now look for the cross-correlations
		                        # what is rigid according to modes with a given threshold:
		                        fig3 = ThisAnalysis.ModeClusters('translations',2e-4)
		                        # how this cross-correlates to rigidy according to pebbles:
		                        P_eig_if_pebble,P_pebble_if_eig,fig5 = ThisAnalysis.RigidModesCorrelate(2e-4)
		                        # These are the conditional probabilities of being rigid by mode while being rigid by pebble and the reverse
		                        print P_eig_if_pebble,P_pebble_if_eig
		                        # if there is a next data set
		                        if ThisConf.Nnext>0:
		                            P_disp_if_pebble,P_pebble_if_disp, fig6 = ThisAnalysis.RigidDisplacementsCorrelate(2e-4)
		                            # Conditional probabilities of being rigid by displacement while being rigid by pebble
		                            print P_disp_if_pebble,P_pebble_if_disp
		                            # D2_min, needs assessment
		                            fig7 = ThisAnalysis.DisplacementCorrelateD2min(True)

		                        AllCorrelations[0].append(float(2.0*ThisConf.ncon/(ThisConf.N-4.0)))
		                        AllCorrelations[1].append(float(2.0*ThisConf.ncon/ThisConf.N_NoRattler))   
		                        AllCorrelations[2].append(P_eig_if_pebble)
		                        AllCorrelations[3].append(P_pebble_if_eig)
		                        AllCorrelations[4].append(P_disp_if_pebble)
		                        AllCorrelations[5].append(P_pebble_if_disp)
		                        AllCorrelations[6].append(direction)
		                        AllCorrelations[7].append(frac)
		                        AllCorrelations[8].append(fracmax)
		                        AllCorrelations[9].append(lenx)
		                        AllCorrelations[10].append(leny)

                                              
                        
#plt.show()       
# print AllCorrelations
np.savetxt('AllCorrelations_of_'+str(len(AllCorrelations[0]))+'_Samples.txt',AllCorrelations)
# data=np.loadtxt('AllCorrelations3_reverse.txt')
# plt.plot(data[:,0],data[:,1],'-o',label='P_eig_if_pebble')
# plt.plot(data[:,0],data[:,2],'-o',label='P_pebble_if_eig')
# plt.plot(data[:,0],data[:,3],'-o',label='P_disp_if_pebble')
# plt.plot(data[:,0],data[:,4],'-o',label='P_pebble_if_disp')
# plt.ylim(0,1)
# plt.legend(loc='upper left')
# plt.show()