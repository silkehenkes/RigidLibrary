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
import Holes as HO


# Location of data to analyse
#topdir=''
topdir=''
# experimental friction coefficient
mu=0.3
# More information about data set
directions=['forward','reverse']
# experiment_nums=['02','03','06','08','09','10','11','12','13','18','19','20','21','22','23','24','25','26','28','29','32','33','34','35']
experiment_nums=['23']
inilabel=5412
# inilabel=3444
dlabel=2
nsteps=2
# nsteps=1700
thresh=2e-4
AllCorrelations=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
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
        NextConf = CF.Configuration(foldername,'experiment',mu,prefix1,prefix2)
        # And finally loop over strain steps for a given experiment 
        for u in range(nsteps-1):
            ########## Reading in the configuration and the next step
            # read in a specific step, careful with format of label
            numlabel0 = inilabel+dlabel*u
            numlabel = "%05d" %numlabel0
            #def ReadExpdata(self,numlabel,scale=False):
            ThisConf.ReadExpdata(numlabel)
            ThisConf.AddBoundaryContacts()
            # also read in the next data at this point: dlabel further along in the numbering 
            #def ReadExpdataNext(self,numlabel,scale=False):
            numlabel2 =  "%05d" %(numlabel0+dlabel)
            ThisConf.ReadExpdataNext(numlabel2)
            ThisConf.AddNextBoundaryContacts()
            if ThisConf.isPosdata==True and ThisConf.isPosdataNext==True:
                # use voronoi tessellation to find areas
                # Area = ThisConf.getArea()
                # filename=str(experiment)+'/'+str(numlabel0)+"_area.txt"
                # np.savetxt(filename,Area)

                #     #evaluate distribution of actual relative displacements
                #     disp2n, disp2t, thresh = ThisConf.Disp2Contacts(minThresh=thresh)
                #     plt.hist(np.log(disp2t+disp2n), bins=40) #-4~3
                #     plt.show()
                if ThisConf.ncon>1:
                    # ThisConf.AddSomeContacts(5) # parameter is percentage of extra bonds.
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
                    # zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx=ThisAnalysis.getStressStat()
                    # cluster statistics
                    frac,fracmax,lenx,leny=ThisAnalysis.clusterStatistics()
                    #def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
                    fig1 = ThisAnalysis.plotStresses(True,False,False,True,False)
                    fig1.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_1.png", bbox_inches='tight')
                    #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):
                    # ThisAnalysis.plotPebbles(True,True,True,False,False)
                    if frac>0:
                        fig2 = ThisAnalysis.plotPebbles(True,True,False,True,False)
                        fig2.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_2.jpg", bbox_inches='tight')
                    
                    
                    ######### continuing with the Hessian now 
                    # constructing the matrix
                    #makeHessian(self,frictional,recomputeFnor,stabilise,verbose=False):
                    ThisHessian.makeHessian(True,False,0,True)
                    # diagonalising the matrix
                    # def getModes(self,debug=False):
                    fig3=ThisHessian.getModes(True)
                    fig3.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_3.jpg", bbox_inches='tight')
                    #### a couple of checks on the modes (optional)
                    # plotModes(self,usepts):
                    # usepts=[3*ThisConf.N-5,3*ThisConf.N-4,3*ThisConf.N-3,3*ThisConf.N-2,3*ThisConf.N-1]
                    # ThisHessian.plotModes(usepts)
                    # plotZeroModes(self,thresh=2e-8,simple=True):
                    fig4=ThisHessian.plotZeroModes()
                    fig4.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_4.jpg", bbox_inches='tight')
                    ########### Now look for the cross-correlations
                    # what is rigid according to modes with a given threshold:
                    fig5 = ThisAnalysis.ModeClusters('translations')
                    fig5.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_5.jpg", bbox_inches='tight')
                    # how this cross-correlates to rigidy according to pebbles:
                    P_eig_if_pebble,P_pebble_if_eig,fig6 = ThisAnalysis.RigidModesCorrelate(thresh)
                    fig6.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_6.jpg", bbox_inches='tight')
                    # These are the conditional probabilities of being rigid by mode while being rigid by pebble and the reverse
                    # print P_eig_if_pebble,P_pebble_if_eig
                    # if there is a next data set
                    if ThisConf.Nnext>0:
                        P_disp_if_pebble,P_pebble_if_disp,P_disp_if_eig,P_eig_if_disp,fig7 = ThisAnalysis.RigidDisplacementsCorrelate(thresh)
                        fig7.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_7.jpg", bbox_inches='tight')
                        # Conditional probabilities of being rigid by displacement while being rigid by pebble
                        # print P_disp_if_pebble,P_pebble_if_disp
                        # D2_min, needs assessment
                        fig8 = ThisAnalysis.DisplacementCorrelateD2min(True)
                        fig8.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_8.jpg", bbox_inches='tight')
                    ThisHo = HO.Holes(ThisConf, ThisPebble, False)        
                    ThisHo.makeHalfEdge()
                    nHole,HoleZ,HolePerim,HoleArea = ThisHo.sortHoles()  
                    # print  nHole
                    ToSave=np.zeros((nHole,3))
                    ToSave[:,0]=HoleZ
                    ToSave[:,1]=HolePerim
                    ToSave[:,2]=HoleArea
                    np.savetxt(str(experiment)+'/'+str(numlabel0)+str(int(10*mu))+"_HoleStat.txt",ToSave)
                    fig9 = ThisHo.plotHoles() 
                    fig9.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_9.jpg", bbox_inches='tight')
                    fig10 = ThisHo.plotFaces() 
                    fig10.savefig(str(experiment)+'/sample'+str(numlabel0)+'_mu='+str(mu)+"_10.jpg", bbox_inches='tight')

                    #create two files containing all information of particles and contacts
                    particle_info=np.zeros((ThisConf.N,3))
                    contact_info=np.zeros((ThisConf.ncon,11))
                    eign, eigt, eigr, eiggear = ThisHessian.ModeContacts()
                    disp2n, disp2t, intThresh = ThisConf.Disp2Contacts(thresh,debug=False)
                    particle_info[:,0]=ThisConf.x
                    particle_info[:,1]=ThisConf.y
                    particle_info[:,2]=ThisConf.rad
                    contact_info[:,0]=np.asarray(ThisConf.I)
                    contact_info[:,1]=np.asarray(ThisConf.J)
                    contact_info[:,2]=ThisConf.fnor
                    contact_info[:,3]=ThisConf.ftan
                    contact_info[:,4]=ThisPebble.cluster[0:ThisConf.ncon]
                    contact_info[:,5]=eign
                    contact_info[:,6]=eigt
                    contact_info[:,7]=eigr
                    contact_info[:,8]=eiggear
                    contact_info[:,9]=disp2n
                    contact_info[:,10]=disp2t

                    
                    filename=str(experiment)+'/'+str(numlabel0)+str(mu)+"_particle_info.txt"
                    np.savetxt(filename,particle_info)
                    filename=str(experiment)+'/'+str(numlabel0)+str(mu)+"_contact_info.txt"
                    np.savetxt(filename,contact_info)
                    filename=str(experiment)+'/'+str(numlabel0)+str(mu)+"_eigval.txt"
                    np.savetxt(filename,ThisHessian.eigval)
                    plt.clf()
                    plt.close()

# plt.show()       

