import random
import sys, os, glob


# Location of the library files: insert into python path
#pebblepath='/home/sh18581/Documents/Friction/RigidLibrary/'
#sys.path.insert(1,pebblepath)

# Load pebble game and dynamical matrix libraries
import Param as PR
import Configuration as CF
import Pebbles as PB
import Hessian as HS
import Analysis as AN
import Holes as HO

import matplotlib.pyplot as plt

# Location of data to analyse
# folder = './DataPrajima/'
folder = './Data_dummy/'
paramfile = 'params.csv'
posfile='particle_positions.txt'
confile='Adjacency_list.txt'
# We will still assume that these are part of a sequence
step=1

ThisPar = PR.Paramreader(folder+paramfile)
# ThisPar.params
print(ThisPar.params)
ThisPar.set_folder(folder)


ThisConf = CF.ConfigurationExpSquare(ThisPar.params)
# ########## Reading in the configuration and the next step
# ReadData(self, posfile,confile,step,verbose=False):    
ThisConf.ReadData(posfile,confile,step,True)
#ThisConf.AddBoundaryContactsSquare()
# # also read in the next data at this point: dlabel further along in the numbering 
# # Revise at a later stage.
# #def ReadExpdataNext(self,numlabel,scale=False):
# #numlabel2 =  "%05d" %(numlabel0+dlabel)
# #ThisConf.ReadExpdataNext(numlabel2)
# #ThisConf.AddNextBoundaryContacts()


# ########## Setting up and playing the pebble game
# # def __init__(self,conf, game10,game20,modifier='nothing',verbose0=False):
ThisPebble = PB.Pebbles(ThisConf,3,3,'nothing',False)
# # play game
ThisPebble.play_game()
# # compute rigid clusters
cidx, clusterall, clusterallBonds, clusteridx, BigCluster=ThisPebble.rigid_cluster()


########### Setting up the dynamical matrix and getting eigenmodes
# This itself does very little, just creates an empty Hessian class
# __init__(self,conf0):
ThisHessian = HS.Hessian(ThisConf)

########## Have a look at some analysis functions of the rigid clusters
#ddef __init__(self,conf0,pebbles0,hessian0,tiling0='skip',fgiven0=0.001,verbose0=False):
ThisAnalysis=AN.Analysis(ThisConf,ThisPebble,ThisHessian,'skip',0.01,False)
# stress statistics
zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx=ThisAnalysis.getStressStat()
# cluster statistics
frac,fracmax,lenx,leny=ThisAnalysis.clusterStatistics()
#def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
fig1 = ThisAnalysis.plotStresses(True,False,False,True,False)
#def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):
#ThisAnalysis.plotPebbles(True,True,True,False,False)
fig2 = ThisAnalysis.plotPebbles(True,True,True,False,False)
fig3 = ThisAnalysis.plotPebbles(True,True,False,True,False)


######### continuing with the Hessian now 
# constructing the matrix
#  makeHessian(self,frictional,recomputeFnor,stabilise,verbose=False):
#ThisHessian.makeHessian(True,False,0,False)
# diagonalising the matrix
# def getModes(self,debug=False):
#ThisHessian.getModes(False)

##### a couple of checks on the modes (optional)
#plotModes(self,usepts):
#usepts=[3*ThisConf.N-5,3*ThisConf.N-4,3*ThisConf.N-3,3*ThisConf.N-2,3*ThisConf.N-1]
#ThisHessian.plotModes(usepts)
#plotZeroModes(self,thresh=2e-8,simple=True):
#ThisHessian.plotZeroModes()

############ Now look for the cross-correlations
# what is rigid according to modes with a given threshold:
#fig3 = ThisAnalysis.ModeClusters('translations',2e-4)
# how this cross-correlates to rigidy according to pebbles:
#P_eig_if_pebble,P_pebble_if_eig,fig5 = ThisAnalysis.RigidModesCorrelate(2e-4)
# These are the conditional probabilities of being rigid by mode while being rigid by pebble and the reverse
#print (P_eig_if_pebble,P_pebble_if_eig)
# if there is a next data set
# Needs revision in any case, don't use for now
#if ThisConf.Nnext>0:
#    P_disp_if_pebble,P_pebble_if_disp, fig6 = ThisAnalysis.RigidDisplacementsCorrelate(2e-4)
#    # Conditional probabilities of being rigid by displacement while being rigid by pebble
#    print (P_disp_if_pebble,P_pebble_if_disp)
#    # D2_min, needs assessment
#    fig7 = ThisAnalysis.DisplacementCorrelateD2min(True)
                                              
