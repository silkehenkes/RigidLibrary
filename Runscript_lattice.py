import random
import sys, os, glob

# Load pebble game and dynamical matrix libraries
import Configuration as CF
import Pebbles as PB
import Hessian as HS
import Analysis as AN
import Tiling as TY

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import csv 

#topdir='/directory/where/experimental/data/is/located/'
topdir='./DataLattice'
#topdir='./KyungLattice'

# experimental friction coefficient
mu=5

# Experiment type, for other types please use a different runscript.
datatype = 'lattice'
bc='open'
#bc='x' 
#bc='y'
#bc='xy'

#nhex1 = 100 ########## please put the lattice size
#nhex2 = 40
nhex1 = 20
nhex2 = 20
# Potential loop over multiple lattices
lattice_nums=['_'+str(nhex1)+'by'+str(nhex2)+'_']
#lattice_nums=['_100by20_']

# Loop over experiment
for lattice in lattice_nums:
        # Extract some information from the data set:
        filename =topdir+'/Adjacency_list' + lattice + bc + '.txt'
        print(filename)
        try:
            #data = np.loadtxt(topdir+'/Adjacency_list' + lattice + bc + '.txt', delimiter=',')
            data = np.loadtxt(filename, delimiter=',')
            print(data[0])
            #nsteps is the maximal frame number
            nsteps = np.max(data[:,0]).astype(int)
        except:
            print('No valid data found for lattice ' + lattice)
            #if no data is found for an experiment number go to the next experiment
            continue
        
        
        #Creating configuration

        ThisConf = CF.Configuration(topdir,datatype, bc, nhex1,nhex2)

        # def __init__(self, folder, datatype, bc='open',nhex1=20, nhex2=20, mu0=0.2, strainstep=0.1):
                
        #Reading in the data
        ThisConf.readLatticedata(lattice, bc,True)

        # Play the pebble game
        ########## Setting up and playing the pebble game
        ThisPebble = PB.Pebbles(ThisConf,3,3,'nothing',False)

        # play game
        ThisPebble.play_game()
        # compute rigid clusters
        cidx, clusterall, clusterallBonds, clusteridx, BigCluster=ThisPebble.rigid_cluster()

        #print(ThisPebble.pcluster)

        maxlabels = ThisPebble.MaxRigidPos()
        
        pdata = np.zeros((len(maxlabels),3))
        pdata[:,0] = maxlabels
        pdata[:,1] = ThisConf.x[maxlabels]
        pdata[:,2] = ThisConf.y[maxlabels]
        filename = topdir + f'/MaxCluster_{nhex1}by{nhex2}_{bc}.dat'
        with open(filename,'w',newline="\n") as f:
            wr=csv.writer(f)
            wr.writerows(pdata)

        contactmat = ThisPebble.MaxRigidContacts()

        ########## Have a look at some analysis functions of the rigid clusters
        #def __init__(self,conf0,pebbles0,hessian0,tiling0='skip',fgiven0=0.001,verbose0=False):
        ThisAnalysis=AN.Analysis(ThisConf,ThisPebble,0)

        # cluster statistics
        frac,fracmax,lenx,leny=ThisAnalysis.clusterStatistics()
        print('We have a system with a total fraction ' + str(frac) + ' and ' + str(fracmax) + ' in the largest rigid cluster ')
        print('with cluster lengths lenx ' + str(lenx) + ' and leny ' + str(leny))

        #def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
        #fig1 = ThisAnalysis.plotStresses(False,False,True,False,False)    

        #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):  
        #fig2 = ThisAnalysis.plotPebbles(True,True,True,False,True)  

        fig3 = ThisAnalysis.plotPebbles(True,False,False,True,False)  
        
        plt.show()