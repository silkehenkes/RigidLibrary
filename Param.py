# Parameters for pebble game dynamics, includes default values
# Location data: not hard-coded, will change every time
# param.folder - location of data
# param.step - assuming it's part of an experimental run / simulation: numerical information about sub-set

# Information about the type of data
# param.setup - type of data: simulation, experiment_square, experiment_annulus, lattice
import csv

class Paramreader:
    def __init__(self,paramfile):
        # store in a dictionary
        self.params = {}
        with open(paramfile, newline='') as csvfile:
            paramreader = csv.reader(csvfile, delimiter=',')
            ## Automate dictionary construction, e.g. file like this
            # datatype, experiment_square,
            # mu ,0.3,
            # periodic , False,
            # has_angles , False,
            # density ,1060, # bulk material density kg/m^3
            # stiffness ,6490, # prefactor of the stiffness coefficient k_n  \frac{\pi}{6} E h  
            # pixel_convert,0.00027,# radius conversion factor from pixel to m.
            # thickness ,0.0031,# thickness of the disks 
            # width ,0.02,  # boundary width
            # prefix1 , DSC      ,# Prefix on the image file names 
            # prefix2 , _solved_Tracks_,
            # experiment , True,
            # position_error,0.1,#not a real value yet
            # force_error,0.1, #not a real value yet
            # pebble_game_m, 3, # pebble game (m,n): m value
            # pebble_game_n, 2, # pebble game (m,n): n value
            for row in paramreader:
                lhs = row[0].strip()
                rhs = row[1]
                # it's either a int, float or a boolean
                try:
                    rhs = eval(row[1])
                except:
                    rhs = row[1].strip()
                self.params[lhs] = rhs
        ## postprocess particular items:
        # Convert main periodic string 'True' / 'False' and make x and y inherit
        # Expect this to be 'checkxy' if x and y are different
        if not self.params['periodic']:
            # if not explicitly set, x and y parts inherit from full one
            self.params['periodicx']=False
            self.params['periodicy']=False
        if self.params['periodic']:
            # if not explicitly set, x and y parts inherit from full one
            self.params['periodicx']=True
            self.params['periodicy']=True
        # Consistency: Reset the global periodic if x and y are different
        if self.params['periodicx'] != self.params['periodicy']:
            self.params['periodic']='checkxy'

    def set_folder(self,folder0):
        self.params["folder"]=folder0
       








