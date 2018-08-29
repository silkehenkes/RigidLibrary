# Silke Henkes, 29.08.18: 
# Created Pebbles class to prepare and execute the pebble game for all situation
# Detangled code parts, contains: 
# - Game set up, including adding frictional or other bonds
# - Pebble game computation
# - Rigid cluster computation
# - Core game functions


# Silke Henkes 11.07.2013: Accelerated version of pebble game and rigid cluster code
#!/usr/bin/python

# Needs a configuration to start with
from Configuration import *

# Copy function as we need to make non-trivial copies of whole pebble configurations
import copy as cp

class Pebbles:
        
        ## ==== Pebble constructor: give it a configuration, the (game10,game20) type of game you want to play, and any modifiers you want to apply
        # (none defined so far, but e.g. adding random bonds, or adding 2nd neighbours)
        def __init__(self,conf, game10,game20,modifier='nothing',verbose0=False):
            self.N=conf.N
            self.game1=game10
            self.game2=game20
            self.verbose=verbose0
            # let's figure out the type of game:
            if modifier=='nothing':
                print "We are running a (" + str(self.game1)+ "," + str(self.game2) + ") game."
                # that's a (x,3) or (x,3) game 
                if (self.game2==3) or (self.game2==2):
                    if (self.game2==2):
                        print "Warning: (2,x) pebble game, using common pebble game code, but results may vary"
                    # if we have a (2,3) game, simply use the contacts as is for the pebble game
                    if (self.game1==2):
                        self.ncon2=conf.ncon
                        self.Ifull=conf.I
                        self.Jfull=conf.J
                    # if we have a (3,3) pebble game, add the frictional bonds as default behaviour
                    elif (self.game1==3):
                        self.AddFrictionalBonds(conf)
                else:
                    print "Atypical pebble game, stopping for now as unclear meaning. If you really must, construct a new modifier type to deal with this game and any extra / deleted contacts."
            else:
                print "Modified pebble games (2nd neighbours, or single contacts, or deleted contacts). Not implemented yet, please construct new functions similar to AddFrictionalBonds(conf) to modify contact lists."
                        
         
        ### ===================== Functions to add extra bonds or otherwise modify the contact list =============
        # Frictional bonds based on mobilisation
        def AddFrictionalBonds(self,conf):
            self.ncon2=conf.ncon+np.sum(conf.fullmobi==0)
            print conf.ncon
            print self.ncon2

            self.Ifull=[0]*self.ncon2
            self.Jfull=[0]*self.ncon2
            jj=0
            # add the fully mobilized contacts at the bottom of the list
            # not concurrently, because numbers of contacts are affected
            self.Ifull[:conf.ncon]=conf.I
            self.Jfull[:conf.ncon]=conf.J
            for k in range (conf.ncon):
                    if (conf.fullmobi[k]==0):
                            self.Ifull[conf.ncon+jj]=conf.I[k]
                            self.Jfull[conf.ncon+jj]=conf.J[k]
                            jj+=1
            print("Sizes of full contact objects:" + str(len(self.Ifull)))
            # Reverse connection object: necessary for marking both rigid clusters and overconstrained regions
            self.conmat=[[] for _ in range(self.N)]
            for k in range(self.ncon2):
                    i=int(self.Ifull[k])
                    j=int(self.Jfull[k])
                    self.conmat[i].append(k)
                    self.conmat[j].append(k)
                    if self.verbose:
                            print 'Contact ' + str(k) + ' with particles ' + str(i) + ' ' + str(j)
        
        # Here go further methods to either add bonds, remove bonds or modify the connectivity
        
        def AddBoundary(self):
                return 0
                
        def AddSecondNeighbour(self):
                return 0
        
	def AddRandombonds(self):
		return 0
         
        ############  Kuang's AddSomeContacts, need to be modified to be compatible with class structure
        #def AddSomeContacts(Contacts,coordinates,size,percentage):
	#new_length=float(len(Contacts[:,0]))*(1+0.01*percentage)
	#Contacts_new=np.zeros((int(new_length),len(Contacts[0,:])))
	#double_bonds=0.0
	#for k in range(len(Contacts[:,0])):
		#if Contacts[k,2]==0:
			#double_bonds+=1.0
	#fraction_double_bonds=double_bonds/float(len(Contacts[:,0]))
	#Num_add_d_bonds=int(fraction_double_bonds*new_length-double_bonds)
	#New_I=[]
	#New_J=[]
	#for i in range(size-1):
		#for j in range(i+1,size):
			#if abs(coordinates[i,0]-coordinates[j,0])<1.2*(coordinates[i,2]+coordinates[j,2]) and abs(coordinates[i,1]-coordinates[j,1])<1.2*(coordinates[i,2]+coordinates[j,2]):
				#for k in range(len(Contacts[:,0])):
					#if Contacts[k,0]==i and Contacts[k,1]==j:
						#break
					#else:
						#if Contacts[k,0]==j and Contacts[k,1]==i:
							#break
						#else:
							#New_I.append(i)
							#New_J.append(j)
							#break
	#for k in range(len(Contacts_new[:,0])):
		#if k<len(Contacts[:,0]):
			#Contacts_new[k,:]=Contacts[k,:]
		#else:
			#x=random.random()
			#temp=int(len(New_I)*x)
			#Contacts_new[k,0]=New_I[temp]
			#Contacts_new[k,1]=New_J[temp]
			#if k<(len(Contacts[:,0])+Num_add_d_bonds):
				#Contacts_new[k,2]=0
				#Contacts_new[k,3]=0.1
				#Contacts_new[k,4]=1.0
			#else:
				#Contacts_new[k,2]=1
				#Contacts_new[k,3]=1.0
				#Contacts_new[k,4]=1.0
	#return Contacts_new
		
	### ================================ The pebble game ==============================================
	def play_game(self):
		# pebbles, defined on particles
		self.pebbles=-1*np.ones((self.N,self.game1))
		self.ptr=-1*np.ones(self.ncon2)
		# tricky: the original [[]]*self.N would create N times the *same* empty list. 
		# overconstrained regions
		self.Overcon=[ False for _ in range(self.ncon2)]
		marked = [False for _ in range(self.N)]
		
		# run pebble game 
		self.fail=0
		for k in range(self.ncon2):
			i=int(self.Ifull[k])
			j=int(self.Jfull[k])
			pebblescopy=cp.copy(self.pebbles)
			for l in range(self.game2+1):
				result=self.enlarge_cover(pebblescopy,i,j)
				#if (self.verbose):
					#if result==1:
						#print "No pebble found in search at contact" + str(k)
			if result==0:
				result=self.enlarge_cover(self.pebbles,i,j)
				self.ptr[k]=0
			else:
				self.enlarge_over(pebblescopy,marked,i,j)
				#print('contact ' + str(k) + ' couldn\'t be added') 
				if (self.verbose):
					print('contact ' + str(k) + ' couldn\'t be added') 
				self.fail+=1
		## print out number of redundant bonds and number of free pebbles
		# number of free pebbles = number of -1 still in pebbles. The other ones contain the contact they're placed on
		# Globally mark the overconstrained regions
		self.over_path(marked)
		self.freepeb=np.sum(self.pebbles==-1)
		print('Statistics:')
		print('Number of free pebbles: ' + str(self.freepeb))
		print('Number of failed contacts: ' + str(self.fail))
		

	# ============================= Rigid clusters =======================================	
	# Problem was the wrong percolation algorithm. Following Jacobs and Hendrickson, only one round is necessary
	# since by construction, every rigid site is found in the pebble search eventually (OR NOT?)
	# Following algorithm on p. 357-358
	def rigid_cluster(self):
		# initialize cluster index and list of lists for particle cluster labels
		self.cidx=-1
		self.maxidx=0
		lenx=0
		leny=0
		frac=0
		fracmax=0
		self.pcluster = [[] for _ in range(self.N)]
		clusterall=[]
		clusterallBonds=[]
		clusteridx=[]
		BigCluster=0
		# Label rattlers as their own cluster
		# nope, no contacts there, and hinges on a technicality. With current algorithm, rattlers will not appear in any cluster
		# nasty to code, too
		self.EMPTY=-1
		self.cluster=self.EMPTY*np.ones(self.ncon2)
		# loop over all contacts which haven't already been marked rigid (or been done)
		# Need the inverse contact matrix, else there is a search necessary every time ...
		
		for k in range(self.ncon2):
		#for k in range(10):
 
			if (self.cluster[k]==-1 and self.ptr[k]==0):
				# Check if k is a double bond
				i=int(self.Ifull[k])
				j=int(self.Jfull[k])
				bonds=[val for val in self.conmat[i] if val in self.conmat[j]]
				stopthis=False
				if len(bonds)>1:
					if self.verbose:
						print "Checking over double bonds" + str(bonds)
					k2=bonds[0]
					if k2==k:
						k2=bonds[1]
					# If the other bond is already part of a rigid cluster ... this one must be too.
					# Abort there
					if self.cluster[k2]!=-1:
						stopthis=True
						self.cluster[k]=self.cluster[k2]
						print "Double bond " + str(k) + " neighbour of " + str(k2) + " was also labeled cluster " + str(self.cluster[k])
						print "Add one to cluster length " + str(k)
						clusterall[self.cluster[k2]]+=1
				if stopthis==False:
					if self.verbose:
						print('contact ' + str(k))
					# We label particles (not bonds!) rigid and floppy based on whether a pebble search fails or succeeds,
					# the whole thing recursively, until either the rigid area is surrounded by floppy sites, or no particles
					# remain accessible through the contact network.
					# A rigid bond is a bond between two particles marked rigid. That includes double bonds where one has not been added
					# during the pebble search.
					# The original bond is rigid if both the original sites are rigid
					# We start the test with the two original sites, and proceed only if neither yields a pebble
					self.cidx+=1
					if self.verbose:
						print('Rigid cluster label # ' + str(self.cidx))
					marked = [False for _ in range(self.N)]
					rigid = [False for _ in range(self.N)]
					# gather 3 pebbles at the site
					# This should always be possible
					i=int(self.Ifull[k])
					j=int(self.Jfull[k])
					pebblescopy=cp.copy(self.pebbles)
					for l in range(self.game2):
						result=self.enlarge_cover(pebblescopy,i,j)
						if result==1:
							print("Error: no " + str(self.game2) + "global DOF pebbles in rigidity search at contact" + str(k))
					
					# First check the two adjoining particles
					done=True
					rig0=self.check_neighbors(i,pebblescopy,marked,rigid)
					rig1=self.check_neighbors(j,pebblescopy,marked,rigid)
					if self.verbose:
						print(rig0)
						print(rig1)
					# only if for both, at least one neighbor is rigid
					if (rig0 and rig1):
						done=False
						rigid[i]=rig0
						rigid[j]=rig1
					# check recursively the neighbors of all rigid sites, until either all those neighbors are floppy
					# or no more sites remain that are unmarked
					while (not done):
						done=True
						if self.verbose:
							print('Entering another rigid loop looking level')
						for i in np.nonzero(rigid)[0]:
							rigtest=self.check_neighbors(i,pebblescopy,marked,rigid)
							if rigtest==True:
								done=False
					cluslen=len(np.nonzero(rigid)[0])
					if (cluslen>0):
						print('Cluster length ' + str(cluslen))
						ctry=self.rig_path(marked,rigid,True)
						if ctry!=self.cidx:
							print "WARNING: Found a (part) duplicate cluster. Merge length here with cluster " + str(ctry)
							cidxstore=self.cidx
							self.cidx=ctry
							ctry=self.rig_path(marked,rigid,False)
							self.cidx=cidxstore
						else:
							if self.verbose:
								print('Found new rigid cluster! # ' + str(self.cidx))	
						cluslenBonds=len([val for val in self.cluster if val==ctry])
						clusterall.append(cluslen)
						clusterallBonds.append(cluslenBonds)
						clusteridx.append(ctry)
						if (cluslen>BigCluster):
							# The biggest cluster and its index
							BigCluster=cluslen
							self.maxidx=ctry
					else:
						self.cidx-=1
                print "Cluster sizes (particles)"
                print clusterall
                print "Cluster sizes (bonds)"
                print clusterallBonds
                print "Cluster labels"
                print clusteridx
                # Note: Further rigid cluster analysis has been moved to Analysis class
                return self.cidx, clusterall, clusterallBonds, clusteridx, BigCluster
		
        
	#========================== The pebble game core functions ========================================	
	#======= These are now generic to be compatible with any (m,n) game as defined above ==============
	#================================ DO NOT TOUCH !!! ================================================
	#
	# Checking neighbors for rigidity
	def check_neighbors(self,i,pebblescopy,marked,rigid):	
		rig=False
		links=self.conmat[i]
		for c in links:
			if (self.Ifull[c]==i):
				u=self.Jfull[c]
			else: 
				u=self.Ifull[c]
			if marked[u]==False:
				seen=-1*np.zeros((self.N,))>0 
				path=-1*np.zeros((self.N,), dtype=np.int)
				found=self.find_pebble2(pebblescopy,int(u),seen,path,marked,rigid)
				# if a pebble is found, mark all those sites as floppy
				# stop looking there
				# else mark them rigid
				if found:
					self.mark_path(marked,rigid,path,u,False)
				else:
					self.mark_path(marked,rigid,path,u,True)
					# for the first round: as long as at least one rigid neighbor has been found
					rig=True
		return rig
					
	# based on rearrange pebbles
	# Mark particles either rigid or floppy
	def mark_path(self,marked,rigid,path,i,marktype):
		if (path[i]==-1):
			marked[i]=True
			rigid[i]=marktype
			if self.verbose:
				if marktype:
					print('Marked particle ' + str(i) + ' ' + ' rigid')
				else:
					print('Marked particle ' + str(i) + ' ' + ' floppy')			
		else:
			while (path[i]!=-1):
				l=path[i]
				marked[i]=True
				rigid[i]=marktype
				if self.verbose:
					if marktype:
						print('Marked particle ' + str(i) + ' ' + ' rigid')
					else:
						print('Marked particle ' + str(i) + ' ' + ' floppy')
				i=l		
	# Close cousin to identify overconstrained regions		
	def mark_path_over(self,marked,path,i):
		if (path[i]==-1):
			marked[i]=True	
		else:
			while (path[i]!=-1):
				l=path[i]
				marked[i]=True
				i=l		
	
	# Define the rigid cluster: Particles are added to cluster cidx if they are rigid
	# contacts are added if they are between two rigid sites (that includes the frictional double bonds)
	# Careful: there are cases where an overlapping cluster has been identified through starting from a bond
	# the game hasn't visited for some reason. Throw an error in that case, return to start
	# This has to be done through the bonds, clusters can have more than one label
	def rig_path(self,marked,rigid,tentative):
		for i in np.nonzero(marked)[0]:
			if self.verbose:
				print('marked particle ' + str(i))
			if (rigid[i]):
				links=self.conmat[i]
				for c in links:
					if (self.Ifull[c]==i):
						u=self.Jfull[c]
					else: 
						u=self.Ifull[c]
					if (marked[u] and rigid[u]):
						if tentative:
							if self.cluster[c]>-1 and self.cluster[c]!=self.cidx:
								print "WARNING!! Attempting to relabel a *bond* from " + str(self.cluster[c])
								return self.cluster[c]
						self.cluster[c]=self.cidx
						if self.verbose:
							print('Added bond ' + str(c)+ ' to cluster ' + str(self.cidx))
				if self.verbose:
					print('rigid')
				self.pcluster[i].append(self.cidx)
				if self.verbose:
					print('Added particle ' + str(i) +' to cluster ' + str(self.cidx))
		return self.cidx
				
							
	# Close cousin to identify overconstrained regions
	def over_path(self,marked):
		for i in np.nonzero(marked)[0]:
			links=self.conmat[i]
			for c in links:
				if (self.Ifull[c]==i):
					u=self.Jfull[c]
				else: 
					u=self.Ifull[c]
				if (marked[u]):
					self.Overcon[c]=True
	
	# Main pebble algorithm	
	def find_pebble(self,pbcopy,i,seen,path):
		seen[i]=True
		path[i]=-1
		found=False
		for k in range(self.game1):
                        if (pbcopy[i,k].astype(int)==-1):
                            found=True
		# if there is not free pebble at vertex v
		if not found:
			k=0
			while (not found) and (k<self.game1):
                            j=pbcopy[i,k].astype(int)
                            if (not seen[j]):
                                    path[i]=pbcopy[i,k].astype(int) # just use the neighbor. rewrite as I and J later if useful
                                    found=self.find_pebble(pbcopy,j,seen,path)
                            k+=1
		return found
		
	# Find pebble in the presence of rigid and floppy sites
	def find_pebble2(self,pbcopy,i,seen,path,marked,rigid):
		## print 'looking for a pebble at vertex ' + str(i)
		seen[i]=True
		path[i]=-1
		found=False
		# if there is a free pebble at vertex v
		for k in range(self.game1):
                        if (pbcopy[i,k].astype(int)==-1):
                            ## print 'found a free pebble here at vertex' + str(i)
                            found=True
                if not found:
			# Look at the contacts covered by pebbles in turn. If they haven't been accessed by the search yet
			# look along those paths, recursively, until you find something (or not)
			# For marked sites: if the site is marked rigid, there is no pebble here, stop the subloop. If it's marked floppy, there must be a pebble
			# somewhere along its path, so stop the whole thing with a positive result
			if marked[int(i)]==True:
				if (rigid[int(i)]==False):
					found=True
			else:
                                k=0
                                while (not found) and (k<self.game1):
                                    j=pbcopy[i,k].astype(int)
                                    if (not seen[j]):
                                            path[i]=pbcopy[i,k].astype(int) # just use the neighbor. rewrite as I and J later if useful
                                            found=self.find_pebble2(pbcopy,j,seen,path,marked,rigid)
                                    k+=1
		return found
		
	# Rearranging search paths after pebble is found
	def rearrange_pebbles(self,pbcopy,i,j,path):
		## print 'entering rearrange pebbles; path of i is' + str(path[i])
		if (path[i]==-1):
			#adding the arrow pointing from i back to j
			## print 'found a pebble right away; adding it in situ: # ' + str(i) + ' to '+ str(j)
			if (pbcopy[i,0]==-1):
				pbcopy[i,0]=j
			else:
				if (pbcopy[i,1]==-1):
					pbcopy[i,1]=j
				else:
					pbcopy[i,2]=j  
		else:
			#flip the arrow from l=pth(i) to j, then head further downstream doing the same thing until the end is reached
			while (path[i]!=-1):
				l=path[i]
				## print 'flipping arrow at ' + str(i) + ' from ' + str(l) + ' to ' + str(j)
				if (l==pbcopy[i,0]):
					pbcopy[i,0]=j # remove pebble and replace on ij contact, i.e. change label from l to i
				else:
					if (l==pbcopy[i,1]):
						pbcopy[i,1]=j
					else:
						pbcopy[i,2]=j;
				j=i
				i=l
			#and now that we've reached the end:    do the same as before
			## print 'finally found a pebble adding arrow # ' + str(i) + ' to '+ str(j)
			if (pbcopy[i,0]==-1):
				pbcopy[i,0]=j
			else:
				if (pbcopy[i,1]==-1):
					pbcopy[i,1]=j
				else:
					pbcopy[i,2]=j


	# attempt adding an edge. The contact info is in pbcopy
	def enlarge_cover(self,pbcopy,i,j):
		seen=-1*np.zeros((self.N,))>0 
		path=-1*np.zeros((self.N,), dtype=np.int)
		found=self.find_pebble(pbcopy,i,seen,path)
		## print 'looked for pebble along i; result is: ' + str(found)
		if found:
			self.rearrange_pebbles(pbcopy,i,j,path)
			## print 'rearranged pebbles'
			return 0 # success
		if (not seen[j]):
			found=self.find_pebble(pbcopy,j,seen,path)
			## print 'looked for pebble along j; result is: ' + str(found)
			if found:
				self.rearrange_pebbles(pbcopy,j,i,path)
				## print 'rearranged pebbles'
				return 0 # success
		## print 'pebble search failed!'
		return 1 # failure

	# Cousin to mark overconstrained
	def enlarge_over(self,pbcopy,marked,i,j):
		seen=-1*np.zeros((self.N,))>0 
		path=-1*np.zeros((self.N,), dtype=np.int)
		found=self.find_pebble(pbcopy,i,seen,path)
		## print 'looked for pebble along i; result is: ' + str(found)
		if (not found):
			self.mark_path_over(marked,path,i)
			
		if (not seen[j]):
			found=self.find_pebble(pbcopy,j,seen,path)
			## print 'looked for pebble along j; result is: ' + str(found)
			if (not found):
				self.mark_path_over(marked,path,j)
				
				
				



