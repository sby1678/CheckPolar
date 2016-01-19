#!/usr/bin/python

import numpy as np
import fractions
import math as m
import os
import decimal as dc
import sys
import json
import time

#
# Material Builder 
# Author: Jakub Kaminski, UCLA 2015
#	Modified by Boya Song, UCLA Winter 2015
#
# Based on Interface Builder
#
if len(sys.argv) <=1:
	print "Usage:\n%s optionsFile"%sys.argv[0]
	exit()

inputFile = sys.argv[1]


def readInput(inputFile):

	file = open(inputFile,'r')
	line = file.readline()
	line = line.split()
	folderPath = line[1]

    line = file.readline()
    line = line.split()
    z_digit = int(line[1])

    line = file.readline()
    line = line.split()
    dis_tol_rate = float(line[1])

    return folderPath, z_digit, dis_tol_rate




def checkPolar(self,atom_coor,atomLabels, z_digit=4, dis_tol_rate=0.01, max_compare = float("Inf")):
    # Description of the input variables:
    # - self: not used here, but required by Python
    # - atom_coor: of the atoms, numpy.array(Natoms,3)
    # - atomLabels: numpy.array with integers denoting atom types, 
    #   i.e [0,1,1,0,2,...,natoms]. The order of the atoms is the same 
    #   as in positions array.
    # - atomTypes: dictionary that allows to decode entries in the atomLabels
    #   in terms of real chemical species. 
    #   Example:
    #    atomTypes = {0: 'Ga', 1: 'As', 2: 'H'}
    #    which means that integer 0 corresponds to "Ga", integer 1 to "As" and 2 to "H"
    #   Usage:
    #    find what is the atom type of the 3rd atom in the structure:
    #    atomLabels = [0,1,1,0,2]
    #    atom = atomLabels[2]  # remeberin Python we count from 0, so 3rd atom is 2nd in the structure
    #    type = atomTypes[atom]
    #    In this case atom will be set to "1", and type to "As"
    #
    # - vecX: lattice vector X, numpy.array[x1, x2, x3]
    # - vecY: lattice vector Y, numpy.array[y1, y2, y3]

    #  Start implementation:
    #  Do lots of cool stuff here

    #  Return values are 
    #  polar = {True,False} - is structure polar or not
    #  periodicity = 0  - double number saying what is the periodicity of the strucure in angstrom
    # minDiffRate - double number saying what is the minimal distance rate found during the chekcing procedure (especially usful for understanding the polar structures)
    
    polar = False
    periodicity = float("nan")
    minDiffRate = -float("Inf")

    vecX = self.vecSubR[0,range(2)];
    vecY = self.vecSubR[1,range(2)];

    # first 180 primes
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069];

    nAtoms = len(atom_type);
    # represent the atom types by prime numbers
    ele_n = copy.deepcopy(atomLabels); 
    for i in xrange(min(atomLabels), max(atomLabels+1)):
        ele_n = np.where(atomLabels == i, primes[i], ele_n);
    ele_n = outer(atomLabels,atomLabels);
    ele_n[range(nAtoms), range(nAtoms)] = 0;
    
    # build the distance matrix (size nAtoms*nAtoms, entry i,j represent for the distance between i-th atom and j-th atom)
    atom_dist = np.zeros([nAtoms, nAtoms], dtype=float);
    for i in xrange(0, nAtoms):
        for j in xrange(i,nAtoms):
            atom_dist[i,j] = compute_min_dist(atom_coor[i, 0:2] - atom_coor[j, 0:2], vecX, vecY);
            atom_dist[i,j] = atom_dist[i,j]+(atom_coor[i, 2] - atom_coor[j, 2])^2;
            atom_dist[j,i] = atom_dist[i,j]

    # round the z-coordinate
    atom_coor = np.floor(atom_coor[:, 2]*(10^z_digit)) /(10**z_digit);
    # atoms with same cut_z coord are considered on the same "surface"
    # return_inverse=False: the resut is ordered from small to large
    surface_z= np.unique(atom_coor[:, 2], return_inverse=False);  
    surface_n=len(surface_z, 1); # number of surfaces
    init_z = surface_z[-1]; # the first surface to consider
    firstRound = False;
    forDelete = np.zeros([nAtoms, 1], dtype=bool);
    while(surface_n>=2 and not((abs(surface_z[0]-surface_z[-1])<abs(init_z - surface_z[-1])) and not(firstRound))): # the first round non-polar must be found within the upper half of the structure
        if max_compare> (surface_n/2): # take advange of integer division
            max_compare=(surface_n/2);
        polar=False;
        doCompare = True;

        if firstRound:
            if abs(init_z -  surface_z[-1]) < perodicity_lower_bound:
                doCompare = False;

        if doCompare:
            for nSurface in xrange(0,max_compare):
            #each row represents for an atom
            # the indices of atoms on the upper surface & lower surface
                u_lidx = atom_coor[:, 2]==surface_z[nSurface];
                u_lidx = u_lidx [~forDelete];
                l_lidx =atom_coor[:, 2]==surface_z[-nSurface-1];
                l_lidx = l_lidx[~forDelete];
            
            
                # if the number of atoms are differnent
                if sum(u_lidx) != sum(l_lidx):
                    polar=True;
                    break;
                   
                nAtomPerLayer = sum(u_lidx);
                data_upper_type= ele_n[u_lidx,:];
                data_lower_type= ele_n[l_lidx,:];
                data_upper_dist= atom_dist[u_lidx,:];
                data_lower_dist= atom_dist[l_lidx,:];
            
                # for each atom, sort the distance from this atom to the others (small to large)
                sort1idx_u = np.argsort(data_upper_dist, axis = 1);
                sort1idx_l = np.argsort(data_lower_dist,axis = 1);
                for i in xrange(nAtomPerLayer):
                    data_upper_type[i,:] = data_upper_type[i,sort1idx_u[i,:]];
                    data_lower_type[i,:] = data_lower_type[i,sort1idx_l[i,:]];
                    data_upper_dist[i,:] = data_upper_dist[i,sort1idx_u[i,:]];
                    data_lower_dist[i,:] = data_lower_dist[i,sort1idx_l[i,:]];
                # for each atom, sort the type from this atom to the others (small to large)
                sort2idx_u = np.argsort(data_upper_type, axis = 1);
                sort2idx_l = np.argsort(data_lower_type,axis = 1);
                for i in xrange(nAtomPerLayer):
                    data_upper_type[i,:] = data_upper_type[i,sort2idx_u[i,:]];
                    data_lower_type[i,:] = data_lower_type[i,sort2idx_l[i,:]];
                    data_upper_dist[i,:] = data_upper_dist[i,sort2idx_u[i,:]];
                    data_lower_dist[i,:] = data_lower_dist[i,sort2idx_l[i,:]];
            
                dist_diff = np.zeros([nAtomPerLayer,nAtomPerLayer], dtype = float); # rate of difference on distance
                type_ok = np.zeros([nAtomPerLayer,nAtomPerLayer], dtype = bool); # true if the type items are matching between a upper-atom and a lower-atom
                for idx_upper in xrange(nAtomPerLayer):
                    for idx_lower in xrange(nAtomPerLayer):
                         dist_diff[idx_upper,idx_lower] = max(abs(np.divide(data_upper_dist[idx_upper,:]-data_lower_dist[idx_lower,:], data_upper_dist[idx_upper,:]+data_lower_dist[idx_lower,:])));
                         type_ok[idx_upper,idx_lower]= all(data_upper_type[idx_upper,:]==data_lower_type[idx_lower,:]);
                minDiffRate = np.min(np.min(dist_diff), minDiffRate)

                g = networkx.to_networkx_graph(dist_diff<=dist_tol_rate & type_ok); # 
                # find the maximal matching of graph g
                if len(networkx.maximal_matching(g))< nAtomPerLayer: 
                    polar=True;
                    break;
    
        if not(polar) and not(firstRound): # first round NonPolar
            firstRound = True
            layer_thickness = 0
            minNP_z = surface_z[-1]
            perodicity_lower_bound = abs(init_z - minNP_z); # the perodicity should be larger than this
        else:
            if not(polar) and firstRound and (abs(minNP_z[0] - surface_z[-1])>perodicity_lower_bound): # the second round non-polar
                np.append(minNP_z, surface_z[-1]);
                layer_thickness = layer_thickness + 1;
                z_thickness = minNP_z[0]-minNP_z[-1];
    #             vz(3) = z_thickness;
    #             isNP = ismember(atom_coor(:,3), minNP_z);
    #             minNPStruc= atom_coor(isNP,:);
    #             minNPStruc(:,3)= minNPStruc(:,3) - minNP_z(end);
    #             atom_type = atom_type(isNP);
    #             atom_cell=[num2cell(minNPStruc) atom_type];
                
    #             # create a new txt file that contains the maximal non-polar structure
    #             isWin = ~isempty(strfind(computer, 'PCWIN'));
    #             [pathstr,name,ext] = fileparts(filepath);
    #             mkdir(pathstr,'minNonPolar5');
    #             if isWin
    #                 pathstr=strcat(pathstr, '\minNonPolar5');
    #             else
    #                 pathstr=strcat(pathstr, '/minNonPolar5');
    #             end
    #             oldpath=cd(pathstr);
    #             fileID=fopen([name '-minNonPolar5' ext], 'w+');
    #             formatSpec1='lattice_vector \t %f \t %f \t %f \n';
    #             fprintf(fileID,formatSpec1,vx);
    #             fprintf(fileID,formatSpec1,vy);
    #             fprintf(fileID,formatSpec1,vz);
    #             formatSpec2='atom \t %f \t %f \t %f \t %s \n';
                
    #             nrows= size(atom_cell,1);
    #             for row = 1:nrows:
    #                 fprintf(fileID,formatSpec2,atom_cell{row, :});
    #             fclose(fileID);
    #             #movefile([name '-nonpolar' ext], pathstr);
    #             cd(oldpath);
    #            
    #            top_z=minNP_z(1);
                return polar, z_thickness, minDiffRate
    
        if polar and firstRound:
            layer_thickness = layer_thickness + 1;
            np.append(minNP_z,surface_z[-1]);

        data_delete = (atom_coor[:, 2]==surface_z[-1]);
        data_delete = data_delete[~forDelete];
        forDelete= (forDelete | (atom_coor[:, 2]==surface_z[-1]));
        atom_data_dist = atom_data_dist[data_delete,:];
        atom_data_dist = atom_data_dist[:,data_delete];
        atom_data_type = atom_data_type[data_delete,:];
        atom_data_type = atom_data_type[:,data_delete];
        surface_n=surface_n-1;
        surface_z = np.delete(surface_z, -1);

    
    z_thickness=float("nan");
    layer_thickness = float("nan");
    #top_z = float("nan");
    #perodicity_lower_bound = float("nan");

    return polar, z_thickness, minDiffRate




#########################################
#										#
#										#
#          Start program 				#
#										#
#										#
#########################################

#start timer
tStart = time.time()
atomTyp = {}
# Read input
folderPath, z_digit, dis_tol_rate = readInput(inputFile)

# # create a list of Miller indices
# if not useMillerList:
# 	MillerList = createMillerList(maxMillerInd)
# else:
# 	MillerList = MillerIndList
# 	print MillerList

#Read CIF file
print
print "********************************************************"
print "Reading structure from .cif file"
idMat,transM,atoms,positions,atomTyp=ReadCIF(subCIF,atomTyp)

# Construt big bulk material that will be reused in all calculations
print "Construction big bulk structure... This might take time."
bigBulk = Surface(transM,positions,atoms,atomTyp,np.array((0,0,0)))
bigBulk.bulkNEW(40)
print "Bulk structure complete"
print "********************************************************"
print 

primFailed = [] # stuctures for which cound find primtive vectors
notExist = [] # stuctures for which given orientation does not exist
# Start the loop on Miller indices 
isPolar = []; # if the structure is polar or not
perodicity  = []; # if polar, report the perodicity in angstrom; otherwise report NaN
minDiffRate = []; # the minimal difference rate
polarFilename=subCIF.split(".")[0]+"-poalrity.txt"
file = open(polarFilename, 'w+')
file.write("Filename\t\tisPolar\t\tperodicity\t\tminDiffRate\n")
file.close()

for subMillerString in MillerList:
	#MATERIAL 1
	print
	print "*****************************************"
	print "Constructing bulk Substrate"
	print "Orientation: %s"%subMillerString
	#idMat,transM,atoms,positions,atomtyp=ReadCIF(subCIF)
	Miller = getMillerFromString(subMillerString)
	#nbulk = max(Miller)*2 # use the max Miller index to ensure non-empty surface

	Sub = Surface(transM,positions,atoms,atomTyp,Miller)
#	Sub.bulk(nbulk)
	Sub.positions = bigBulk.positions.copy()
	Sub.atoms = bigBulk.atoms.copy()
	Sub.positionsSmall = bigBulk.positionsSmall.copy()
	Sub.atomsSmall = bigBulk.atomsSmall.copy()

	Sub.construct()
	Sub.plane()
	if not Sub.exists: 
		print "!!! Orientation does not exists!!!"
		print "!!! Proceeding to next one !!!"
		print "*****************************************"
		notExist.append(subMillerString)
		continue # if given plane does not exists, continue to next one

	Sub.initpvecNEW2()
	if not Sub.exists: 
		print "!!! Failed to find primitive vectors !!!"
		print "!!! Proceeding to next orientation  !!!"
		print "*****************************************"
		primFailed.append(subMillerString)
		continue # if given plane does not exists, continue to next one

	Sub.primitivecell()

	vecsS = Superlattice(Sub.a,Sub.b,1) # multiplication of unit cells
	# vecsS has set of vecors that span on the lattice. They give the same area,
	# but the two vectors might be of different lengths. Find such a pair, for
	# which length is similar, so we have the most "square" lattice
	maxRatio = 0
	ii = 0 #counter
	for vec in vecsS:
		normA = np.linalg.norm(vec[0])
		normB = np.linalg.norm(vec[1])
		ratio = normA/normB
		if normA >= normB : ratio = normB/normA
		if ratio > maxRatio:
			maxRatio = ratio
			vecIndex = ii
		ii += 1

	#Construct big planes for Substrate and Deposit
	#print "CREATING BULK - this may take time"
	#nbulk2 = max(Miller)*2
	#Sub.bulk(nbulk2)
	Sub.positions = bigBulk.positions.copy()
	Sub.atoms = bigBulk.atoms
	Sub.construct()
	Sub.plane()

	print "CREATING STRUCTURE"
	iface = Interface(vecsS[vecIndex],Sub,nL)

	print "OUTPUTTING COORDINATES"
	# Output in cartesian coordiantes
	millerLab= subMillerString.split()
	strufilename=subCIF.split(".")[0]+"-%s%s%s"%(millerLab[0],
			                            millerLab[1],
						    millerLab[2],)

	"""
	file = open(strufilename+".xyz",'w')
	natoms = len(iface.SSurfPos)
	file.write("%i\n\n"%natoms)
	ii = 0
	for i in iface.SSurfPos:
		file.write("%5s  %12.6f  %12.6f  %12.6f\n"%(iface.SSurfAtm[ii], i[0],i[1],i[2]))
		ii+= 1
	file.close()
	"""

	# Output in FHI-AIMS .in format
	file = open(strufilename+".in",'w')
	file.write("lattice_vector   %12.6f   %12.6f   %12.6f\n"%(iface.vecSubR[0][0],iface.vecSubR[0][1],iface.vecSubR[0][2]))
	file.write("lattice_vector   %12.6f   %12.6f   %12.6f\n"%(iface.vecSubR[1][0],iface.vecSubR[1][1],iface.vecSubR[1][2]))
	file.write("lattice_vector   %12.6f   %12.6f   %12.6f\n"%(0.0, 0.0, nL))
	ii = 0
	atom_coor = np.zeros([len(iface.SSurfPos), 3])
    atomLabels = np.zeros([len(iface.SSurfPos), 1])
    for i in iface.SSurfPos:
        file.write("atom %12.6f  %12.6f  %12.6f  %5s\n"%(i[0],i[1],i[2],atomTyp[iface.SSurfAtm[ii]]))
        atom_coor[ii,:] = i[range(3)];
        atomLabels[ii] = iface.SSurfAtm[ii];
        ii+= 1
    file.close()
    # check the polarity here
    if checkPolarity:
        print "Check the poalrity"
        polar, ped, minDiff = iface.checkPolar(atom_coor,atomLabels, z_digit, dis_tol_rate);
        isPolar.append(polar);
        perodicity.append(ped);
        minDiffRate.append(minDiff);
        if polar:
            print "The structure is polar with perodicity %s angstrom."%ped  
        else:
            print "The structure is non-polar, with a minimal distance difference rate %s."%minDiff 
        file = open(polarFilename, 'a')
        file.write(strufilename+"\t\t%d\t\t%f\t\t%f\n"%(polar,ped,minDiff))
        file.close()

	# end of the loop on Miller indices 

# nlayers? #
	print "*****************************************"

# Stop timer
tStop = time.time()

runtime = calcRuntime(tStart,tStop)

#Output statistics
nStruc = len(MillerList)
nPrimFailed = len(primFailed)
nNotExist = len(notExist)
nCreated = nStruc - nPrimFailed - nNotExist
file = open('stats.out','w')
file.write("Max. Miller Index: %i\n"%maxMillerInd)
file.write("Total number of possible orientations: %i\n"%nStruc)
file.write("Number created orientations: %i\n"%nCreated)
file.write("Number of not existing orientations: %i\n"%nNotExist)
for i in notExist:
	file.write("%s\n"%i)

file.write("\nNumber of orientations where primitive vectors failed: %i\n"%nPrimFailed)
for i in primFailed:
	file.write("%s\n"%i)
file.write("\nRuntime: %i min %i sec\n"%(runtime[0],runtime[1]))
file.close()

# End of program