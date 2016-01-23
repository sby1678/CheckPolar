#!/usr/bin/python

import numpy as np
import fractions
import math as m
import os
#import decimal as dc
import time
import networkx
import platform
import sys
import copy
import pdb # for debug 

# Check polarity of .in files
# Author: Jakub Kaminski, UCLA 2015
#    Modified by Boya Song, UCLA Winter 2016
#
# Based on MatBuilder-V1


if len(sys.argv) <=1:
    print "Usage:\n%s optionsFile"%sys.argv[0]
    exit()

inputFile = sys.argv[1]
systemInfo = platform.system();
if systemInfo == "Windows":
    isWin = True
elif systemInfo == "Linux":
    isLinux = True
else:
    print "The program is designed to run on Windows or Linux system"
    exit();




def readInput(inputFile):

    file = open(inputFile,'r')
    line = file.readline()
    line = line.split("\t");
    line = filter(None, line)
    folderPath = line[1]

    line = file.readline()
    line = line.split()
    z_digit = int(line[1])

    line = file.readline()
    line = line.split()
    dis_tol_rate = float(line[1])

    line = file.readline()
    line = line.split()
    checkAllLayer = bool(line[1])

    return folderPath, z_digit, dis_tol_rate, checkAllLayer

def compute_min_dist(vec_v, vec_a, vec_b):
# algorithm provided by Craig Schroeder
# return: best_dist, best_x, best_y
    a = vec_a[0]*vec_a[0]+ vec_a[1]*vec_a[1];
    c = vec_b[0]*vec_b[0] + vec_b[1]*vec_b[1];
    swp = False;
    if a < c:
        temp = vec_a; 
        vec_a = vec_b; 
        vec_b = temp;
        temp = a; 
        a=c; 
        c=temp;
        swp = True;

    f = vec_v[0]*vec_v[0] + vec_v[1]*vec_v[1];
    b = 2 * (vec_a[0]*vec_b[0] + vec_a[1]*vec_b[1]);
    d = 2* (vec_a[0]* vec_v[0] + vec_a[1]*vec_v[1]);
    e  = 2* (vec_b[0]* vec_v[0]+ vec_b[1]*vec_v[1]);
    x_c = (b*e-2*c*d)/(-b*b+4*c*a);
    x0 = np.round(x_c); 
    best_x=x0;
    y0 = np.floor((b*x0+e)/(-2.*c));
    y1 = y0+1;
    f0 = a*x0*x0 + b*x0*y0 + c*y0*y0 + d*x0 + e*y0 + f;
    f1 = a*x0*x0 + b*x0*y1 + c*y1*y1 + d*x0 + e*y1 + f;
    if f0<f1:
        best_y = y0;
        best_dist = f0;
    else:
        best_y = y1;
        best_dist = f1;

    C=4*a*c-b*b;
    D=4*d*c-2*b*e;
    A=C*x0*x0+D*x0-e*e+4*f*c;
    B=2*C*x0+D;

    j=1;
    while A+B*j+C*j*j<=4*best_dist*c:
        y0=np.floor((b*(x0+j)+e)/(-2*c)); 
        y1=y0+1;
        if(2*c*y0+b*x0+b*j+e+c>0):
            y1=y0;
        f1=a*(x0+j)*(x0+j) + b*(x0+j)*y1 + c*y1*y1 + d*(x0+j) + e*y1 + f;
        if f1 < best_dist:
            best_dist = f1;
            best_y = y1;
            best_x = x0+j;
        j = j+1;


    j= 1;
    while A-B*j+C*j*j<=4*best_dist*c:
        y0=np.floor((b*(x0-j)+e)/(-2.*c)); 
        y1=y0+1;
        if(2*c*y0+b*x0-b*j+e+c>0):
             y1=y0
        f1=a*(x0-j)*(x0-j) + b*(x0-j)*y1 + c*y1*y1 + d*(x0-j) + e*y1 + f;
        if(f1 < best_dist):
            best_dist = f1;
            best_y = y1;
            best_x = x0-j;
        j = j+1;
    if swp:
        temp = best_x; 
        best_x = best_y; 
        best_y = temp;
        
    # return best_dist, best_x, best_y
    return best_dist



def checkPolar(atom_coor,atomLabels, vecX, vecY, z_digit=4, dist_tol_rate=0.01, max_compare = float("Inf"), inverse_struc = False):
    # Description of the input variables:
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
    # - vec_x: lattice vector X, numpy.array[x1, x2, x3]
    # - vec_y: lattice vector Y, numpy.array[y1, y2, y3]

    #  Return values are 
    #  polar = {True,False} - is structure polar or not
    #  periodicity = 0  - double number saying what is the periodicity of the structure in angstrom
    #  minDiffRate - double number saying what is the minimal distance rate found during the checking procedure (especially useful for understanding the polar structures)
    
    polar = False
    periodicity = float("nan")
    minDiffRate = 1
    typeDiff = False

    # first 180 primes
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069];

    nAtoms = len(atomLabels);
    # represent the atom types by prime numbers
    ele_n = np.empty([nAtoms, 1]); 
    ii=0;
    atomLabels = np.array(atomLabels);
    for i in xrange(min(atomLabels), max(atomLabels)+1):
        ele_n[atomLabels == i]= primes[ii];
        ii += 1;
    ele_n = np.outer(ele_n,ele_n);
    ele_n[range(nAtoms), range(nAtoms)] = 0;

    # build the distance matrix (size nAtoms*nAtoms, entry i,j represent for the distance between i-th atom and j-th atom)
    atom_dist = np.empty([nAtoms, nAtoms], dtype=float);
    for i in xrange(0, nAtoms):
        for j in xrange(i,nAtoms):
            atom_dist[i,j] = compute_min_dist(atom_coor[i, 0:2] - atom_coor[j, 0:2], vecX, vecY);
            atom_dist[i,j] = atom_dist[i,j]+(atom_coor[i, 2] - atom_coor[j, 2])**2;
            atom_dist[j,i] = atom_dist[i,j]
    atom_dist[range(nAtoms), range(nAtoms)] = -1; # avoid zero-division in later steps

    # round the z-coordinate
    atom_coor[:, 2] = np.around(atom_coor[:, 2], decimals = z_digit);
    # atoms with same cut_z coord are considered on the same "surface"
    # inverse_struc=False: the resut is ordered from small to large
    if inverse_struc:
        surface_z= np.squeeze(np.unique(atom_coor[:, 2], return_inverse=True))[0]; 
    else:
        surface_z= np.squeeze(np.unique(atom_coor[:, 2], return_inverse=False)); 

    surface_n=len(surface_z); # number of surfaces
    init_z = surface_z[-1]; # the first surface to consider
    firstRound = False;
    forDelete = np.zeros([nAtoms], dtype=bool);

    while(surface_n>=2 and not(abs(surface_z[0]-surface_z[-1])<abs(init_z-surface_z[-1]) and not(firstRound))): # the first round non-polar must be found within the upper half of the structure
        if max_compare> (surface_n/2): # take advange of integer division
            max_compare=(surface_n/2);
        polar=False;
        doCompare = True;

        if firstRound:
            if abs(init_z-surface_z[-1]) < perodicity_lower_bound:
                doCompare = False;

        if doCompare:
            for nSurface in xrange(0,max_compare):
            #each row represents for an atom
            # the indices of atoms on the upper surface & lower surface
                u_lidx = atom_coor[:, 2]==surface_z[nSurface];
                u_lidx = u_lidx [~forDelete];
                l_lidx =atom_coor[:, 2]==surface_z[-nSurface-1];
                l_lidx = l_lidx[~forDelete];
                # if the number of atoms are different
                if sum(u_lidx) != sum(l_lidx):
                    polar=True;
                    break;
                   
                nAtomPerLayer = sum(u_lidx);

                data_upper_dist= atom_dist[u_lidx,:];
                data_lower_dist= atom_dist[l_lidx,:];

                # # can speed up the code if we only want to get polar/non-polar 
                # if np.setdiff1d(data_upper_dist,data_lower_dist):
                #     polar = True;
                #     break;


                data_upper_type= ele_n[u_lidx,:];
                data_lower_type= ele_n[l_lidx,:];

                # # can speed up the code if we only want to get polar/non-polar 
                # if np.setdiff1d(data_upper_type, data_lower_type).size:
                #     polar = True;
                #     break;

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

                match_matrix = (dist_diff<=dist_tol_rate) & type_ok;

                g = networkx.to_networkx_graph(match_matrix); # 
                # find the maximal matching of graph g
                if len(networkx.maximal_matching(g))< nAtomPerLayer: 
                    polar=True;

                    g = networkx.to_networkx_graph(type_ok);
                    if len(networkx.maximal_matching(g))< nAtomPerLayer: 
                        typeDiff = True;
                        minDiffRate = float("nan");
                    else:
                        minDiffRate = np.min([np.min(dist_diff), minDiffRate]);

                    break;
                # END of nSurface loop

            if not(polar) and not(firstRound): # first round NonPolar
                firstRound = True
                layer_thickness = 0
                minNP_z = surface_z[-1]
                perodicity_lower_bound = abs(init_z - minNP_z); # the periodicity should be larger than this
            else:
                if not(polar) and firstRound and (abs(minNP_z[0] - surface_z[-1])>perodicity_lower_bound): # the second round non-polar
                    minNP_z = np.append(minNP_z, surface_z[-1]);
                    layer_thickness = layer_thickness + 1;
                    z_thickness = minNP_z[0]-minNP_z[-1];
                    minDiffRate = float("nan")
                    #  vz(3) = z_thickness;
                    #  isNP = ismember(atom_coor(:,3), minNP_z);
                    #  minNPStruc= atom_coor(isNP,:);
                    #  minNPStruc(:,3)= minNPStruc(:,3) - minNP_z(end);
                    #  atom_type = atom_type(isNP);
                    #  atom_cell=[num2cell(minNPStruc) atom_type];
                    
                     #  # create a new txt file that contains the maximal non-polar structure
                    # isWin = ~isempty(strfind(computer, 'PCWIN'));
                    # [pathstr,name,ext] = fileparts(filepath);
                    # mkdir(pathstr,'minNonPolar5');
                    # if isWin
                    #   pathstr=strcat(pathstr, '\minNonPolar5');
                    # else
                    #    pathstr=strcat(pathstr, '/minNonPolar5');
                    # end
                    # oldpath=cd(pathstr);
                    # fileID=fopen([name '-minNonPolar5' ext], 'w+');
                    # formatSpec1='lattice_vector \t %f \t %f \t %f \n';
                    # fprintf(fileID,formatSpec1,vx);
                    # fprintf(fileID,formatSpec1,vy);
                    # fprintf(fileID,formatSpec1,vz);
                    # formatSpec2='atom \t %f \t %f \t %f \t %s \n';
                    
                    # nrows= size(atom_cell,1);
                    # for row = 1:nrows:
                    #     fprintf(fileID,formatSpec2,atom_cell{row, :});
                    # fclose(fileID);
                    # #movefile([name '-nonpolar' ext], pathstr);
                    # cd(oldpath);
                    
                    # top_z=minNP_z(1);
                    return polar, z_thickness, minDiffRate, typeDiff
        
            if polar and firstRound:
                layer_thickness = layer_thickness + 1;
                minNP_z = np.append(minNP_z,surface_z[-1]);

        data_delete = (atom_coor[:, 2]==surface_z[-1]);
        data_delete = data_delete[~forDelete];
        forDelete= np.squeeze((forDelete | (atom_coor[:, 2]==surface_z[-1])));
        atom_dist = atom_dist[~data_delete,:];
        atom_dist = atom_dist[:,~data_delete];
        ele_n = ele_n[~data_delete,:];
        ele_n = ele_n[:,~data_delete];
        surface_n=surface_n-1;
        surface_z = np.delete(surface_z, -1);

    
    z_thickness=float("nan");
    layer_thickness = float("nan");
    #top_z = float("nan");
    #perodicity_lower_bound = float("nan");

    return polar, z_thickness, minDiffRate, typeDiff

def checkPolar_all(atom_coor,atomLabels, vecX, vecY, z_digit=4, dist_tol_rate=0.01, max_compare = float("Inf"), inverse_struc = False):
    # Description of the input variables:
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
    # - vec_x: lattice vector X, numpy.array[x1, x2, x3]
    # - vec_y: lattice vector Y, numpy.array[y1, y2, y3]

    # first 180 primes
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069];

    nAtoms = len(atomLabels);
    # represent the atom types by prime numbers
    ele_n = np.empty([nAtoms, 1]); 
    ii=0;
    atomLabels = np.array(atomLabels);
    for i in xrange(min(atomLabels), max(atomLabels)+1):
        ele_n[atomLabels == i]= primes[ii];
        ii += 1;
    ele_n = np.outer(ele_n,ele_n);
    ele_n[range(nAtoms), range(nAtoms)] = 0;

    # build the distance matrix (size nAtoms*nAtoms, entry i,j represent for the distance between i-th atom and j-th atom)
    atom_dist = np.empty([nAtoms, nAtoms], dtype=float);
    for i in xrange(0, nAtoms):
        for j in xrange(i,nAtoms):
            atom_dist[i,j] = compute_min_dist(atom_coor[i, 0:2] - atom_coor[j, 0:2], vecX, vecY);
            atom_dist[i,j] = atom_dist[i,j]+(atom_coor[i, 2] - atom_coor[j, 2])**2;
            atom_dist[j,i] = atom_dist[i,j]
    atom_dist[range(nAtoms), range(nAtoms)] = -1; # avoid zero-division in later steps

    # round the z-coordinate
    atom_coor[:, 2] = np.around(atom_coor[:, 2], decimals = z_digit);
    # atoms with same cut_z coord are considered on the same "surface"
    # inverse_struc=False: the resut is ordered from small to large
    if inverse_struc:
        surface_z= np.squeeze(np.unique(atom_coor[:, 2], return_inverse=True))[0]; 
    else:
        surface_z= np.squeeze(np.unique(atom_coor[:, 2], return_inverse=False)); 

    surface_n=len(surface_z); # number of surfaces
    surface_thickness = surface_n;
    init_z = surface_z[-1]; # the first surface to consider
    firstRound = False;
    forDelete = np.zeros([nAtoms], dtype=bool);
    is_Polar = np.zeros([surface_n], dtype=bool);
    z_coords = surface_z;
    minDiffRate = np.empty([surface_n])
    typeDiff = np.zeros([surface_n], dtype=bool);

    for current_surface in xrange(surface_thickness-1): 
        if max_compare> (surface_n/2): # take advange of integer division
            max_compare=(surface_n/2);

        for nSurface in xrange(0,max_compare):
        #each row represents for an atom
        # the indices of atoms on the upper surface & lower surface
            u_lidx = atom_coor[:, 2]==surface_z[nSurface];
            l_lidx =atom_coor[:, 2]==surface_z[-nSurface-1];
            # if the number of atoms are differnent
            if sum(u_lidx) != sum(l_lidx):
                is_Polar[current_surface]=True;
                break;
               
            nAtomPerLayer = sum(u_lidx);

            data_upper_dist= atom_dist[u_lidx,:];
            data_lower_dist= atom_dist[l_lidx,:];

            # # can speed up the code if we only want to get polar/non-polar 
            # if np.setdiff1d(data_upper_dist,data_lower_dist):
            #     polar = True;
            #     break;


            data_upper_type= ele_n[u_lidx,:];
            data_lower_type= ele_n[l_lidx,:];

            # # can speed up the code if we only want to get polar/non-polar 
            # if np.setdiff1d(data_upper_type, data_lower_type).size:
            #     polar = True;
            #     break;

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

            match_matrix = (dist_diff<=dist_tol_rate) & type_ok;

            g = networkx.to_networkx_graph(match_matrix); # 
            # find the maximal matching of graph g
            if len(networkx.maximal_matching(g))< nAtomPerLayer: 
                is_Polar[current_surface]=True;

                g = networkx.to_networkx_graph(type_ok);
                if len(networkx.maximal_matching(g))< nAtomPerLayer: 
                    typeDiff[current_surface] = True;
                    minDiffRate[current_surface] = float("nan");
                else:
                    minDiffRate[current_surface] = np.min([np.min(dist_diff), minDiffRate[current_surface]]);

                break;
            else:
                minDiffRate[current_surface] = float("nan");


            # END of nSurface loop


        data_delete = (atom_coor[:, 2]==surface_z[-1]);
        atom_dist = atom_dist[~data_delete,:];
        atom_dist = atom_dist[:,~data_delete];
        ele_n = ele_n[~data_delete,:];
        ele_n = ele_n[:,~data_delete];
        atom_coor = atom_coor[~data_delete,:];
        surface_n=surface_n-1;
        surface_z = np.delete(surface_z, -1);

    minDiffRate[-1] = float("nan")
    is_Polar = is_Polar[::-1];
    minDiffRate = minDiffRate[::-1]
    typeDiff = typeDiff[::-1]

    return is_Polar, z_coords, minDiffRate, typeDiff


def calcRuntime(tStart,tStop):
    # Converts time between tStop and tStart (given in seconds)
    # to minues:seconds
    # 
    # Outputs: list where first element is number of minutes,
    #          second element number of seconds

    delta = tStop - tStart

    minutes = int(m.floor(delta/60))
    seconds = delta%60
    
    return [minutes,seconds]


#########################################
#                                       #
#                                       #
#          Start program                #
#                                       #
#                                       #
#########################################

#start timer
tStart = time.time()
nStruc = 0;
# Read input
folderPath, z_digit, dis_tol_rate, checkAllLayer = readInput(inputFile)
if isWin:
    typeName = folderPath.split("\\")[-1]
elif isLinux:
    typeName = folderPath.split("/")[-1]

if isWin:
    fPath = folderPath + "\\";
elif isLinux:
    fPath = folderPath + "/";
# # create a list of Miller indices
# if not useMillerList:
#     MillerList = createMillerList(maxMillerInd)
# else:
#     MillerList = MillerIndList
#     print MillerList

isPolar  = []; 
if not(checkAllLayer):
    polarFilename=typeName+"-poalrity.txt"
    file = open(fPath+polarFilename, 'w+')
    file.write("z_digit: %d \n"%z_digit)
    file.write("relative tolerance of distance difference: %.2e \n \n"%dis_tol_rate)
    file.write("Filename\t\tisPolar\tperodicity\tTypeDiff?\tminDiffRate\n")
    file.close()
else:
    rstPath = fPath + "allNonPolarRst"
    if not os.path.exists(rstPath):
        os.makedirs(rstPath)
    if isWin:
        rstPath = rstPath + "\\";
    elif isLinux:
        rstPath = rstPath + "/";



print "********************************************************"
for fn in os.listdir(folderPath):
    if fn.endswith(".in"):
        print "Reading structure from file "+fn
        atom_coor = [];
        atom_type = [];
        atom_label = [];
        nline = 0;
        nEle = 0;


        with open(fPath+fn, 'r') as openfileobject:
            for line in openfileobject:
                nline +=1;
                line  = line.split();
                if nline == 1:
                    vec_x = np.array([float(line[1]),float(line[2]), float(line[3])])
                elif nline ==2:
                    vec_y = np.array([float(line[1]),float(line[2]), float(line[3])])
                elif nline ==3:
                    vec_z = np.array([float(line[1]),float(line[2]), float(line[3])])
                else:
                    atom_coor.append([float(line[1]),float(line[2]), float(line[3])]);
                    if not(line[-1] in atom_type):
                        atom_type.append(line[-1]);
                    atom_label.append(atom_type.index(line[-1]));

        openfileobject.close()

        atom_coor = np.array(atom_coor);
        atom_coor = np.reshape(atom_coor, [nline-3, 3])

        # check polarity
        if not(checkAllLayer):
            print "Check polarity"
            polar,ped,minDiff, typeDiff = checkPolar(atom_coor,atom_label, vec_x, vec_y, z_digit, dis_tol_rate, inverse_struc = True)
            isPolar.append(polar)

            if not(polar):
                print "The structure is non-polar with periodicity %s angstrom."%ped  
            elif typeDiff:
                print "The structure is polar, with a mismatch on atom type."   
            else:
                print "The structure is polar, with a minimal distance difference rate %.3e."%minDiff 

            file = open(fPath+polarFilename, 'a')
            file.write(fn+"\t\t%d\t\t%f\t\t%d\t\t%.3e\n"%(polar,ped,typeDiff, minDiff))
            file.close()

            nStruc+=1;
        else:
            print "Check polarity on all layers"
            is_Polar, z_coords, minDiffRate, typeDiff = checkPolar_all(atom_coor,atom_label, 
                vec_x, vec_y, z_digit, dis_tol_rate, inverse_struc = True)
            file = open(rstPath+fn.split(".")[0]+".txt", 'w+')
            file.write("z_digit: %d \n"%z_digit)
            file.write("relative tolerance of distance difference: %.2e \n \n"%dis_tol_rate)
            file.write("z_coords\tis_Polar\ttypeDiff\tminDiffRate \n")
                      
            for i in xrange(0, len(is_Polar)):
                file.write("%f\t\t%d\t\t%d\t\t%.3e\n"%(z_coords[i],is_Polar[i],typeDiff[i], minDiffRate[i]))
            file.close()
            nStruc+=1;

        print "********************************************************"
# Stop timer
tStop = time.time()
runtime = calcRuntime(tStart,tStop)
#Output statistics
if not(checkAllLayer):
    nNonPolar = nStruc - sum(isPolar)
    file = open(fPath+'stats.out','w+')
    file.write("Total number of structures: %i\n"%nStruc)
    file.write("Total number of polar structure: %i\n"%nNonPolar)
    file.write("\nRuntime: %i min %i sec\n"%(runtime[0],runtime[1]))
    file.close()
else:
    file = open(fPath+'stats.out','w+')
    file.write("Total number of structures: %i\n"%nStruc)
    file.write("\nRuntime: %i min %i sec\n"%(runtime[0],runtime[1]))
    file.close()

# # End of program
