#!/usr/bin/env python
'''
Created on Sep 16, 2016

@author: fu
'''

import numpy as np
import itertools
import os
from pymatgen.io.vaspio import Poscar
from progressbar import *

class EffectiveBandStructure(object):
    '''
    generate the structures for effective band calculation.
    
    Args:
        path: path of origin POSCAR excluding self-POSCAR
    '''


    def __init__(self, path=None):
        if path == None:
            self.path='./'
        else:
            self.path=path
            if not(os.path.exists(self.path)):
                print 'Warning! path is wrong.'
            elif not(os.path.isdir(self.path)):
                print "Warning! not a directory."
                
        self.rstructures=[] # irreducible structure
        self.weight=[] # weight of irreducible structures
        self.irrlist=[] # array of irreducible structures
                
    def read(self, filename='POSCAR'):
        poscar=Poscar.from_file(self.path+'/'+filename, True, True)
        return poscar.structure
    
    def getEffectiveBandStructure(self, filename='POSCAR', size=[3,3,3], numreplace=1, elementType='Se', flength=3):
        '''
        generate the irreducible structures for effective band calculation
        
        Args:
            filename: POSCAR
            size: size of supercell
            numreplace: number of replacing atom
            flength: length of X in POSCAR-XXX
        '''
        structure=self.read(filename)
        structure.make_supercell(size)

        atoms=structure.num_sites
        rslist=np.array([x for x in itertools.combinations(range(atoms), numreplace)]) # save list of replacing site of atom
        
        pbar = ProgressBar(widgets=['Read: ', Percentage(), ' ', Bar(),
                                    ' ', ETA()], maxval=len(rslist)).start()
        
        counter=0 # total number of structures
        icounter=0 # number of irreducible structures
        for rs in rslist:
            # replacing atom
            tmpstr=structure.copy()
            tmpfilename='POSCAR-'+(flength-len(str(icounter)))*'0'+str(icounter)
            for site in rs:
                tmpstr.replace(site, elementType)
            #tmpstr.sort()
            
            # category of structure
            if self.rstructures == []:
                self.rstructures.append(tmpstr)
                self.weight.append(1)
                self.irrlist.append(rs)
                icounter+=1
                tmpstr=tmpstr.get_primitive_structure(tolerance=0.01)
                tmpstr.to(fmt='poscar', filename=self.path+'/'+tmpfilename)
            else:
                if not(self.exist(self.rstructures, tmpstr, rs)):
                    self.rstructures.append(tmpstr)
                    icounter+=1
                    tmpstr=tmpstr.get_primitive_structure(tolerance=0.01)
                    tmpstr.to(fmt='poscar', filename=self.path+'/'+tmpfilename)
                
            pbar.update(counter)
            counter+=1
        pbar.finish()
        
        # output weight and irrlist
        self.outputWeightIrrlist(numreplace)
        
    def exist(self, structureSet, structure, rs):
        '''
        1. judge whether the structure exist in the set of structures
        2. modify the weight and irrlist array according to the result of comparison
        Args:
            structureSet: set of irreducible structures
            structure: comparing structure
            rs: replacing site of this structure
            
        Return:
            boolean value (True/False) of comparison
        '''
        counter=0
        flag=False
        while counter < len(structureSet) and not(flag):   
            flag=structure.matches(structureSet[counter])
            if flag == True:
                self.weight[counter]+=1
                self.irrlist[counter]=np.vstack((self.irrlist[counter],rs))
                flag=True
            elif flag == None:
                print "Warning! "+str(rs)
            counter+=1
            
        if flag == False:
            self.weight.append(1)
            self.irrlist.append(rs)
        #print '>> '+str(counter), str(len(structureSet))
        return flag
    
    def exist2(self, structureSet, structure, rs):
        '''
        1. judge whether the structure exist in the set of structures
        2. modify the weight and irrlist array according to the result of comparison
        Args:
            structureSet: set of irreducible structures
            structure: comparing structure
            rs: replacing site of this structure
            
        Return:
            boolean value (True/False) of comparison
        '''
        counter=0
        flag=False
        for s0 in structureSet:
            flag=structure.matches(s0)
            if flag == True:
                self.weight[counter]+=1
                self.irrlist[counter]=np.vstack((self.irrlist[counter],rs))
                flag=True
                break
            elif flag == None:
                print "Warning! "+str(rs)
            counter+=1
            
        if flag == False:
            self.weight.append(1)
            self.irrlist.append(rs)
        #print '>> '+str(counter), str(len(structureSet))
        return flag
    
    def outputWeightIrrlist(self, numreplace, filename='irreducible'):
        '''
        output the weight and list of irreducible replacing
        '''
        with open(self.path+'/'+filename, 'w') as out:
            if numreplace == 1:
                for i in xrange(0, len(self.weight)): # irreducible
                    out.write('%d %6d' %(i, self.weight[i]))
                    for j in xrange(0, len(self.irrlist[i])): # replacing
                        out.write('\t[%d]' %(self.irrlist[i][j]))
                    out.write('\n')
            else:           
                for i in xrange(0, len(self.weight)): # irreducible
                    out.write('%d %6d' %(i, self.weight[i]))
                    for j in xrange(0, len(self.irrlist[i])): # replacing
                            out.write('\t[')
                            for k in xrange(0, len(self.irrlist[i][j])): # atom
                                if k == len(self.irrlist[i][j])-1:
                                    out.write('%d]' %self.irrlist[i][j][k])
                                else:
                                    out.write('%d,' %self.irrlist[i][j][k])
                    out.write('\n')
                    

# -------------------- test --------------------
e=EffectiveBandStructure('/home/fu/workspace/SeTe/alloy/333/test/2')
e.getEffectiveBandStructure(numreplace=2)                   

        
