#!/usr/bin/env python
'''
Created on Sep 9, 2016

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
        rlist=np.array([x for x in itertools.combinations(range(atoms), numreplace)]) # save list of replacing atom
        rstructures=[]
        weight=[] # weight of irreducible structures
        irrlist=[] # array of irreducible structures
        
        pbar = ProgressBar(widgets=['Read: ', Percentage(), ' ', Bar(),
                                    ' ', ETA()], maxval=len(rlist)).start()
        counter=0 # number of irreducible structures
        for i in xrange(0, len(rlist)):
            tmpstr=structure.copy()
            tmpname='POSCAR-'+(flength-len(str(counter)))*'0'+str(counter)
            for atom in rlist[i]:
                tmpstr.replace(atom, elementType)
                #tmpname=tmpname+'-'+str(atom) # index of replacing atom
            tmpstr.sort()
            
            # category of structure
            if rstructures == []:
                rstructures.append(tmpstr)
                weight.append(1)
                irrlist.append(rlist[i])
                counter+=1
                tmpstr.to(fmt='poscar', filename=self.path+'/'+tmpname)
            else:
                isExist=False
                for j in xrange(0, len(rstructures)):
                    tmpstate=tmpstr.matches(rstructures[j])
                    if tmpstate == True:
                        isExist=True
                        weight[j]=weight[j]+1
                        irrlist[j]=np.vstack((irrlist[j], rlist[i]))
                        break
                    elif tmpstate == None:
                        print 'Warning! '+str(rlist[i])
                if not(isExist):
                    rstructures.append(tmpstr)
                    weight.append(1)
                    irrlist.append(rlist[i])
                    counter+=1
                    tmpstr.to(fmt='poscar', filename=self.path+'/'+tmpname)
                    
            pbar.update(i)
        irrlist=np.array(irrlist)
        pbar.finish()
        
        # out weight and irrlist
        with open(self.path+'/irreducible', 'w') as out:
            if numreplace == 1:
                for i in xrange(0, len(weight)): # irreducible
                    out.write('%d %6d' %(i, weight[i]))
                    for j in xrange(0, len(irrlist[i])): # replacing
                        out.write('\t[%d]' %(irrlist[i][j]))
                    out.write('\n')
            else:           
                for i in xrange(0, len(weight)): # irreducible
                    out.write('%d %6d' %(i, weight[i]))
                    for j in xrange(0, len(irrlist[i])): # ieplacing
                            out.write('\t[')
                            for k in xrange(0, len(irrlist[i][j])): # atom
                                if k == len(irrlist[i][j])-1:
                                    out.write('%d]' %irrlist[i][j][k])
                                else:
                                    out.write('%d,' %irrlist[i][j][k])
                    out.write('\n')
                    
# -------------------- test --------------------
e=EffectiveBandStructure('/home/fu/workspace/SeTe/alloy/333/3')
e.getEffectiveBandStructure(numreplace=3)
