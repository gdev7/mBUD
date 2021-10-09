import os
import pandas as pd
import numpy as np
from numpy import sin,cos
import math
import argparse
import cif2cell

import networkx as nx

import datetime

class LIBRARY:
    
    def __init__ (self):
        
        # Library dir
        self.dir = os.getcwd() + "/lib/"
        
        # General infor
        self.dim = 3
        self.a2r = np.pi / 180.0
        self.cifextension = '.cif'
        self.xyzextension = '.xyz'
        
        # Metal info
        self.METALcsv = self.dir + "metal.csv"
        self.dfmetal = pd.read_csv(self.METALcsv)
        self.metalsymbol = self.dfmetal.iloc[:,1]

        # Bond info
        self.ATRcsv = self.dir + "atr.csv"
        self.dfatr = pd.read_csv(self.ATRcsv)
        
class ATOM(LIBRARY):
    
    dim = 3
    lib = LIBRARY()
    
    def __init__ (self, data, index):
        #self.label = data[0]
        self.label = data[1] + str(index + 1)
        self.symbol = data[1]
        self.index = index
        
        self.neighborlist = []
        self.nneighbor = 0
        
        self.gridindex = 0
        
        self.ismetal = False
        
        self.x = []
        self.writex = []
        
        for idim in range(self.dim):
            self.x.append(data[idim+3])
            self.writex.append(data[idim+3])
            
        self.atr = self.get_atr()
        
    def destroy(self, array):
        del array
        array = []
        return array
    
    def destroy_all(self):
        del self.label
        del self.symbol

        del self.neighborlist
        del self.x
        del self.writex
            
    def get_atr(self):
        try:
            atr = self.lib.dfatr.iloc[self.lib.dfatr[(self.lib.dfatr['Symbol'] == self.symbol)].index[0]]['atr']
        except:
            raise ValueError('ATR not found')
        return atr
    
    def print_info(self):
        print(self.index, end = ' ')
        print(self.label, end = ' ')
        print(self.symbol, end = ' ')
        print(self.x, end = ' ')
        print(self.atr)
    
    def setnneighborzero(self):
        self.nneighbor = 0
    
    def check_neighbor(self, iindex):
        
        check = False
        for ineighbor in self.neighborlist:
            if ineighbor == iindex:
                check = True
                
        return check
    
    def add_neighbor(self, iindex):
        check = False
        
        if len(self.neighborlist) > 0:
            check = self.check_neighbor(iindex)
            
        if check == False:
            self.neighborlist.append(iindex)
        
    def remove_neighbor(self, iindex):
        currentneighborlist = []
        currentneighborlist = self.neighborlist
        self.neighborlist = self.destroy(self.neighborlist)

        for ineighbor in currentneighborlist:
            if ineighbor != iindex:
                self.add_neighbor(ineighbor)
                
        del currentneighborlist
        
class MOF(ATOM,LIBRARY):
    
    lib = LIBRARY()
    dim = 3
    
    def __init__ (self, name):
        self.ciffile = name
        self.G = nx.Graph(name = 'mbud')
        self.fragG = nx.Graph(name = 'mbud')
        self.compG = nx.Graph(name = 'mbud')
        
        self.detH = 1.0
        
        self.lx = []
        self.ar = []
        self.atom = []
        self.h = []
        self.hinv = []
        self.atomtypelist = []
        self.metaltypelist = []
        
        #SBU lists (frag)
        self.linkerlist = []
        self.funcgrouplist = []
        self.solventlist = []
        self.metalnodelist = []
        self.capairlist = []
        
        #SBU lists (comp)
        self.complinkerlist = []
        self.compmetalnodelist = []
        self.compcapairlist = []
        
        #Grid Informations
        self.skin = 0.18
        self.gridlxmax = self.skin
        #self.minbondlength2 = 0.63*0.63
        self.totgrid = 1
        
        self.gridlx = []
        self.ngrid = []
        self.neighgrid = []
        self.gridatomlist = []
        
    def destroy(self, array):
        del array
        array = []
        return array
    
    def destroy_all(self):
        del self.ciffile
        
        self.G.clear()
        self.fragG.clear()
        self.compG.clear()
        
        del self.lx
        del self.ar
        
        for iatom in self.atom:
            iatom.destroy_all()
        del self.atom
        
        del self.h
        del self.hinv
        del self.atomtypelist
        del self.metaltypelist
        
        #SBU lists (frag)
        del self.linkerlist
        del self.funcgrouplist
        del self.solventlist
        del self.metalnodelist
        del self.capairlist
        
        #SBU lists (comp)
        del self.complinkerlist
        del self.compmetalnodelist
        del self.compcapairlist
        
        del self.gridlx
        del self.ngrid
        del self.neighgrid
        del self.gridatomlist
        
    def get_boxinfo(self):
        
        self.loop = False
        self.atom = self.destroy(self.atom)
        self.lx = self.destroy(self.lx)
        self.ar = self.destroy(self.ar)
        
        index = 0
        
        f = open (self.ciffile,'r')
        lines = f.readlines()
        f.close()
        
        for line in lines:
            list = line.split()

            if len(list) > 0:
                if list[0] == '_cell_length_a':
                    self.lx.append(float(list[1]))
                if list[0] == '_cell_length_b':
                    self.lx.append(float(list[1]))
                if list[0] == '_cell_length_c':
                    self.lx.append(float(list[1]))

                if list[0] == '_cell_angle_alpha':
                    self.ar.append(self.lib.a2r*float(list[1]))
                if list[0] == '_cell_angle_beta':
                    self.ar.append(self.lib.a2r*float(list[1]))
                if list[0] == '_cell_angle_gamma':
                    self.ar.append(self.lib.a2r*float(list[1]))

                if list[0] == 'loop_':
                    self.loop = True

                if self.loop == True and len(list) > 5:
                    data = []
                    #label, symbol, symmetry, ux, uy, uz, occupancy
                    data.append(list[0])
                    data.append(list[1])
                    data.append(int(list[2]))
                    data.append(float(list[3]))
                    data.append(float(list[4]))
                    data.append(float(list[5]))
                    data.append(float(list[6]))
                    self.atom.append(ATOM(data,index))
                    del data
                    index += 1   
        del lines

    def get_hmatrix(self):
        self.h = self.destroy(self.h)
            
        for idim in range(self.dim):
            data = []
            for jdim in range(self.dim):
                data.append(0.0)
            self.h.append(data)
            del data
            
        #get h-matrix
        self.h[0][0] = self.lx[0]
        self.h[0][1] = self.lx[1]*cos(self.ar[2])
        self.h[0][2] = self.lx[2]*cos(self.ar[1])
        self.h[1][0] = 0.0
        self.h[1][1] = self.lx[1]*sin(self.ar[2])
        local = cos(self.ar[0]) - cos(self.ar[1])*cos(self.ar[2])
        local /= sin(self.ar[2])
        self.h[1][2] = self.lx[2]*local
        self.h[2][0] = 0.0
        self.h[2][1] = 0.0
        self.h[2][2] = self.lx[2]*np.sqrt(1 - cos(self.ar[1])*cos(self.ar[1]) - local*local)
        
        for idim in range(self.dim):
            self.detH *= self.h[idim][idim]
            
        self.hinv = self.destroy(self.hinv)
        ht = []
        
        for idim in range(self.dim):
            data = []
            for jdim in range(self.dim):
                data.append(0.0)
            self.hinv.append(data)
            ht.append(data)
            del data
            
        #get hinv-matrix
        ht[0][0] = self.h[1][1]*self.h[2][2] - self.h[1][2]*self.h[2][1]
        ht[0][1] = 0.0
        ht[0][2] = 0.0

        ht[1][0] = self.h[0][2]*self.h[2][1] - self.h[0][1]*self.h[2][2]
        ht[1][1] = self.h[0][0]*self.h[2][2] - self.h[0][2]*self.h[2][0]
        ht[1][2] = self.h[0][1]*self.h[2][0] - self.h[0][0]*self.h[2][1]

        ht[2][0] = self.h[0][1]*self.h[1][2] - self.h[0][2]*self.h[1][1]
        ht[2][1] = self.h[0][2]*self.h[1][0] - self.h[0][0]*self.h[1][2]
        ht[2][2] = self.h[0][0]*self.h[1][1] - self.h[0][1]*self.h[1][0]

        for idim in range(self.dim):
            for jdim in range(self.dim):
                self.hinv[idim][jdim] = ht[jdim][idim]

        for idim in range(self.dim):
            for jdim in range(self.dim):
                self.hinv[idim][jdim] /= self.detH
                
        del ht
                
    def get_gridinfo(self):

        self.gridlxmax = self.get_maxcovbl() + self.skin
        #print('Minimum grid length (maximum covalent radius + skin/2):', self.gridlxmax)
        
        if self.gridlxmax < self.skin:
            raise NameError('Grid length is not right')
            
        self.gridlx = self.destroy(self.gridlx)
        self.ngrid = self.destroy(self.ngrid)
        
        for idim in range(self.dim):
            self.gridlx.append(0.0)
            self.ngrid.append(1)
            
        for idim in range(self.dim):
            for jdim in range(self.dim):
                if jdim == idim:
                    self.gridlx[jdim] = 1.0 / self.ngrid[idim]
                else:
                    self.gridlx[jdim] = 0.0

            self.gridlx = self.fractional_to_cartesian(self.gridlx)

            while self.gridlx[idim] > self.gridlxmax:
                self.ngrid[idim] += 1
                
                for jdim in range(self.dim):
                    if jdim == idim:
                        self.gridlx[jdim] = 1.0 / self.ngrid[idim]
                    else:
                        self.gridlx[jdim] = 0.0
                self.gridlx = self.fractional_to_cartesian(self.gridlx)
            self.ngrid[idim] -= 1
            
        for idim in range(self.dim):
            if self.ngrid[idim] == 0:
                self.ngrid[idim] = 1

        for idim in range(self.dim):
            self.gridlx[idim] = 1.0 / self.ngrid[idim]
        
        #print('Grid dimension (updated fractional):', self.gridlx)
        #print('Grid dimension (updated cartesian):', self.fractional_to_cartesian(self.gridlx))
        #print('Number of grids in axis:', self.ngrid)

        self.totgrid = 1
        for idim in range(self.dim):
            self.totgrid *= self.ngrid[idim]
        #print('Total number of grids:', self.totgrid)
        
        self.neighgrid = self.destroy(self.neighgrid)

        for igrid in range(self.totgrid):
            igz = math.floor(igrid/(self.ngrid[0]*self.ngrid[1]))
            igy = math.floor((igrid - igz*self.ngrid[0]*self.ngrid[1])/self.ngrid[0])
            igx = igrid - igz*self.ngrid[0]*self.ngrid[1] - igy*self.ngrid[0]
            #print('\nCurrent grid:', igrid)
            data = []

            for iz in range(-1,2):
                iggz = igz + iz

                if iggz < 0:
                    iggz += self.ngrid[2]
                elif iggz >= self.ngrid[2]:
                    iggz -= self.ngrid[2]

                for iy in range(-1,2):
                    iggy = igy + iy

                    if iggy < 0:
                        iggy += self.ngrid[1]
                    elif iggy >= self.ngrid[1]:
                        iggy -= self.ngrid[1]

                    for ix in range(-1,2):
                        iggx = igx + ix

                        if iggx < 0:
                            iggx += self.ngrid[0]
                        elif iggx >= self.ngrid[0]:
                            iggx -= self.ngrid[0]

                        currgrid = iggx + iggy * self.ngrid[0] + iggz * self.ngrid[0] * self.ngrid[1]
                        if currgrid not in data:
                            data.append(currgrid)
                        #print(currgrid, end = ' ')
            self.neighgrid.append(data)
            del data
        #print('Neighgrids', self.neighgrid)

    def get_gridindex(self,x):

        igrid = []
        for i in range(self.dim):
            if x[i] < 0.0:
                x[i] += 1.0
            if x[i] > 1.0:
                x[i] -= 1.0

            igrid.append(math.ceil(x[i] / self.gridlx[i]))
            
        for idim in range(self.dim):
            if igrid[idim] == 0:
                igrid[idim] += 1

        gridindex = 0
        gridindex += self.ngrid[0] * self.ngrid[1] * (igrid[2] - 1)
        gridindex += self.ngrid[0] * (igrid[1] - 1)
        gridindex += igrid[0] - 1
        #print(igrid, self.ngrid)
        del igrid
        
        found = True
        if gridindex < 0:
            found = False

        return gridindex, found
        
    def get_atomgridinfo(self):
        
        self.get_gridinfo()
        
        for iatom in self.atom:
            iatom.gridindex, found = self.get_gridindex(iatom.x)
            if found == False:
                print(iatom.label,iatom.x, iatom.gridindex)
                raise NameError('Wrong gridindex')
                
        self.gridatomlist = self.destroy(self.gridatomlist)
        
        for i in range(self.totgrid):
            self.gridatomlist.append([])
        #print('Length of gridatomlist:', len(self.gridatomlist))

        for iatom in self.atom:
            self.gridatomlist[iatom.gridindex].append(iatom.index)
        #print(self.gridatomlist)
                    
    def get_atomtypelist(self):
        self.atomtypelist = self.destroy(self.atomtypelist)

        tempatomtypelist = []

        for iatom in self.atom:
            if iatom.symbol not in tempatomtypelist:
                tempatomtypelist.append(iatom.symbol)
                data = []
                data.append(iatom.symbol)
                data.append(self.lib.dfatr.iloc[self.lib.dfatr[(self.lib.dfatr['Symbol'] == iatom.symbol)].index[0]]['atr'])
                self.atomtypelist.append(data)
                del data
        del tempatomtypelist
                    
    def get_metaltypelist(self):
        self.metaltypelist = self.destroy(self.metaltypelist)
        
        if len(self.atomtypelist) == 0:
            raise NameError("Get atomtypelist first")
        else:
            for iatomtype in self.atomtypelist:
                for imetalsymbol in self.lib.metalsymbol:
                    if iatomtype[0] == imetalsymbol:
                        data1 = []
                        for iatom in self.atom:
                            if iatom.symbol == imetalsymbol:
                                data1.append(iatom.index)
                                iatom.ismetal = True
                        data2 = []
                        data2.append(imetalsymbol)
                        data2.append(data1)
                        del data1
                        self.metaltypelist.append(data2)
                        del data2
    
    def fractional_to_cartesian(self, dux):
        dx = []
        dx = self.destroy(dx)
        
        for i in range(self.dim):
            val = 0.0
            
            for j in range(self.dim):
                val += self.h[i][j] * dux[j]
                
            dx.append(val)            
        return dx
    
    def cartesian_to_fractional(self, dx):
        dux = []
        dux = self.destroy(dux)
        
        for i in range(self.dim):
            val = 0.0
            
            for j in range(self.dim):
                val += self.hinv[i][j] * dx[j]
                
            dux.append(val)
        return dux
    
    def get_maxcovbl(self):
        
        maxcovbl = self.skin
        for iatomtype in self.atomtypelist:
            for jatomtype in self.atomtypelist:
                covbl = iatomtype[1] + jatomtype[1]
                if covbl > maxcovbl:
                    maxcovbl = covbl
                        
        return maxcovbl
    
    def check_bond(self, iatom, jatom):
        '''1. Check this website: https://www.slideshare.net/NextMoveSoftware/rdkit-gems
            Slide number: 20+
           2. r_cov values are taken from openbabel elementtable.h
        '''
        maxdx2 = iatom.atr + jatom.atr + self.skin
        maxdx2 *= maxdx2
        
        x1 = iatom.x
        x2 = jatom.x

        dx = []
        for idim in range(self.dim):
            val = x1[idim] - x2[idim]
            if val > 0.5:
                val -= 1.0
            if val < -0.5:
                val += 1.0
            dx.append(val)

        dx = self.fractional_to_cartesian(dx)
        val = 0.0
        for idim in range(self.dim):
            val += dx[idim]*dx[idim]
            if val > maxdx2:
                del dx
                return False
        
        del dx
        
        minbondlength2 = (iatom.atr + jatom.atr)/2.0
        minbondlength2 *= minbondlength2
        
        if val < minbondlength2:
            raise ValueError('Atom overlap detected')
        return True

    def get_distance(self, x1, x2):
        """ Calculate the distance between two atoms, takes self.atom[index] as input, return distance in float"""

        dx = []
        for idim in range(self.dim):
            val = x1[idim] - x2[idim]
            if val > 0.5:
                val -= 1.0
            if val < -0.5:
                val += 1.0
            dx.append(val)

        dx = self.fractional_to_cartesian(dx)
        val = 0.0
        for idim in range(self.dim):
            val += dx[idim]*dx[idim]
        val = np.sqrt(val)
        del dx
        return val
            
    def check_repeat_atomlabel(self):
        
        repeat = False
        atomcount = len(self.atom)
        
        for atomcount1 in range(0, atomcount - 1):
            label1 = self.atom[atomcount1].label
            
            for atomcount2 in range(atomcount1 + 1, atomcount):
                label2 = self.atom[atomcount2].label
                if label1 == label2:
                    repeat = True
        
        return repeat         
        
    def get_neighborlist_without_grid(self):
        
        self.G.clear()
        self.G.add_nodes_from([i for i in range(len(self.atom))])
            
        for iatom in self.atom:
            #print(iatom.index, end = ' ')
            for jatom in self.atom:
                if iatom.index < jatom.index:
                    if self.check_bond(iatom,jatom) == True:
                        iatom.add_neighbor(jatom.index)
                        jatom.add_neighbor(iatom.index)
                            
                        if self.G.has_edge(iatom.index, jatom.index) == False:
                            self.G.add_edge(iatom.index, jatom.index)

        #print('Bond Calculation Time (s):\t', (end-start))
        
    def get_neighborlist_with_grid(self):
        
        self.G.clear()
        self.G.add_nodes_from([i for i in range(len(self.atom))])
        
        for iatom in self.atom:
            #print(iatom.index, end = ' ')
            for igrid in self.neighgrid[iatom.gridindex]:
                for jindex in self.gridatomlist[igrid]:
                    if iatom.index < jindex:
                        jatom = self.atom[jindex]
                        if self.check_bond(iatom,jatom) == True:
                            iatom.add_neighbor(jatom.index)
                            jatom.add_neighbor(iatom.index)
                            
                            if self.G.has_edge(iatom.index, jatom.index) == False:
                                self.G.add_edge(iatom.index, jatom.index)
                                
        #print('Bond Calculation Time (s):\t', (end-start))
        
    def get_neighborlist(self, grid):
        
        if grid:
            self.get_neighborlist_with_grid()
        else:
            self.get_neighborlist_without_grid()
        
    def clear_neighborlist(self):
        
        for iatom in self.atom:
            for ineighbor in iatom.neighborlist:
                iatom.remove_neighbor(ineighbor)
                
    def remove_neighbor(self, iindex, jindex):
        
        self.atom[iindex].remove_neighbor(jindex)
        self.atom[jindex].remove_neighbor(iindex)
        
    def get_graph(self):
        
        self.G.clear()
        self.G.add_nodes_from([i for i in range(len(self.atom))])
        
        for iatom in self.atom:
            for ineighbor in iatom.neighborlist:
                if self.G.has_edge(iatom.index, ineighbor) == False:
                    self.G.add_edge(iatom.index, ineighbor)
            
    def get_solvent(self):
        
        fragmentlist = nx.connected_components(self.G)
        
        for ifragment in fragmentlist:
            metal = False
            for iindex in ifragment:
                if self.atom[iindex].ismetal:
                    metal = True
                    break
                    
            if metal == False:
                self.solventlist.append(list(ifragment))
                
        #Remove solvent from graph
        if len(self.solventlist) > 0:
            print('Solvent presence:\t\t','True')
            for isolvent in self.solventlist:
                self.G.remove_nodes_from(isolvent)
        '''
        else:
            print('Solvent presence:\t\t','False')
        '''
    
    def check_multimetal(self, iindex, imetalindex):
    
        for ineighbor in self.atom[iindex].neighborlist:
            if self.atom[ineighbor].ismetal == True:
                if ineighbor != imetalindex:
                    return True

        return False

    def break_mof(self):
        
        self.capairlist = self.destroy(self.capairlist)

        #get capairlist
        for imetaltype in self.metaltypelist:
            for imetalindex in imetaltype[1]:
                for ineighbor in self.atom[imetalindex].neighborlist:
                    #print('All Neigh',self.atom[imetalindex].label,self.atom[ineighbor].label)
                    if self.atom[ineighbor].ismetal == False:
                        count = 0
                        for jneighbor in self.atom[ineighbor].neighborlist:
                            if self.atom[jneighbor].ismetal == True:
                                count += 1
                            elif self.atom[jneighbor].symbol == 'H':
                                count += 1
                        if count != len(self.atom[ineighbor].neighborlist):
                            self.capairlist.append([imetalindex, ineighbor])
                            #print('Only CA',self.atom[imetalindex].label,self.atom[ineighbor].label)

        ############################################
        #This part is for actual building block
        ############################################
        self.fragG = self.G.copy()

        #break capairlist
        for inode, jnode in self.capairlist:
            if self.fragG.has_edge(inode, jnode):
                #print('Breaking path between', self.atom[inode].label, self.atom[jnode].label)
                self.fragG.remove_edge(inode, jnode)

        self.metalnodelist = self.destroy(self.metalnodelist)
        self.funcgrouplist = self.destroy(self.funcgrouplist)
        self.linkerlist = self.destroy(self.linkerlist)

        fragmentlist = nx.connected_components(self.fragG)

        for ifragment in fragmentlist:
            found = False
            count = 0
            for iindex in ifragment:
                if self.atom[iindex].ismetal:
                    found = True

                if found:
                    break
                else:
                    for inode, jnode in self.capairlist:
                        if iindex == jnode:
                            #print(iindex,inode,jnode)
                            count += 1

            #print(len(ifragment),count)
            if found:
                self.metalnodelist.append(list(ifragment))
            else:
                if count == 1:
                    self.funcgrouplist.append(list(ifragment))
                elif count > 1:
                    self.linkerlist.append(list(ifragment))
                else:
                    raise ValueError('Wrong Breaking')

        """
        Breaking of Building Block Critical Info:

        Ref to: Figure 3 of (Cryst. Growth Des. 2017, 17, 5801−5810)
        The buliding blocks can be differentiated into to types based on the breaking point.
        We will take COO group for the example.

        1. Building blocks consistent with chemical synthesis. In this case the breaking point
        will be between metal atom and O atoms.

        2. Building blocks useful for computational construction. In this case the breaking point
        will be between C of COO group and the atom, part of the linker.

        Update:
            1. Currently this scenario is observered for COO group only.
            2. We are searching for more cases.
        """
        ############################################
        #This part is for computational construction
        ############################################

        self.compG = self.G.copy()
        if self.funcgrouplist:
            for ifuncgroup in self.funcgrouplist:
                self.compG.remove_nodes_from(ifuncgroup)
                
        self.compcapairlist = self.destroy(self.compcapairlist)
        self.compmetalnodelist = self.destroy(self.compmetalnodelist)
        self.complinkerlist = self.destroy(self.complinkerlist)
        
        allcalist = []
        for inode, jnode in self.capairlist:
            allcalist.append(jnode)

        for imetalnode in self.metalnodelist:
            calist = []
            for iindex in imetalnode:
                #print(self.atom[iindex].label)
                for inode, jnode in self.capairlist:
                    if iindex == inode:
                        calist.append(jnode)

            capairlist = []
            removecalist = []
            for k, ica in enumerate(calist):
                for ineighbor in self.atom[ica].neighborlist:
                    for jca in calist[k+1:]:
                        for jneighbor in self.atom[jca].neighborlist:
                            if ineighbor == jneighbor:
                                if ineighbor not in imetalnode:
                                    update = True
                                    for kneighbor in self.atom[ica].neighborlist:
                                        if kneighbor not in imetalnode:
                                            if kneighbor != ineighbor:
                                                if self.atom[kneighbor].symbol != 'H':
                                                    update = False
                                    for kneighbor in self.atom[jca].neighborlist:
                                        if kneighbor not in imetalnode:
                                            if kneighbor != ineighbor:
                                                if self.atom[kneighbor].symbol != 'H':
                                                    update = False
                                    
                                    if update:
                                        for kneighbor in self.atom[jneighbor].neighborlist:
                                            if kneighbor not in allcalist:
                                                if kneighbor not in imetalnode:
                                                    capairlist.append([jneighbor,kneighbor])
                                                else:
                                                    removecalist.append(jneighbor)
                                                if ica not in removecalist:
                                                    removecalist.append(ica)
                                                if jca not in removecalist:
                                                    removecalist.append(jca)

            for inode, jnode in self.capairlist:
                if jnode in calist:
                    if inode in imetalnode:
                        if jnode not in removecalist:
                            self.compcapairlist.append([inode,jnode])

            if capairlist:
                for inode, jnode in capairlist:
                    self.compcapairlist.append([inode,jnode])
                    
            del calist
            del capairlist
            del removecalist

        del allcalist
        
        #break compcapairlist
        for inode, jnode in self.compcapairlist:
            if self.compG.has_edge(inode, jnode):
                #print('Breaking path between', self.atom[inode].label, self.atom[jnode].label)
                self.compG.remove_edge(inode, jnode)

        fragmentlist = nx.connected_components(self.compG)

        for ifragment in fragmentlist:
            metal = False
            for iindex in ifragment:
                if self.atom[iindex].ismetal:
                    metal = True
                    break

            if metal:
                self.compmetalnodelist.append(list(ifragment))
            else:
                self.complinkerlist.append(list(ifragment))

        #print(self.metalnodelist)
        #print(self.compmetalnodelist)
        
        count = 0
        for imetalnode in self.compmetalnodelist:
            addindex = []
            for iindex in imetalnode:
                for inode, jnode in self.compcapairlist:
                    if inode == iindex:
                        addindex.append(jnode)

            for iindex in addindex:
                self.compmetalnodelist[count].append(iindex)
            del addindex
            count += 1

        #print(self.compmetalnodelist)
        #print(self.linkerlist)

        count = 0
        for ilinker in self.complinkerlist:
            addindex = []
            for iindex in ilinker:
                for inode, jnode in self.compcapairlist:
                    if jnode == iindex:
                        if inode not in addindex:
                            addindex.append(inode)

            for iindex in addindex:
                self.complinkerlist[count].append(iindex)
            del addindex
            count += 1       
        #print('Metal Node Count:\t\t', len(self.metalnodelist))
        #print('Functional Group Count:\t\t', len(self.funcgrouplist))
        #print('Organic Linker Counts:\t\t', len(self.linkerlist))
        
    def get_fragment_data(self, fragment, case):
        
        fragmentatomtypelist = []
        fragmentatomtypelist = self.destroy(fragmentatomtypelist)
        for iindex in fragment:
            if self.atom[iindex].symbol not in fragmentatomtypelist:
                fragmentatomtypelist.append(self.atom[iindex].symbol)
        fragmentatomtypelist.sort()
        
        fragmentatomtypecountlist = []
        fragmentatomtypecountlist = self.destroy(fragmentatomtypecountlist)
        for fragmentatomtype in fragmentatomtypelist:
            fragmentatomtypecountlist.append(0)
            
        for iindex in fragment:
            count = 0
            for fragmentatomtype in fragmentatomtypelist:
                if fragmentatomtype != self.atom[iindex].symbol:
                    count += 1
                else:
                    break
            fragmentatomtypecountlist[count] += 1
            
        fragmentnca = 0
        if case == 0:
            for iindex in fragment:
                for inode, jnode in self.capairlist:
                    if inode == iindex:
                        fragmentnca += 1
        elif case == 1:
            for iindex in fragment:
                for inode, jnode in self.capairlist:
                    if jnode == iindex:
                        fragmentnca += 1
            
        return fragmentatomtypelist, fragmentatomtypecountlist, fragmentnca
        
    def check_uniq_fragmentlist(self, fragment, uniqfragmentlist, case):
        
        if len(uniqfragmentlist) == 0:
            uniqfragmentlist.append(fragment)
            return uniqfragmentlist
        else:
            uniq = True
            for uniqfragment in uniqfragmentlist:
                if len(uniqfragment) == len(fragment):
                    uniqfragmentatomtypelist, uniqfragmentatomtypecountlist, uniqfragmentnca = self.get_fragment_data(uniqfragment, case)
                    fragmentatomtypelist, fragmentatomtypecountlist, fragmentnca = self.get_fragment_data(fragment, case)
                    if uniqfragmentatomtypelist == fragmentatomtypelist:
                        if uniqfragmentatomtypecountlist == fragmentatomtypecountlist:
                            if case == 0:
                                if uniqfragmentnca == fragmentnca:
                                    if fragmentnca != 0:
                                        uniq = False
                            elif case == 1:
                                if uniqfragmentnca == fragmentnca:
                                    if fragmentnca != 0:
                                        uniq = False
                            else:
                                uniq = False
                            
            if uniq:
                uniqfragmentlist.append(fragment)
        
        return uniqfragmentlist
        
    def get_uniq_fragmentlist(self, fragmentlist, case):
        
        '''
        Cases:
        
        0 - ca connected to metalnode
        1 - ca connected to linker
        # - others
        '''
        
        uniqfragmentlist = []
        uniqfragmentlist = self.destroy(uniqfragmentlist)
        
        for fragment in fragmentlist:
            uniqfragmentlist = self.check_uniq_fragmentlist(fragment, uniqfragmentlist, case)
            
        return uniqfragmentlist
    
    def wrap_fragment(self, fragment):
        
        for iindex in fragment:
            for idim in range(self.dim):
                self.atom[iindex].writex[idim] = self.atom[iindex].x[idim]
                
        atomshift = []
        for i, iindex in enumerate(fragment):
            atomshift.append([iindex,[False,False,False]])

        done = False
        while done == False:
            shift = False
            for iindex in fragment:
                iatom = self.atom[iindex]
                for jindex in iatom.neighborlist:
                    if iindex > jindex:
                        jatom = self.atom[jindex]
                        for idim in range(self.dim):
                            val = abs(iatom.writex[idim] - jatom.writex[idim])
                            if val > 0.5:
                                #print('before:', iatom.label, iatom.x, jatom.label, jatom.x)
                                #print(idim, iatom.label, iatom.x, jatom.label, jatom.x)
                                if iatom.writex[idim] > jatom.writex[idim]:
                                    for i, dim in atomshift:
                                        if i == iindex:
                                            if dim[idim] == False:
                                                iatom.writex[idim] -= 1.0
                                                shift = True
                                                dim[idim] = True
                                                break
                                else:
                                    for i, dim in atomshift:
                                        if i == jindex:
                                            if dim[idim] == False:
                                                jatom.writex[idim] -= 1.0
                                                shift = True
                                                dim[idim] = True
                                                break
                                #print('shift:', iindex, jindex)
                                #print('after:', iatom.label, iatom.x, jatom.label, jatom.x, '\n')
                                break
                    if shift:
                        break
                if shift:
                    break

            if shift == False:
                done = True
                
        del atomshift
        return fragment
    
    def update_label(self, label):
        
        num = ''
        for ichar in label:
            if ichar.isdigit():
                num += ichar
        return 'X' + num
    
    def get_label(self, iindex, case):
        
        status = False
        label = ''
        
        '''
        Cases:
        
        0 - compmetalnodelist
        1 - complinkerlist
        # - others
        '''
        
        if case == 0:
            for inode, jnode in self.compcapairlist:
                if status == False:
                    if iindex == jnode:
                        status = True
                        label = self.update_label(self.atom[iindex].label)
            return status, label
        elif case == 1:
            for inode, jnode in self.compcapairlist:
                if status == False:
                    if iindex == inode:
                        status = True
                        label = self.update_label(self.atom[iindex].label)
            return status, label
        else:
            return status, label

    def write_xyz(self, fragment, xyzfile, case):
        
        fragment = self.wrap_fragment(fragment)
        f = open(xyzfile,'w')
        f.write("%-4d\n\n" %len(fragment))
    
        for iindex in fragment:
            iatom = self.atom[iindex]
            status, label = self.get_label(iindex, case)
            if status:
                f.write('Ar')
            else:
                f.write(iatom.symbol)
                
            xc = self.fractional_to_cartesian(iatom.writex)    
            for idim in range(self.dim):
                f.write('\t%-10.6f' %xc[idim])
            f.write('\n')
                
        f.close()
        
    def write_cif(self, fragment, ciffile, case):

        f = open(ciffile,'w')
        f.write(self.ciffile + '\n')
        today = datetime.date.today()
        date = today.strftime("%Y-%m-%d")
        f.write('_audit_creation_date\t\t\t%s\n' %date)
        del date
        f.write('_audit_creation_method\t\t\t\'mfg\'\n')
        f.write('_symmetry_space_group_name_H-M\t\t\'P1\'\n')
        f.write('_symmetry_Int_Tables_number\t\t1\n')
        f.write('_symmetry_cell_setting\t\t\ttriclinic\n')
        f.write('loop_\n')
        f.write('_symmetry_equiv_pos_as_xyz\n')
        f.write('\tx,y,z\n')
        f.write('_cell_length_a\t\t\t\t%-10.6f\n' %self.lx[0])
        f.write('_cell_length_b\t\t\t\t%-10.6f\n' %self.lx[1])
        f.write('_cell_length_c\t\t\t\t%-10.6f\n' %self.lx[2])
        ad = self.ar[0]/self.lib.a2r
        f.write('_cell_angle_alpha\t\t\t%-10.6f\n' %ad)
        ad = self.ar[1]/self.lib.a2r
        f.write('_cell_angle_beta\t\t\t%-10.6f\n' %ad)
        ad = self.ar[2]/self.lib.a2r
        f.write('_cell_angle_gamma\t\t\t%-10.6f\n' %ad)
        f.write('loop_\n')
        f.write('_atom_site_label\n')
        f.write('_atom_site_type_symbol\n')
        f.write('_atom_site_fract_x\n')
        f.write('_atom_site_fract_y\n')
        f.write('_atom_site_fract_z\n')
        f.write('_atom_site_U_iso_or_equiv\n')
        f.write('_atom_site_adp_type\n')
        f.write('_atom_site_occupancy\n')
        f.write('_atom_site_charge\n')

        for iindex in fragment:
            iatom = self.atom[iindex]
            status, label = self.get_label(iindex, case)

            if status:
                f.write(label)
            else:
                f.write(iatom.label)
            f.write('\t')
            if status:
                f.write('Ar')
            else:
                f.write(iatom.symbol)
            for idim in range(self.dim):
                f.write('\t%-10.6f' %iatom.writex[idim])
            f.write('\t0.00000\tUsio\t1.00\t0.00\n')
        f.write('loop_\n')
        f.write('_geom_bond_atom_site_label_1\n')
        f.write('_geom_bond_atom_site_label_2\n')
        f.write('_geom_bond_distance\n')
        f.write('_geom_bond_site_symmetry_2\n')
        f.write('_ccdc_bond_type\n')

        bondpairlist = []
        for iindex in fragment:
            status1, label = self.get_label(iindex, case)
            if status1 == False:
                for ineighbor in self.atom[iindex].neighborlist:
                    if iindex > ineighbor:
                        ibondpair = [ineighbor, iindex]
                    else:
                        ibondpair = [iindex, ineighbor]
                    if ibondpair not in bondpairlist:
                        bondlength = self.get_distance(self.atom[iindex].writex,self.atom[ineighbor].writex)
                        f.write(self.atom[iindex].label)
                        f.write('\t')
                        status2, label = self.get_label(ineighbor, case)
                        if status2:
                            f.write(label)
                        else:
                            f.write(self.atom[ineighbor].label)
                        f.write('\t%-10.6f\t.\tS\n' %bondlength)
                        bondpairlist.append(ibondpair)
        del bondpairlist
        f.close()
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--ciffile', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--cif2cell', type=str)
    
    args = parser.parse_args()
    
    mofname = args.ciffile
    if args.outdir[len(args.outdir)-1] == '/':
        outdir = args.outdir
    else:
        outdir = args.outdir + '/'
    
    if args.cif2cell:
        cif2cell = args.cif2cell
    else:
        cif2cell = "cif2cell"
        
    os.system(f"cp {mofname} mbud_{mofname}")
    os.system(f"{cif2cell} -f mbud_{mofname} -p cif -o {mofname}")
    os.system(f"rm mbud_{mofname}")

    iMOF = MOF(mofname)
    iMOF.get_boxinfo()
    iMOF.get_hmatrix()
    iMOF.get_atomtypelist()
    iMOF.get_metaltypelist()
    iMOF.get_atomgridinfo()
    iMOF.clear_neighborlist()
    iMOF.get_neighborlist(True)
    iMOF.get_solvent()
    iMOF.break_mof()

    if len(iMOF.compmetalnodelist) > 0:
        uniqmetalnodelist = iMOF.get_uniq_fragmentlist(iMOF.compmetalnodelist, 2)
        count = 0
        for uniq in uniqmetalnodelist:
            xyzfile = outdir + 'comp-node-' + str(count) + '.xyz'
            iMOF.write_xyz(uniq, xyzfile, 0)

            ciffile = outdir + 'comp-node-' + str(count) + '.cif'
            iMOF.write_cif(uniq, ciffile, 0)
            count += 1

    if len(iMOF.complinkerlist) > 0:
        uniqlinkerlist = iMOF.get_uniq_fragmentlist(iMOF.complinkerlist, 2)
        count = 0
        for uniq in uniqlinkerlist:
            xyzfile = outdir + 'comp-linker-' + str(count) + '.xyz'
            iMOF.write_xyz(uniq, xyzfile, 1)

            ciffile = outdir + 'comp-linker-' + str(count) + '.cif'
            iMOF.write_cif(uniq, ciffile, 1)
            count += 1

    if len(iMOF.metalnodelist) > 0:
        uniqmetalnodelist = iMOF.get_uniq_fragmentlist(iMOF.metalnodelist, 0)
        count = 0
        for uniq in uniqmetalnodelist:
            xyzfile = outdir + 'node-' + str(count) + '.xyz'
            iMOF.write_xyz(uniq, xyzfile, 2)

            ciffile = outdir + 'node-' + str(count) + '.cif'
            iMOF.write_cif(uniq, ciffile, 2)
            count += 1

    if len(iMOF.linkerlist) > 0:
        uniqlinkerlist = iMOF.get_uniq_fragmentlist(iMOF.linkerlist, 1)
        count = 0
        for uniq in uniqlinkerlist:
            xyzfile = outdir + 'linker-' + str(count) + '.xyz'
            iMOF.write_xyz(uniq, xyzfile, 2)

            ciffile = outdir + 'linker-' + str(count) + '.cif'
            iMOF.write_cif(uniq, ciffile, 2)
            count += 1

    if len(iMOF.funcgrouplist) > 0:
        uniqfuncgrouplist = iMOF.get_uniq_fragmentlist(iMOF.funcgrouplist, 1)
        count = 0
        for uniq in uniqfuncgrouplist:
            xyzfile = outdir + 'func-' + str(count) + '.xyz'
            iMOF.write_xyz(uniq, xyzfile, 2)

            ciffile = outdir + 'funcgroup-' + str(count) + '.cif'
            iMOF.write_cif(uniq, ciffile, 2)
            count += 1

    if len(iMOF.solventlist) > 0:
        uniqsolventlist = iMOF.get_uniq_fragmentlist(iMOF.solventlist, 2)
        count = 0
        for uniq in uniqsolventlist:
            xyzfile = outdir + 'solvent-' + str(count) + '.xyz'
            iMOF.write_xyz(uniq, xyzfile, 2)

            ciffile = outdir + 'solvent-' + str(count) + '.cif'
            iMOF.write_cif(uniq, ciffile, 2)
            count += 1
    iMOF.destroy_all()
