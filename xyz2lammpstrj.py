import numpy as np
import os

class TrjGenerator():

    def __init__(self,working_dir,title,notes='',index='06'):
        self.title = title 
        self.wdir = working_dir
        self.notes = notes 
        self.index = index

        if self.title=='QC':
            self.xyzfile = os.path.join(self.wdir,'02/coord_array.npy')
        elif self.title=='ML':
            self.xyzfile = os.path.join(self.wdir,'04/MLcoord.npy')
        self.atomfile= os.path.join(self.wdir,'00/atom.dat')
        self.xyz = np.load(self.xyzfile)
        self.atom = np.loadtxt(self.atomfile,dtype=str)
        self.atom = [each for each in self.atom]


    def main(self):
        minvalue = []
        for iatom in range(len(self.atom)):
            for icoord in range(3):
                minvalue.append(min([self.xyz[i][iatom][icoord] for i in range(len(self.xyz))]))

        self.xyz = self.xyz - min(minvalue)


        with open(os.path.join(self.wdir,self.index+'/%s.lammpstrj'%(self.title)),'w') as trjf:
            for istep in range(len(self.xyz)):
                trjf.write('ITEM: TIMESTEP\n')
                trjf.write('%d\n'%istep)
                trjf.write('ITEM: NUMBER OF ATOMS\n')
                trjf.write('%d\n'%len(self.atom))
                trjf.write('ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
                trjf.write('1000.000   0.000   0.000\n')
                trjf.write('0.000   1000.000   0.000\n')
                trjf.write('0.000   0.000   1000.000\n')
                trjf.write('ITEM: ATOMS id type x y z\n')
                for iatom in range(len(self.atom)):
                    trjf.write('%d %d %.7f %.7f %.7f\n'%(iatom+1,self.atom.index(self.atom[iatom])+1,self.xyz[istep][iatom][0],self.xyz[istep][iatom][1],self.xyz[istep][iatom][2]))
                
