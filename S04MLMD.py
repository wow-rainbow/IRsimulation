# 04 Machine Learning Molecular Dynamics
# This script intends to run molecular dynamics with machine learning PES

import numpy as np 
import os, sys

class MLMD():
    def __init__(self,dt,nsteps,working_dir,mlatom,nproc,notes='',index='04'):
        # Input Options
        # Specify working directory (current directory by default)
        self.working_dir = working_dir

        # MLatom.py 
        self.mlatom = mlatom

        # Time step (Unit: fs)
        self.dt = dt

        # Total number of steps 
        self.nsteps = nsteps

        self.nproc=nproc

        self.notes=notes 
        self.index=index

    # Define functions
    def ele2mass(self,ele):
        ele2mass_dic = {'H':1.008,
                    'He':4.0026,
                    'Li':6.94,
                    'Be':9.0122,
                    'B':10.81,
                    'C':12.011,
                    'N':14.007,
                    'O':15.999,
                    'F':18.998,
                    'Ne':20.18,
                    'Na':22.99,
                    'Mg':24.305,
                    'Al':26.982,
                    'Si':28.085,
                    'P':30.974,
                    'S':32.06,
                    'Cl':35.45,
                    'Ar':39.948,
                    'K':39.098,
                    'Ca':40.078,}
        return ele2mass_dic.get(ele)

    def calcDipole(self,charges,coord):
        dipolex, dipoley, dipolez = 0.0, 0.0, 0.0
        Tcoord = np.array(coord).T
        
        for i in range(len(charges)):
            dipolex += charges[i] * Tcoord[0][i]
            dipoley += charges[i] * Tcoord[1][i]
            dipolez += charges[i] * Tcoord[2][i]
            
            dipole = np.sqrt(dipolex**2+dipoley**2+dipolez**2)
            
        return dipole

    def writeXYZ(self,wdir,atom,coord):
        with open(wdir+'xyz.dat','w') as xyzf:
            xyzf.write('%d\n\n'%(len(atom)))
            for iatom in range(len(atom)):
                xyzf.write('%s\t%f\t%f\t%f\n'%(atom[iatom],coord[iatom][0],coord[iatom][1],coord[iatom][2]))

    def useMLmodelInp(self,wdir,file_name,MLpred):
        if MLpred == 'grad':
            with open(wdir+file_name+'.inp','w') as gradf:
                gradf.write('nthreads=%d\n'%self.nproc)
                gradf.write('useMLmodel\n')
                gradf.write('MLmodelIn=%s\n'%(wdir+file_name+'.unf'))
                gradf.write('XYZfile=%sxyz.dat\n'%wdir)
                gradf.write('YgradXYZestfile=%sest.dat\n'%(wdir+file_name))
                gradf.write('kernel=Gaussian\n')
                gradf.write('molDescriptor=ID\n')
                gradf.write('molDescrType=unsorted\n')
        if MLpred == 'en':
            with open(wdir+file_name+'.inp','w') as gradf:
                gradf.write('nthreads=%d\n'%self.nproc)
                gradf.write('useMLmodel\n')
                gradf.write('MLmodelIn=%s\n'%(wdir+file_name+'.unf'))
                gradf.write('XYZfile=%sxyz.dat\n'%wdir)
                gradf.write('Yestfile=%sest.dat\n'%(wdir+file_name))
                gradf.write('kernel=Gaussian\n')
                gradf.write('molDescriptor=ID\n')
                gradf.write('molDescrType=unsorted\n')
                
                
    def useMLmodel(self,wdir,atom):
        os.system('rm %s*est.dat'%(wdir))
        # Predict forces
        #os.system('python %s %s.inp > %s.out'%(mlatom,wdir+'engrad',wdir+'engrad'))
        os.system('%s useMLmodel MLmodelIn=%sengrad.unf XYZfile=%sxyz.dat YgradXYZestFile=%sengradest.dat YestFile=%senest.dat kernel=Gaussian molDescriptor=ID molDescrType=unsorted > %sengrad.out'%(self.mlatom,wdir,wdir,wdir,wdir,wdir))
        with open(wdir+'enest.dat','r') as enf:
            energy = eval(enf.readline()[:-1])
        forces = np.genfromtxt(wdir+'engradest.dat',dtype=float,skip_header=2)
        forces = -forces # Forces are negative gradients
        
        # Predict charges
        charges = []
        for iatom in range(1,1+len(atom)):
            #os.system('python %s %s.inp > %s.out'%(mlatom,wdir+'charge'+str(iatom).zfill(2),wdir+'charge'+str(iatom).zfill(2)))
            os.system('%s useMLmodel MLmodelIn=%s.unf XYZfile=%sxyz.dat YestFile=%sest.dat kernel=Gaussian molDescriptor=ID molDescrType=unsorted > %s.out'%(self.mlatom,wdir+'charge'+str(iatom).zfill(2),wdir,wdir+'charge'+str(iatom).zfill(2),wdir+'charge'+str(iatom).zfill(2)))
            with open(wdir+'charge%sest.dat'%(str(iatom).zfill(2)),'r') as chf:
                charges.append(eval(chf.readline()[:-1]))
        charges = np.array(charges)
        
        return forces, charges, energy

        
    def MLdynamics(self,working_dir,dt,nsteps,atom,init_coord,init_velocity):
        coord_array = []
        dipole_array = [] 
        
        atom_mass = [self.ele2mass(each) for each in atom]
        atom_mass = np.array(atom_mass).reshape(len(atom),1)
        
        coord = init_coord
        velocity = init_velocity
        
        # Initialize forces (Hartree/Angstrom)
        forces = np.zeros([len(atom),3])
        # Acceleration(Anstrom/fs^2)
        acceleration = forces / atom_mass / 1822.0 * (0.52917706**2) * (100.0/2.4189)**2
        
        with open(working_dir+self.index+'/'+self.index+'.log','a') as logf:
            for istep in range(1,nsteps+1):
                if (istep%100 == 0):
                    print('Step %d ...'%(istep))
                    sys.stdout.flush() 

                self.writeXYZ(working_dir+self.index+'/',atom,coord)
                forces, charges, energy= self.useMLmodel(working_dir+self.index+'/',atom)

                

                acceleration = forces / atom_mass / 1822.0 * (0.52917706**2) * (100.0/2.4189)**2
                
                # index.log
                logf.write('    Step %-6d : t=%f fs \n'%(istep,istep*dt))
                logf.write('    Geometry (Angstrom): \n')
                for iatom in range(len(atom)):
                    logf.write('         %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],coord[iatom][0],coord[iatom][1],coord[iatom][2]))
                logf.write('    Forces (Hartree/Angstrom): \n')
                for iatom in range(len(atom)):
                    logf.write('         %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],forces[iatom][0],forces[iatom][1],forces[iatom][2]))
                logf.write('    Mulliken charges: \n')
                for iatom in range(len(atom)):
                    logf.write('         %s\t%.7f\n'%(atom[iatom],charges[iatom]))
                

                dipole = self.calcDipole(charges,coord)
                dipole_array.append(dipole)
                coord_array.append(coord)

                # Verlet algorithm
                if istep==1:
                    coord_old = coord - velocity * dt 
                coord_new = 2*coord - coord_old + acceleration*(dt**2)
                velocity = (coord_new-coord_old) / (2*dt)
                coord_old = coord 
                coord = coord_new

                kin_en = np.sum(velocity**2 * atom_mass)
                kin_en = kin_en * 1822.0 / 2.0 / (0.52917706**2) / (100.0/2.4189)**2

                logf.write('    Velocity (Angstrom/fs): \n')
                for iatom in range(len(atom)):
                    logf.write('         %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],velocity[iatom][0],velocity[iatom][1],velocity[iatom][2]))
                logf.write('\n')
                logf.write('\n')
                logf.write('='*60+'\n\n')
                logf.write('    Kinetic energies   : %.7f \n'%(kin_en))
                logf.write('    Electronic energies: %.7f \n'%(energy))
                logf.write('    Total energies     : %.7f \n'%(kin_en+energy))
                logf.write('\n'+'='*60+'\n\n')
            
        dipole_array = np.array(dipole_array)
        coord_array = np.array(coord_array)
        
        return coord_array, dipole_array
        
    def main(self):   
        # Main Code
        if self.working_dir == '':
            self.working_dir = os.getcwd()
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index+'/*')
        os.system('cp '+self.working_dir+'03/*.unf '+self.working_dir+self.index+'/')

        # Load data
        atom = np.loadtxt(self.working_dir+'00/atom.dat',dtype=str)
        ML_init_coord = np.load(self.working_dir+'01/ML_init_coord.npy')
        ML_init_velocity = np.load(self.working_dir+'01/ML_init_velocity.npy')

        # Generate MLatom input files
        self.useMLmodelInp(self.working_dir+self.index+'/','engrad','grad')
        for i in range(len(atom)):
            self.useMLmodelInp(self.working_dir+self.index+'/','charge'+str(i+1).zfill(2),'en')

        # 02.log
        with open(self.working_dir+self.index+'/'+self.index+'.log','w') as logf:
            logf.write('Input Options:\n')
            logf.write('    working_dir=%s\n'%self.working_dir)
            logf.write('    mlatom=%s\n'%self.mlatom)
            logf.write('    Time step: %f fs\n'%self.dt)
            logf.write('    Number of steps: %d\n'%self.nsteps)
            logf.write('\n')
            logf.write('Main:\n')
            
        # ML dynamics
        coords, dipoles = self.MLdynamics(self.working_dir,self.dt,self.nsteps,atom,ML_init_coord,ML_init_velocity)
        np.save(self.working_dir+self.index+'/MLdipoles.npy',dipoles)
        np.save(self.working_dir+self.index+'/MLcoord.npy',coords)

        # index.log 
        with open(self.working_dir+self.index+'/'+self.index+'.log','a') as logf:
            logf.write('Output:\n')
            logf.write('    MLMD geometries saved in %s\n'%(self.working_dir+self.index+'/MLcoord.npy'))
            logf.write('    MLMD dipoles saved in %s\n'%(self.working_dir+self.index+'/MLdipoles.npy'))
