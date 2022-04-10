# 02 Ab intio Molecular Dynamics
# This script intends to run molecular dynamics with the help of Gaussian

import numpy as np 
import os, sys


class QCMD():

    def __init__(self,dt,nsteps,scaling,Ntotal,working_dir,gauss_path,nproc,keywords,notes='',index='02'):
        # Input options
        # Specify working directory (current directory by default)
        self.working_dir = working_dir

        # Path to Gaussian
        self.gauss_path = gauss_path

        # Number of CPU
        self.nproc = nproc

        # Keywords for single point energy calculation, must include 'force' to calculate forces
        self.keywords = keywords

        # Time step (Unit: fs)
        self.dt = dt

        # Total number of steps 
        self.nsteps = nsteps

        if Ntotal < nsteps:
            self.Ntotal = self.nsteps 
        else:
            self.Ntotal = Ntotal

        self.scaling = scaling

        self.notes = notes
        self.index = index

    # Define functions
    def readCom(self,input_file):
        file_name = input_file.split('/')[-1][:-4]
        rawdata = np.genfromtxt(input_file,usecols=[0,1,2,3],dtype=str,skip_header=6)
        #print(rawdata)
        atom = [rawdata[iatom][0] for iatom in range(len(rawdata))]
        coord = [[eval(rawdata[iatom][i]) for i in range(1,4)] for iatom in range(len(rawdata))]
        
        return file_name, atom, coord

    def writeCom(self,nproc,param,output_file,atom,coord):
        param_line = param
        title_line = output_file
        charge_line = '0 1'
        with open(output_file+'.com','w') as fcom:
            fcom.write('%'+'nproc=%d\n'%nproc)
            fcom.write(param_line+'\n')
            fcom.write('\n')
            fcom.write(title_line+'\n')
            fcom.write('\n')
            fcom.write(charge_line+'\n')
            for iatom in range(len(atom)):
                fcom.write(' %s\t%f\t%f\t%f\n'%(atom[iatom],coord[iatom][0],coord[iatom][1],coord[iatom][2]))
            fcom.write('\n')

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

    def readGaussOut(self,input_file,Natom):
        with open(input_file,'r') as fgauss:
            lines = fgauss.readlines()
            for line in lines:
                # Read forces (Hartree/Bohr)
                if 'Forces (Hartree' in line:
                    f_index = lines.index(line)+3
                    with open('temp','w') as temp:
                        for i in range(Natom):
                            temp.write(lines[f_index+i])
                    forces = np.genfromtxt('temp',usecols=[2,3,4],dtype=float)
                    #print(forces)
                # Read charges
                if 'Mulliken charges:' in line:
                    c_index = lines.index(line)+2
                    with open('temp','w') as temp:
                        for i in range(Natom):
                            temp.write(lines[c_index+i])
                    charges = np.genfromtxt('temp',usecols=[2],dtype=float)
                    #print(charges)
                # Read single point energy & dipole moment
                if '1\\1\\' in line:
                    archive_index = lines.index(line)
                    tmpline = line
                    archive = ''
                    while(tmpline!='\n'):
                        archive += tmpline[1:-1]
                        archive_index += 1
                        tmpline = lines[archive_index]
                    archive_split = archive.split('\\')
                    for each in archive_split:
                        # Single point energy 
                        if 'HF=' in each:
                            energy = eval(each[3:])
                            #print(energy)
                        if 'Dipole=' in each:
                            dipole = each[7:].split(',')
                            dipole = [eval(dipole[i]) for i in range(3)]
                            #print(dipole)
        return energy, forces, charges, dipole

    def QCdynamics(self,wdir,file_name,dt,nsteps,atom,init_coord,init_velocity): # Verlet algorithm
        coord_array = []
        energy_array = []
        forces_array = []
        charges_array = [] 
        dipole_array = []
        
        coord = init_coord # (Angstrom)
        velocity = init_velocity # (Angstrom/fs)
        
        atom_mass = [self.ele2mass(each) for each in atom]
        
        
        atom_mass = np.array(atom_mass).reshape(len(atom),1)
        coord = np.array(coord)
        velocity = np.array(velocity)
        
        # Initialize forces (Hartree/Angstrom)
        forces = np.zeros([len(atom),3])
        # Acceleration(Anstrom/fs^2)
        acceleration = forces / atom_mass / 1822.0 * (0.52917706**2) * (100.0/2.4189)**2    
        
        with open(wdir+self.index+'.log','a') as logf:
            for istep in range(1,nsteps+1):
                
                comfile = wdir+file_name+str(istep).zfill(5)
                self.writeCom(self.nproc,self.keywords,comfile,atom,coord)
                print('Submit %s.com to g09'%(comfile))
                sys.stdout.flush()
                os.system('$GAUSS_EXEDIR/g09 %s.com'%(comfile))
                energy, forces, charges, dipole = self.readGaussOut(comfile+'.log',len(atom))
            
                forces = np.array(forces) / 0.52917706  #Convert from Hartree/Bohr to Hartree/Angstrom
                
                #acceleration = forces / atom_mass / 1822.0 * (0.52917706**2) * (100.0/2.4189)**2
                acceleration = forces / atom_mass * 0.2626744409
                
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
                


                coord_array.append(coord)
                energy_array.append(energy)
                charges_array.append(charges)
                dipole_array.append(dipole)    
                forces_array.append(forces)
                
                # Verlet algorithm
                if istep==1:
                    coord_old = coord - velocity * dt
                coord_new = 2*coord - coord_old + acceleration*(dt**2)
                velocity = (coord_new-coord_old) / (2*dt)
                coord_old=coord 
                coord=coord_new   

                kin_en = np.sum(velocity**2 * atom_mass)
                kin_en = kin_en * 1822.0 / 2.0 / (0.52917706**2) / (100.0/2.4189)**2
                

                logf.write('    Velocity (Angstrom/fs): \n')
                for iatom in range(len(atom)):
                    logf.write('         %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],velocity[iatom][0],velocity[iatom][1],velocity[iatom][2]))
                logf.write('\n')
                logf.write('='*60+'\n\n')
                logf.write('    Kinetic energies   : %.7f \n'%(kin_en))
                logf.write('    Electronic energies: %.7f \n'%(energy))
                logf.write('    Total energies     : %.7f \n'%(kin_en+energy))
                logf.write('\n'+'='*60+'\n\n')

                

        
        print('MD calculations terminated normally')
        
        return coord_array, energy_array, forces_array, charges_array, dipole_array

    def perturbation(self,wdir,file_name,coord_array,Nsample,scaling,atom):
        coord_sample = [] 
        energy_sample = [] 
        force_sample = [] 
        charges_sample = [] 
        dipole_sample = []
        
        np.random.seed()
        
        atom_mass = [self.ele2mass(each) for each in atom]
        atom_mass = np.array(atom_mass).reshape(len(atom),1)
        
        Nbase = len(coord_array)
        
        print('Nbase = %d'%Nbase)
        
        Ntimes = Nsample // Nbase 
        Nrest = Nsample - Ntimes*Nbase
        
        for i in range(Nbase):
            for j in range(Ntimes):
                randnum = np.random.randn(len(atom),3) * scaling / np.sqrt(atom_mass)
                coord = np.array(coord_array[i]) + randnum
                
                comfile = wdir+file_name+str(Nbase+j+1+Ntimes*i).zfill(5)
                self.writeCom(self.nproc,self.keywords,comfile,atom,coord)
                print('Submit %s.com to g09'%(comfile))
                os.system('$GAUSS_EXEDIR/g09 %s.com'%(comfile))
                energy, forces, charges, dipole = self.readGaussOut(comfile+'.log',len(atom))
                
                coord_sample.append(coord)
                energy_sample.append(energy)
                force_sample.append(forces)
                charges_sample.append(charges)
                dipole_sample.append(dipole)
                
        for i in range(Nrest):
            randnum = np.random.randn(len(atom),3) * scaling / np.sqrt(atom_mass)
            coord = np.array(coord_array[i]) + randnum
                
            comfile = wdir+file_name+str(Nbase*(Ntimes+1)+i+1).zfill(5)
            self.writeCom(self.nproc,self.keywords,comfile,atom,coord)
            print('Submit %s.com to g09'%(comfile))
            os.system('$GAUSS_EXEDIR/g09 %s.com'%(comfile))
            energy, forces, charges, dipole = self.readGaussOut(comfile+'.log',len(atom))
                
            coord_sample.append(coord)
            energy_sample.append(energy)
            force_sample.append(forces)
            charges_sample.append(charges)
            dipole_sample.append(dipole)
        
        print('Perturbation sample completed')
        return coord_sample,energy_sample,force_sample,charges_sample,dipole_sample

    def main(self):
        # Main Code
        if self.working_dir == '':
            self.working_dir = os.getcwd()
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index + '/*')

        with open(self.working_dir+'00/file_name.dat','r') as namef:
            file_name = namef.readline()[:-1]
        atom = np.loadtxt(self.working_dir+'01/atom.dat',dtype=str)
        init_coord = np.load(self.working_dir+'01/QC_init_coord.npy')
        init_velocity = np.load(self.working_dir+'01/QC_init_velocity.npy')

        # index.log
        with open(self.working_dir+self.index+'/'+self.index+'.log','w') as logf:
            logf.write('Input Options:\n')
            logf.write('    working_dir=%s\n'%self.working_dir)
            logf.write('    gauss_path=%s\n'%self.gauss_path)
            logf.write('    nproc=%d\n'%self.nproc)
            logf.write('    keywords=%s\n'%self.keywords)
            logf.write('    Time step: %f fs\n'%self.dt)
            logf.write('    Number of steps: %d'%self.nsteps)
            logf.write('\n')
            logf.write('Main:\n')
            
        coord_array,energy_array,forces_array,charges_array,dipole_array = self.QCdynamics(self.working_dir+self.index+'/',file_name,self.dt,self.nsteps,atom,init_coord,init_velocity)
        coord_sample,energy_sample,forces_sample,charges_sample,dipole_sample = self.perturbation(self.working_dir+self.index+'/',file_name,coord_array,self.Ntotal-self.nsteps,self.scaling,atom)

        coord_array = np.array(coord_array+coord_sample)
        energy_array = np.array(energy_array+energy_sample)
        forces_array = np.array(forces_array+forces_sample)
        charges_array = np.array(charges_array+charges_sample)
        dipole_array = np.array(dipole_array+dipole_sample)


        np.save(self.working_dir+self.index+'/coord_array.npy',coord_array)
        np.save(self.working_dir+self.index+'/energy_array.npy',energy_array)
        np.save(self.working_dir+self.index+'/forces_array.npy',forces_array)
        np.save(self.working_dir+self.index+'/charges_array.npy',charges_array)
        np.save(self.working_dir+self.index+'/dipole_array.npy',dipole_array)


        # 02.log 
        with open(self.working_dir+self.index+'/'+self.index+'.log','a') as logf:
            logf.write('Output:\n')
            logf.write('    QCMD geometries saved in %s\n'%(self.working_dir+self.index+'/coord_array.npy'))
            logf.write('    QCMD energies saved in %s\n'%(self.working_dir+self.index+'/energy_array.npy'))
            logf.write('    QCMD forces saved in %s\n'%(self.working_dir+self.index+'/forces_array.npy'))
            logf.write('    QCMD charges saved in %s\n'%(self.working_dir+self.index+'/charges_array.npy'))
            logf.write('    QCMD dipoles saved in %s\n'%(self.working_dir+self.index+'/dipole_array.npy'))
