# 01 Initial Condition Generation
# This scripts intend to generate initial conditions for QMMD (should also for MLMD)

import numpy as np
import os 


class InitialConditionGeneration():

    def __init__(self,QC_temp,ML_temp,gauss_path,working_dir,notes='',index='01'):
        # Input options
        # Specify working directory (current directory by default)
        self.working_dir = working_dir

        # Path to Gaussian
        self.gauss_path = gauss_path

        # Temperature for QC molecular dynamics (Unit: K)
        self.QC_temp = QC_temp

        # Temperature for ML molecular dynamics (Unit: K)
        self.ML_temp = ML_temp

        self.notes=notes
        self.index=index

    # Define functions
    def readCom(self,input_file):
        file_name = input_file.split('/')[-1][:-4]
        rawdata = np.genfromtxt(input_file,usecols=[0,1,2,3],dtype=str,skip_header=6)
        #print(rawdata)
        atom = [rawdata[iatom][0] for iatom in range(len(rawdata))]
        coord = [[eval(rawdata[iatom][i]) for i in range(1,4)] for iatom in range(len(rawdata))]
        
        return file_name, atom, coord

    def readOptXYZ(self,input_file,Natom):
        with open(input_file,'r') as fgauss:
            lines=fgauss.readlines()
            for line in lines:
                if 'Standard orientation:' in line:
                    index = lines.index(line)+5
                    with open('temp','w') as tempf:
                        for i in range(Natom):
                            tempf.write(lines[index+i])
                    coord = np.genfromtxt('temp',usecols=[3,4,5],dtype=float)
                    os.system('rm temp')
        return coord

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

    def randVelocity(self,temp,atom):
        np.random.seed()
        randnum = np.random.randn(len(atom),3)
        mass = np.array([self.ele2mass(each) for each in atom])
        mass = np.array(mass).reshape(len(atom),1)
        kb = 1.380649E-23  # Unit: J/K
        hartree = 27.21070 * 1.602176565E-19 # Unit: J/Hartree
        kb = kb / hartree # Unit: Hartree/K
        
        total_energy = 3.0/2.0 * kb * temp * len(atom)
        rand_velocity = randnum / np.sqrt(mass*1822)

        rand_energy = np.sum((rand_velocity**2) * (mass*1822)) / 2.0

        ratio = rand_energy / total_energy

        velocity = rand_velocity / np.sqrt(ratio) # Unit: a.u.

        bohr = 0.5291772083 # 1 Bohr = 0.5291772083 Angstrom
        time = 2.4189E-2 # Time Unit: 1 a.u. = 2.4189E-17 s = 2.4189E-2 fs

        velocity = velocity * bohr / time # Unit: Angstrom / fs

        return velocity


    def main(self):
        # Main Code
        if self.working_dir == '':
            self.working_dir = os.getcwd()
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index + '/*')

        listdir00 = os.listdir(self.working_dir+'00/')
        for eachfile in listdir00:
            if ('.log' in eachfile) and (eachfile!='00.log'):
                gauss_out = self.working_dir+'00/'+eachfile

        # Load data
        os.system('cp ' + self.working_dir+'00/atom.dat' + ' ' + self.working_dir + self.index+'/')
        atom = np.loadtxt(self.working_dir+self.index+'/atom.dat',dtype=str)
        init_coord = self.readOptXYZ(gauss_out,len(atom))

        # Generate velocities
        QC_init_velocity = self.randVelocity(self.QC_temp,atom)
        if not self.notes=='aimd_only':
            ML_init_velocity = self.randVelocity(self.ML_temp,atom)

        np.save(self.working_dir+self.index+'/QC_init_coord.npy',init_coord)
        np.save(self.working_dir+self.index+'/QC_init_velocity.npy',QC_init_velocity)

        if not self.notes=='aimd_only':
            np.save(self.working_dir+self.index+'/ML_init_coord.npy',init_coord)
            np.save(self.working_dir+self.index+'/ML_init_velocity.npy',ML_init_velocity)

        # index.log
        with open(self.working_dir+self.index+'/'+self.index+'.log','w') as logf:
            logf.write('Input Options:\n')
            logf.write('    working_dir=%s\n'%self.working_dir)
            logf.write('    Temperature for QC molecular dynamics:  %d K\n'%self.QC_temp)
            if not self.notes=='aimd_only':
                logf.write('    Temperature for ML molecular dynamics:  %d K\n'%self.ML_temp)
            logf.write('    Initial geometry for QC molecular dynamics saved in %s\n'%(self.working_dir+self.index+'/QC_init_coord.npy'))
            logf.write('    Initial velocity for QC molecular dynamics saved in %s\n'%(self.working_dir+self.index+'/QC_init_velocity.npy'))
            if not self.notes=='aimd_only':
                logf.write('    Initial geometry for ML molecular dynamics saved in %s\n'%(self.working_dir+self.index+'/ML_init_coord.npy'))
                logf.write('    Initial velocity for ML molecular dynamics saved in %s\n'%(self.working_dir+self.index+'/ML_init_velocity.npy'))
            logf.write('\n')
            logf.write('Main:\n')
            logf.write('    Initial geometry for QC molecular dynamics (Angstrom):\n')
            for iatom in range(len(atom)):
                logf.write('        %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],init_coord[iatom][0],init_coord[iatom][1],init_coord[iatom][2]))
            logf.write('\n')
            logf.write('    Initial velocity for QC molecular dynamics (Angstrom/fs):\n')
            for iatom in range(len(atom)):
                logf.write('        %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],QC_init_velocity[iatom][0],QC_init_velocity[iatom][1],QC_init_velocity[iatom][2]))
            logf.write('\n')
            if not self.notes=='aimd_only':
                logf.write('    Initial geometry for ML molecular dynamics (Angstrom):\n')
                for iatom in range(len(atom)):
                    logf.write('        %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],init_coord[iatom][0],init_coord[iatom][1],init_coord[iatom][2]))
                logf.write('\n')
                logf.write('    Initial velocity for ML molecular dynamics (Angstrom/fs):\n')
                for iatom in range(len(atom)):
                    logf.write('        %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],ML_init_velocity[iatom][0],ML_init_velocity[iatom][1],ML_init_velocity[iatom][2]))
