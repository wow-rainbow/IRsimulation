# 00DataPreparation
#     This script intends to assist you in data preparation, including molecule geometry optimization, 
# using QC calculation software like Gaussian

import numpy as np
import os

class DataPreparation():
    # Initialization
    def __init__(self,gauss_inp,working_dir,gauss_path,nproc,keywords,notes='',index='00'):
        # Input options
        # Path to Gaussian input file (Supported input format: .com or .gjf)
        self.gauss_inp = gauss_inp

        # Specify working directory (current directory by default)
        self.working_dir = working_dir

        # Path to Gaussian
        self.gauss_path = gauss_path

        # Number of CPU
        self.nproc = nproc

        # Keywords for optimization ,must include 'opt' 
        self.keywords = keywords

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
            fcom.write('%'+'nproc=%d\n'%(nproc))
            fcom.write(param_line+'\n')
            fcom.write('\n')
            fcom.write(title_line+'\n')
            fcom.write('\n')
            fcom.write(charge_line+'\n')
            for iatom in range(len(atom)):
                fcom.write(' %s\t%f\t%f\t%f\n'%(atom[iatom],coord[iatom][0],coord[iatom][1],coord[iatom][2]))
            fcom.write('\n')
            
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

    def main(self):
        if self.working_dir == '':
            self.working_dir = os.getcwd()
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index +'/*')
        os.system('cp ' + self.gauss_inp + ' ' + self.working_dir + self.index + '/')
        file_name, atom, coord = self.readCom(self.gauss_inp)
        self.writeCom(self.nproc,self.keywords,self.working_dir+self.index+'/'+file_name+'_opt',atom,coord)
        os.system(self.gauss_path+' '+self.working_dir+self.index+'/'+file_name+'_opt.com')
        with open(self.working_dir+self.index+'/atom.dat','w') as atomf:
            for each in atom:
                atomf.write(each+'\n')
        with open(self.working_dir+self.index+'/file_name.dat','w') as namef:
            namef.write('%s\n'%file_name)

        # index.log
        with open(self.working_dir+self.index+'/'+self.index+'.log','w') as logf:
            logf.write('Input Options:\n')
            logf.write('    working_dir=%s\n'%self.working_dir)
            logf.write('    gauss_path=%s\n'%self.gauss_path)
            logf.write('    nproc=%d\n'%self.nproc)
            logf.write('    keywords=%s\n'%self.keywords)
            logf.write('    gauss_inp=%s\n'%self.gauss_inp)
            logf.write('\n')
            logf.write('Main:\n')
            logf.write('    Gaussian input file for optimization is saved in %s\n'%(self.working_dir+self.index+'/'+file_name+'_opt.com'))
            logf.write('        Gaussian log file: %s\n'%(self.working_dir+self.index+'/'+file_name+'_opt.log'))
            logf.write('\n')
            logf.write('    Input geometry (Angstrom): \n')
            for iatom in range(len(atom)):
                logf.write('        %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],coord[iatom][0],coord[iatom][1],coord[iatom][2]))
            logf.write('\n')
            logf.write('    Optimized geometry (Angstrom): \n')
            opt_coord = self.readOptXYZ(self.working_dir+self.index+'/'+file_name+'_opt.log',len(atom))
            for iatom in range(len(atom)):
                logf.write('        %s\t%.7f\t%.7f\t%.7f\n'%(atom[iatom],opt_coord[iatom][0],opt_coord[iatom][1],opt_coord[iatom][2]))
            logf.write('\n')
