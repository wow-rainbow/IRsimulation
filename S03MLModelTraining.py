# 03 ML Model Training 
# This script intends to train ML model using dataset generated from QCMD

import numpy as np
import os

class MLModelTraining():
    def __init__(self,mlatom,working_dir,notes='',index='03'):
        # Input Options 
        # Specify working directory (current directory by default)
        self.working_dir = working_dir

        # MLatom.py 
        self.mlatom = mlatom

        self.notes=notes
        self.index=index

    # Define Functions
    def prepareMLatomData(self,wdir,atom,coord,energy,forces,charges):
        # XYZ file
        with open(wdir+'xyz.dat','w') as xyzf:
            for xyz in coord:
                xyzf.write(str(len(atom))+'\n\n')
                for iatom in range(len(xyz)):
                    xyzf.write('%s\t%f\t%f\t%f\n'%(atom[iatom],xyz[iatom][0],xyz[iatom][1],xyz[iatom][2]))

        # Energy file
        with open(wdir+'en.dat','w') as enf:
            for en in energy:
                enf.write('%.10f\n'%(en))

        # Grad file (Gradients should be negative forces)
        with open(wdir+'grad.dat','w') as gradf:
            for grad in forces:
                gradf.write(str(len(atom))+'\n\n')
                for iatom in range(len(atom)):
                    gradf.write('\t%.10f\t%.10f\t%.10f\n'%(-grad[iatom][0],-grad[iatom][1],-grad[iatom][2]))

        # Charge file
        for icharge in range(len(atom)):
            with open(wdir+'charge'+str(icharge+1).zfill(2)+'.dat','w') as chargef:
                for each in charges.T[icharge]:
                    chargef.write('%.8f\n'%(each))
                    
        return len(charges)
                    
    def trainMLmodel(self,wdir,Ndata,file_name, MLtask):
        if MLtask == 'learnValGradXYZ':
            with open(wdir+file_name+'.inp','w') as engradf:
                engradf.write('createMLmodel\n')
                engradf.write('KRRtask=%s\n'%MLtask)
                engradf.write('MLmodelOut=%s.unf\n'%(wdir+file_name))
                engradf.write('XYZfile=%sxyz.dat\n'%wdir)
                engradf.write('Yfile=%sen.dat\n'%wdir)
                engradf.write('YgradXYZfile=%sgrad.dat\n'%wdir)
                engradf.write('sampling=random\n')
                engradf.write('Ntrain=%d Nsubtrain=%d\n'%(Ndata,int(Ndata*0.8)))
                engradf.write('kernel=Gaussian\n')
                engradf.write('molDescriptor=ID\n')
                engradf.write('molDescrType=unsorted\n')
                engradf.write('sigma=opt lambda=opt\n')
            os.system('python %s %s.inp > %s.out'%(self.mlatom,(wdir+file_name),(wdir+file_name)))
        elif MLtask == 'learnGradXYZ':
            with open(wdir+file_name+'.inp','w') as engradf:
                engradf.write('createMLmodel\n')
                engradf.write('KRRtask=%s\n'%MLtask)
                engradf.write('MLmodelOut=%s.unf\n'%(wdir+file_name))
                engradf.write('XYZfile=%sxyz.dat\n'%wdir)
                engradf.write('YgradXYZfile=%sgrad.dat\n'%wdir)
                engradf.write('sampling=random\n')
                engradf.write('Ntrain=%d Nsubtrain=%d\n'%(Ndata,int(Ndata*0.8)))
                engradf.write('kernel=Gaussian\n')
                engradf.write('molDescriptor=ID\n')
                engradf.write('molDescrType=unsorted\n')
                engradf.write('sigma=opt lambda=opt\n')
            os.system('python %s %s.inp > %s.out'%(self.mlatom,(wdir+file_name),(wdir+file_name)))
        elif MLtask == 'learnVal':
            with open(wdir+file_name+'.inp','w') as engradf:
                engradf.write('createMLmodel\n')
                engradf.write('KRRtask=%s\n'%MLtask)
                engradf.write('MLmodelOut=%s.unf\n'%(wdir+file_name))
                engradf.write('XYZfile=%sxyz.dat\n'%wdir)
                engradf.write('Yfile=%s.dat\n'%(wdir+file_name))
                engradf.write('sampling=random\n')
                engradf.write('Ntrain=%d Nsubtrain=%d\n'%(Ndata,int(Ndata*0.8)))
                engradf.write('kernel=Gaussian\n')
                engradf.write('molDescriptor=ID\n')
                engradf.write('molDescrType=unsorted\n')
                engradf.write('sigma=opt lambda=opt\n')
            os.system('python %s %s.inp > %s.out'%(self.mlatom,(wdir+file_name),(wdir+file_name)))
            

    def main(self):
        # Main Code
        if self.working_dir == '':
            self.working_dir = os.getcwd()
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index+'/*')

        # Load data
        coord_array = np.load(self.working_dir+'02/coord_array.npy')
        energy_array = np.load(self.working_dir+'02/energy_array.npy')
        forces_array = np.load(self.working_dir+'02/forces_array.npy')
        charges_array = np.load(self.working_dir+'02/charges_array.npy')
        dipole_array = np.load(self.working_dir+'02/dipole_array.npy')
        atom = np.loadtxt(self.working_dir+'00/atom.dat',dtype=str)

        # Prepare MLatom-type data file 
        Ndata = self.prepareMLatomData(self.working_dir+self.index+'/',atom,coord_array,energy_array,forces_array,charges_array)

        # Train ML Model
        # ML PES
        self.trainMLmodel(self.working_dir+self.index+'/',Ndata,'engrad','learnValGradXYZ')

        # ML charges
        for i in range(len(atom)):
            self.trainMLmodel(self.working_dir+self.index+'/',Ndata,'charge'+str(i+1).zfill(2),'learnVal')
            
        # index.log
        with open(self.working_dir+self.index+'/'+self.index+'.log','w') as logf:
            logf.write('Input Options:\n')
            logf.write('    working_dir=%s\n'%self.working_dir)
            logf.write('    mlatom=%s\n'%self.mlatom)
            logf.write('    additional input=')
            logf.write('Output:\n')
            logf.write('    ML PES model saved in %s\n'%(self.working_dir+self.index+'/engrad.unf'))
            for iatom in range(len(atom)):
                logf.write('    ML charge%s model saved in %s\n'%(str(iatom).zfill(2),self.working_dir+self.index+'/charges'+str(iatom).zfill(2)+'.unf'))
            logf.write('\n')
            logf.write('See %s'%(self.working_dir)+self.index+'/*.out for MLatom output files\n')
            
