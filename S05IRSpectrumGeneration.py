# 05 IR Spectrum Generation
# This script intends to generate IR spectrum using MD generated dipoles 
# Autocorrelation function of dipoles with respect to time is used

import os
import numpy as np
import matplotlib.pyplot as plt 
import numpy.fft as nf

class IRSpectrumGenerator():

    def __init__(self,dt,nsteps,working_dir,notes='',index='05'):
        # Input Options
        self.dt = dt
        self.nsteps = nsteps
        self.working_dir = working_dir

        self.notes=notes 
        self.index=index

    def main(self):

        # Main code
        if self.working_dir == '':
            self.working_dir = './'
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index+'/*')
        if self.notes=='aimd_only':
            os.system('cp '+self.working_dir+'02/dipole_array.npy '+self.working_dir+self.index+'/')
            dipoles = np.load(os.path.join(self.working_dir+self.index+'/dipole_array.npy'))
            dipoles = [np.sqrt(each[0]**2+each[1]**2+each[2]**2) for each in dipoles]
        else:
            os.system('cp '+self.working_dir+'04/MLdipoles.npy '+self.working_dir+self.index+'/')
            dipoles = np.load(os.path.join(self.working_dir+self.index+'/MLdipoles.npy'))
        #dipoles = np.load('MLdipoles.npy')
        dipole_prime = [(dipoles[i+1]-dipoles[i])/self.dt for i in range(self.nsteps-1)]
        dipole_prime = np.array(dipole_prime)
        time = np.linspace(0.0,(self.nsteps-2)*self.dt,self.nsteps-1)

        comp_arr = nf.fft(dipole_prime)

        freqs = nf.fftfreq(comp_arr.size,self.dt)

        pows = np.array(abs(comp_arr))
        max_pow = max(pows)
        pows = 1.0-(pows / max_pow)

        pows = pows[freqs>0]
        freqs = freqs[freqs>0]

        freqs = freqs*1.0e15/3e10


        fig, ax = plt.subplots(1,1)
        ax.plot(freqs[freqs<4000], pows[freqs<4000])
        ax.invert_xaxis()
        plt.xlabel('$cm^{-1}$')
        plt.ylabel('$Transmition$')

        plt.savefig(os.path.join(self.working_dir,self.index+'/ir.png'))
