import os
from xyz2lammpstrj import TrjGenerator

class TrajGeneration():
    def __init__(self,qc,ml,working_dir,notes='',index='06'):
        self.qc = qc 
        self.ml = ml 
        self.working_dir = working_dir

        self.notes=notes 
        self.index=index
    def main(self):
        # Main Code
        if self.working_dir == '':
            self.working_dir = os.getcwd()
        if self.working_dir[-1] != '/':
            self.working_dir += '/'
        if not os.path.exists(self.working_dir+self.index):
            os.mkdir(self.working_dir+self.index)
        os.system('rm -rf '+ self.working_dir + self.index+'/*')

        if self.qc:
            if self.notes=='aimd_only':
                TrjGenerator(self.working_dir,'QC','aimd_only','04').main()
            else:
                TrjGenerator(self.working_dir,'QC').main()
        if self.ml:
            TrjGenerator(self.working_dir,'ML').main()
            
