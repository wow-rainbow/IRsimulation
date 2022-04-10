import stopper

class Doc():
    @classmethod 
    def printDoc(self):
        doc = """
    ____________________________________________________________________________

        Input options:

            gauss_path=S            Path to Gaussian
            mlatompy=S              Path to MLatom.py
            mlatomf=S               Path to MLatomF
            gauss_inp               Path to Gaussian input file 
            working_dir=S           Path to working directory
                '' [default]            Current directory
                <user-defined>
            nproc=N                 Number of CPU
                4 [default]
                <user-defined>
            optwords=S              Keywords for optimization
                # opt b3lyp/6-311++g(d,p) [default]
                <user-defined>          Must include '# opt'
            qc_temp=R               Temperature for AI molecular dynamics
                300.0 [default]         Unit: K
                <user-defined>
            ml_temp=R               Temperature for ML molecular dynamics
                300.0 [default]         Unit: K
                <user-defined>
            spwords=S               Keywords for single point energy calculation
                '# b3lyp/6-311++g(d,p) force' [default]
                <user-defined>          Must include '# force'
            qc_dt=R                 Time step of AIMD
                1.0 [default]           Unit: fs
                <user-defined>
            qc_nsteps=N             Number of steps of AIMD
                500 [default]
                <user-defined>
            ml_dt=R                 Time step of MLMD
                0.1 [default]           Unit: fs
                <user-defined>
            ml_nsteps=N             Number of steps of MLMD
                30000 [default]
                <user-defined>
            qc_traj=B               Generate AIMD trajectory
                True [default]
                <user-defined>
            ml_traj=B               Generate MLMD trajectory
                True [default]
                <user-defined>

    ____________________________________________________________________________    
        """
        print(doc)
        

