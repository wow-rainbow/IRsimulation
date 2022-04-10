import time, sys, os, parse
import numpy as np
import stopper
from doc import Doc

from S00DataPreparation import DataPreparation
from S01InitialConditionGeneration import InitialConditionGeneration
from S02QCMD import QCMD
from S03MLModelTraining import MLModelTraining
from S04MLMD import MLMD
from S05IRSpectrumGeneration import IRSpectrumGenerator
from S06TrajGeneration import TrajGeneration


class IRSS_ML():
    def __init__(self,args):
        starttime = time.time()
        print(time.strftime(" IRSS started on %d.%m.%Y at %H:%M:%S",time.localtime()))
        print()

        # Info 
        if args['qc_traj']:
            info_qc_traj = 'Generate AIMD trajectory'
        else: 
            info_qc_traj = 'Do not generate AIMD trajectory' 
        
        if args['ml_traj']:
            info_ml_traj = 'Generate MLMD trajectory'
        else: 
            info_ml_traj = 'Do not generate MLMD trajectory'

        info = """
    __________________________________________________________________

    Simulations with the following options:

    Interface options:
        Gaussian program: %s
        MLatom.py       : %s 
        MLatomF         : %s
            
    Number of CPU used:   %d

    $ Step 0 $ Data preparation options:
        Initial structure             : %s
        Gaussian optimization keywords: %s

    $ Step 1 $ Generate initial conditions:
        Temperature for ab initio molecular dynamics       :  %f K
        Temperature for machine learning molecular dynamics:  %f K

    $ Step 2 $ Perform ab initio molecular dynamics 
        Gaussian single point energy calculation keywords: %s
        MD time step                                     : %f fs
        Number of steps                                  : %d

    $ Step 3 $ Train machine learning models

    $ Step 4 $ Perform machine learning molecular dynamics
        MD time step                                     : %f fs
        Number of steps                                  : %d

    $ Step 5 $ Generate IR spectrum

    $ Step 6 $ Generate MD trajectories [optional]
        %s
        %s

    __________________________________________________________________


        """%(args['gauss_path'],args['mlatompy'],args['mlatomf'],args['nproc'],args['gauss_inp'],args['optwords'],args['qc_temp'],args['ml_temp'],args['spwords'],args['qc_dt'],args['qc_nsteps'],args['ml_dt'],args['ml_nsteps'],info_qc_traj,info_ml_traj)

        print(info)
        sys.stdout.flush()


        kwargs00 = {
            'gauss_inp'   : args['gauss_inp'],
            'working_dir' : args['working_dir'],
            'gauss_path'  : args['gauss_path'],
            'nproc'       : args['nproc'],
            'keywords'    : args['optwords'] 
        }

        kwargs01 = {
            'working_dir' : args['working_dir'],
            'gauss_path'  : args['gauss_path'],
            'QC_temp'     : args['qc_temp'],
            'ML_temp'     : args['ml_temp']
        }

        kwargs02 = {
            'working_dir' : args['working_dir'],
            'gauss_path'  : args['gauss_path'],
            'nproc'       : args['nproc'],
            'keywords'    : args['spwords'],
            'dt'          : args['qc_dt'],
            'nsteps'      : args['qc_nsteps'],
            'Ntotal'      : args['qc_ntotal'],
            'scaling'     : args['qc_scaling']
        }

        kwargs03 = {
            'working_dir' : args['working_dir'],
            'mlatom'      : args['mlatompy']
            #'nproc'       : args['nproc']
        }

        kwargs04 = {
            'working_dir' : args['working_dir'],
            'mlatom'      : args['mlatomf'],
            'dt'          : args['ml_dt'],
            'nsteps'      : args['ml_nsteps'],
            'nproc'       : args['nproc']
        }

        kwargs05 = {
            'working_dir' : args['working_dir'],
            'dt'          : args['ml_dt'],
            'nsteps'      : args['ml_nsteps']
        }

        kwargs06 = {
            'working_dir' : args['working_dir'],
            'qc'          : args['qc_traj'],
            'ml'          : args['ml_traj']
        }

        wdir = args['working_dir']
        if not os.path.exists(os.path.join(wdir,'record')):
            os.mknod(os.path.join(wdir,'record'))

        with open(os.path.join(wdir,'record'),'r') as recordf:
            raw_record = recordf.readlines()
            record = []
            for each in raw_record:
                each_check = each.rstrip()
                if each_check != '':
                    record.append(each_check)

        if not '00' in record: 
            DataPreparation(**kwargs00).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('00\n')
            endtime0 = time.time()
            wallclock = endtime0 - starttime
            print(' Step 0 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime0 = time.time()
            print(' Skip step 0 ')

        sys.stdout.flush()

        if not '01' in record:
            InitialConditionGeneration(**kwargs01).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('01\n')
            endtime1 = time.time()
            wallclock = endtime1 - endtime0
            print(' Step 1 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime1 = time.time()
            print(' Skip step 1 ')

        sys.stdout.flush()

        if not '02' in record:
            QCMD(**kwargs02).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('02\n')
            endtime2 = time.time()
            wallclock = endtime2 - endtime1
            print(' Step 2 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime2 = time.time()
            print(' Skip step 2 ')

        sys.stdout.flush()

        if not '03' in record:
            MLModelTraining(**kwargs03).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('03\n')
            endtime3 = time.time()
            wallclock = endtime3 - endtime2
            print(' Step 3 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime3 = time.time()
            print(' Skip step 3 ')

        sys.stdout.flush()

        if not '04' in record:
            MLMD(**kwargs04).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('04\n')
            endtime4 = time.time()
            wallclock = endtime4 - endtime3
            print(' Step 4 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime4 = time.time()
            print(' Skip step 4 ')

        sys.stdout.flush()

        if not '05' in record:
            IRSpectrumGenerator(**kwargs05).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('05\n')
            endtime5 = time.time()
            wallclock = endtime5 - endtime4
            print(' Step 5 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime5 = time.time()
            print(' Skip step 5')

        sys.stdout.flush()

        if not '06' in record:
            TrajGeneration(**kwargs06).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('06\n')
            endtime6 = time.time()
            wallclock = endtime6 - endtime5
            print(' Step 6 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime5 = time.time()
            print(' Skip step 6')

        sys.stdout.flush()

        endtime = time.time()
        wallclock = endtime - starttime
        print('='*78)
        print(' Wall-clock time: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        print(time.strftime(" IRSS terminated on %d.%m.%Y at %H:%M:%S", time.localtime()))
        print('='*78)

        sys.stdout.flush()


class IRSS_AIMD():
    def __init__(self,args):
        starttime = time.time()
        print(time.strftime(" IRSS started on %d.%m.%Y at %H:%M:%S",time.localtime()))
        print()

        # Info
        
        if args['qc_traj']:
            info_qc_traj = 'Generate AIMD trajectory'
        else: 
            info_qc_traj = 'Do not generate AIMD trajectory' 

        info = """
    __________________________________________________________________

    !!! Generate IR spectrum using AIMD only !!!

    Simulations with the following options:

    Interface options:
        Gaussian program: %s
            
    Number of CPU used:   %d

    $ Step 0 $ Data preparation options:
        Initial structure             : %s
        Gaussian optimization keywords: %s

    $ Step 1 $ Generate initial conditions:
        Temperature for ab initio molecular dynamics       :  %f K

    $ Step 2 $ Perform ab initio molecular dynamics 
        Gaussian single point energy calculation keywords: %s
        MD time step                                     : %f fs
        Number of steps                                  : %d

    $ Step 3 $ Generate IR spectrum

    $ Step 4 $ Generate MD trajectories [optional]
        %s

    __________________________________________________________________


        """%(args['gauss_path'],args['nproc'],args['gauss_inp'],args['optwords'],args['qc_temp'],args['spwords'],args['qc_dt'],args['qc_nsteps'],info_qc_traj)

        print(info)
        sys.stdout.flush()

        kwargs00 = {
            'gauss_inp'   : args['gauss_inp'],
            'working_dir' : args['working_dir'],
            'gauss_path'  : args['gauss_path'],
            'nproc'       : args['nproc'],
            'keywords'    : args['optwords'] 
        }

        kwargs01 = {
            'working_dir' : args['working_dir'],
            'gauss_path'  : args['gauss_path'],
            'QC_temp'     : args['qc_temp'],
            'ML_temp'     : args['ml_temp'],
            'notes'       : 'aimd_only'
        }

        kwargs02 = {
            'working_dir' : args['working_dir'],
            'gauss_path'  : args['gauss_path'],
            'nproc'       : args['nproc'],
            'keywords'    : args['spwords'],
            'dt'          : args['qc_dt'],
            'nsteps'      : args['qc_nsteps'],
            'Ntotal'      : args['qc_ntotal'],
            'scaling'     : args['qc_scaling']
        }

        kwargs03 = {
            'working_dir' : args['working_dir'],
            'dt'          : args['qc_dt'],
            'nsteps'      : args['qc_nsteps'],
            'notes'       : 'aimd_only',
            'index'       : '03'
        }

        kwargs04 = {
            'working_dir' : args['working_dir'],
            'qc'          : args['qc_traj'],
            'ml'          : False,
            'notes'       : 'aimd_only',
            'index'       : '04'
        }

        wdir = args['working_dir']
        if not os.path.exists(os.path.join(wdir,'record')):
            os.mknod(os.path.join(wdir,'record'))
        
        with open(os.path.join(wdir,'record'),'r') as recordf:
            raw_record = recordf.readlines()
            record = []
            for each in raw_record:
                each_check = each.rstrip()
                if each_check != '':
                    record.append(each_check)

        if not '00' in record: 
            DataPreparation(**kwargs00).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('00\n')
            endtime0 = time.time()
            wallclock = endtime0 - starttime
            print(' Step 0 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime0 = time.time()
            print(' Skip step 0 ')

        sys.stdout.flush()

        if not '01' in record:
            InitialConditionGeneration(**kwargs01).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('01\n')
            endtime1 = time.time()
            wallclock = endtime1 - endtime0
            print(' Step 1 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime1 = time.time()
            print(' Skip step 1 ')

        sys.stdout.flush()

        if not '02' in record:
            QCMD(**kwargs02).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('02\n')
            endtime2 = time.time()
            wallclock = endtime2 - endtime1
            print(' Step 2 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime2 = time.time()
            print(' Skip step 2 ')

        sys.stdout.flush()

        if not '03' in record:
            IRSpectrumGenerator(**kwargs03).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('03\n')
            endtime3 = time.time()
            wallclock = endtime3 - endtime2
            print(' Step 3 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime3 = time.time()
            print(' Skip step 3')

        sys.stdout.flush()

        if not '04' in record:
            TrajGeneration(**kwargs04).main()
            with open(os.path.join(wdir,'record'),'a') as recordf: recordf.write('04\n')
            endtime4 = time.time()
            wallclock = endtime4 - endtime3
            print(' Step 4 finished in: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        else:
            endtime4 = time.time()
            print(' Skip step 4')

        sys.stdout.flush()

        endtime = time.time()
        wallclock = endtime - starttime
        print('='*78)
        print(' Wall-clock time: %.2f s (%.2f min, %.2f hours)\n' % (wallclock, wallclock / 60.0, wallclock / 3600.0))
        print(time.strftime(" IRSS terminated on %d.%m.%Y at %H:%M:%S", time.localtime()))
        print('='*78)

        sys.stdout.flush()

class IRSS_ML_AL():
    def __init__(self,args):
        starttime = time.time()
        print(time.strftime(" IRSS started on %d.%m.%Y at %H:%M:%S",time.localtime()))
        print()

        # Info 
        if args['qc_traj']:
            info_qc_traj = 'Generate AIMD trajectory'
        else: 
            info_qc_traj = 'Do not generate AIMD trajectory' 
        
        if args['ml_traj']:
            info_ml_traj = 'Generate MLMD trajectory'
        else: 
            info_ml_traj = 'Do not generate MLMD trajectory'

        info = """
    __________________________________________________________________

    !!! Training ML model with active learning !!!

    Simulations with the following options:

    Interface options:
        Gaussian program: %s
        MLatom.py       : %s 
        MLatomF         : %s
            
    Number of CPU used:   %d

    $ Step 0 $ Data preparation options:
        Initial structure             : %s
        Gaussian optimization keywords: %s

    $ Step 1 $ Generate initial conditions:
        Temperature for ab initio molecular dynamics       :  %f K
        Temperature for machine learning molecular dynamics:  %f K

    $ Step 2 $ Perform ab initio molecular dynamics 
        Gaussian single point energy calculation keywords: %s
        MD time step                                     : %f fs
        Number of steps                                  : %d

    $ Step 3 $ Train machine learning models

    $ Step 4 $ Perform machine learning molecular dynamics
        MD time step                                     : %f fs
        Number of steps                                  : %d

    $ Step 5 $ Generate IR spectrum

    $ Step 6 $ Generate MD trajectories [optional]
        %s
        %s

    __________________________________________________________________


        """%(args['gauss_path'],args['mlatompy'],args['mlatomf'],args['nproc'],args['gauss_inp'],args['optwords'],args['qc_temp'],args['ml_temp'],args['spwords'],args['qc_dt'],args['qc_nsteps'],args['ml_dt'],args['ml_nsteps'],info_qc_traj,info_ml_traj)

        print(info)
        sys.stdout.flush()

