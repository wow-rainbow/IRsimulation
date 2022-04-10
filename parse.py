from curses import raw
import os, sys
import stopper

def checkArgs(args):
    keys = args.keys()

    all_args = ['working_dir','gauss_inp','gauss_path','nproc','optwords','qc_temp','ml_temp','spwords','qc_dt','qc_nsteps','qc_ntotal','qc_scaling','mlatompy','ml_dt','ml_nsteps','mlatomf','qc_traj','ml_traj','aimd_only','active_learning']

    float_args = ['qc_temp','ml_temp','qc_scaling','qc_dt','ml_dt']
    int_args = ['nproc','qc_nsteps','qc_ntotal','ml_nsteps']
    bool_args = ['qc_traj','ml_traj','aimd_only','active_learning']
    


    # Check arguments validity
    for key in keys:
        if not key in all_args:
            stopper.stopIRSS('Unknown argument %s'%(key))
        if key in float_args:
            try:
                args[key]=float(args[key])
            except:
                stopper.stopIRSS('Argument type of %s should be float'%(key))
        if key in int_args:
            try:
                args[key]=int(args[key])
            except:
                stopper.stopIRSS('Argument type of %s should be int'%(key))
        if key in bool_args:
            if args[key] in ['0','f','false','.false.','.f.']:
                args[key]=False 
            elif args[key] in ['1','t','true','.true.','.t.']:
                args[key]=True
            else: 
                stopper.stopIRSS('Argument type of %s should be bool'%(key))

    # Check necessary arguments
    if 'aimd_only' in keys and args['aimd_only']:
        necessary_args = ['gauss_inp','gauss_path']
    else:
        necessary_args = ['gauss_inp','gauss_path','mlatompy','mlatomf']
    for each in necessary_args: 
        if not each in keys:
            stopper.stopIRSS("Argument %s should be specified"%(each))



def updateArgs(args):
    default_args = {}
    default_args['working_dir'] = ''
    default_args['nproc'] = 4
    default_args['optwords'] = '# opt b3lyp/6-311++g(d,p)'
    default_args['qc_temp'] = 300
    default_args['ml_temp'] = 300
    default_args['spwords'] = '# b3lyp/6-311++g(d,p) force'
    default_args['qc_dt'] = 1.0
    default_args['qc_nsteps'] = 500
    default_args['qc_ntotal'] = 0
    default_args['qc_scaling'] = 0.05
    default_args['ml_dt'] = 0.1
    default_args['ml_nsteps'] = 30000
    default_args['qc_traj'] = True 
    default_args['ml_traj'] = True
    default_args['aimd_only'] = False

    # Update
    for each in args.keys():
        default_args[each] = args[each]

    return default_args        

def parseArgs(args):
    rawargs = {}
    for each in args:
        if each.rstrip() != '':
            each_split = each.split('=')
            #print(each_split)
            rawargs[each_split[0].rstrip()] = each_split[1].rstrip()

    print('654321')
    checkArgs(rawargs)

    print('123456')
    fineargs = updateArgs(rawargs)

    return fineargs
