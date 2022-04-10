"""

    !=====================================================!
    !                                                     !
    !          Infrared Spectra Simulation (IRSS)         !
    !                     Yi-Fan Hou                      !
    !                                                     !
    ! State Key Laboratory of Physical Chemistry of Solid !
    ! Surfaces, Fujian Provincial Key Laboratory of Theo- !
    ! retical and Computational Chemistry, Department of  !
    ! Chemistry, and College of Chemistry and Chemical    !
    !           Engineering, Xiamen University,           !
    !                Xiamen 361005, China                 !
    !                                                     !
    !=====================================================!


"""


import time, sys, os, parse
import numpy as np
import stopper
from doc import Doc

from options import IRSS_ML, IRSS_AIMD

class IRSS(object):
    def __init__(self):
        
        print(__doc__)
        
        if len(sys.argv[1:]) == 1:
            if os.path.exists(sys.argv[1]):
                print("\t\tInput file:")
                with open(sys.argv[1],'r') as inpf:
                    rawargs = inpf.readlines()
                    inpf.seek(0)
                    for line in inpf:
                        print(line.rstrip())
                print()
                args = parse.parseArgs(rawargs)
            elif sys.argv[1].lower() in ['help','-h','-help','--help']:
                Doc.printDoc()
            else: 
                stopper.stopIRSS("Input file %s not found"%(sys.argv[1]))
        
        if args['aimd_only']:
            IRSS_AIMD(args)
        else:
            IRSS_ML(args)
if __name__ == '__main__':
    IRSS()