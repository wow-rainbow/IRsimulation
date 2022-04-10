import sys
from doc import Doc

def stopIRSS(errorMsg):
    if errorMsg != '':
        Doc.printDoc()
        print('<!> %s <!>'%(errorMsg))
    sys.exit(0)
