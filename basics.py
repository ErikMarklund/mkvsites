from sys import stderr, stdout

# Some basic functions

class boolBase(object):
    def __nonzero__(self):
        V = vars(self)
        for v in V.values():
            if v:
                return True
        return False

def warn(s, bEndline=True, bWarn=True, bError=False, ostream=stderr, bPrint=True):
    if not bPrint:
        return

    if bEndline:
        s += '\n'

    if bError:
        s = 'ERROR: '+s
    elif bWarn:
        s = 'WARNING: '+s
    else:
        pass

    stderr.write(s)
    
    
def output(s, bEndline=True, ostream=stdout, bPrint=True):
    if not bPrint:
        return

    if bEndline:
        s += '\n'

    ostream.write(s)
