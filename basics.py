from sys import stderr, stdout

# Some basic functions

class boolBase():
    def __nonzero__(self):
        V = vars(self)
        for v in V.values():
            if v:
                return True
        return False

def warn(s, bEndline=True, bWarn=True, bError=False):
    if bEndline:
        s += '\n'

    if bError:
        s = 'ERROR: '+s
    elif bWarn:
        s = 'WARNING: '+s
    else:
        pass

    stderr.write(s)
    
    
def output(s, bEndline=True):
    if bEndline:
        s += '\n'

    stdout.write(s)
