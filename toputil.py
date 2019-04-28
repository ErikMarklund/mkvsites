from basics import *
import subprocess as sp

class TopologyError(Exception):
    pass
class RtpError(Exception):
    pass
class FFError(Exception):
    pass


class topBase(boolBase):
    """Base class for itp and rtp"""

    def __init__(self, name='', atoms=[], bonds=[], angles=[]):
        self.name = ''
        self.atoms = atoms[:]
        self.bonds = bonds[:]
        self.angles = angles[:]

        
class atom(object):
    """Atom class"""
    
    def __init__(self, name='', type=''):
        self.name = name
        self.type = type

    def isSameType(self, atom):
        return self.type == atom.type

    def dump(self, ostream=stdout):
        output("{:8s}{:8s}".format(self.name, self.type), ostream=ostream)
    

class bond(object):
    """Bond class"""

    def __init__(self, atoms=[]):
        if len(atoms) != 0 and len(atoms) != 2:
            raise(TopologyError)
        
        else:
            self.atoms = atoms


class angle(object):
    """Angle class"""

    def __init__(self, atoms=[]):
        if len(atoms) !=0 and len(atoms) != 3:
            raise(TopologyError)
        
        else:
            self.atoms = atoms


class Cpp():
    """C-preprocessor wrapper class"""
    def __init__(self, defs=[], exe='cpp', flags='-E -P'):
        self.exe = exe[:]
        self.flags = flags[:]

    def parse(self, fname):
        """Returns a list of lines that can be read like an open file."""
        cmd = [self.exe] + self.flags.split() + [fname]
        p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        so, se = p.communicate()

        if so:
            return so.splitlines()
        else:
            return []
        

def readDirective(s):
    """Reads a directive line from forcefield, itp, or rtp file."""
    # Remove comments, then split normally
    sline = s.split(';')[0].split()
    
    if len(sline) == 3 and sline[0] == '[' and sline[2] == ']':
        return sline[1]
    else:
        return ''
