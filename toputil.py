from basics import *

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

def readDirective(s):
    # Remove comments, then split normally
    sline = s.split(';')[0].split()
    
    if len(sline) == 3 and sline[0] == '[' and sline[2] == ']':
        return sline[1]
    else:
        return ''

