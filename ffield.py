from basics import *
from toputil import *
from os import path

class ffAtom(atom):
    def __init__(self, name='', type='', element=''):
        atom.__init__(self, name=name, type=type)
        self.element = ''

    def setElement(self, element):
        self.element = element

    def isSameType(self, a):
        return self.type == a.type

class ffAtomType(boolBase):
    def __init__(self, type='', atnum=0, mass=0, charge=0.0, ptype='', sigma=0.0, epsilon=0.0):
        self.type    = type
        self.atnum   = atnum
        self.mass    = mass
        self.charge  = charge
        self.ptype   = ptype
        self.sigma   = sigma
        self.epsilon = epsilon
        
    def readFFAtomTypeLine(self, line):
        if readDirective(line):
            raise(FFError)
        
        sline = line.split(';')[0].split()

        if not sline or len(sline) != 7:
            raise(FFError)

        self.type    = sline[0]
        self.atnum   = int(sline[1])
        self.mass    = float(sline[2])
        self.charge  = float(sline[3])
        self.ptype   = sline[4]
        self.sigma   = float(sline[5])
        self.epsilon = float(sline[6])
    
    def getElement(self):
        if self.atnum == 1:
            return 'H'
        elif self.atnum == 6:
            return 'C'
        elif self.atnum == 7:
            return 'N'
        elif self.atnum == 8:
            return 'O'
        elif self.atnum == 16:
            return 'S'
        elif self.atnum == 15:
            return 'P'
        else:
            warn('Unknown element (number {:d})'.format(self.atnum))
            return ''

    def getMass(self):
        return self.mass

class ffBondType(boolBase):
    def __init__(self, atom=[], func=0, b0=0.0, kb=0.0):
        self.atom = atom
        self.func = func
        self.b0 = b0
        self.kb = kb

    def readFFBondTypeLine(self, line):
        if readDirective(line):
            raise(FFError)

        sline = line.split(';')[0].split()
        if not sline or len(sline) != 5:
            raise(FFError)

        self.atom = [sline[0], sline[1]]
        self.func = int(sline[2])
        self.b0   = float(sline[3])
        self.kb   = float(sline[4])

class ffConstraintType(boolBase):
    def __init__(self, atom=[], func=0, b0=0.0):
        self.atom = atom
        self.func = func
        self.b0 = b0

    def readFFConstraintTypeLine(self, line):
        if readDirective(line):
            raise(FFError)
        
        sline = line.split(';')[0].split()
        if not sline or len(sline) != 4:
            raise(FFError)
        
        self.atom = [sline[0], sline[1]]
        self.func = int(sline[2])
        self.b0   = float(sline[3])
        
class ffAngleType(boolBase):
    def __init__(self, atom=[], func=0, theta0=0.0, ktheta=0.0, ub0=0.0, kub=0.0):
        self.atom = atom
        self.func = func
        self.theta0 = theta0
        self.ktheta = ktheta
        self.ub0 = ub0
        self.kub = kub
        
    def readFFAngleTypeLine(self, line):
        if readDirective(line):
            raise(FFError)
        
        sline = line.split(';')[0].split()
        if not sline or len(sline) != 8:
            raise(FFError)

        self.atom   = [sline[0], sline[1], sline[2]]
        self.func   = int(sline[3])
        self.theta0 = float(sline[4])
        self.ktheta = float(sline[5])
        self.ub0    = float(sline[6])
        self.kub    = float(sline[7])


class forceField(boolBase):
    def __init__(self, ff=''):
        self.ff = ff
        self.atoms = []
        self.bonds = []
        self.constraints = []
        self.angles = []
        self.vsites = []

    def setDirectory(self, d):
        self.ff = d
        
    def read(self):
        # Read in atomtypes from ffnonbonded.itp
        fname = path.join(self.ff, 'ffnonbonded.itp')
        preprocessor = Cpp()
        bSkip = False
        f = preprocessor.parse(fname)
        reading = 'None'
        for line in f:
            if not line.split(';')[0].split():
                continue

            if line and line.lstrip()[0] == '#':
                continue

            directive = readDirective(line)
            if directive:
                reading = directive
                continue

            if reading == 'atomtypes':
                ffat = ffAtomType()

                try:
                    ffat.readFFAtomTypeLine(line)
                    self.atoms.append(ffat)
                except FFError:
                    warn('Unexpected line under atomtypes:')
                    warn(line, bWarn=False)

            else:
                continue
                
        # Read in bonded types from ffbonded.itp
        fname = path.join(self.ff, 'ffbonded.itp')
        f = preprocessor.parse(fname)
        reading = 'None'
        for line in f:
            if not line.split(';')[0].split():
                continue

            if line and line.lstrip()[0] == '#':
                continue

            directive = readDirective(line)
            if directive:
                reading = directive
                continue

            elif reading == 'bondtypes':
                ffbo = ffBondType()

                try:
                    ffbo.readFFBondTypeLine(line)
                    self.bonds.append(ffbo)
                except FFError:
                    warn('Unexpected line under bondtypes:')
                    warn(line, bWarn=False)

            elif reading == 'constrainttypes':
                ffco = ffConstraintType()

                try:
                    ffco.readFFConstraintTypeLine(line)
                    self.constraints.append(ffco)
                except FFError:
                    warn('Unexpected line under constrainttypes:')
                    warn(line, bWarn=False)

            elif reading == 'angletypes':
                ffan = ffAngleType()

                try:
                    ffan.readFFAngleTypeLine(line)
                    self.angles.append(ffan)
                except FFError:
                    warn('Unexpected line under angletypes:')
                    warn(line, bWarn=False)
            else:
                continue


    def getAtom(self, n):
        for at in self.atoms:
            if at.type == n:
                return at
        # If not found by now, return False
        return ffAtomType()
        
    def getElement(self, a):
        for at in self.atoms:
            if at.type == a.type:
                return at.getElement()
        # If not found by now, return False
        return ''

    def getMass(self, a):
        for at in self.atoms:
            if at.type == a.type:
                return at.getMass()
        # If not found by now, return False
        return 0.0

    def getBond(self, i, j):
        for b in self.bonds:
            if \
              (b.atom[0] == i and b.atom[1] == j) or \
              (b.atom[0] == j and b.atom[1] == i):
              return b.b0

        # If not found by now, return False
        return 0.0

    def getConstraint(self, i, j):
        for c in self.constraints:
            if \
              (c.atom[0] == i and c.atom[1] == j) or \
              (c.atom[0] == j and c.atom[1] == k):
              return c.b0

        # If not found by now, return False
        return 0.0
    
    def getAngle(self, i, j, k):
        for a in self.angles:
            if a.atom[1] == j and \
              (
                  # [i, j, k] and [k, j, i] are equivalent 
                  (a.atom[0] == i and a.atom[2] == k) or \
                  (a.atom[0] == k and a.atom[2] == i) \
              ):
                return a.theta0

        # If not found by now, return False
        return 0.0
    
    #def gatherVsites(self):

        

