# MkVsites

# Copyright (c) 2020 Erik Marklund

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


from basics import *
from toputil import *
from os import path
from sys import stderr

class ffAtom(atom):
    """Forcefield atom class"""
    def __init__(self, name='', type='', element=''):
        atom.__init__(self, name=name, type=type)
        self.element = ''

    def setElement(self, element):
        self.element = element

    def isSameType(self, a):
        return self.type == a.type

class ffAtomType(boolBase):
    """Forcefield atom type class"""
    def __init__(self, type='', atnum=0, mass=0, charge=0.0, ptype='', sigma=0.0, epsilon=0.0):
        self.type    = type
        self.atnum   = atnum
        self.mass    = mass
        self.charge  = charge
        self.ptype   = ptype
        self.sigma   = sigma
        self.epsilon = epsilon
        
    def readFFAtomTypeLine(self, line):
        """Read forcefield atom type"""
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
        """Derive element from atom number"""
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
        """Return mass"""
        return self.mass

class ffBondType(boolBase):
    """Class for bond types"""
    def __init__(self, atom=[], func=0, b0=0.0, kb=0.0):
        self.atom = atom[:]
        self.func = func
        self.b0 = b0
        self.kb = kb

    def readFFBondTypeLine(self, line):
        """Read forcefield bond type"""
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
    """Class for constraint types"""
    def __init__(self, atom=[], func=0, b0=0.0):
        self.atom = atom[:]
        self.func = func
        self.b0 = b0

    def readFFConstraintTypeLine(self, line):
        """Read forcefield constraint type"""
        if readDirective(line):
            raise(FFError)
        
        sline = line.split(';')[0].split()
        if not sline or len(sline) != 4:
            raise(FFError)
        
        self.atom = [sline[0], sline[1]]
        self.func = int(sline[2])
        self.b0   = float(sline[3])
        
class ffAngleType(boolBase):
    """Class for angle types"""
    def __init__(self, atom=[], func=0, theta0=0.0, ktheta=0.0, ub0=0.0, kub=0.0):
        self.atom = atom[:]
        self.func = func
        self.theta0 = theta0
        self.ktheta = ktheta
        self.ub0 = ub0
        self.kub = kub
        
    def readFFAngleTypeLine(self, line):
        """Read forcefield angle type"""
        if readDirective(line):
            raise(FFError)
        
        sline = line.split(';')[0].split()
        if not sline or len(sline) < 6:
            raise(FFError)

        self.atom   = [sline[0], sline[1], sline[2]]
        self.func   = int(sline[3])
        self.theta0 = float(sline[4])
        self.ktheta = float(sline[5])
        if len(sline) > 6:
            self.ub0    = float(sline[6])
        if len(sline) > 7:
            self.kub    = float(sline[7])

        
class forceField(boolBase):
    """Class for forcefield data"""
    def __init__(self, ff='', bVerbose=False):
        self.ff = ff
        self.atoms = []
        self.bonds = []
        self.constraints = []
        self.angles = []
        self.vsites = []
        self.bVerbose = bVerbose
        self.uniqueVsites = []

    def setPath(self, d):
        """Set path to forcefield directory"""
        self.ff = d

    def getPath(self):
        """Returns path to forcefield"""
        return self.ff
        
    def read(self):
        """Read forcefield data"""
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
        """Get atom of given type"""
        for at in self.atoms:
            if at.type == n:
                return at
        # If not found by now, return False
        # TODO: RAISE ERROR INSTEAD
        return ffAtomType()
        
    def getElement(self, a):
        """Return element for an atom"""
        for at in self.atoms:
            if at.type == a.type:
                return at.getElement()
        # If not found by now, return False
        return ''

    def getMass(self, a):
        """Get mass of an atom"""
        for at in self.atoms:
            if at.type == a.type:
                return at.getMass()
        # If not found by now, return False
        return 0.0

    def getBond(self, i, j):
        """Get bond length for two atoms"""
        for b in self.bonds:
            if \
              (b.atom[0] == i and b.atom[1] == j) or \
              (b.atom[0] == j and b.atom[1] == i):
              return b.b0

        # If not found by now, return False
        return 0.0

    def getConstraint(self, i, j):
        """Get constraint length for two atoms"""
        for c in self.constraints:
            if \
              (c.atom[0] == i and c.atom[1] == j) or \
              (c.atom[0] == j and c.atom[1] == k):
              return c.b0

        # If not found by now, return False
        return 0.0
    
    def getAngle(self, i, j, k):
        """Get angle for three atoms"""
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
    
    def gatherVsites(self):
        """Gather vsites. EXPLAIN MORE!"""
        for a in self.atoms:
            if not isDummy(a):
                continue

            v = ffVsite(DummyName=a.type)

            # find constraints
            for c in self.constraints:
                if a.type in c.atom:
                    # AnchorConstraint or DummyConstraint?
                    if c.atom[0] == c.atom[1]:
                        v.setDummyConstraint(c.b0)
                    else:
                        try:
                            if a.type == c.atom[0]:
                                anchor = c.atom[1]
                            else:
                                anchor = c.atom[0]
                                
                            v.addAnchorConstraint(anchor, c.b0)
                            
                        except FFError:
                            v.dump(o=stderr)

            if v.isComplete():
                self.vsites.append(v)
            else:
                if self.bVerbose:
                    warn('Discarding incomplete vsite:', bWarn=False)
                    v.dump(stderr)

    def getVsiteType(self, v):
        """Returns the ffVsite if match, otherwise None"""
        for vf in self.vsites:
            vt = vt.getVsiteType()
            if vt:
                return vt

        return None


class ffVsiteType(boolBase):
    """Class for vsite types"""
    def __init__(self, DummyName='', DummyConstraint=0.0, Anchor='', AnchorConstraint=0.0, Center=''):
        self.DummyName = DummyName
        self.DummyConstraint = DummyConstraint
        self.Anchor = Anchor
        self.AnchorConstraint = AnchorConstraint
        self.Center = Center

    def getDummyName(self):
        return self.DummyName

    def dump(self, ostream=stdout):
        output('DummyName        {:10s}'.format(self.DummyName), ostream=ostream)
        output('DummyConstraint  {:<10.6f}'.format(self.DummyConstraint), ostream=ostream)
        output('Anchor           {:10s}'.format(self.Anchor), ostream=ostream)
        output('AnchorConstraint {:<10.6f}'.format(self.AnchorConstraint), ostream=ostream)
        output('Center           {:10s}'.format(self.Center), ostream=ostream)

class ffVsite(boolBase):
    """Class for forcefield vsites"""
    def __init__(self, DummyName='', DummyConstraint=0.0, Anchors=[], AnchorConstraints=[]):
        self.DummyName = DummyName
        self.DummyConstraint = DummyConstraint
        self.Anchors = Anchors[:]
        self.AnchorConstraints = AnchorConstraints[:]

    def setDummyName(self, name):
        """Set the name of dummies"""
        self.DummyName = name

    def getDummyName(self):
        """Get the name of the dummies"""
        return self.DummyName
        
    def setDummyConstraint(self, DummyConstraint):
        """Set inter-dummy constraint"""
        if self.DummyConstraint:
            warn('Dummy constraint already set.')
            warn('  was: {:f}'.format(self.DummyConstraint), bWarn=False)
            warn('   is: {:f}'.format(DummyConstraint), bWarn=False)
            
        self.DummyConstraint = DummyConstraint

    def addAnchorConstraint(self, Anchor, AnchorConstraint):
        """Set dymmy-anchor constraints"""
        if Anchor not in self.Anchors:
            self.Anchors.append(Anchor)
            self.AnchorConstraints.append(AnchorConstraint)
        else:
            warn('Ignoring duplicate Anchor: {:s}, {:f}'.format(Anchor, AnchorConstraint))
            raise FFError

    def isComplete(self):
        """Checks if a vsite is complete"""
        return self.DummyName and self.DummyConstraint and self.AnchorConstraints

    def dump(self, o=stdout):
        """Prints a vsite to stdout (or other stream)"""
        output('Name:            {:s}'.format(self.DummyName), ostream=o)
        output('DummyConstraint: {:f}'.format(self.DummyConstraint), ostream=o)
        for a, c in zip(self.Anchors, self.AnchorConstraints):
            output('Anchor {:s}: {:f}'.format(a, c), ostream=o)

    def getVsiteType(self, v):
        """Takes a topology.vsite and compares with self."""
        # Allow 1% discrepancy (Order of 1/100 A)
        if (abs(v.dDD-self.DummyConstraint)/self.DummyConstraint > 0.01 or \
            v.Center != self.Center or \
            v.Anchor != self.Anchor):
            return None

        for a, ac in zip(self.Anchors, self.AnchorConstraints):
            if abs(v.dDHeavy-ac)/ac < 0.01:
                return ffVsiteType(DummyName=self.DummyName,
                                       DummyConstraint=self.DummyConstraint,
                                       Anchor=a, AnchorConstraint=ac)

        return None

        
def isDummy(a):
    """True if atom is a dummy atom. False otherwise."""
    return (a.mass == 0 and a.atnum == 0 and a.charge == 0.0 and a.sigma == 0.0 and a.epsilon == 0.0)
