from sys import exit, argv
from os import environ, path
from basics import *
from toputil import *
from math import *
import ffield
import rtp

# Topology class and functions for reading rtp files and itp/top files.

# Class for storing the angle constraints.
class constraint(boolBase):
    def __init__(self, atoms=[], b0=0.0):
        self.atoms = atoms
        self.b0 = b0

    def setConstraint(self, atoms, b0):
        self.atoms = atoms
        self.b0 = b0

    # isSame() uses a ~0.01% tolerance because floating point arithmetic.
    def isSame(self, c):
        return (
                (self.i == c.i and self.j == c.j) or \
                (self.i == c.j and self.j == c.i) \
            ) and abs(self.b0-c.b0)*2/(self.b0+c.b0) < 0.0001

# Class for storing vsites.
class vsite(boolBase):
    def __init__(self, center='', heavy='', nH=0, mDummy=0.0, dDummy=0.0, dDummyHeavy=0.0, name=''):
        self.center = center
        self.heavy  = heavy
        self.nH     = nH
        # Mass of a single dummy mass
        self.mDummy = mDummy
        # Distance between dummy masses
        self.dDummy = dDummy
        # Distance between dummy and heavy anchor atom
        self.dDummyHeavy = dDummyHeavy
        self.name   = name

    def set_mDummy(self, m):
        self.mDummy = m

    def set_dDummy(self, d):
        self.dDummy = d

    def set_dDummyHeavy(self, d):
        self.dDummyHeavy = d

    def setCenter(self, c):
        self.center = c

    def setHeavy(self, h):
        self.heavy = h

    def set_nH(self, nH):
        self.nH = nH

    def setName(self, name):
        self.name = name

    def deriveName(self, ff):
        c = ff.getAtom(self.center)
        if not c:
            self.name = ''
            return

        # We give a name based on center and nH, with trailing underscore to distinguish from existing vsites.
        self.name = 'M{:s}H{:d}_'.format(c.getElement(), self.nH)

    # Does not consider nH or name.
    def isSame(self, v):
        return (self.center == v.center and \
                self.heavy == v.heavy and \
                self.mDummy == v.mDummy and \
                abs(self.dDummy-v.dDummy) < 0.00001 and \
                abs(self.dDummyHeavy-v.dDummyHeavy) < 0.00001)
        
        
# Atom with neighbours
class node(boolBase):
    def __init__(self, center='', neighbours=[]):
        self.center = center
        self.neighbours = neighbours

    def addCenter(self, center):
        if type(center) != atom:
            warn('Adding center that is not an atom to node')
            raise(TopologyError)
        else:
            self.center = center
            
    def addNeighbours(self, neighbours):
        if neighbours and type(neighbours) == list:
            if type(neighbours[0]):
                self.neighbours.extend(neighbours)
            else:
                warn('Adding empty list of neighbours to node')
                raise(TopologyError)
        else:
            self.neighbours.append(neighbours)
            
    # Returns empty constraint in case of any error
    def makeAngleConstraint(self, ff):
        """Returns constraint(i,j,d), where i and are the atoms
        and d is the constraint distance."""
        
        # We only do OH angle constraints for now.
        if ff.getElement(self.center) != 'O':
            return constraint()
        
        nH = 0
        H = -1
        nHeavy = 0
        Heavy = -1

        for i, n in enumerate(self.neighbours):
            if ff.getElement(n) == 'H':
                H = i
                nH += 1
            else:
                Heavy = i
                nHeavy += 1

        if nH != 1 or nHeavy != 1:
            return constraint()

        a = ff.getAngle(self.neighbours[Heavy].type, \
                self.center.type, \
                self.neighbours[H].type)

        b1 = ff.getBond(self.neighbours[Heavy].type, self.center.type)
        b2 = ff.getBond(self.neighbours[H].type,     self.center.type)

        x = b1 + b2 * -cos(a * pi/180)
        y = b2 * sin(a * pi/180)

        d = sqrt(x*x + y*y)

        return constraint([self.neighbours[Heavy].type, self.neighbours[H].type], d)

    
    # Returns empty vsite in case of any error
    def makeVsite(self, ff):

        v = vsite()
        
        if len(self.neighbours) not in [3, 4]:
            return v

        nH = 0
        H = []
        nHeavy = 0
        Heavy = -1

        for i,n in enumerate(self.neighbours):
            if ff.getElement(n) == 'H':
                H.append(i)
                nH += 1
            else:
                Heavy = i
                nHeavy += 1

        if nHeavy != 1 or nH != len(self.neighbours)-1:
            return v

        mTot = ff.getMass(self.center)
        I = 0.0
        # H-angle
        aH = ff.getAngle(self.neighbours[H[0]].type, self.center.type, self.neighbours[Heavy].type)
        dH = ff.getBond(self.neighbours[H[0]].type, self.center.type)
        # H-distance orthogonal to the anchor-bond.
        dHy = dH * sin(aH * pi/180)
        
        for i in H:
            m = ff.getMass(self.neighbours[i])
            I += m * dHy*dHy
            mTot += m

        v.set_mDummy(mTot/2.0)

        v.set_dDummy(2 * sqrt(I/(mTot)))

        # Distance between center and heavy
        dHeavyC = ff.getBond(self.center.type, self.neighbours[Heavy].type)
        # H-distance parallel to the anchor-bond
        dHx = dH * -cos(aH * pi/180)

        mH = ff.getMass(self.neighbours[H[0]])
        mC = ff.getMass(self.center)
        dX = (mC*dHeavyC + mH*nH*(dHeavyC+dHx)) / mTot

        d = sqrt(dX*dX + (v.dDummy/2)**2)
            
        v.set_dDummyHeavy(d)

        v.setCenter(self.center.type)
        v.setHeavy(self.neighbours[Heavy].type)
        v.set_nH(nH)
        v.deriveName(ff)
        
        return v

class topology:
    """Topology class"""
    
    def __init__(self):
        self.directories = ['./']
        try:
            self.gmxPath = environ['GMXDATA']
            if not self.gmxPath:
                raise(KeyError)
        except KeyError:
            warn('GMXDATA not set.')
            
        self.ffieldName = 'amber99sb-ildn.ff'

        self.ffield = None

        self.atoms  = []
        self.bonds  = []
        self.angles = []
        self.nodes  = []
        self.angleConstraints = []
        self.vsites = []

    def addDirectory(self, d):
        if d not in self.directories:
            self.directories.append(d)

    def FFread(self):
        for d in self.directories:
            try:
                self.ffield.read()
            except FFError:
                warn('No success reading forcefield files in '+d, bWarn=False)

    def setFF(self, ff):
        self.ffieldName = ff

    def finalise(self):
        self.directories.append(path.join(self.gmxPath, 'top', self.ffieldName))
        self.ffield = ffield.forceField()
        self.ffield.setDirectory(self.directories[-1])

    def itpRead(self, fname='topol.top'):
        warn('itpRead not implemented yet.')
        pass

    def rtpRead(self, resName, fileName='aminoacids.rtp'):

        residue = rtp.rtp()
        
        for d in self.directories:
            f = path.join(d, fileName)
            try:
                residue.readResidue(resName, fileName=f)
            except IOError:
                warn('Could not read '+f, bWarn=False)
                continue
                
            if residue.name:
                warn('Successfully read '+f, bWarn=False)
                break

        # Copy to self
        self.atoms = residue.atoms
        self.bonds = residue.bonds
        self.angles = residue.angles

        self.makeNodes()
        
    def atomExists(self, aname):
        for a in self.atoms:
            if aname == a.name:
                return True
        return False
    
    def addAtom(self, name, type):
        self.atoms.append(atom(name=name,type=type))

    def addBond(self, atoms):
        """Takes list of two atom names as argument"""
        
        # Do the atoms exist?
        self.bonds.append(bond(atoms=atoms))
        
    def addAngle(self, atoms):
        """Takes list of three atom names as argument"""
        self.angles.apend(angle(atoms=atoms))

    def getAtom(self, an):
        for a in self.atoms:
            if a.name == an:
                return a
        return None
        
    def getAtomType(self, an):
        for a in self.atoms:
            if a.name == an:
                return a.type
        return ''
    
    def makeNodes(self):
        for a in self.atoms:
            neighbours = []
            for b in self.bonds:

                a1 = self.getAtom(b.atoms[0].lstrip('+-'))
                a2 = self.getAtom(b.atoms[1].lstrip('+-'))
                
                if not (a1 and a2):
                    warn('Atom not found', bError=True)
                    warn('Bond {:s} -- {:s}'.format(b.atoms[0], b.atoms[1]), bWarn=False)
                    raise TopologyError
                
                if (a.type == a1.type):
                    neighbours.append(a2)
                elif (a.type == a2.type):
                    neighbours.append(a1)
                else:
                    continue

            self.nodes.append(node(center=a, neighbours=neighbours))
        
    def makeAngleConstraints(self):
        for a in self.nodes:
            c = a.makeAngleConstraint(self.ffield)
            if c:
                self.angleConstraints.append(c)
                
    def makeVsites(self):
        for a in self.nodes:
            v = a.makeVsite(self.ffield)
            if v:
                self.vsites.append(v)
        
    def dumpAtoms(self):
        output('--- ATOMS ---')
        for a in self.atoms:
            output('{:10s} : {:10s}'.format(a.name, a.type))

    def dumpBonds(self):
        output('--- BONDS ---')
        for b in self.bonds:
            output('{:10s} - {:10s}'.format(b.atoms[0], b.atoms[1]))

    def dumpAngles(self):
        output('--- ANGLES ---')
        for a in self.angles:
            output('{:10s} - {:10s} - {:10s}'.format(a.atoms[0], a.atoms[1], a.atoms[2]))

    def dumpAngleConstraints(self):
        output('--- ANGLE CONSTRAINTS ---')
        output('{:10s}   {:10s}   {:10s}'.format('i', 'j', 'b0'))
        for a in self.angleConstraints:
            output('{:10s} - {:10s}   {:<10f}'.format(a.atoms[0], a.atoms[1], a.b0))

    def dumpVsites(self):
        output('--- VSITES ---')
        output('{:10s}   {:10s}   {:10s}   {:10s}   {:10s}   {:10s}'.format(\
            'Center', 'Heavy', 'mDummy', 'dDD', 'dDHeavy', 'name'))
        for v in self.vsites:
            output('{:10s}   {:10s}   {:<10f}   {:<10f}   {:<10f}   {:10s}'.format(\
            v.center, v.heavy, v.mDummy, v.dDummy, v.dDummyHeavy, v.name))
def test():
    top = topology()

    top.setFF('charmm36.ff')
    top.finalise()
    top.rtpRead('THR', fileName='merged.rtp')
    top.FFread()

    top.makeAngleConstraints()
    top.makeVsites()
    
    top.dumpAtoms()
    top.dumpBonds()
    top.dumpAngles()
    top.dumpAngleConstraints()
    top.dumpVsites()

if __name__ == '__main__':
    test()
