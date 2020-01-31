from sys import exit, argv, stdout
from os import environ, path
from basics import *
from toputil import *
from math import *
import ffield
import rtp
import itp
# Topology class and functions for reading rtp files and itp/top files.

# Class for storing the angle constraints.
class constraint(boolBase):
    """Class fopr storing the angle constraints"""
    def __init__(self, atoms=[], b0=0.0):
        self.atoms = atoms
        self.b0 = b0

    def setConstraint(self, atoms, b0):
        """Set atom indices and distance"""
        self.atoms = atoms
        self.b0 = b0

    # isSame() uses a ~0.01% tolerance because floating point arithmetic.
    def isSame(self, c):
        """Returns True if and only if indices and distances both match (within 0.01%)"""
        return (
                (self.i == c.i and self.j == c.j) or \
                (self.i == c.j and self.j == c.i) \
            ) and abs(self.b0-c.b0)*2/(self.b0+c.b0) < 0.0001

# Class for storing vsites.
class vsite(boolBase):
    """Class for storing and making vsites with dummy atoms"""
    def __init__(self, center='', heavy='', cName='', hName='', nH=0, mDummy=0.0, dDummy=0.0, dDummyHeavy=0.0, name=''):
        self.center = center
        self.heavy  = heavy
        self.centerName = cName
        self.heavyName = hName
        self.nH     = nH
        # Mass of a single dummy mass
        self.mDummy = mDummy
        # Distance between dummy masses
        self.dDummy = dDummy
        # Distance between dummy and heavy anchor atom
        self.dDummyHeavy = dDummyHeavy
        self.name   = name

    def set_mDummy(self, m):
        """Set the mass of dummy atoms"""
        self.mDummy = m

    def set_dDummy(self, d):
        """Set distance between dummy atoms"""
        self.dDummy = d

    def set_dDummyHeavy(self, d):
        """Set distance between dummies and heavy"""
        self.dDummyHeavy = d

    def setCenter(self, c):
        """Set center atom"""
        self.center = c

    def setHeavy(self, h):
        """Set heavy atom"""
        self.heavy = h

    def set_nH(self, nH):
        """Set number of H"""
        self.nH = nH

    def setName(self, name):
        """Set name of vsite"""
        self.name = name

    def deriveName(self, ff):
        """Derive name from center atom and number of H"""
        c = ff.getAtom(self.center)
        if not c:
            self.name = ''
            return

        # We give a name based on center and nH, with trailing underscore to distinguish from existing vsites.
        self.name = 'M{:s}H{:d}_'.format(c.getElement(), self.nH)

    # Does not consider nH or name.
    def isSame(self, v):
        """Are two vsites the same? Based on atoms and distances."""
        return (self.center == v.center and \
                self.heavy == v.heavy and \
                self.mDummy == v.mDummy and \
                abs(self.dDummy-v.dDummy) < 0.00001 and \
                abs(self.dDummyHeavy-v.dDummyHeavy) < 0.00001)

    def dump(self, ostream=stdout):
        """Prints the vsite to stream"""
        output('+ + + + {:8s} vsite + + + +'.format(self.name))
        output('center       {:10s}'.format(self.center), ostream=ostream)
        output('centerName   {:10s}'.format(self.centerName), ostream=ostream)
        output('heavy        {:10s}'.format(self.heavy), ostream=ostream)
        output('heavyName    {:10s}'.format(self.heavyName), ostream=ostream)
        output('nH           {:<10d}'.format(self.nH), ostream=ostream)
        output('mDummy       {:<10.6}'.format(self.mDummy), ostream=ostream)
        output('dDummy       {:<10.6}'.format(self.dDummy), ostream=ostream)
        output('dDummyHeavy  {:<10.6}'.format(self.dDummyHeavy), ostream=ostream)

# Atom with neighbours
class node(boolBase):
    """Class for handling nodes. A node is a heavy atom and its immediate neighbours"""
    def __init__(self, center='', neighbours=[]):
        self.center = center
        self.neighbours = neighbours

    def addCenter(self, center):
        """Add center atom to node"""
        if type(center) != atom:
            warn('Adding center that is not an atom to node')
            raise(TopologyError)
        else:
            self.center = center
            
    def addNeighbours(self, neighbours):
        """Add neighbours to node"""
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
        """Creates and returns constraint(i,j,d), where i and are the atoms
        and d is the constraint distance.
        """
        
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

    
    # Returns empty vsite in case of any error.
    def makeVsite(self, ff):
        """Makes and returns a vsite if the node meets the criteria.
        Returns empty vsite if node fails to meet any criterion
        """
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
        self.ffield = ffield.forceField()

        self.atoms  = []
        self.bonds  = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.nodes  = []
        self.angleConstraints = []
        self.vsites = []


    def FFread(self):
        """Reads the forcefield"""
        try:
            self.ffield.read()
        except FFError:
            fdir = self.ffield.getPath()
            if not fdir:
                fdir = 'current directory'
            warn('No success reading forcefield files in '+fdir, bWarn=False)

    def setFF(self, ff):
        """Sets the name of the forcefield"""
        self.ffield.setPath(ff)

    def getFF(self):
        """Returns the name of the forcefield"""
        return self.ffield.getPath()

    def mol2top(self, r):
        """Stores the atoms, bonds, and angles from a molecule/residue."""
        self.atoms = r.atoms
        self.bonds = r.bonds
        self.angles = r.angles
        if hasattr(r, 'dihedrals'):
            self.dihedrals = r.dihedrals
        if hasattr(r, 'impropers'):
            self.impropers = r.impropers

        self.makeNodes()
        
    def itpRead(self, fileName='topol.top', bVerbose=False):
        """Reads itp for a molecule and stores the atoms, bonds, and angles."""
        mol = itp.Itp()

        f = fileName

        try:
            mol.read(f)
        except IOError:
            warn('Could not read '+f, bWarn=False, bPrint=bVerbose)

        if mol.name:
            output('Successfully read '+f, bPrint=bVerbose)
        else:
            raise TopologyError

        self.mol2top(mol)

            
    def rtpRead(self, resName, fileName='aminoacids.rtp', bVerbose=False):
        """Reads specific entry from rtp file and stores the atoms, bonds, and angles."""
        residue = rtp.Rtp()
        
        f = fileName
        try:
            residue.readResidue(resName, fileName=f)
        except IOError:
            warn('Could not read '+f, bWarn=False, bPrint=bVerbose)

        if residue.name:
            output('Successfully read '+f, bPrint=bVerbose)
        else:
            raise TopologyError

        self.mol2top(residue)

    def rtpWriteFile(self, fileName, resName='RES', bVerbose=False):
        """Writes residue to a new rtp file."""
        with open(fileName,'w') as f:
            self.rtpWrite(resName=resName, ostream=f, bVerbose=bVerbose)

    def rtpWrite(self, resName='RES', ostream=stdout, bVerbose=False):
        """Writes residue to a new rtp entry."""
        f = ostream
        
        try:
            f.write('; Generated by mkvsites.py.\n; Please cite:\n')
            c = cite_mkvsites()
            for line in c.split('\n'):
                f.write(';   '+line+'\n')
            f.write('\n')

            f.write('[ bondedtypes ]\n')
            f.write('; bonds  angles  dihedrals  impropers\n')
            bt = 0
            at = 0
            dt = 0
            it = 0
            if self.bonds:
                bt = self.bonds[0].ftype
            if self.angles:
                at = self.angles[0].ftype
            if self.dihedrals:
                dt = self.dihedrals[0].ftype
            if self.impropers:
                it = self.impropers[0].ftype
            f.write("{:6d}{:6d}{:6d}{:6d} ; Derived from first parameters found. Check!\n\n".format(bt,at,dt,it))

            f.write('[ {:s} ]\n'.format(resName))
            f.write(' [ atoms ]\n')

            for a in self.atoms:
                f.write('{:>6s}  {:>6s}         {:9.5f}{:6d}\n'.format(a.name, a.type, a.q, a.cgnr))
            if len(self.bonds) != 0:
                f.write('\n [ bonds ]\n')
                f.write(';atom1 atom2\tparams\n')
                for b in self.bonds:
                    b1 = b.atoms[0]
                    b2 = b.atoms[1]
                    f.write('{:>6s}{:>6s}'.format(b1,b2))
                    for p in b.params:
                        f.write('\t{:s}'.format(p))
                    f.write('\n')

            if len(self.angles) != 0:
                f.write('\n [ angles ]\n')
                f.write(';atom1 atom2 atom3\tparams\n')
                for a in self.angles:
                    at = [a.atoms[i] for i in range(3)]
                    f.write('{:>6s}{:>6s}{:>6s}'.format(at[0],at[1],at[2]))
                    for p in a.params:
                        f.write('\t{:s}'.format(p))
                    f.write('\n')

            if len(self.dihedrals) != 0:
                f.write('\n [ dihedrals ]\n')
                f.write(';atom1 atom2 atom3 atom4\tparams\n')
                for d in self.dihedrals:
                    at = [d.atoms[i] for i in range(4)]
                    f.write('{:>6s}{:>6s}{:>6s}{:>6s}'.format(at[0],at[1],at[2],at[3]))
                    for p in d.params[1:]:
                        f.write('\t{:s}'.format(p))
                    f.write('\n')

            if len(self.impropers) != 0:
                f.write('\n [ impropers ]\n')
                f.write(';atom1 atom2 atom3 atom4\tparams\n')
                for d in self.impropers:
                    at = [self.atoms[d.atoms[i]].name for i in range(4)]
                    f.write('{:>6s}{:>6s}{:>6s}{:>6s}'.format(at[0],at[1],at[2],at[3]))
                    for p in d.params[1:]:
                        f.write('\t{:10.5f}'.format(p))
                    f.write('\n')

        except IOError:
            warn('Could not write rtp file.', bWarn=False, bPrint=bVerbose)
            raise

    def atomExists(self, aname):
        """True if atom exists in self.atoms. False if not."""
        for a in self.atoms:
            if aname == a.name:
                return True
        return False
    
    def addAtom(self, name, type):
        """Adds the atom to self.atoms:"""
        self.atoms.append(atom(name=name,type=type))

    def addBond(self, atoms):
        """Takes list of two atom names as argument"""
        # Do the atoms exist?
        self.bonds.append(bond(atoms=atoms))
        
    def addAngle(self, atoms):
        """Takes list of three atom names as argument"""
        self.angles.apend(angle(atoms=atoms))

    def getAtom(self, an):
        """Returns the atom if it exists, None if not"""
        for a in self.atoms:
            if a.name == an:
                return a
        return None
        
    def getAtomType(self, an):
        """Gives the atom type for the atom"""
        for a in self.atoms:
            if a.name == an:
                return a.type
        return ''
    
    def makeNodes(self):
        """Makes nodes based on the topology data"""
        for a in self.atoms:
            neighbours = []
            for b in self.bonds:
                # Takes atom names, returns atom types. To test if atoms exist in atoms list
                # Strip '+-' to handle rtp inter-residue bonds.
                a1 = self.getAtom(b.atoms[0].lstrip('+-'))
                a2 = self.getAtom(b.atoms[1].lstrip('+-'))

                if not (a1 and a2):
                    warn('Atom not found', bError=True)
                    warn('Bond {:s} -- {:s}'.format(b.atoms[0], b.atoms[1]), bWarn=False)
                    raise TopologyError

                if (a.name == a1.name):
                    neighbours.append(a2)
                elif (a.name == a2.name):
                    neighbours.append(a1)
                else:
                    continue

            self.nodes.append(node(center=a, neighbours=neighbours))

    def spotAlienAtoms(self):
        """Checks if all atoms are described by the forcefield"""
        aliens = []
        for n in self.nodes:
            A = [n.center] + n.neighbours
            for a in A:
                if not self.ffield.getAtom(a.type):
                    aliens.append(a)

        return aliens

    def makeAngleConstraints(self):
        """Derives angle constraints"""

        for a in self.nodes:
            c = a.makeAngleConstraint(self.ffield)
            if c:
                self.angleConstraints.append(c)
                
    def makeVsites(self):
        """Creates vsites for all atoms where applicable"""
        for a in self.nodes:
            v = a.makeVsite(self.ffield)
            if v:
                self.vsites.append(v)
        
    def dumpAtoms(self):
        """Print all atoms"""
        output('--- ATOMS ---')
        for a in self.atoms:
            output('{:10s} : {:10s}'.format(a.name, a.type))

    def dumpBonds(self):
        """Print all bonds"""
        output('--- BONDS ---')
        for b in self.bonds:
            output('{:10s} - {:10s}'.format(b.atoms[0], b.atoms[1]))

    def dumpAngles(self):
        """Print all angles"""
        output('--- ANGLES ---')
        for a in self.angles:
            output('{:10s} - {:10s} - {:10s}'.format(a.atoms[0], a.atoms[1], a.atoms[2]))

    def dumpAngleConstraints(self):
        """Print all angle constraints"""
        output('--- ANGLE CONSTRAINTS ---')
        output('{:10s}  {:10s}  {:10s}'.format('i', 'j', 'b0'))
        for a in self.angleConstraints:
            output('{:10s}  {:10s}  {:<10f}'.format(a.atoms[0], a.atoms[1], a.b0))

    def dumpVsites(self):
        """Print all vsites"""
        output('--- VSITES ---')
        output('{:10s}   {:10s}   {:10s}   {:10s}   {:10s}   {:10s}'.format(\
            'Center', 'Heavy', 'mDummy', 'dDD', 'dDHeavy', 'name'))
        for v in self.vsites:
            output('{:10s}   {:10s}   {:<10f}   {:<10f}   {:<10f}   {:10s}'.format(\
            v.center, v.heavy, v.mDummy, v.dDummy, v.dDummyHeavy, v.name))

    def addNewVsiteType(self, v):
        """Add a new vsite type if it does not already exist"""
        # Does it already exist?
        for u in self.uniqueVsites:
            if u.DummyConstraint == v.DummyConstraint \
              and u.Anchor == v.Anchor \
              and u.AnchorConstraint == v.AnchorConstraint:
                return

        # Ok, so a new one.
        if v.DummyName[-1] == '_':
            Mnames = [u.DummyName for u in self.uniqueVsites]
            i = 1
            while v.DummyName+str(i) in Mnames:
                i += 1

            v.DummyName = v.DummyName+str(i)

        output('Adding vsite type {:8s}'.format(v.DummyName))
        self.uniqueVsites.append(v)

    def getVsiteType(self, v):
        """ """
        for uv in self.uniqueVsites:
            # Allow 1% discrepancy (Order of 1/100 A)
            if (abs(v.dDummy-uv.DummyConstraint)/uv.DummyConstraint > 0.01 or \
                    abs(v.dDummyHeavy-uv.AnchorConstraint)/uv.AnchorConstraint > 0.01 or \
                v.center != uv.Center or \
                v.heavy != uv.Anchor):
                continue

            return uv

        return None

    def identifyVsites(self):
        """Identify vsite types and rename vsites in topology"""
        self.uniqueVsites = []

        for v in self.vsites:
            vt = self.getVsiteType(v)
            if not vt:
                vt = ffield.ffVsiteType(DummyName=v.name,
                                        DummyConstraint=v.dDummy,
                                        Anchor=v.heavy,
                                        AnchorConstraint=v.dDummyHeavy,
                                        Center=v.center)

            self.addNewVsiteType(vt)

        for v in self.vsites:
            vt = self.getVsiteType(v)
            if not vt:
                warn('Encountered unknown vsite type', bError=True, ostream=stderr)
                v.dump()
                output('{:8d} Known vsite types:'.format(len(self.uniqueVsites)), ostream=stderr)
                for uv in self.uniqueVsites:
                    uv.dump()
                raise TopologyError

            v.setName(vt.getDummyName())

    def dumpVsiteTypes(self):
        """Print all vsite types"""
        output('### Vsite data to add to forcefield')
        output('# Constraints (ffbonded.itp)')
        for u in self.uniqueVsites:
            output('{:8s}{:8s} 2    {:8.6f}'.format(u.DummyName, u.DummyName, u.DummyConstraint))
            output('{:8s}{:8s} 2    {:8.6f}'.format(u.DummyName, u.Anchor, u.AnchorConstraint))

        output('# Dummy masses (ffnonbonded.itp)')
        for u in self.uniqueVsites:
            output('{:8s}0   0.000000    0.00    A   0.0 0.0'.format(u.DummyName))

        output('# Dummy masses (atomtypes.itp)')
        for u in self.uniqueVsites:
            output('{:8s}  0.00000'.format(u.DummyName))

        output('# Add to the .vsd file under the appropriate group')
        for v in self.vsites:
            output('{:8s}{:8s}{:8s}'.format(v.center, v.heavy, v.name))


def testRtp():
    """Function for testing the Rtp class"""
    top = topology()
    ffpath='/Users/erikmarklund/src/mkvsites/testdata/BDDM_FF_CHARMM36/charmm36-jun2015.ff'
    top.setFF(ffpath)
    rtppath=path.join(ffpath,'merged.rtp')
    top.rtpRead('THR', fileName=rtppath)

    top.FFread()

    top.makeAngleConstraints()
    top.makeVsites()
    
    top.dumpAtoms()
    top.dumpBonds()
    top.dumpAngles()
    top.dumpAngleConstraints()
    top.dumpVsites()

def testItp():
    """Function for testing the Itp class"""
    top = topology()
    ffpath='/Users/erikmarklund/src/mkvsites/testdata/BDDM_FF_CHARMM36/charmm36-jun2015.ff'
    top.setFF(ffpath)
    itppath='/Users/erikmarklund/src/mkvsites/testdata/BDDM_FF_CHARMM36/1-bDM_CHARMM36.itp'
    output(itppath)
    top.itpRead(fileName=itppath)

    top.FFread()

    top.makeAngleConstraints()
    top.makeVsites()
    
    top.dumpAtoms()
    top.dumpBonds()
    top.dumpAngles()
    top.dumpAngleConstraints()
    top.dumpVsites()

    output('### ffVsites ###')
    top.ffield.gatherVsites()
    for i, v in enumerate(top.ffield.vsites):
        output('Vsite {:d}:'.format(i))
        v.dump()

def testBoth():
    """Wrapper for testing both Rtp and Itp classes"""
    output('##### Testing rtp reader #####')
    testRtp()
    output('\n##### Testing itp reader #####')
    testItp()

    
if __name__ == '__main__':
    testBoth()
