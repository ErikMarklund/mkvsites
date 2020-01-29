from basics import *
from toputil import *

class Itp(topBase):
    """ITP (molecule type) reader"""

    def __init__(self, *args, **kwargs):
       super(Itp, self).__init__(*args, **kwargs)

    def read(self, fileName):
        """Read itp file"""
        preprocessor = Cpp()
        f = preprocessor.parse(fileName)
        if not f:
            raise IOError
        
        lineNo = 0
        readType = ''
        try:
            for line in f:
                lineNo += 1

                if not line.strip() or (line[0] in ';#') or (line.split(';')[0].strip()[0] in ';#'):
                    # Comments, empty lines, ...
                    continue

                sline = line.split(';')[0].split()

                directive = readDirective(line)

                if directive:
                    readType = directive
                    continue

                if readType == 'moleculetype':
                    self.name = sline[0]

                elif readType == 'atoms':
                    if len(sline) < 8:
                        warn('Ill-formated atom line in '+fileName+':', bError=True)
                        warn(line, bWarn=False)
                        raise(ItpError)

                    atype = sline[1]
                    aname = sline[4]
                    cgnr  = int(sline[5])
                    q     = float(sline[6])

                    if not (aname and atype and cgnr and q!=None):
                        warn('Could not derive name/type/cgnr/q of atom:', bError=True)
                        warn(line, bWarn=False)
                        warn('aname: '+aname, bWarn=False)
                        warn('cgnr: {:d}'.format(cgnr), bWarn=False)
                        warn('q: {:f}'.format(q), bWarn=False)
                        raise(ItpError)

                    self.atoms.append(atom(name=aname, type=atype, q=q, cgnr=cgnr))

                elif readType == 'bonds':
                    if len(sline) < 3:
                        warn('Ill-formated bond line in '+fileName+':', bError=True)
                        warn(line, bWarn=False)
                        raise(ItpError)

                    atomIndices = [int(i)-1 for i in sline[0:2]]

                    atoms = [self.atoms[i].name for i in atomIndices]
                    ftype = int(sline[2])
                    params = sline[3:]
                    self.bonds.append(bond(atoms=atoms, ftype=ftype, params=params))

                elif readType == 'angles':
                    if len(sline) < 4:
                        warn('Ill-formated angle line in '+fileName+':', bError=True)
                        warn(line, bWarn=False)
                        raise(ItpError)

                    atomIndices = [int(i)-1 for i in sline[0:3]]

                    atoms = [self.atoms[i].name for i in atomIndices]
                    ftype = int(sline[3])
                    params = sline[4:]
                    self.angles.append(angle(atoms=atoms, ftype=ftype, params=params))

                elif readType == 'dihedrals':
                    if len(sline) < 5:
                        warn('Ill-formated dihedral line in '+fileName+':', bError=True)
                        warn(line, bWarn=False)
                        raise(ItpError)

                    atomIndices = [int(i)-1 for i in sline[0:4]]

                    atoms = [self.atoms[i].name for i in atomIndices]
                    ftype = int(sline[4])
                    params = sline[5:]
                    if ftype in (2,4):
                        self.impropers.append(dihedral(atoms=atoms, ftype=ftype, params=params))
                    else:
                        self.dihedrals.append(dihedral(atoms=atoms, ftype=ftype, params=params))

                else:
                    # Irrelevant directive
                    pass


        except IOError:
            warn('Error while reading from file '+f.__name__, bError=True)
            raise
        except IndexError:
            warn('Choked on this line ({:d}):'.format(lineNo), bError=True)
            warn(line, bWarn=False)
