from basics import *
from toputil import *

class Rtp(topBase):
    """RTP reader"""

    def __init__(self, *args, **kwargs):
        super(Rtp, self).__init__(*args, **kwargs)
        
    def readResidue(self, resName, fileName='aminoacids.rtp'):
        """Read residue from rtp file and store it"""
        preprocessor = Cpp()
        f = preprocessor.parse(fileName)
        if not f:
            raise IOError
        
        bMatchRes = False
        readType = 'res'
        lineNo = 0
        try:
            for line in f:
                lineNo += 1

                if len(line.strip()) == 0 or line.strip()[0] in ';#':
                # Comments & empty lines
                    continue

                sline = (line.split(';')[0]).split()

                # Try to read as directive
                directive = readDirective(line)
                if directive:
                    # Find start of residue
                    if readType == 'res':
                        if directive == resName:
                            readType = directive

                        continue

                    readType = directive
                    continue
                elif readType == 'res':
                    # Only interested in scanning directives when we haven't found our residue
                    continue                            

                # Find [ atoms ]
                elif readType == 'atomhead':
                    if directive:
                        if directive == 'atoms':
                            readType = directive
                            continue
                        else:
                            warn('Unexpected directive', bError=True)
                            raise(RtpError)
                        readType = 'atom'
                        continue

                # Read atoms
                elif readType == 'atoms':                        
                    if len(sline) != 4:
                        if directive:
                            if directive == 'bonds':
                                readType = directive
                                continue
                            else:
                                warn('Unexpected directive', bError=True)
                                raise(RtpError)

                        warn('Ill-formatted atom on line {:d}:'.format(lineNo))
                        warn(line, bWarn=False)
                        continue

                    aname = sline[0]
                    atype = sline[1]
                    self.atoms.append(atom(name=aname, type=atype))

                # Read bonds
                elif readType == 'bonds':
                    if len(sline) != 2:
                        if directive:
                            if directive == 'angles':
                                readType = directive
                                continue
                            else:
                                warn('Unexpected directive:', bError=True)
                                warn(line, bWarn=False)
                                raise(RtpError)

                        warn('Ill-formatted bond on line {:d}:'.format(lineNo))
                        warn(line, bWarn=False)

                    anames = sline[:]
                    self.bonds.append(bond(atoms=anames))

                elif readType == 'angles':
                    if len(sline) != 3:
                        if directive:
                            if directive:
                                readType = 'Done' # or whatever
                                continue
                            else:
                                warn('Unexpected directive', bError=True)
                                raise(RtpError)

                        warn('Ill-formatted angle on line {:d}:'.format(lineNo))
                        warn(line, bWarn=False)

                    anames = sline[:]
                    self.angles.append(angle(atoms=anames))

                else:
                    # We appear to be done.
                    break

        except IOError:
            warn('Error while reading from file '+f.__name__, bError=True)
            raise(IOError)

        
        if len(self.atoms) < 3:
            warn('{:d} atoms. No extra params needed'.format(len(self.atoms)), bWarn=False)

        elif self.bonds != [] and self.angles != []:
            warn('{:d} bonds and {:d} angles. Is that correct?'.format(len(self.bonds), len(self.angles)))

        else:
            # evrything ok
            pass

        self.name = resName
