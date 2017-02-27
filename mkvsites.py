from math import *
from sys import exit
import topology as top
from os import path
from basics import *
import argparse



# Read topology/rtp and make vsite parameters and angle-constraints for OH groups.
#
# Reads topology/rtp, finds geometries in need of dummy masses or angle constraints,
# and generates the necessary parameters using bond lengths and angles from forcefield.

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate stuff")
    parser.add_argument('-ff', '--force-field', action='store', nargs=1, type=str, required=True, metavar='FF', dest='ff', help='Path to force field.')
    parser.add_argument('-res', '--residue', action='store', nargs=1, type=str, metavar='RES', dest='res', help='Residue name, rtp files only.')
    parser.add_argument('file', action='store', nargs=1, type=str, help='itp or rtp file')

    args = parser.parse_args()

    t = top.topology()

    ff = args.ff[0]
    res = args.res[0]
    topfile = args.file[0]

    t.setFF(ff)

    t.finalise()

    
    if topfile:
        topdir, topbase = path.split(topfile)
        t.addDirectory(topdir)
        ext = path.splitext(topbase)[-1]

        if ext == '.itp':
            output('Input is itp file')
            t.itpRead(fileName=topbase)

        elif ext == '.rtp':
            output('Input is rtp file')
            t.rtpRead(res, fileName=topbase)

            if not res:
                warn('Need to provide residue name with rtp files (-res=resname).', bError=True)
                exit(1)

        else:
            warn('Unsupported file type: {:s}'.format(topbase), bError=True)
            exit(1)

        output('Will find vsites and angle constraints for {:s} dressed in the forcefield {:s}'.format(topfile, ff))
        t.FFread()
        t.makeAngleConstraints()
        t.makeVsites()

        t.dumpAngleConstraints()
        t.dumpVsites()

        t.identifyVsites()
        t.dumpVsiteTypes()
