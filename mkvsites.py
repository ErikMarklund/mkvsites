from math import *
from sys import exit, argv
import topology as top
from os import path
from basics import *

def dump_help():
    output('  USAGE: python mkvsites.py <options> rtp/itp-file\n')
    output('  mkvsites.py reads an rtp entry or an itp file, and derives')
    output('  the angle constraints and vsite parameters needed for')
    output('  simulating with a longer time step.\n')
    output('  OPTIONS:')
    output('  -h        : Display this help text and quit')
    output('  -ff=FF    : Employ force field FF')
    output('  -res=RES  : With rtp files, make vsites for residue RES')


# Read topology/rtp and make vsite parameters and angle-constraints for OH groups.
#
# Reads topology/rtp, finds geometries in need of dummy masses or angle constraints,
# and generates the necessary parameters using bond lengths and angles from forcefield.

if __name__ == '__main__':
    t = top.topology()

    ff = ''
    res = ''

    for a in argv[1:]:
        if a[0] == '-':
            # Flag/option

            sarg = a.split('=')
            flag = sarg[0]

            if flag == '-h':
                dump_help()
                exit(0)

            elif flag == '-ff':
                if ff:
                    warn('Forcefield can only be chosen once.', bError=True)

                try:
                    ff = sarg[1]
                except IndexError:
                    warn('Expected -ff=forcefield.ff, got {:a}'.format(a), bError=True)

            elif flag == '-res':
                res = sarg[1]

        else:
            # Topology file to analyse
            topfile = a


    if ff:
        t.setFF(ff)
    else:
        output('No force field specified. Will use {:s}.'.format(t.getFF()))

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
