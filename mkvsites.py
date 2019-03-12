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
    __version__ = '0.1'

    parser = argparse.ArgumentParser(description="Calculates missing constraints for virtual-site constructions.")
    parser.add_argument('-ff', '--force-field', action='store', nargs=1, type=str, required=True, metavar='FF', dest='ff', help='Path to force field.')
    parser.add_argument('-res', '--residue', action='store', nargs=1, type=str, metavar='RES', dest='res', help='Residue name, rtp files only.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Be loud and noisy')
    parser.add_argument('-ver', '--version', action='store_true', help='Print version')
    parser.add_argument('file', action='store', nargs=1, type=str, help='itp or rtp file')

    args = parser.parse_args()

    bVersion = vars(args)['version']
    bVerbose = vars(args)['verbose']

    if bVersion:
        output('Version {:s}'.format(__version__))

    t = top.topology()

    ff = args.ff[0]
    if args.res:
        res = args.res[0]
    else:
        res = None
    topfile = args.file[0]

    t.setFF(ff)

    t.finalise()

    
    if topfile:
        topdir, topbase = path.split(topfile)
        t.addDirectory(topdir)
        ext = path.splitext(topbase)[-1]

        if ext in [ '.itp', '.top' ] :
            output('Input is {:s} file'.format(ext.lstrip('.')), bPrint=bVerbose)
            t.itpRead(fileName=topbase, bVerbose=bVerbose)

        elif ext == '.rtp':
            output('Input is rtp file', bPrint=bVerbose)
            if not res:
                parser.error('Need to provide residue name with rtp files (-res).')
            t.rtpRead(res, fileName=topbase, bVerbose=bVerbose)

        else:
            parser.error('Unsupported file type: {:s}'.format(topbase))

        output('Will find vsites and angle constraints for {:s} dressed in the forcefield {:s}'.format(topfile, ff))
        t.FFread()
        t.makeAngleConstraints()
        t.makeVsites()

        if bVerbose:
            t.dumpAngleConstraints()
            t.dumpVsites()

        t.identifyVsites()
        t.dumpVsiteTypes()
