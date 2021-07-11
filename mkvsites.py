#!/usr/bin/env python

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


from math import *
from sys import exit
import topology as top
from toputil import *
from os import path
from basics import *
import argparse



# Read topology/rtp and make vsite parameters and angle-constraints for OH groups.
#
# Reads topology/rtp, finds geometries in need of dummy masses or angle constraints,
# and generates the necessary parameters using bond lengths and angles from forcefield.

if __name__ == '__main__':
    __version__ = '1.1'

    parser = argparse.ArgumentParser(description="Calculates missing constraints for virtual-site constructions.")
    parser.add_argument('-ff', '--force-field', action='store', nargs=1, type=str, required=True, metavar='FF', dest='ff', help='Path to force field.')
    parser.add_argument('-res', '--residue', action='store', nargs=1, type=str, metavar='RES', dest='res', help='Residue name, rtp files only.')
    parser.add_argument('-ro', '--rtp-out', action='store', nargs=1, type=str, metavar='mymol.rtp', dest='ro', help='RTP output (for itp/top input). Enables vsite and dummy construction using pdb2gmx.')
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

    if args.ro:
        rofile = args.ro[0]
    else:
        rofile = None

    t.setFF(ff)
    
    if topfile:
        #topdir, topbase = path.split(topfile)
        #t.addDirectory(topdir)
        ext = path.splitext(topfile)[-1]

        if ext in [ '.itp', '.top' ] :
            output('Input is {:s} file'.format(ext.lstrip('.')), bPrint=bVerbose)
            try:
                t.itpRead(fileName=topfile, bVerbose=bVerbose)
            except TopologyError:
                warn('Incomplete reading of {:s}. Faulty topology?'.format(topfile), bError=True)
                raise TopologyError

        elif ext == '.rtp':
            output('Input is rtp file', bPrint=bVerbose)
            if not res:
                parser.error('Need to provide residue name with rtp files (-res).')
            if rofile:
                parser.error('RTP output (-ro) only makes sense for itp/top input.')
            try:
                t.rtpRead(res, fileName=topfile, bVerbose=bVerbose)
            except TopologyError:
                warn('Incomplete reading of {:s}. Faulty rtp-file?'.format(res), bError=True)

        else:
            parser.error('Unsupported file type: {:s}'.format(topbase))

        output('Will find vsites and angle constraints for {:s} dressed in the forcefield {:s}'.format(topfile, ff))
        t.FFread()
        aliens = t.spotAlienAtoms()
        if aliens:
            output("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =", ostream=stderr)
            warn("The following atoms are not found in the force field,\n" \
                     "and will thus not be considered for vsites or angle constraints.\n" \
                     "Perhaps derive_params.py can help.\n")
            output("{:8s}{:8s}".format("Name", "Type"), ostream=stderr)
            output("------------", ostream=stderr)
            for a in aliens:
                a.dump(ostream=stderr)
            output("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n", ostream=stderr)

        t.makeAngleConstraints()
        t.makeVsites()
        t.identifyVsites()
        output('***********************')
        t.ffield.gatherVsites()
        output('*** Dummy-atom constructs in force field')
        for i, v in enumerate(t.ffield.vsites):
            output('Vsite {:d}:'.format(i))
            v.dump()
        output('***********************')

        if bVerbose:
            #t.dumpAngleConstraints()
            t.dumpVsites()

        t.dumpVsiteTypes()


        if rofile:
            # Dump RTP file too
            t.rtpWriteFile(rofile, bVerbose=bVerbose)
