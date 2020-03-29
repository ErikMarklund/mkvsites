#!/usr/bin/env python

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


from sys import exit
import pdb
import argparse
import dihedral_convert as dconv
from sys import stdout, stderr
from basics import *
from os import path

class interaction:
    def __init__(self):
        self.a = []
        self.func = '0'
        self.r = '0.0'
        self.k = '0.0'

    def issame(self, I):
        if len(self.a) != len(I.a):
            return False

        if type(self.func) != type(I.func) or type(self.r) != type(I.r) or type(self.k) != type(I.k):
            return False

        if self.a != I.a or self.func != I.func or self.r != self.r or self.k != self.k:
            return False

        return True

    def set_a(self, a):
        for aa in a:
            self.a.append(aa)

    def set_func(self, f):
        self.func = f

    def set_r(self, r):
        self.r = r

    def set_k(self, k):
        if type(k) == type([]):
            self.k = []
            for kk in k:
                self.k.append(kk)
        else:
            self.k = k

    def get_atoms(self):
        return self.a

    def get_func(self):
        return self.func

    def get_r(self):
        return self.r

    def get_k(self):
        return self.k

    def get_type(self):
        l = len(self.a)

        if l==2:
            return 'bond'
        elif l==3:
            return 'angle'
        elif l==4:
            return 'dihedral'
        else:
            raise ValueError

    def dump(self, bSkip=False, macro=None):
        if bSkip:
            s = ";   "
        else:
            s = "    "
        for a in self.a:
            s = s + "{:4s}".format(a)

        s = s + "{:4s}".format(self.func)

        if self.get_type() == 'dihedral':
            if macro:
                s = s + "    {:s}".format(macro.name)
            else:
                for p in self.k:
                    s = s + "    {:s}".format(p)
        else:
            s = s + "    {:s}".format(self.r)
            s = s + "    {:s}".format(self.k)

        output(s)

    def rtp_dump_interaction(self, ostream=stdout, bDumpParams=False, macro=None):
        s = "  "
        for a in self.a:
            s = s + "{:6s}".format(a)

        #Function type given at the start of the file for all bonded interaction classes.
        #s = s + "{:4s}".format(self.func)

        if bDumpParams:
            if macro:
                s = s + "    {:s}".format(macro.get_name())

            else:
                if self.get_type() == 'dihedral':
                    if int(self.get_func()) in (1,9):
                        s = s + "    {: 8.6e}".format(float(self.k[0]))
                        s = s + "    {: 8.6e}".format(float(self.k[1]))
                        s = s + "    {:3d}".format(int(self.k[2]))
                    else:
                        for p in self.k:
                            s = s + "    {: 8.6e}".format(float(p))
                else:
                    s = s + "    {:s}".format(self.r)
                    for p in self.k.split():
                        s = s + "    {: 8.6e}".format(float(p))

        s = s + "\n"
        ostream.write(s)

    def make_macro(self, name='MACRO'):
        r = self.get_r()
        k = self.get_k()
        if type(k) == type([]):
            macrodef = ' '.join([r]+k)
        else:
            macrodef = ' '.join([r, k])

        return parameter_macro(name=name, params=macrodef)


class parameter_macro():
    def __init__(self, name='MACRO', params=''):
        self.name = name
        self.params = params

    def get_name(self):
        return self.name

    def get_params(self):
        return self.params

    def dump(self, f):
        f.write('#define  {:20s}  {:s}\n'.format(self.name, self. params))


def make_unique(I):
    UI = []
    for i in I:
        bFound = False
        for ui in UI:
            if i.issame(ui):
                bFound = True
                continue
        if not bFound:
            UI.append(i)

    return UI

class atom_dict:
    def __init__(self):
        self.d = dict()

    def store(self, a, at):
        self.d[a] = at

    def abstract(self, a):
        return self.d[a]


    def make_abstract_interaction(self, I):
        ai = interaction()
        ai.set_a([self.abstract(a) for a in I.get_atoms()])
        ai.set_func(I.func)
        ai.set_r(I.r)
        ai.set_k(I.k)
        return ai


class atom:
    def __init__(self, nr=0, atype='', rnr=0, rname='', aname='', cgnr=0, q=0.0, m=0.0):
        self.nr = nr
        self.atype = atype
        self.rnr = rnr
        self.rname = rname
        self.aname = aname
        self.cgnr = cgnr
        self.q = q
        self.m = m


    def from_itp_line(self, itp):
        sline = itp.split()

        self.nr    = int(sline[0])
        self.atype = sline[1]
        self.rnr   = int(sline[2])
        self.rname = sline[3]
        self.aname = sline[4]
        self.cgnr  = int(sline[5])
        self.q     = float(sline[6])
        self.m     = float(sline[7])


    def set_nr(self, nr):
        self.nr = nr

    def set_atype(self, atype):
        self.atype = atype

    def set_rnr(self, rnr):
        self.rnr = rnr

    def set_rname(self, rname):
        self.rname = rname

    def set_aname(self, aname):
        self.aname = aname

    def set_cgnr(self, cgnr):
        self.cgnr = cgnr

    def set_q(self, q):
        self.q = q

    def set_m(self, m):
        self.m = m

    def rtp_dump(self, ostream=stdout):
        s = '{:>6s} {:>6s}   {: 8.5f}  {:>4d}\n'.format(self.aname, self.atype, self.q, self.cgnr)
        ostream.write(s)

if __name__=='__main__':

    convertable_dihedrals = {1:'proper', 3:'Ryckaert-Bellemans', 5:'Fourier', 9:'proper (multiple)'}

    parser = argparse.ArgumentParser(description="Extracts force field parameters for e.g. GAFF topologies.")
    parser.add_argument('file', action='store', nargs=1, type=str, help='itp or rtp file')
    parser.add_argument('-cv', '--convert-dihedrals', action='store', nargs=1, type=str, metavar='"X Y ..."', dest='conv_str', help='Convert dihedrals of type X to Y. Multiple pairs can be provided.')
    parser.add_argument('-minf', '--min-force-constant', action='store', nargs=1, type=str, metavar='Fmin', dest='fmin', help='Remove converted dihedrals if force constant is below this value (default: 1e-5 kJ/mol).')
    parser.add_argument('-r', '--rtp-name', action='store', nargs=1, type=str, metavar='mymol.rtp', dest='rtp', help='Create rtp file.')
    parser.add_argument('-rb', '--rtp-bond-parameters', action='store_true', dest='rb', help='Write bond parameters in rtp file.')
    parser.add_argument('-ra', '--rtp-angle-parameters', action='store_true', dest='ra', help='Write angle parameters in rtp file.')
    parser.add_argument('-rd', '--rtp-dihedral-parameters', action='store_true', dest='rd', help='Write dihedral parameters in rtp file.')
    parser.add_argument('-m', '--rtp-dihedral-macros', action='store_true', dest='pm', help='Define and use macros for dihedral parameter values in rtp.')
    args = parser.parse_args()

    #fname='/Users/erikmarklund/Projects/Interface_evolution/Structures/topologies/SR1_cofactor_new/FF_CONV/BCC/GMX_SR1_cofactor_new/SR1_cofactor_new.top'

    fname = args.file[0]

    conv = []
    if args.conv_str:
        conv_split = args.conv_str[0].split()

        if len(conv_split)%2:
            warn('Conversion requires pairs fo dihedral types. Odd number given.', bError=True)
            exit(1)

        conv_ = [ int(c) for c in conv_split ]

        conv = dict(zip(conv_[::2], conv_[1::2]))

    if args.fmin:
        kmin = float(args.fmin[0])
    else:
        kmin = 0.0

    if args.rtp:
        rtp_name = args.rtp[0]
    else:
        rtp_name = ''

    bRtpB = args.rb
    bRtpA = args.ra
    bRtpD = args.rd
    bMacros = args.pm

    adict = atom_dict()
    idict = atom_dict()

    atomtypes = []
    atoms     = []
    bonds     = []
    angles    = []
    dihedrals = []

    with open(fname, 'r') as f:
        section = None
        try:
            for line in f:
                sline = (line.split(';')[0]).strip()

                if not sline:
                    continue

                if sline[0] == '[':
                    # Directive

                    directive = sline.strip('[ ]')
                    section = directive

                    continue

                ssline=sline.split()

                if section == 'atomtypes':
                    atomtypes.append(sline)

                if section == 'atoms':
                    adict.store(ssline[4], ssline[1])
                    idict.store(ssline[0], ssline[4])
                    a = atom()
                    a.from_itp_line(line)
                    atoms.append(a)
                    continue

                if section == 'bonds':
                    #pdb.set_trace()

                    b = interaction()
                    b.set_a([ssline[0], ssline[1]])
                    b.set_func(ssline[2])
                    b.set_r(ssline[3])
                    b.set_k(ssline[4])
                    bonds.append(b)
                    continue

                if section == 'angles':
                    a = interaction()
                    a.set_a([ssline[0], ssline[1], ssline[2]])
                    a.set_func(ssline[3])
                    a.set_r(ssline[4])
                    a.set_k(ssline[5])
                    angles.append(a)
                    continue

                if section == 'dihedrals':
                    d = interaction()
                    d.set_a([ssline[0], ssline[1], ssline[2], ssline[3]])
                    d.set_func(ssline[4])
                    d.set_r('0.0')
                    d.set_k(ssline[5:11])
                    dihedrals.append(d)
                    continue

        except IOError:
            warn("Failed to read from file "+fname, bError=True)
            exit(1)


    abonds = [idict.make_abstract_interaction(b) for b in bonds]
    tbonds = [adict.make_abstract_interaction(b) for b in abonds]

    aangles = [idict.make_abstract_interaction(a) for a in angles]
    tangles = [adict.make_abstract_interaction(a) for a in aangles]

    adihedrals = [idict.make_abstract_interaction(d) for d in dihedrals]
    tdihedrals = [adict.make_abstract_interaction(d) for d in adihedrals]

    ubonds = make_unique(tbonds)
    uangles = make_unique(tangles)
    udihedrals = make_unique(tdihedrals)


    output("=== ATOM TYPES ===")
    for a in atomtypes:
        output(a)

    if rtp_name:
        try:
            rtp_file = open(rtp_name, 'w')

            rtp_file.write('; Generated by derive_params.py, part of mkvsites.\n; Please cite:\n')
            c = cite_mkvsites()
            for line in c.split('\n'):
                rtp_file.write(';   '+line+'\n')
            rtp_file.write('\n')

            rtp_file.write('[ bondedtypes ]\n')
            rtp_file.write('; bonds  angles  dihedrals  impropers\n')
            bt = 0
            at = 0
            dt = 0
            it = 0
            if bonds:
                bt = bonds[0].get_func()
            if angles:
                at = angles[0].get_func()
            if udihedrals:
                dt = dihedrals[0].get_func()
                if conv:
                    idt = int(dt)
                    output('converting default {:d} to {:d}'.format(idt, conv[idt]))
                    if idt in conv.keys():
                        dt = str(conv[idt])
                it = dt
            rtp_file.write("{:>6d}{:>6d}{:>6d}{:>6d} ; Derived from first parameters found. Check!\n\n".format(int(bt),int(at),int(dt),int(it)))

            rtp_file.write("[ {:s} ]\n".format(atoms[0].rname))
            rtp_file.write(" [ atoms ]\n")
            for a in atoms:
                a.rtp_dump(ostream=rtp_file)

        except IOError:
            warn('Failed to write atoms to rtp file '+rtp_name, bError=True)
            raise

    output("=== INTERACTION TYPES ===")

    output("--- bonds ---")
    for b in ubonds:
        b.dump()

    if rtp_name:
        rtp_file.write("\n [ bonds ]\n")
        for b in abonds:
            b.rtp_dump_interaction(ostream=rtp_file, bDumpParams=bRtpB)

    output("--- angles ---")
    for a in uangles:
        a.dump()

    if rtp_name:
        rtp_file.write("\n [ angles ]\n")
        for a in aangles:
            a.rtp_dump_interaction(ostream=rtp_file, bDumpParams=bRtpA)

    output("--- dihedrals ---")
    if conv:
        output(' - Converting dihedrals - ')
        for c in conv.keys():
            output('   type {:d} to {:d}'.format(c,conv[c]))

    bKok = True

    if rtp_name:
        rtp_file.write("\n [ dihedrals ]\n")


    # Functionalise this since it is reused below
    for d in udihedrals:
        if conv:
            t = int(d.get_func())

            if t not in conv.keys():
                warn('Could not convert type {:d}, because it is not in list of conversions:'.format(t), bError=True)
                output(conv, ostream=stderr)
                exit(1)

            if t not in convertable_dihedrals.keys() or \
              conv[t] not in convertable_dihedrals.keys():
                warn('Dihedral type {:d} (or {:d}) cannot be converted (to). Convertable dihedrals are:'.format(t, conv[t]), bError=True)
                for k,v in convertable_dihedrals.items():
                    output('{:4d} - {:s}'.format(k, v))
                exit(1)

            # Make general dihedral, then convert.
            if type(d.k) == type([]):
                params = ' '.join([str(s) for s in d.k])
            else:
                if type(d.k) == type(''):
                    params = d.k
                else:
                    warn('d.k is of unsupported type {:s}'.format(type(d.k)), bError=True)
                    raise TypeError

            cd = dconv.cdihedral(atoms=d.a, ft=t, params=params)

            dnew = cd.make_specific()

            if type(dnew) != type([]):
                dnew = [dnew]

            for dn in dnew:
                # This 'switch' must match the supported conversions
                if (conv[t] == 1):
                    # dc = dn.make_proper()
                    raise ValueError('Most dihedrals cannot be converted to proper. Try proper multiple (ft=9)')
                elif (conv[t] == 3):
                    dc = dn.make_rb()
                elif (conv[t] == 5):
                    dc = dn.make_fourier()
                elif (conv[t] == 9):
                    dc = dn.make_propermultiple()
                else:
                    raise ValueError('Unsupported dihedral type')

                if type(dc) != type([]):
                    dc = [dc]
                for dcc in dc:
                    ft = dcc.get_ft()
                    dni = interaction()
                    dni.set_a(dcc.atoms)
                    dni.set_func(str(ft))
                    pp = dcc.get_params()
                    if type(pp) == type([]):
                        dni.set_k([str(p) for p in dcc.get_params()])
                    else:
                        dni.set_k(str(pp))

                    if conv[t] in (1,9):
                        bSkipit = False
                        pparam = dcc.get_params()
                        if abs(pparam[1]) < kmin:
                            bKok = False
                            bSkipit = True

                    dni.dump(bSkip=bSkipit)

        else:
            d.dump()


    if rtp_name:
        d_macros = []
        for pnr,d in enumerate(adihedrals):
            if conv:
                # Skip all checks. They were passed above.
                t = int(d.get_func())

                # Make general dihedral, then convert.
                if type(d.k) == type([]):
                    params = ' '.join([str(s) for s in d.k])
                else:
                    params = d.k

                cd = dconv.cdihedral(atoms=d.a, ft=t, params=params)
                dnew = cd.make_specific()

                if type(dnew) != type([]):
                    dnew = [dnew]

                for dn in dnew:
                    # This 'switch' must match the supported conversions
                    if (conv[t] == 1):
                        # dc = dn.make_proper()
                        raise ValueError('Most dihedrals cannot be converted to proper. Try proper multiple (ft=9)')
                    elif (conv[t] == 3):
                        dc = dn.make_rb()
                    elif (conv[t] == 5):
                        dc = dn.make_fourier()
                    elif (conv[t] == 9):
                        dc = dn.make_propermultiple()
                    else:
                        raise ValueError('Unsupported dihedral type')

                    if type(dc) != type([]):
                        dc = [dc]
                    for dcc in dc:
                        ft = dcc.get_ft()
                        dni = interaction()
                        dni.set_a(dcc.atoms)
                        dni.set_func(str(ft))
                        pp = dcc.get_params()
                        if type(pp) == type([]):
                            dni.set_k([str(p) for p in dcc.get_params()])
                        else:
                            dni.set_k(str(pp))

                        if conv[t] in (1,9):
                            bSkipit = False
                            pparam = dcc.get_params()
                            if abs(pparam[1]) < kmin:
                                bKok = False
                                bSkipit = True

                        if not bSkipit:
                            if bMacros:
                                macroname = 'MKV_{:s}_{:03d}'.format(atoms[0].rname, pnr)
                                pmacro = dni.make_macro(name=macroname)
                                d_macros.append(pmacro)
                            else:
                                pmacro = None

                            dni.rtp_dump_interaction(ostream=rtp_file, bDumpParams=bRtpD, macro=pmacro)

            else:
                if bMacros:
                    macroname = 'MKV_{:s}_{:03d}'.format(atoms[0].rname, pnr)
                    pmacro = d.make_macro(name=macroname)
                    d_macros.append(pmacro)
                else:
                    pmacro = None

                d.rtp_dump_interaction(ostream=rtp_file, bDumpParams=bRtpD, macro=pmacro)

        if d_macros:
            mfname, mfext = path.splitext(rtp_file.name)
            mfname = mfname + '_macros.txt'
            with open(mfname, 'w') as mf:
                try:
                    mf.write('; Add these dihedral parameters to ffbonded.itp\n;\n')
                    for m in d_macros:
                        m.dump(mf)
                except IOError:
                    warn('Failed to write to file {:s}.'.format(mfname), bError=True)
                    raise

    if conv:
        output('NOTE: Not all dihedrals contain constant energy offsets,')
        output('      hence the dihedral energy might differ, but that')
        output('      matters little for the force and hence the dynamics.')
        output('      If the dihedrals are used for, e.g. lambda simulations,')
        output('      this might be a problem however, since the offsets')
        output('      might differ for different lambda values.')
    if not bKok:
        output('NOTE: Some dihedrals have small force coefficients for ')
        output('      numerical reasons, lower than the stipulated cut-off.')
        output('      These have been commented out above.')

    output('################################################')
    output(' Note that new atomtypes have nonsense masses.')
    output(' We cannot safely derive elements/masses in all')
    output(' cases, so that is left to the user to sort out')
