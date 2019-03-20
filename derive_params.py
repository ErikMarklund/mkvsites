#!/usr/bin/env python

from sys import exit
import pdb

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

    def dump(self):
        s = "    "
        for a in self.a:
            s = s + "{:4s}".format(a)

        s = s + "{:4s}".format(self.func)

        if self.get_type() == 'dihedral':
            for p in self.k:
                s = s + "    {:s}".format(p)
        else:
            s = s + "    {:s}".format(self.r)
            s = s + "    {:s}".format(self.k)

        print s

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

if __name__=='__main__':
    fname = 'SR1_cofactor.top'

    adict = atom_dict()
    idict = atom_dict()

    atomtypes = []
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
            print "Failed to read from file "+fname
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

    print "=== ATOM TYPES ==="
    for a in atomtypes:
        print a

    print "=== INTERACTION TYPES ==="

    print "--- bonds ---"
    for b in ubonds:
        b.dump()

    print "--- angles ---"
    for a in uangles:
        a.dump()

    print "--- dihedrals ---"
    for d in udihedrals:
        d.dump()
