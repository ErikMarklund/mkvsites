#!/usr/bin/env python

# MkVsites

# Copyright © 2020 Erik Marklund

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
from matplotlib import pyplot as plt
import numpy as np
from sys import stdout

class DihedralError(Exception):
    pass


def plot_dihedral(D, ax=None):
    if not ax:
        ax = plt.gca()

    npoints = 100
    ran = [-pi, pi]

    #x = [ ran[0]+(xx/100.)*(ran[1]-ran[0]) for xx in range(101) ]
    x = np.arange(ran[0], ran[1], (ran[1]-ran[0])/npoints)
    y = np.zeros(npoints)
    if type(D) == type([]):
        for d in D:
            y += np.array([d.calcpot(xx) for xx in x ])
    else:
        y = np.array([D.calcpot(xx) for xx in x ])

    ax.plot(x,y)


# base class for dihedrals
class cdihedral:
    def __init__(self, atoms=[], ft=0, params=''):
        self.atoms = atoms
        self.ft = ft
        self.params = params

    def __make_proper(self, bMultiple=False):
        if bMultiple:
            pdtype = pmdihedral
        else:
            pdtype = pdihedral

        pd = pdtype(atoms=self.atoms)
        pspl = self.params.strip().split()
        pd.set_phi(float(pspl[0]))
        pd.set_k(float(pspl[1]))
        pd.set_m(float(pspl[2]))
        return pd

    def __make_fourier(self):
        pspl = self.params.strip().split()
        F = []
        for p in self.params.strip().split():
            F.append(float(p))

        pf = fdihedral(atoms=self.atoms, F=F)
        return pf

    def __make_RB(self):
        pspl = self.params.strip().split()
        C = []
        for p in self.params.strip().split():
            C.append(float(p))

        pd = rbdihedral(atoms=self.atoms, C=C)
        return pd

    def make_specific(self):
        ft = self.ft
        if (ft == 1):
            d = self.__make_proper(bMultiple=False)
        elif(ft == 3):
            d = self.__make_RB()
        elif(ft == 5):
            d = self.__make_fourier()
        elif(ft == 9):
            d = self.__make_proper(bMultiple=True)
        else:
            raise DihedralError("Unsupported dihedral function type")
        return d

    def get_ft(self):
        return self.ft

    def get_params(self):
        return self.params

    def calcpot(self, PHI):
        return 0

    def plot(self, ax=None):
        plot_dihedral(self)

    def dump(self, ostream=stdout):
        for a in self.atoms:
            ostream.write('{:5s}'.format(str(a)))
        ostream.write('{:>5d}'.format(self.ft))
        self.dump_specific(ostream=ostream)
        ostream.write('\n')

    def dump_specific(self, ostream=stdout):
        ostream.write(self.params)


# Regular proper dihedral
# k(1+cos(n*PHIi-phi))
# PHI is actual angle
class pdihedral(cdihedral):
    def __init__(self, atoms=[], phi=0.0, k=0.0, m=0):
        self.atoms = atoms
        self.ft = 1
        self.phi = phi
        self.k = k
        self.m = m

    def set_phi(self, phi):
        self.phi = phi

    def set_k(self, k):
        self.k = k

    def set_m(self, m):
        self.m = m

    def get_params(self):
        return [self.phi, self.k, self.m]

    def make_fourier(self):
        if self.phi not in (0.0, 180.0):
            raise DihedralError("Unsupported angle phi={:f}".format(self.phi))

        F1 = 0.0
        F2 = 0.0
        F3 = 0.0
        F4 = 0.0

        m = self.m
        k = self.k
        if (m == 1):
            F1 = -2*k
        elif (m == 2):
            F2 = 2*k
        elif (m == 3):
            F3 = -2*k
        elif (m == 4):
            F4 = 2*k
        else:
            raise DihedralError('Unsupported "multiplicity" for Fourier dihedral: {:d}'.format(m))

        F = [F1, F2, F3, F4]
        fd = fdihedral(atoms=self.atoms, F=F)
        return fd

    def make_rb(self):
        if self.phi not in (0.0, 180.0):
            raise DihedralError("Unsupported angle phi={:f}".format(self.phi))

        df = self.make_fourier()
        drb = df.make_rb()
        return drb

    def calcpot(self, PHI):
        try:
            V = self.k * (1+cos(self.m * PHI - self.phi*pi/180.))
        except:
            print self.k
            print self.m
            print self.phi
            raise
        return V

    def dump_specific(self, ostream=stdout):
        ostream.write('  {: 6.4e}'.format(self.phi))
        ostream.write('  {: 6.4e}'.format(self.k))
        ostream.write('  {:6d}'.format(self.m))


# Proper dihedral multiple
class pmdihedral(pdihedral):
    def __init__(self, atoms=[], phi=0.0, k=0.0, m=0):
        pdihedral.__init__(self, atoms=atoms, phi=phi, k=k, m=m)
        self.ft = 9


# Fourier dihedral
# 1/2 * [ F1 (1+cos(PHI)) +
#         F2 (1+cos(2*PHI)) +
#         F3 (1+cos(2*PHI)) +]
#         F4 (1+cos(2*PHI)) ]
class fdihedral(cdihedral):
    def __init__(self, atoms=[],  F=[.0, .0, .0, .0]):
        self.atoms = atoms
        self.ft = 5
        self.F = F

    def set_F(self, F):
        self.F = F

    def get_params(self):
        return self.F

    def make_rb(self):
        F = self.F
        C = [ \
                  F[1] + 0.5*(F[0] + F[2]), \
                  0.5*(-1*F[0] + 3*F[2]), \
                  -1*F[1] + 4*F[3], \
                  -2*F[2], \
                  -4*F[3], \
                  0.0 \
            ]

        rb = rbdihedral(atoms=self.atoms, C=C)
        return rb

    def make_propermultiple(self):
        PD = []

        # If all coefficients are zero, just output a single flat dihedral to make grompp happy
        bFlat = True
        for f in self.F:
            if f != 0.0:
                bFlat = False
                break

        if bFlat:
            pd = pmdihedral(atoms=self.atoms, k=0.0, phi=0.0, m=1)
            return [pd]

        for i,f in enumerate(self.F):
            if f == 0.0:
                continue

            m = i+1
            if m%2 == 0:
                # even
                inv = -1
            else:
                # odd
                inv = 1

            k = 0.5*f*inv
            phi = 0.0
            pd = pmdihedral(atoms=self.atoms, k=k, phi=phi, m=m)
            PD.append(pd)

        return PD

    def calcpot(self, PHI):
        F = self.F
        F1 = F[0]
        F2 = F[1]
        F3 = F[2]
        F4 = F[3]

        V = 0.5 * ( \
                        F1 * (1+cos(PHI)) + \
                        F2 * (1-cos(2*PHI)) + \
                        F3 * (1+cos(3*PHI)) + \
                        F4 * (1-cos(4*PHI)) \
                    )

        return V

    def dump_specific(self, ostream=stdout):
        for f in self.F:
            ostream.write('  {: 6.4e}'.format(f))

# Ryckaert-Belleman dihedral
class rbdihedral(cdihedral):
    def __init__(self, atoms=[], ft=3, C=[.0, .0, .0, .0, .0, .0]):
        self.atoms = atoms
        self.ft = ft
        self.C = C

    def set_C(self, C):
        self.C = C

    def get_params(self):
        return self.C

    def make_fourier(self):
        C = self.C

        F4 = -.25*C[4]
        F3 = -.5*C[3]
        F2 = -C[2]+4*F4
        F1 = -2*C[1] + 3*F3
        F = [F1,F2,F3,F4]

        df = fdihedral(atoms=self.atoms, F=F)
        return df

    def make_propermultiple(self):
        df = self.make_fourier()
        DP = df.make_propermultiple()
        return DP

    def calcpot(self, PHI):
        psi = PHI-pi
        #psi = PHI
        C = self.C
        V = 0
        for i,c in enumerate(C):
            V += c * (cos(psi)**i)

        return V

    def dump_specific(self, ostream=stdout):
        for c in self.C:
            ostream.write('  {: 6.4e}'.format(c))


def test():
    dc = cdihedral(atoms=[1,2,3,4], ft=1, params="180. 1000 1")

    dc.plot()

    dp = dc.make_specific()
    dp.plot()
    plt.savefig('dc_dp.pdf')
    plt.close()

    df = dp.make_fourier()
    df.plot()
    plt.savefig('dp_df.pdf')
    plt.close()

    drb = df.make_rb()
    drb.plot()
    plt.savefig('df_drb.pdf')
    plt.close()
    print drb.C

    drb.set_C([9.28, 12.16, -13.12, -3.06, 26.24, 0.0])
    drb.plot()
    drb.dump()
    plt.savefig('drb.pdf')
    plt.close()

    df = drb.make_fourier()
    df.dump()
    df.plot()
    plt.savefig('drb_df.pdf')
    plt.close()

    dp = df.make_propermultiple()
    plot_dihedral(dp)
    for d in dp:
        d.dump()
    plt.savefig('df_dp.pdf')
    plt.close()


if __name__ == '__main__':
    test()
