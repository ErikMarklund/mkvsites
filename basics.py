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


from sys import stderr, stdout

# Some basic functions

class boolBase(object):
    """Generic base class with method for assessing True/False.
    If none of its values evaluates to True, the object is false.
    """
    def __nonzero__(self):
        V = vars(self)
        for v in V.values():
            if v:
                return True
        return False

def warn(s, bEndline=True, bWarn=True, bError=False, ostream=stderr, bPrint=True):
    """Issues a warning or error to stderr or other stream."""
    if not bPrint:
        return

    if bEndline:
        s += '\n'

    if bError:
        s = 'ERROR: '+s
    elif bWarn:
        s = 'WARNING: '+s
    else:
        pass

    stderr.write(s)
    
    
def output(s, bEndline=True, ostream=stdout, bPrint=True):
    """Print to stdout or other stream."""
    if not bPrint:
        return

    if bEndline:
        s += '\n'

    ostream.write(s)


def cite_mkvsites():
    """Returns a string with citataion info"""
    s = 'Per Larsson, Rosita Carolina Kneisz, Erik G Marklund\n' \
      + 'MkVsites: A general-purpose python tool to increase\n' \
      + 'performance  in  all-atom  molecular  dynamics simulations\n' \
      + 'Journal X(XX):XXX--XXXX.'
    return s
