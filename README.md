# README #

### What is this repository for? ###

Derives Gromacs virtual site parameters (angle constraint, dummy constraints) for an itp/rtp file based on force field data.

[Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

Install Gromacs if needed, set up paths with GMXRC. Install python if needed.

* Configuration

* Dependencies

    + Gromacs
    + Python 3

* How to run tests

`python topology.py`

* Deployment instructions

    + For itp files: `python mkvsites.py -ff <path-to>/<forcefield>.ff topology.itp`
    + For rtp files: `python mkvsites.py -ff <path-to>/<forcefield>.ff -res=POPE residues.rtp`

See mkvsites.py -h and topinspect.py -h for more information.

### Contribution guidelines ###

### License ###

MkVsites is distributed under an "MIT license":

MkVsites

Copyright © 2020 Erik Marklund

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### Who do I talk to? ###

* egmarklund
