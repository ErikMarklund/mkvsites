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
    + Python

* How to run tests

`python topology.py`

* Deployment instructions

    + For itp files: `python mkvsites.py -ff <path-to>/<forcefield>.ff topology.itp`
    + For rtp files: `python mkvsites.py -ff <path-to>/<forcefield>.ff -res=POPE residues.rtp`

See mkvsites.py -h and topinspect.py -h for more information.

### Contribution guidelines ###

### Who do I talk to? ###

* egmarklund
