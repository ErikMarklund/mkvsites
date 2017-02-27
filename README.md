# README #

### What is this repository for? ###

Derives Gromacs virtual site parameters (angle constraint, dummy constraints) for an itp/rtp file based on force field data.

[Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

Install Gromacs if needed, set up paths with GMXRC. Install python if needed.

* Configuration

issue `source <path-to-gromacs-installation>/bin/GMXRC.bash` (assuming bash)

* Dependencies

    + Gromacs
    + Python

* How to run tests

`python topology.py`

* Deployment instructions

    The amber99sb-ildn force field is used by default. -ff can be used to specify another force field.

    + For itp files: `python mkvsites.py -ff=charmm36.ff <path-to>/topology.itp`
    + For rtp files: `python mkvsites.py -ff=charmm36.ff -res=POPE $GMXDATA/top/charmm36.ff/merged.rtp`

### Contribution guidelines ###

### Who do I talk to? ###

* egmarklund