# bash script for fast testing. Issue 'source test_DT.sh'

if [[ -z ${GMXDATA} ]]; then
    echo 'Need access to a working copy of Gromacs.'
    echo 'issue "source path/to/gromacs/bin/GMXRC.yourshell"'
else

    gmx pdb2gmx -f CATG.pdb -o p2g -i p2g -p p2g -vsite hydrogens -ignh
    gmx editconf -f p2g.gro -o inbox.pdb -d 1 -bt dodecahedron

    cp p2g.top solvated.top

    gmx solvate -cp inbox.pdb -o solvated.pdb -p solvated.top

    gmx grompp -f genion.mdp -o genion.tpr -p solvated.top -c solvated.pdb -maxwarn 2
fi
