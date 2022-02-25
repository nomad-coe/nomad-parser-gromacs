#!/bin/bash

#####################################################################

#System Specific
sys="2awx"
######################################################################

#Directories
bin="/sw/linux/gromacs/5.1.2/bin"

######################################################################

#Executables
gmx="$bin/gmx"
#grompp="$bin/grompp"
#mdrun="$bin/mdrun"
#tpbconv="$bin/tpbconv"

######################################################################

#File Names
mdp="${sys}.NPT.mdp"
mdout="mdout-${sys}.mdp"
gro="${sys}.NPT-out.gro"
top="${sys}.top"
edr="${sys}.edr"
tpr="${sys}.tpr"
oldtrr="${sys}.PAR.trr"
trr="${sys}.trr"
xtc="${sys}.xtc"
ndx="index_types.ndx"
gro_out="${sys}.confout.gro"
md_log="mdlog-${sys}.log"
mpirun="mpirun"
gr_log="grompp.log"
md_log="mdrun.log"
cpi="state_old.cpt"
time="1000"
frame="4641"
table="table.xvg"

#####################################################################

# Determine processors for mpirun
NODELIST=$( cat $PBS_NODEFILE )
echo Processors used
count=0
for i in $NODELIST
do
    let count=$count+1
    echo Processor: $count $i
done
nnodes=$count

nnodes=16
#####################################################################

#grompp
#Normal
$gmx grompp -f $mdp -c $gro -p $top -o $tpr -po $mdout\
        >& $gr_log

#Extend
#$grompp -f $mdp -c $oldtpr -p $top -o $newtpr -po $mdout\
#        >& $gr_log

#$tpbconv -s $newtpr -extend $time -o $new2tpr 

#use trr file for starting conf and velocities
#$grompp -f $mdp -c $gro -t $oldtrr -n index.ndx -time $frame -p $top -o $tpr -po $mdout\
#        >& $gr_log
#Tables
#$grompp -f $mdp -c $gro -p $top -o $tpr -po $mdout\
#        -n $ndx >& $gr_log

#Suffle
#$grompp -f $mdp -c $gro -p $top -o $tpr -po $mdout\
#        -shuffle  -deshuf $gr_des >& $gr_log

#####################################################################

#mdrun
#normal
#$MPIEXEC 
$gmx mdrun -s $tpr -o $trr -c $gro_out -nt $nnodes -e $edr -g $md_log >& $md_log
#$mdrun -s $tpr -o $trr -x $xtc -c $gro_out -e $edr -g $md_log >& $md_log

#$mdrun -s $sys -o $sys -x $sys -c $gro_out -e $sys -g $md_log >& $md_log

#normal-old-CGmpirun
#$mpirun $mdrun -pd -s $tpr -o $trr -c $gro_out -table $table -tableb $table -n $nnodes -e $edr -g $md_log >& $md_log --mca pml ob1

#extend
#$mpirun $mdrun -s $new2tpr -cpi $cpi -o $trr -c $gro_out -n $nnodes -e $edr -g $md_log >& $md_log

#continue
#$mdrun -s $tpr -o $trr -cpi $cpi -noappend -nt $nnodes -c $gro_out -e $edr -g $md_log >& $md_log




