#!/bin/bash
#PBS -N intronic_read_overlap
#PBS -q hotel
#PBS -o /projects/ps-jcvi/projects/Lasken_Single_Cell/scripts/intronic_read_overlap.run.oe
#PBS -j oe
#PBS -M kunalbhutani@gmail.com
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -A schork-group
#PBS -m abe

cd /projects/ps-jcvi/projects/Lasken_Single_Cell/scripts
echo "<startTime>"`date`"</startTime>"
echo "<output>"
sh exonic_all.sh
echo "</output>"
echo "<exitStatus>"$?"</exitStatus>"
echo "<stopTime>"`date`"</stopTime>"
qstat -f $PBS_JOBID | grep Job
qstat -f $PBS_JOBID | grep Resource
