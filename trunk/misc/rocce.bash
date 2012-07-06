#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -M huskeypm@mccammon.ucsd.edu
#$ -e smol.err
#$ -o smol.out


echo "Execute host is "
/bin/hostname


export CASE=tnc_isolated
export SMOLPATH=/home/huskeypm/sources/dolfin_smol/
source $SMOLPATH/config.bash 

# conversion 
export SMOLTODOLFIN=$SMOLPATH/smoltodolfin.py
export SMOL=$SMOLPATH/smol.py
export VALID=$SMOLPATH/validation.py

export B=`ls -t potential-*dx | perl -ane 'chomp $_; print "-p $_ "'`
export PYTHON=/share/apps/python2.7/bin/python
##$PYTHON $SMOLTODOLFIN -mcsf $CASE.m $B > smoltodolfin.out


# simulation 
$PYTHON $VALID run > validation.out

