#!/bin/sh

# Your job name 
#$ -N single-tt-run-trpn-koff-0.02_0

# Resources
#$ -l h_cpu=48:00:00

# Use current working directory
#$ -cwd

# Join stdout and stderr
#$ -j y

# Execute in current environment
#$ -V

# Set queue
#$ -q all.q

# Set your number of processors here. 
# Requests mpich environment although we actually are using openmpi
#$ -pe orte 2

# add /usr/local/bin in path
PATH=$PATH:/usr/local/bin

# Run job through bash shell
#$ -S /bin/bash

# enter the job-directory
cd /home/huskeypm/single-tt

# setup environment
source ~/.bashrc
source /share/apps/si2011/mesoscale/mesoscale.conf
export DOLFIN_NOPLOT=1
export SUBCELL_NOPLOT=1




