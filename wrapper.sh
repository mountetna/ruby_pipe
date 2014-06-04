#!/bin/bash
#$ -S /bin/bash

if [ -f ~/.bash_profile ] ; then
    . ~/.bash_profile
fi

cd $SGE_O_WORKDIR

eval $*
c=$?
exit $c
