#!/bin/bash

# $1 = cel dir

cwd=`pwd`

#for genotype in M #D
#do
#    for tissue in TES #KID SPN
#    do
#	for data in $1/$genotype[1-9]*$tissue.CEL
#	do
#	    for method in ITER_PLIER #PLIER GCRMA
#	    do
#		for level in core #mps
#		do
#		    for treatment in random SNP
#		    do
#	    	        base=`	basename $data`
#		        echo -n "Submitting $base; $method; $treatment: "
#		        qsub -l nodes=1:ppn=1:brno:x86_64,mem=8gb,walltime=1:16:00:00,scratch=2gb -N $level-$base-$method-$treatment -m a -j oe -d $cwd/jobs -v cwd=$cwd,genotype=$genotype,tissue=$tissue,method=$method,data=$base,level=$level,treatment=$treatment $cwd/mcmc.sh
#		    done
#		done
#	    done
#	done
#    done
#done

for genotype in D M
do
    for tissue in SPN KID TES
    do
	for data in $1/$genotype[1-9]*$tissue.CEL
	do
	    for method in PLIER ITER_PLIER #GCRMA
	    do
		for level in mps #core
		do
		    for treatment in SNP random
		    do
	    	        base=`basename $data`
		        echo -n "Submitting $base; $method; $treatment: "
		        qsub -l nodes=1:ppn=1:brno:x86_64,mem=4gb,walltime=1:16:00:00,scratch=2gb -N $level-$base-$method-$treatment -m a -j oe -d $cwd/jobs -v cwd=$cwd,genotype=$genotype,tissue=$tissue,method=$method,data=$base,level=$level,treatment=$treatment $cwd/mcmc.sh
		    done
		done
	    done
	done
    done
done
