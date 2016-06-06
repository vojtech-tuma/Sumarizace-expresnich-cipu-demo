#!/bin/bash

# Run in this folder after all summarization is completed

scripts=`pwd`/scripts

for genotype in D M
do
    for tissue in KID SPN TES
    do
	for method in DABG PLIER ITER_PLIER GCRMA
	do
	    # Core with SNP
	    pushd out/core_snp/$method/${genotype}_${tissue}_${method} > /dev/null	    
	    bash $scripts/join_summary
	    popd > /dev/null

	    # Core without SNP
	    pushd out/core_nosnp/$method/${genotype}_${tissue}_${method} > /dev/null	    
	    bash $scripts/join_summary
	    popd > /dev/null

	    # Control with ranom kill list
	    pushd out/core_random/$method/${genotype}_${tissue}_${method} > /dev/null	    
	    bash $scripts/join_summary
	    popd > /dev/null

	    # MPS with SNP
	    pushd out/core_mps_snp/$method/${genotype}_${tissue}_${method} > /dev/null	    
	    bash $scripts/join_summary
	    popd > /dev/null

	    # MPS without SNP
	    pushd out/core_mps_nosnp/$method/${genotype}_${tissue}_${method} > /dev/null	    
	    bash $scripts/join_summary
	    popd > /dev/null

	    # MPS with ranom kill list
	    pushd out/core_mps_random/$method/${genotype}_${tissue}_${method} > /dev/null	    
	    bash $scripts/join_summary
	    popd > /dev/null


	done
    done
done
