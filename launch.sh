#!/bin/bash

#THIS IS INTENDED TO BE RUN IN THIS FOLDER! Expects "scripts" folder here.
./scripts/prepare_data ~/moex/affy/annot/dt/gff/ ~/moex/snp_data/*final*annots.vcf 0

# Core with SNP
for genotype in D M
do
    for tissue in KID SPN TES
    do
	for method in DABG PLIER ITER_PLIER GCRMA
	do
	    ./scripts/submit_job $genotype $tissue $method ~/moex/moex2010/lists/ps/core ~/moex/moex2010/cel_genotype ~/moex/moex2010/release/out/empty.kill ~/moex/affy/annot/lib ~/moex/moex2010/release/out/core_snp
	    ./scripts/submit_job $genotype $tissue $method ~/moex/moex2010/lists/mps/core ~/moex/moex2010/cel_genotype ~/moex/moex2010/release/out/empty.kill ~/moex/affy/annot/lib ~/moex/moex2010/release/out/core_mps_snp
	done
    done
done

# Core without SNP
for genotype in D M
do
    for tissue in KID SPN TES
    do
	for method in DABG PLIER ITER_PLIER GCRMA
	do
	    ./scripts/submit_job $genotype $tissue $method ~/moex/moex2010/lists/ps/core ~/moex/moex2010/cel_genotype ~/moex/moex2010/release/out/probes.kill ~/moex/affy/annot/lib ~/moex/moex2010/release/out/core_nosnp
	    ./scripts/submit_job $genotype $tissue $method ~/moex/moex2010/lists/mps/core ~/moex/moex2010/cel_genotype ~/moex/moex2010/release/out/probes.kill ~/moex/affy/annot/lib ~/moex/moex2010/release/out/core_mps_nosnp
	done
    done
done

# Control with random kill list
for genotype in D M
do
    for tissue in KID SPN TES
    do
	for method in DABG PLIER ITER_PLIER GCRMA
	do
	    ./scripts/submit_job $genotype $tissue $method ~/moex/moex2010/lists/ps/core ~/moex/moex2010/cel_genotype ~/moex/moex2010/release/out/random.kill ~/moex/affy/annot/lib ~/moex/moex2010/release/out/core_random
	    ./scripts/submit_job $genotype $tissue $method ~/moex/moex2010/lists/mps/core ~/moex/moex2010/cel_genotype ~/moex/moex2010/release/out/random.kill ~/moex/affy/annot/lib ~/moex/moex2010/release/out/core_mps_random
	done
    done
done

		