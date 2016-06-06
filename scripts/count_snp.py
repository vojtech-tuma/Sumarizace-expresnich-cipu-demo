### COUNT SNP
import csv
import sys

try:
    quality_treshold = sys.argv[1]
except Exception:
    quality_treshold = 0
    pass

with open("probes.snp", 'wb') as probes, open("probesets.snp", 'wb') as probesets, open("wsb_snphits", 'wb') as wsb_snps, open("pwk_snphits", 'wb') as pwk_snps, open("vcf.error", 'wb') as vcf_error, open("vcf.snp", 'rb') as vcf, open("probesets", 'rb') as gff:
    #probes
    fieldnames = ["chrom","probe_id","probeset_id","PWK_SNP","WSB_SNP"]
    probe_writer = csv.DictWriter(probes, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    probe_writer.writeheader()

    #probesets
    fieldnames = ["chrom","probeset_id","PWK_SNP","WSB_SNP"]
    probeset_writer = csv.DictWriter(probesets, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    probeset_writer.writeheader()
    
    #wsb_snps
    fieldnames = ["chrom", "probe_id", "probeset_id", "REF", "ALT", "SNP", "QUAL"]
    wsb_snps_writer = csv.DictWriter(wsb_snps, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    wsb_snps_writer.writeheader()
    
    #pwk_snps
    fieldnames = ["chrom", "probe_id", "probeset_id", "REF", "ALT", "SNP", "QUAL"]
    pwk_snps_writer = csv.DictWriter(pwk_snps, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    pwk_snps_writer.writeheader()
    
    #vcf_error
    fieldnames = ["chrom", "probe_id", "WSB SNP", "WSB QUAL", "PWK SNP", "PWK QUAL", "vcf REF", "gff probe base", "snp in probe", "probe sequence"]
    vcf_error_writer = csv.DictWriter(vcf_error, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    vcf_error_writer.writeheader()
    
    #Write SNP info in pseudo-vcf file
    #Format: chrom    POS    REF    ALT    PWK    WSB (tab separated; sorted by chrom and by pos for every chrom)
    vcf_reader = csv.DictReader(vcf, delimiter='\t', quotechar='"')
    snp = next(vcf_reader, "END")
    
    #Write probe info in pseudo-gff file
    #Format: chrom    start    end    probe_id    probeset_id    (tab separated; sorted by chrom and by start for every chrom)
    gff_reader = csv.DictReader(gff, delimiter='\t', quotechar='"')
    probe = next(gff_reader, "END")
    
    sys.stderr.write('chromozome ' + probe["chrom"] + '\n')                
    probeset_wsb_count = 0
    probeset_pwk_count = 0
    
    while ((snp != "END") and (probe != "END")):
        if(snp["CHROM"] != probe["chrom"]):
                sys.stderr.write("1 " + '\t'.join(probe.values()) + '\t' + "  ---  " + '\t'.join(snp.values()) + '\n')       
        wsb_snp_count = 0
        pwk_snp_count = 0
        
        start = int(probe["start"])
        end = int(probe["end"])
        pos = int(snp["POS"])
        
        #skip snps at positions < start
        while (snp != "END" and 
                   ((pos < start and probe["chrom"] == snp["CHROM"]) or 
                    (int(snp["CHROM"])<int(probe["chrom"]) 
                        if (snp["CHROM"].isdigit() and probe["chrom"].isdigit()) 
                        else (True 
                            if (snp["CHROM"].isdigit() and not probe["chrom"].isdigit())
                            else snp["CHROM"] < probe["chrom"]
                            )
                         )
                    )
                ): ## POZOR!! neosetruje pripad, kdyz SNP prekroci chromozom, ale existuje vice prob ve starem chromozomu
            if(snp["CHROM"] != probe["chrom"]):
                sys.stderr.write("skip " + '\t'.join(probe.values()) + '\t' + "  ---  " + '\t'.join(snp.values()) + '\n')       
            snp = next(vcf_reader, "END")
            pos = int(snp["POS"]) if snp != "END" else 0  
            
        
        while (snp != "END" and (pos >= start and pos <= end and (snp["CHROM"] == probe["chrom"]))):
            pwk_snp = snp["PWK"].split(':')
            if (pwk_snp[0] != "0/0"):
                if(pwk_snp[1] >= quality_treshold):  # only count SNPs with quality above treshold
                   pwk_snp_count += 1
                pwk_snps_writer.writerow({"chrom":probe["chrom"], "probe_id":probe["probe_id"],"probeset_id":probe["probeset_id"],"REF":snp["REF"],"ALT":snp["ALT"],"SNP":pwk_snp[0], "QUAL":pwk_snp[1]})
                
            wsb_snp = snp["WSB"].split(':')
            if (wsb_snp[0] != "0/0"):
                if(wsb_snp[1] >= quality_treshold):  # only count SNPs with quality above treshold
                    wsb_snp_count += 1
                wsb_snps_writer.writerow({"chrom":probe["chrom"], "probe_id":probe["probe_id"],"probeset_id":probe["probeset_id"],"REF":snp["REF"],"ALT":snp["ALT"],"SNP":wsb_snp[0], "QUAL":wsb_snp[1]})
    
            if( probe["sequence"][pos-start] != snp["REF"].lower() ):
                vcf_error_writer.writerow({"chrom":probe["chrom"], "probe_id":probe["probe_id"], "WSB SNP":wsb_snp[0], "WSB QUAL":wsb_snp[1], "PWK SNP":pwk_snp[0], "PWK QUAL":pwk_snp[1], "vcf REF":snp["REF"], "gff probe base":probe["sequence"][pos-start].upper(), "snp in probe":str(pos-start), "probe sequence": probe["sequence"][0:pos-start] + probe["sequence"][pos-start].upper() + probe["sequence"][pos-start:]})
            snp = next(vcf_reader, "END")
            pos = int(snp["POS"]) if snp != "END" else 0 
            if(snp["CHROM"] != probe["chrom"]):
                sys.stderr.write("2 " + '\t'.join(probe.values()) + '\t' + "  ---  " + '\t'.join(snp.values()) + '\n')       
        probeset_pwk_count += pwk_snp_count
        probeset_wsb_count += wsb_snp_count
        
        chrom = False
        while(not chrom):
            
            probe_writer.writerow({"chrom":probe["chrom"],"probe_id":probe["probe_id"],"probeset_id":probe["probeset_id"],"PWK_SNP":pwk_snp_count,"WSB_SNP":wsb_snp_count})
            new_probe = next(gff_reader, "END")
            if(new_probe != "END"):
                if(new_probe["chrom"] != probe["chrom"]):
                    sys.stderr.write("# " + '\t'.join(probe.values()) + '\t' + "  ---  " + '\t'.join(new_probe.values()) + '\n')       
                    #print "chromozome " + new_probe["chrom"]
                    chrom = True
                    sys.stderr.write('chromozome ' + new_probe["chrom"] + '\n')                
                if((new_probe["probeset_id"] != probe["probeset_id"]) or new_probe == "END" ):# or snp == "END"):
                    probeset_writer.writerow({"chrom":probe["chrom"], "probeset_id":probe["probeset_id"],"PWK_SNP":probeset_pwk_count, "WSB_SNP":probeset_wsb_count})
                    probeset_pwk_count = 0
                    probeset_wsb_count = 0
                    
            probe = new_probe
            if(snp != "END" and probe != "END" and snp["CHROM"] != probe["chrom"]):
                sys.stderr.write("3 " + '\t'.join(probe.values()) + '\t' + "  ---  " + '\t'.join(snp.values()) + '\n')       
            if(snp == "END" or probe == "END" or (probe["chrom"] == snp["CHROM"])):
                chrom = True
            
    while(probe != "END"):
        probe_writer.writerow({"chrom":probe["chrom"],"probe_id":probe["probe_id"],"probeset_id":probe["probeset_id"],"PWK_SNP":0,"WSB_SNP":0})
        new_probe = next(gff_reader, "END")
        if(new_probe != "END"):      
            if((new_probe["probeset_id"] != probe["probeset_id"]) or new_probe == "END"):
                probeset_writer.writerow({"chrom":probe["chrom"], "probeset_id":probe["probeset_id"],"PWK_SNP":probeset_pwk_count, "WSB_SNP":probeset_wsb_count})
                probeset_pwk_count = 0
                probeset_wsb_count = 0
                
        probe = new_probe
