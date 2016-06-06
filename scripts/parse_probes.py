import csv

with open("probes","rb") as probes, open("probesets","wb") as probesets, open("probeset_exon_transcript","wb") as clusters:
    fieldnames = ["chrom","start","end","probe_id","probeset_id","sequence"]
    probeset_writer = csv.DictWriter(probesets, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    probeset_writer.writeheader()
    
    fieldnames = ["probe_id","probeset_id","exon_cluster_id","transcript_cluster_id"]
    cluster_writer = csv.DictWriter(clusters, fieldnames, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    cluster_writer.writeheader()
    
    for line in probes:
        line = line.replace(";","")
        line = line.replace("\r","")
        line = line.replace("\n","")
        line = line.replace(" ","\t")
        line = line.replace("chr","") #or delete with sed beforehand
        line = line.replace("_random","") #or delete with sed beforehand
        line_parse_tab = line.split('\t')
        
        probeset_writer.writerow({"chrom":line_parse_tab[0],"start":line_parse_tab[1],"end":line_parse_tab[2],"probe_id":line_parse_tab[line_parse_tab.index("probe_id")+1],"probeset_id":line_parse_tab[line_parse_tab.index("probeset_id")+1], "sequence":line_parse_tab[line_parse_tab.index("probe")+1].strip('"')})
        
        cluster_writer.writerow({"probe_id":line_parse_tab[line_parse_tab.index("probe_id")+1],"probeset_id":line_parse_tab[line_parse_tab.index("probeset_id")+1],"exon_cluster_id":line_parse_tab[line_parse_tab.index("exon_cluster_id")+1],"transcript_cluster_id":line_parse_tab[line_parse_tab.index("transcript_cluster_id")+1]})