import sys
import vcfpy
import json

def vcf_reader(path):
    print(f"Loading {path}") 
    
    reader = vcfpy.Reader.from_path(path)
    
    csq_header_str = reader.header.get_info_field_info("CSQ").description.split(":")[-1].strip()
    csq_headers = csq_header_str.split("|")
    
    max_csq = 0
    processedRecs = []
    for rec in reader:
        
        pr = {}
        
        #csq_values = rec.INFO["CSQ"][0].split("|")
        #csq = {csq_headers[x]:csq_values[x] for x in range(0,len(csq_headers))} 
        csq_entries = []
        
        if len(rec.INFO["CSQ"]) > max_csq: 
            max_csq = len(rec.INFO["CSQ"])
            print(max_csq) 
        for str_csq in rec.INFO["CSQ"]: 
            csq_values = str_csq.split("|")
            csq_entries.append({csq_headers[x]:csq_values[x] for x in range(0,len(csq_headers))})
        
        #Merge csqs
        alt_alleles = [csq["Allele"] for csq in csq_entries]
        if len(alt_alleles) == 0: 
            print("NO CSQ!!!!!!!!!!!!!!!!!!")
        
        csq_facts = {}
        
        #"Consequence": "upstream_gene_variant",
        #"IMPACT": "MODIFIER",
        #"SYMBOL": "FAM138F",
        #"Gene": "ENSG00000282591",
        #"Feature_type": "Transcript",
        #"Feature": "ENST00000631376.1",
        #"BIOTYPE": "lncRNA",
        #"Existing_variation": "rs868831437",
        
        
        for csq in csq_entries:
            allele = csq["Allele"]
            if allele not in csq_facts:
                csq_facts[allele] = {
                    "freq":[],
                    "features":{},
                    "variant":None,
                    "symbol":[],
                    "impact":[],
                }
            
            if csq["AF"]:
                csq_facts[allele]["freq"].append(csq["AF"])
            else:
                csq_facts[allele]["freq"].append("unknown")
                
            csq_facts[allele]["variant"] = csq["Existing_variation"] # Go with ID column - will be . or user
            strand = "forward strand"
            if csq["STRAND"] == "-1": 
                strand = "reverse strand"
            csq_facts[allele]["strand"] = strand
            
            if csq["Feature"]: 
                feature_id = csq["Feature"]
                if feature_id not in csq_facts[allele]["features"]:  
                     csq_facts[allele]["features"][feature_id] = {
                         "type":csq["Feature_type"],
                         "consequence":[],
                         "gene":csq["Gene"],
                         "impact":csq["IMPACT"],
                         "biotype":csq['BIOTYPE'],
                     }
                    
                cons = csq["Consequence"].split("&")
                for c in cons: 
                    if c not in csq_facts[allele]["features"][feature_id]["consequence"]:
                        csq_facts[allele]["features"][feature_id]["consequence"].append(c)
        
        #print(csq)
        #break 
        
        pr["ref"] = rec.REF
        pr["location"] = f"{rec.CHROM}:{rec.POS}"
        if len(alt_alleles) > 0:
            pr["alt"] = alt_alleles
        else:
            pr["alt"] = str(rec.ALT[0].value)
        pr["feature"] = csq_facts
        
        processedRecs.append(pr)
    
    dump_json = True
    dump_html = True
    
    
    if dump_json:
        with open("dump.json","w") as dump:
            json.dump(processedRecs[0:10],dump)
        
    if dump_html:
        with open("dump.html","w") as dump:
            dump.write("<HTML><head><title>vep test dump</title>")
            dump.write('<link rel="stylesheet" href="http://getskeleton.com/dist/css/normalize.css">')
            dump.write('<link rel="stylesheet" href="http://getskeleton.com/dist/css/skeleton.css">')
            dump.write("</head><body>")
            dump.write('<table style="width:90%;margin-left:auto;margin-right:auto"><tr>')
            dump.write("<th>Varient</th>")
            dump.write("<th>Ref</th>")
            dump.write("<th>Location</th>")
            dump.write("<th>Alt Allele</th>")
            dump.write("<th>Alt Freq</th>")
            dump.write("<th>Features</th>")
            dump.write("<th>Transcripts</th>")
            dump.write("<th>Consequences</th>")
            dump.write("</tr>\n")
            
            for rec in processedRecs[0:200]:
                
                if len(rec["feature"]) > 0: 
                    print(f"Feature length - {len(rec['feature'])}")
                    for alt_allele, facts in rec["feature"].items():
                        dump.write("<tr>")
                        dump.write(f"<td>{facts['variant']}</td>")
                        dump.write(f"<td>{rec['ref']}</td>")
                        dump.write(f"<td>{rec['location']}</td>")
                        dump.write(f"<td>{alt_allele}</td>")
                        dump.write(f"<td>{facts['freq'][0]}</td>")
                        dump.write(f"<td>")
                        for f_id,f in facts["features"].items():
                            if f['type'] != "Transcript":
                                dump.write(f"{f_id} - {f['type']} - {f['gene']} - {f['impact']} - {facts['strand']}<br>")
                        dump.write("</td>")
                        
                        dump.write(f"<td>")
                        
                        for f_id,f in facts["features"].items():
                            if f['type'] == "Transcript":
                                dump.write(f"{f_id} {f['biotype']}<br>")
                        dump.write("</td>")
                        
                        dump.write("<td>")
                        for f_id,f in facts["features"].items():
                            for con in f['consequence']:
                                dump.write(f"{con}<br>")
                        dump.write("</td>")
                        dump.write("</tr>\n")
                        
                else:
                    dump.write("<tr>")
                    dump.write(f"<td>unknown</td><td>{rec['ref']}</td>")
                    dump.write(f"<td>{rec['location']}</td>")
                    dump.write(f"<td>{rec['alt']}</td>")
                    dump.write(f"<td>unknown</td>")
                    dump.write(f"<td></td>")
                    dump.write(f"<td></td>")
                    dump.write(f"<td></td>")
                    dump.write("</tr>\n")
                    
                    
                
            
            dump.write("</table>")
            dump.write("</body></html>")
    
        
    
        
        #print(pr)
        #print("--------------")
        #for f in rec.INFO["CSQ"]:
        #    print(f)
        #print(reader.header)
        
        #break
        
    #Convert to VCF 
    #
    #writer = vcfpy.Writer.from_path("output.vcf",reader.header)
    #for rec in reader:
    #    writer.write_record(rec)
    
    print(f"Max CSQ entries per record { max_csq}")
    
    
        
        



if __name__ == "__main__": 
    vcf_reader(sys.argv[1])