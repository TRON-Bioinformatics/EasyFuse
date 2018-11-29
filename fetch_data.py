#!/usr/bin/env python

import sys,shutil,os,stat,subprocess,datetime,csv

from argparse import ArgumentParser

from config import Config
from fusion_expression import FusionExpression
import misc.logger as logger

class ResultParser(object):
    
    def __init__(self,cfg,working_dir,prio_cutoff,stringency):
        self.cfg = cfg
        self.working_dir = working_dir
        self.prio_cutoff = prio_cutoff
        self.stringency = stringency
        self.fetchdata_dir = os.path.join(self.working_dir,"fetchdata_" + str(self.prio_cutoff) + "tool")
        self.sample_id = self.working_dir.rstrip("/").rstrip("/scratch").rsplit("/",1)[1]
        self.sub_dict = {"orf":"ORF","GCN1L1":"GCN1","PTPLB":"HACD2"}
        self.chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"]

        logfile = os.path.join(self.working_dir,"fusion.log")
        self.logger = logger.Logger(logfile)
        self.logger.setLevel(logger.DEBUG)
        self.logger.debug("Starting Result Parser")

        self.create_folder(self.fetchdata_dir)


    def create_folder(self,folder):
        '''This function creates a specified folder, if not already existing and grants the respective permission.'''
        if not os.path.exists(folder):
            self.logger.debug("Creating folder " + folder)
            os.makedirs(folder)
            self.grant_permissions(folder, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
        else:
            self.logger.debug("Folder " + folder + " already exists")
        
    def grant_permissions(self,path, mode):
        '''This function grants specific permissions for a given folder or file.'''
        for root, dirs, files in os.walk(path, topdown=False):
            for dir in dirs:
                os.chmod(dir, mode)
            for file in files:
                os.chmod(file, mode)

    def generate_matrix(self,res):
        keys = {}
        for tool in res:
            if len(res[tool]) == 0:
                keys[tool] = []

            for key in res[tool]:
                try:
                    keys[tool].append(key)
                except:
                    keys[tool] = [key]

        total_fusion_genes = []
        for tool in keys:
            total_fusion_genes += keys[tool]
        fusion_genes = list(set(total_fusion_genes))

        matrix = {}
        for key in fusion_genes:
            key_split = key.split("_")

            t = []
            for tool in res:
                found = False
                for ele in keys[tool]:
                    ele_split = ele.split("_")
                    if key_split[0] in ele_split and key_split[2] in ele_split:
                        found = True

                t.append(found)

            sum_tools = sum(t)
            if sum_tools >= self.prio_cutoff:
                matrix[key] = t
        return matrix

    def generate_results(self,tool_state_path):
        outf = open(tool_state_path,"a")
        tools = self.cfg.get('general','tools').split(",")
        res = {}
        err = {}

        self.logger.info("Processing " + self.sample_id)
        output_folder_path = os.path.join(self.fetchdata_dir,"tool_res")
        self.create_folder(output_folder_path)
        sample_string = self.sample_id
        for tool in tools:
            res[tool],err[tool] = self.get_tool_results(output_folder_path,tool)

            if err[tool]:
                sample_string += ";0"
            else:
                sample_string += ";1"
        sample_string.rstrip(";")
        outf.write(sample_string+"\n")
        outf.close()
        return (res,err)

    def get_tool_results(self,output_folder_path,tool):
        res = {}
        err = False
        self.logger.info("Parsing results for " + tool)
        if tool == "fusioncatcher":
            try:
                res = self.get_fusioncatcher_results()
            except:
                self.logger.error(tool + " results incomplete. Please check!")
                err = True
                return (res,err)
        elif tool == "starfusion":
            try:
                res = self.get_starfusion_results()
            except:
                self.logger.error(tool + " results incomplete. Please check!")
                err = True
                return (res,err)
        elif tool == "soapfuse":
            try:
                res = self.get_soapfuse_results()
            except:
                self.logger.error(tool + " results incomplete. Please check!")
                err = True
                return (res,err)
        elif tool == "mapsplice":
            try:
                res = self.get_mapsplice_results()
            except:
                self.logger.error(tool + " results incomplete. Please check!")
                err = True
                return (res,err)
        elif tool == "jaffa":
            try:
                res = self.get_jaffa_results()
            except:
                self.logger.error(tool + " results incomplete. Please check!")
                err = True
                return (res,err)
        elif tool == "infusion":
            try:
                res = self.get_infusion_results()
            except:
                self.logger.error(tool + " results incomplete. Please check!")
                err = True
                return (res,err)
        tool_res_file = os.path.join(output_folder_path,tool + "_res.csv")
        tool_outf = open(tool_res_file,"w")
        tool_outf.write("fgid;fusion_gene;breakpoint1;breakpoint2;junc_reads;span_reads;sample_id;tool\n")
        for key in res:
            tool_outf.write(key + ";" + ";".join(res[key]) + "\n")
        tool_outf.close()
        return (res,err)

    def get_fusioncatcher_results(self):
        summary_file = os.path.join(self.working_dir,"fusion","fusioncatcher","summary_candidate_fusions.txt")
        final_list_cand_fus_genes = os.path.join(self.working_dir,"fusion","fusioncatcher","final-list_candidate-fusion-genes.txt")
        reciprocal_fusions = []
        with open(summary_file) as summary:
            start = False
            for line in summary:
                if line.strip().startswith("*"):
                    gene_line = line.rstrip().replace("(",",(").replace(" ","").replace("*","").replace("--","_")

                    if "reciprocal" in line:
                        reciprocal_fusions.append((gene_line.split(",")[0]).upper())
        fusion_map = {}
        with open(final_list_cand_fus_genes) as final:
            next(final)
            for i,line in enumerate(final):
                elements = line.rstrip().split("\t")

                fusion_gene = (elements[0] + "_" + elements[1]).upper()

                if fusion_gene in reciprocal_fusions:
                    fusion_gene = (elements[1] + "_" + elements[0]).upper()
                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = elements[4]
                span_reads_num = elements[5]

                up_gene_bp = "chr" + elements[8]
                dn_gene_bp = "chr" + elements[9]
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp

                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"fusioncatcher"]
        return fusion_map

    def get_starfusion_results(self):
        full_tab = os.path.join(self.working_dir,"fusion","starfusion","star-fusion.fusion_candidates.final")
        fusion_map = {}
        with open(full_tab) as tab:
            next(tab)
            for i,line in enumerate(tab):
                elements = line.rstrip().split("\t")
                fusion_gene = elements[0].replace("--","_").upper()

                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = elements[1]
                span_reads_num = elements[2]
                up_gene_bp = elements[5]
                dn_gene_bp = elements[7]
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp

                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"starfusion"]
        return fusion_map


    def get_soapfuse_results(self):
        summary_file = ""
        folder_to_scan = os.path.join(self.working_dir,"fusion","soapfuse","final_fusion_genes")
        for file in os.listdir(folder_to_scan):
            folder_path = os.path.join(folder_to_scan,file)
            if os.path.isdir(folder_path):
                for res in os.listdir(folder_path):
                    if res.endswith(".final.Fusion.specific.for.genes"):
                        summary_file = os.path.join(folder_path,res)
        if not summary_file:
            summary_file = os.path.join(self.working_dir,"fusion","soapfuse","final_fusion_genes",self.sample_id,self.sample_id+".final.Fusion.specific.for.genes")
        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for i,line in enumerate(summary):
                elements = line.rstrip().split("\t")
                fusion_gene = (elements[0] + "_" + elements[5]).upper()

                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = elements[11]
                span_reads_num = elements[10]
                up_gene_bp = elements[1] + ":" + elements[3] + ":" + elements[2]
                dn_gene_bp = elements[6] + ":" + elements[8] + ":" + elements[7]
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp
                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"soapfuse"]
        return fusion_map

    def get_mapsplice_results(self):
        summary_file = os.path.join(self.working_dir,"fusion","mapsplice","fusions_well_annotated.txt")

        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for i,line in enumerate(summary):
                elements = line.rstrip().split("\t")

                fusion_gene = (elements[60] + elements[61]).replace(",","_").rstrip("_").upper()

                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = elements[4]
                strands = elements[5]
                strand_up = strands[0]
                strand_dn = strands[1]
                span_reads_num = elements[19]
                chrs = elements[0].split("~")
                up_gene_bp = "chr" + chrs[0] + ":" + elements[1] + ":" + strand_up
                dn_gene_bp = "chr" + chrs[1] + ":" + elements[2] + ":" + strand_dn
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp

                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"mapsplice"]
        return fusion_map

    def get_jaffa_results(self):
        summary_file = os.path.join(self.working_dir,"fusion","jaffa","jaffa_results.csv")
        if not os.path.exists(summary_file):
            summary_file = os.path.join(self.working_dir,"fusion","jaffa","results_JAFFA.csv")
        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for i,line in enumerate(summary):
                elements = line.rstrip().split(",")
                elements = [ele.strip("\"") for ele in elements]
                fusion_gene = elements[1].replace(":","_").upper()

                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = elements[9]
                span_reads_num = elements[10]
                up_gene_bp = elements[2] + ":" + elements[3] + ":" + elements[4]
                dn_gene_bp = elements[5] + ":" + elements[6] + ":" + elements[7]
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp

                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"jaffa"]
        return fusion_map

    def get_infusion_results(self):
        summary_file = os.path.join(self.working_dir,"fusion","infusion","fusions.detailed.txt")

        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for i,line in enumerate(summary):
                elements = line.rstrip().split("\t")

                fusion_gene = (elements[21] + "_" + elements[27]).upper().replace(";","-")

                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = elements[7]
                span_reads_num = elements[8]

                # To-Do: correct none-type strand ("." -> "+"/"-")
                up_gene_bp = "chr" + elements[1] + ":" + elements[2] + ":" + elements[23]
                dn_gene_bp = "chr" + elements[4] + ":" + elements[5] + ":" + elements[29]
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp

                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"infusion"]
        return fusion_map

    def get_chimpipe_results(self):
        summary_file = os.path.join(self.working_dir,"fusion","chimpipe","blablabla.txt")
        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for i,line in enumerate(summary):
                elements = line.rstrip().split("\t")
                fusion_gene = ""
                for key in self.sub_dict:
                    fusion_gene = fusion_gene.replace(key,self.sub_dict[key])

                junc_reads_num = 0
                span_reads_num = 0
                up_gene_bp = ""
                dn_gene_bp = ""
                up_gene_bp_chr = up_gene_bp.split(":")[0]
                dn_gene_bp_chr = dn_gene_bp.split(":")[0]
                
                if up_gene_bp_chr not in self.chr_list or dn_gene_bp_chr not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = fg_split[0] + "_" + up_gene_bp + "_" + fg_split[1] + "_" + dn_gene_bp

                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"chimpipe"]
        return fusion_map

    def csv_to_fasta(self,fusions_table,fasta):
        outf = open(fasta,"w")
        with open(fusions_table) as csvfile:
            inf = csv.reader(csvfile, delimiter=';')
            header = inf.next()

            col = {}
            for i,colname in enumerate(header):
                col[colname] = i
            for c,row in enumerate(inf):
                ftid = row[col["FTID"]]
#                fgid = row[col["FGID"]]
#                context_seq = row[col["context_sequence"]]
#                md5_hash = row[col["context_sequence_id"]]
#                id = ftid + "_fg_" + str(c)
#                context_id = fgid + "_" + md5_hash
                
#                breakpoint1 = row[col["Breakpoint1"]]
#                breakpoint2 = row[col["Breakpoint2"]]
                transcript_sequence = row[col["transcript_sequence"]].rstrip("NA").upper()
                wt_transcript1_sequence = row[col["original_transcript1"]].rstrip("NA").upper()
                wt_transcript2_sequence = row[col["original_transcript2"]].rstrip("NA").upper()
#                exbnd1 = row[col["exon_boundary1"]]
#                exbnd2 = row[col["exon_boundary2"]]
#                rel_bp_1 = row[col["relative_breakpoint1"]]
#                rel_bp_2 = row[col["relative_breakpoint2"]]

                outf.write(">" + ftid + "_ft" + "\n")
                for i in range(0,len(transcript_sequence),60):
                    outf.write(transcript_sequence[i:i+60]+"\n")

                outf.write(">" + ftid + "_wt1" + "\n")
                for i in range(0,len(wt_transcript1_sequence),60):
                    outf.write(wt_transcript1_sequence[i:i+60]+"\n")

                outf.write(">" + ftid + "_wt2" + "\n")
                for i in range(0,len(wt_transcript2_sequence),60):
                    outf.write(wt_transcript2_sequence[i:i+60]+"\n")
        outf.close()

    def get_gene_mapping(self,mapping_table):
        gene_map = {}
        with open(mapping_table) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(",")
                id = elements[0]
                ensembl_id = elements[1]
                try:
                    gene_map[id].append(ensembl_id)
                except:
                    gene_map[id] = [ensembl_id]
        return gene_map

    def add_chromosomal_dist(self,detected_fusions):
        chrom_state = {}
        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                id = elements[0] + "_" + elements[6] + "_" + elements[7]
                br1_split = elements[6].strip("chr").split(":")
                br2_split = elements[7].strip("chr").split(":")
                
                chr_dist = "NA"
                chr_dist_abs = ""
                chr_rel = "Interchromosomal"
                if br1_split[0] == br2_split[0]:
                    chr_dist = int(br1_split[1]) - int(br2_split[1])
                    if chr_dist < 0:
                        self.logger.debug("ALARM: Negative value detected!")
                    chr_dist_abs = str(abs(chr_dist))
                    chr_rel = "Intrachromosomal"

                chrom_state[id] = [chr_dist_abs,chr_rel]
        tools = self.cfg.get('general','tools')
        outf = open(detected_fusions+".tmp","w")
        outf.write("fusion_gene;context_id;ensembl_id_1;ensembl_id_2;junc_reads;span_reads;breakpoint1;breakpoint2;chrom_dist;chrom_rel;sequence;sample;known;tool;"+";".join(tools.split(","))+";gene1_rpkm;gene2_rpkm;fusion_gene_counts_kallisto;fusion_gene_tpm_kallisto;fusion_gene_counts_salmon;fusion_gene_tpm_salmon;wt_gene1_counts_kallisto;wt_gene1_tpm_kallisto;wt_gene1_counts_salmon;wt_gene1_tpm_salmon;wt_gene2_counts_kallisto;wt_gene2_tpm_kallisto;wt_gene2_counts_salmon;wt_gene2_tpm_salmon;frequency;exon_boundary_1;exon_boundary_2;fusion_transcript;unique_id\n")

        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                print elements
                id = elements[0] + "_" + elements[6] + "_" + elements[7]
                dist,rel = chrom_state[id]
                outf.write(";".join(elements[0:8])+";"+dist+";"+rel+";"+";".join(elements[8:])+"\n")
        outf.close()
        shutil.move(detected_fusions+".tmp",detected_fusions)

    def add_tool_state(self,detected_fusions):
        tool_state_map = {}

        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                print elements
                id = elements[0]
                tool = elements[13]

                try:
                    tool_state_map[id].append(tool)
                except:
                    tool_state_map[id] = [tool]
        tools = self.cfg.get('general','tools')
        outf = open(detected_fusions + ".state","w")
        outf.write("fusion_gene;context_id;ensembl_id_1;ensembl_id_2;junc_reads;span_reads;breakpoint1;breakpoint2;chrom_dist;chrom_rel;sequence;sample;known;tool;"+";".join(tools.split(","))+";gene1_rpkm;gene2_rpkm;fusion_gene_counts_kallisto;fusion_gene_tpm_kallisto;fusion_gene_counts_salmon;fusion_gene_tpm_salmon;wt_gene1_counts_kallisto;wt_gene1_tpm_kallisto;wt_gene1_counts_salmon;wt_gene1_tpm_salmon;wt_gene2_counts_kallisto;wt_gene2_tpm_kallisto;wt_gene2_counts_salmon;wt_gene2_tpm_salmon;frequency;exon_boundary_1;exon_boundary_2;fusion_transcript;unique_id\n")

        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                id = elements[0]
                tool_string = ""
                for tool in tools.split(","):
                    if tool in tool_state_map[id]:
                        tool_string += ";1"
                    else:
                        tool_string += ";0"
                outf.write(";".join(elements[0:14])+tool_string+";".join(elements[14:])+"\n")
        outf.close()        

    def add_gene_expression(self,detected_fusions):
        sample_expr_map = {}

        trsp6_in_rpkm = os.path.join(self.working_dir,"expression","gene_RPKM.csv")
        trsp6_in_counts = os.path.join(self.working_dir,"expression","gene_counts.csv")
        shutil.copy(trsp6_in_rpkm,self.fetchdata_dir)
        shutil.copy(trsp6_in_counts,self.fetchdata_dir)

        trsp6_res_map = {}
        try:
            with open(trsp6_in_rpkm) as res:
                for line in res:
                    elements = line.rstrip().split("\t")
                    id = elements[0]
                    rpkm_val = elements[1]
                    trsp6_res_map[id] = rpkm_val
            sample_expr_map = trsp6_res_map
        except:
            sample_expr_map = ""

        tools = self.cfg.get('general','tools')
        outf = open(detected_fusions + ".expression","w")
        outf.write("fusion_gene;context_id;ensembl_id_1;ensembl_id_2;junc_reads;span_reads;breakpoint1;breakpoint2;chrom_dist;chrom_rel;sequence;sample;known;tool;"+";".join(tools.split(","))+";gene1_rpkm;gene2_rpkm;fusion_gene_counts_kallisto;fusion_gene_tpm_kallisto;fusion_gene_counts_salmon;fusion_gene_tpm_salmon;wt_gene1_counts_kallisto;wt_gene1_tpm_kallisto;wt_gene1_counts_salmon;wt_gene1_tpm_salmon;wt_gene2_counts_kallisto;wt_gene2_tpm_kallisto;wt_gene2_counts_salmon;wt_gene2_tpm_salmon;frequency;exon_boundary_1;exon_boundary_2;fusion_transcript;unique_id\n")

        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                genes = elements[0].split("_")
                gene_1 = genes[0]
                gene_2 = genes[1]
                sample_id = elements[11]
                exprs = sample_expr_map
                try:
                    gene_1_expr = exprs[gene_1]
                    gene_2_expr = exprs[gene_2]
                    outf.write(";".join(elements)+";"+gene_1_expr+";"+gene_2_expr+"\n")
                except:
                    outf.write(";".join(elements)+";na;na\n")
        outf.close()

    def add_fusion_expression(self,context_seqs,quant_file):
        transcripts_fasta = self.cfg.get('references','ensembl_genes_hg38_fasta')

        sample_expr_map = {}
        sequences_fasta = os.path.join(self.fetchdata_dir,"fusion_genes.fa")
        all_fasta = os.path.join(self.fetchdata_dir,"all.fa")
        self.csv_to_fasta(context_seqs,sequences_fasta)
        
        with open(all_fasta, 'w') as outfile:
            with open(transcripts_fasta) as infile:
                for line in infile:
                    outfile.write(line)
            with open(sequences_fasta) as infile:
                for line in infile:
                    outfile.write(line)

        kallisto_index = os.path.join(self.fetchdata_dir,"kallisto.idx")
        salmon_index = os.path.join(self.fetchdata_dir,"salmon_idx")
        
        self.logger.debug("Parsing context sequence table")
        
        context_seqs_map = {}
        with open(context_seqs) as con:
            inf = csv.reader(con, delimiter=";")
            header = inf.next()
            col = {}
            for i,colname in enumerate(header):
                col[colname] = i
            for c,row in enumerate(inf):
                ftid = row[col["FTID"]]
                fgid = row[col["FGID"]]
                context_sequence = row[col["context_sequence"]]
                breakpoint1 = row[col["Breakpoint1"]]
                breakpoint2 = row[col["Breakpoint2"]]
                fusion_transcript = row[col["transcript"]].rstrip("NA").upper()
                orig_trans1 = row[col["original_transcript1"]]
                orig_trans2 = row[col["original_transcript2"]]
                exbnd1 = row[col["exon_boundary1"]]
                exbnd2 = row[col["exon_boundary2"]]

                context_seqs_map[ftid] = [fgid,context_sequence,fusion_transcript,exbnd1,exbnd2]
        self.logger.debug(str(len(context_seqs_map.keys())) + " unique FTIDs in table")

        fusexpr = FusionExpression(self.cfg,all_fasta,kallisto_index,salmon_index)
        calc_expr = False
        if not os.path.exists(kallisto_index) or not os.path.exists(salmon_index):
            calc_expr = True

        if calc_expr:
            fusexpr.write_index_kallisto()
            fusexpr.write_index_salmon()

        left = os.path.join(self.working_dir,"skewer","out_file-trimmed-pair1.fastq.gz")
        right = os.path.join(self.working_dir,"skewer","out_file-trimmed-pair2.fastq.gz")
        if not os.path.exists(left) and not os.path.exists(right):
            left = os.path.join(self.working_dir,"out_file-trimmed-pair1.fastq.gz")
            right = os.path.join(self.working_dir,"out_file-trimmed-pair2.fastq.gz")

        kallisto_output = os.path.join(self.fetchdata_dir,"kallisto")
        salmon_output = os.path.join(self.fetchdata_dir,"salmon")

        self.create_folder(kallisto_output)
        self.create_folder(salmon_output)

        kallisto_res = os.path.join(kallisto_output,"abundance.tsv")
        salmon_res = os.path.join(salmon_output,"quant.sf")

        if not os.path.exists(kallisto_res):
            self.logger.info("Quantifying fusion genes using kallisto (" + self.sample_id + ")")
            fusexpr.quantification_kallisto(left,right,kallisto_output)
        if not os.path.exists(salmon_res):
            self.logger.info("Quantifying fusion genes using salmon (" + self.sample_id + ")")
            fusexpr.quantification_salmon(left,right,salmon_output)
        kallisto_res_map = {}
        with open(kallisto_res) as res:
            next(res)
            for line in res:
                elements = line.rstrip().split("\t")
                ftid = elements[0]
                length = elements[1]
                eff_length = elements[2]
                est_counts = elements[3]
                tpm_val = elements[4]

                kallisto_res_map[ftid] = [length,eff_length,est_counts,tpm_val]

        salmon_res_map = {}
        with open(salmon_res) as res:
            next(res)
            for line in res:
                elements = line.rstrip().split("\t")
                ftid = elements[0]
                length = elements[1]
                eff_length = elements[2]
                est_counts = elements[4]
                tpm_val = elements[3]

                salmon_res_map[ftid] = [length,eff_length,est_counts,tpm_val]

        sample_expr_map = [kallisto_res_map,salmon_res_map]

        tools = self.cfg.get('general','tools')
        outf = open(quant_file,"w")

        outf.write("FTID;FGID;Context_Sequence;Fusion_Transcript;Exon_Boundary1;Exon_Boundary2;fg_kallisto_counts;fg_kallisto_tpm;fg_salmon_counts;fg_salmon_tpm;wt1_kallisto_counts;wt1_kallisto_tpm;wt1_salmon_counts;wt1_salmon_tpm;wt2_kallisto_counts;wt2_kallisto_tpm;wt2_salmon_counts;wt2_salmon_tpm\n")

        self.logger.debug("Generating Quantification table")
        c = 0
        for ftid in context_seqs_map:
            elements = context_seqs_map[ftid]
            outf.write(ftid+";"+elements[0] + ";" + elements[1] + ";" + elements[2] + ";" + elements[3] + ";" + elements[4])

            try:
                kallisto_res_fg = sample_expr_map[0][ftid + "_ft"]
            except:
                kallisto_res_fg = [0,0,0,0]
            try:
                kallisto_res_wt1 = sample_expr_map[0][ftid + "_wt1"]
            except:
                kallisto_res_wt1 = [0,0,0,0]
            try:
                kallisto_res_wt2 = sample_expr_map[0][ftid + "_wt2"]
            except:
                kallisto_res_wt2 = [0,0,0,0]

            try:
                salmon_res_fg = sample_expr_map[1][ftid + "_ft"]
            except:
                salmon_res_fg = [0,0,0,0]
            try:
                salmon_res_wt1 = sample_expr_map[1][ftid + "_wt1"]
            except:
                salmon_res_wt1 = [0,0,0,0]
            try:
                salmon_res_wt2 = sample_expr_map[1][ftid + "_wt2"]
            except:
                salmon_res_wt2 = [0,0,0,0]


            outf.write(";"+str(kallisto_res_fg[2])+";"+str(kallisto_res_fg[3])+";"+str(salmon_res_fg[2])+";"+str(salmon_res_fg[3])+";"+str(kallisto_res_wt1[2])+";"+str(kallisto_res_wt1[3])+";"+str(salmon_res_wt1[2])+";"+str(salmon_res_wt1[3])+";"+str(kallisto_res_wt2[2])+";"+str(kallisto_res_wt2[3])+";"+str(salmon_res_wt2[2])+";"+str(salmon_res_wt2[3])+"\n")

            c += 1

        outf.close()
        self.logger.debug("Wrote " + str(c) + " lines for " + str(len(context_seqs_map.keys())) + " unique FTIDs.")
        

    def add_frequency(self,input,detected_fusions):
        num_samples = len(input)
        sample_frequency_map = {}
        frequency_map = {}
        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                fusion_gene = elements[0]
                sample_id = elements[11]
                try:
                    frequency_map[fusion_gene].append(sample_id)
                except:
                    frequency_map[fusion_gene] = [sample_id]

        for gene in frequency_map.iterkeys():
            freq = float(len(list(set(frequency_map[gene]))))/num_samples
            frequency_map[gene] = freq            
        
        tools = self.cfg.get('general','tools')
        outf = open(detected_fusions+".freq","w")

        outf.write("fusion_gene;context_id;ensembl_id_1;ensembl_id_2;junc_reads;span_reads;breakpoint1;breakpoint2;chrom_dist;chrom_rel;sequence;sample;known;tool;"+";".join(tools.split(","))+";gene1_rpkm;gene2_rpkm;fusion_gene_counts_kallisto;fusion_gene_tpm_kallisto;fusion_gene_counts_salmon;fusion_gene_tpm_salmon;wt_gene1_counts_kallisto;wt_gene1_tpm_kallisto;wt_gene1_counts_salmon;wt_gene1_tpm_salmon;wt_gene2_counts_kallisto;wt_gene2_tpm_kallisto;wt_gene2_counts_salmon;wt_gene2_tpm_salmon;frequency;exon_boundary_1;exon_boundary_2;fusion_transcript;unique_id\n")

        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                fusion_gene = elements[0]
                
                outf.write(";".join(elements)+";"+str(frequency_map[fusion_gene])+"\n")
        outf.close()

    def add_exon_boundaries(self,context_seqs,detected_fusions):
        context_seqs_map = {}
        with open(context_seqs) as con:
            inf = csv.reader(con, delimiter=";")
            header = inf.next()

            col = {}
            for i,colname in enumerate(header):
                col[colname] = i
            for c,row in enumerate(con):

                context_id = row[col["context_sequence_id"]]
                sequence = row[col["context_sequence"]]
                breakpoint1 = row[col["Breakpoint1"]]
                breakpoint2 = row[col["Breakpoint2"]]
                fusion_gene = row[col["FGID"]]
                ex_bnd_1 = row[col["exon_boundary1"]]
                ex_bnd_2 = row[col["exon_boundary2"]]
                fusion_transcript = row[col["transcript"]].rstrip("NA").upper()

                context_split = context_id.split("_")

                unique_id = ""
                if len(context_split) > 4:
                    unique_id = fusion_gene + "_" + breakpoint1 + "_" + breakpoint2 + "_" + context_split[1] + "_" + context_split[4]
                else:
                    unique_id = fusion_gene + "_" + breakpoint1 + "_" + breakpoint2 + "_" + context_split[1] + "_" + context_split[3]
                id = context_id + "_" + breakpoint1 + "_" + breakpoint2
                context_seqs_map[id] = [ex_bnd_1,ex_bnd_2,fusion_transcript,unique_id]

        tools = self.cfg.get('general','tools')
        outf = open(detected_fusions+".exbnd.csv","w")

        outf.write("fusion_gene;context_id;ensembl_id_1;ensembl_id_2;junc_reads;span_reads;breakpoint1;breakpoint2;chrom_dist;chrom_rel;sequence;sample;known;tool;"+";".join(tools.split(","))+";gene1_rpkm;gene2_rpkm;fusion_gene_counts_kallisto;fusion_gene_tpm_kallisto;fusion_gene_counts_salmon;fusion_gene_tpm_salmon;wt_gene1_counts_kallisto;wt_gene1_tpm_kallisto;wt_gene1_counts_salmon;wt_gene1_tpm_salmon;wt_gene2_counts_kallisto;wt_gene2_tpm_kallisto;wt_gene2_counts_salmon;wt_gene2_tpm_salmon;frequency;exon_boundary_1;exon_boundary_2;fusion_transcript;unique_id\n")

        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                context_id = elements[1]
                breakpoint1 = elements[6]
                breakpoint2 = elements[7]
                
                id = context_id + "_" + breakpoint1 + "_" + breakpoint2
                
                outf.write(";".join(elements)+";"+";".join(context_seqs_map[id])+"\n")
        outf.close()


    def filter_final_table(self,detected_fusions):
        filtered_map = {}
        with open(detected_fusions) as csv:
            next(csv)
            for line in csv:
                elements = line.rstrip().split(";")
                fusion_gene = elements[0]
                breakpoint1 = elements[6]
                breakpoint2 = elements[7]
                sequence = elements[10]
                sample = elements[11]
                tool = elements[13]

                id = fusion_gene + "_" + breakpoint1 + "_" + breakpoint2 + "_" + sequence
                if id in filtered_map:

                    if not sample in filtered_map[id][11]:
                        filtered_map[id][11] += ":" + sample
                    if not tool in filtered_map[id][13]:
                        filtered_map[id][13] += ":" + tool

                else:
                    filtered_map[id] = elements

        self.logger.info(str(len(filtered_map.keys())) + " unique fusions remaining")
        tools = self.cfg.get('general','tools')
        outf = open(detected_fusions+".filtered.csv","w")

        outf.write("fusion_gene;context_id;ensembl_id_1;ensembl_id_2;junc_reads;span_reads;breakpoint1;breakpoint2;chrom_dist;chrom_rel;sequence;sample;known;tool;"+";".join(tools.split(","))+";gene1_rpkm;gene2_rpkm;fusion_gene_counts_kallisto;fusion_gene_tpm_kallisto;fusion_gene_counts_salmon;fusion_gene_tpm_salmon;wt_gene1_counts_kallisto;wt_gene1_tpm_kallisto;wt_gene1_counts_salmon;wt_gene1_tpm_salmon;wt_gene2_counts_kallisto;wt_gene2_tpm_kallisto;wt_gene2_counts_salmon;wt_gene2_tpm_salmon;frequency;exon_boundary_1;exon_boundary_2;fusion_transcript;unique_id\n")

        for key in filtered_map:
            elements = filtered_map[key]
            outf.write(";".join(elements)+"\n")
        outf.close()

    def run(self):
        tools = self.cfg.get('general','tools').split(",")
        tool_state_path = os.path.join(self.fetchdata_dir,"tool_state.csv")
        tool_state = open(tool_state_path,"w")
        tool_state.write("Sample ID")
        for tool in tools:
            tool_state.write(","+tool)
        tool_state.write("\n")
        tool_state.close()
        incomplete = False

        detected_file = os.path.join(self.fetchdata_dir,"Detected_Fusions.csv")
        quant_file = os.path.join(self.fetchdata_dir,"Quantification.csv")
        outf = open(detected_file,"w")

        outf.write("FGID;Fusion_Gene;Breakpoint1;Breakpoint2;Junction_Reads;Spanning_Reads;Sample;Tool;Chromosomal_Distance;Chrom_Rel;Frequency\n")
        c = 0
        self.logger.debug("Generating Detected Fusions table")

        res,err = self.generate_results(tool_state_path)

        if sum(err.values()) != 0:
            incomplete = True
            if not self.stringency:
                sys.exit(1)

        matrix = self.generate_matrix(res)
        for key in matrix:
            for i,tool in enumerate(res):
                if matrix[key][i]:
                    for f_key in res[tool]:
                        if f_key == key:
                            c += 1
                            outf.write(key + ";" + ";".join(res[tool][f_key]) + "\n")

        outf.close()

        self.logger.debug("Wrote " + str(c) + " lines.")
        self.logger.info("Detected fusion genes table created.")

#    print "Adding chromosomal distance to gene fusions"
#    res_parse.add_chromosomal_dist(outfile)

#    print "Adding tool state to gene fusions"
#    res_parse.add_tool_state(outfile)

        if incomplete and not self.stringency:
            self.logger.error("Results incomplete. Aborting. Please make sure that all tools have run completely on dataset.")
            sys.exit(0)

        context_seqs = os.path.join(self.fetchdata_dir,"Context_seqs.csv")

        expr_file = os.path.join(self.working_dir,"expression","transcript_BPKM.csv")
#    cmd = "%s %s %s %s" % (cfg.get('commands','fetch_context_cmd'),detected_file,fetchdata_dir,organism)
        cmd = "%s -i %s -o %s -f %s -e %s --context_seq_len %s" % (self.cfg.get('commands','fetch_context_cmd'),detected_file,context_seqs,self.cfg.get('bowtie_indexes','bowtie_index_hg38_fasta'),self.cfg.get('references','ensembl_genes_hg38_csv'), self.cfg.get('general', 'context_seq_len'))

        if os.path.exists(detected_file) and not os.path.exists(context_seqs):

            self.logger.info("Generating context sequences out of detected fusions.")
            self.logger.debug("Command was: " + cmd)
            proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            (stdoutdata,stderrdata) = proc.communicate()
            if proc.returncode != 0:
                self.logger.error(stderrdata)

#        seq2hla_classI_res = os.path.join(self.working_dir,"seq2hla","run-ClassI-class.HLAgenotype4digits")
#        seq2hla_classII_res = os.path.join(self.working_dir,"seq2hla","run-ClassII.HLAgenotype4digits")

#        cmd = "%s %s %s %s %s" % (self.cfg.get('commands','mhc_prediction_cmd'),context_seqs,seq2hla_classI_res,seq2hla_classII_res,self.fetchdata_dir)

#        classI_pred = os.path.join(self.fetchdata_dir,"hla_classI.csv")
#        classII_pred = os.path.join(self.fetchdata_dir,"hla_classII.csv")

#        if not os.path.exists(classI_pred) and not os.path.exists(classII_pred) and os.path.exists(context_seqs):
#            self.logger.info("Predicting MHC eptitopes for peptides.")
#            self.logger.debug("Command was: " + cmd)
#            proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#            (stdoutdata,stderrdata) = proc.communicate()
#            if proc.returncode != 0:
#                self.logger.error(stderrdata)

        file_list = []

        file_r1 = os.path.join(self.working_dir,"skewer","out_file-trimmed-pair1.fastq.gz")
        file_r2 = os.path.join(self.working_dir,"skewer","out_file-trimmed-pair2.fastq.gz")
        if not os.path.exists(file_r1) and not os.path.exists(file_r2):
            file_r1 = os.path.join(self.working_dir,"out_file-trimmed-pair1.fastq.gz")
            file_r2 = os.path.join(self.working_dir,"out_file-trimmed-pair2.fastq.gz")
        file_list.append(file_r1)
        file_list.append(file_r2)


#    requantification_dir = os.path.join(working_dir,"requantification")
        cmd = "%s -i %s -f %s -o %s -p prod" % (self.cfg.get('commands','custom_trans_cmd')," ".join(file_list),context_seqs,self.fetchdata_dir)
    
        classification_file = os.path.join(self.fetchdata_dir,"Classification.csv")
#    classification_file_extend = os.path.join(working_dir,"Classification_extend.csv")

        if not os.path.exists(classification_file) and os.path.exists(context_seqs):
            self.logger.info("Creating custom transcriptome and quantifying samples.")
            self.logger.debug("Command was: " + cmd)
            proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            (stdoutdata,stderrdata) = proc.communicate()
            if proc.returncode != 0:
                self.logger.error(stderrdata)

        if os.path.exists(context_seqs):
            self.logger.info("Adding fusion expression to table")
            self.add_fusion_expression(context_seqs,quant_file)

####
#    res_parse.add_gene_expression(args.input,args.output,quant_file)

#    print "Adding fusion frequency to table"
#    res_parse.add_frequency(args.input,outfile+".state.expression.fusion")
    
#    print "Adding exon boundaries to table"
#    res_parse.add_exon_boundaries(context_seqs,outfile+".state.expression.fusion.freq")

#    print "Filtering final table"
#    res_parse.filter_final_table(outfile+".state.expression.fusion.freq.exbnd.csv")
    
#    primerblat_results = os.path.join(args.output,"fusion_primerblat.csv")
#    primerblat_err = os.path.join(args.output,"fusion_primerblat.err")

#    cmd = "python /kitty/code/primerblat/primerblat.py --primer3_target=80,40 --primer3_size=100-150  --column_values=1 --inputdelim=";"  --debug_p3 --delim=";" --dec="," %s /kitty/data/human/hg38.2bit > %s 2> %s" % (context_seqs,primerblat_results,primerblat_err)
#    log.info("Generating primers out of context sequences.")
#    log.debug("Command was: " + cmd)
#    proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#    (stdoutdata,stderrdata) = proc.communicate()
#    if proc.returncode != 0:
#        log.error(stderrdata)        
                

def main():
    parser = ArgumentParser(description='Extracts information on fusion genes')
    parser.add_argument('-i', '--input', dest='input', help='Specify sample working directory.',required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.',default="")
#    parser.add_argument('-g', '--gene-symbols', dest='gene_symbols', help='Specify hugo to enembl gene symbol mapping file',default='/kitty/data/human/ensembl/GRCh38.86/Homo_sapiens.GRCh38.86.HGNC2ENS.csv')
    parser.add_argument('-s', '--species', dest='species', choices=['human','mouse'], help='Specify species to be processed',default='human')
    parser.add_argument('-v', '--stringent', dest='stringency', action='store_true', help='Special case where not all tools have results')

    args = parser.parse_args()

    if args.config == "":
        cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.ini'))
    else:
        cfg = Config(args.config)
    working_dir = args.input
    organism = args.species

    if not os.path.exists(working_dir) or not os.path.isdir(working_dir):
        sys.exit(1)
    
    script_call = "python " + os.path.realpath(__file__) + " " + " ".join(sys.argv[1:])

    outf = open(os.path.join(working_dir,"fetch_data.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

    
    print "Running with 1 tool cutoff"
    try:
        res_parse_1 = ResultParser(cfg,working_dir,1,args.stringency)
        res_parse_1.run()
    except:
        pass

    print "Running with 2 tool cutoff"
    try:
        res_parse_2 = ResultParser(cfg,working_dir,2,args.stringency)
        res_parse_2.run()
    except:
        pass

            
                
if __name__ == '__main__':
    main()
