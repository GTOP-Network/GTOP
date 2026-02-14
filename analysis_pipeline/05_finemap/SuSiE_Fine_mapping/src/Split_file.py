
import sys
import os
from collections import defaultdict


def split_qtl(tissue_name):
	path="/lustre/home/xdzou/2024-10-21-GTBMap/2025-06-15-joint-split-gt-QTL/output/QTL_mapping/all_joint/cis_QTL/text_format/"
	qtl_file = path+tissue_name+".cis_eQTL.all_pairs.add_tstat.txt"
#	snp_file = "/lustre/home/xdzou/data/GTEx_SNP_info.txt"
#	snp_info = {}
#	fh = open(snp_file,'r')
#	for line in fh.readlines():
#		line = line.strip()
#		w = line.split("\t")
#		snp_info[w[3]] = (w[4],w[5])
#	fh.close()

	gene_dict={}
	gene_dict_gcta=defaultdict(list)
	gene_dict_caviar=defaultdict(list)
	snp_dict={}
	fh = open(qtl_file,'r')
	for qtl_line in fh.readlines()[1:]:
		qtl_lines = qtl_line.rstrip().split("\t")
	#	gene_name = qtl_lines[1]
		gene_name = qtl_lines[0]
	#	snp_infos = qtl_lines[0]
		snp_id = qtl_lines[1]
		chr = qtl_lines[8]
		pos = qtl_lines[9]
	#	allel1 = snp_infos.split("_")[2]
	#	allel2 = snp_infos.split("_")[3]
	#	allele1 = snp_info[snp_id][0]
	#	allele2 = snp_info[snp_id][1]
		if gene_name not in gene_dict:
		#	gcta_id=chr+"_"+pos+"_"+allel1+"_"+allel2+"_b38"
			gcta_id=snp_id
			gene_dict_gcta[gene_name]=[gcta_id]
		#	gene_dict_caviar[gene_name]=[[gcta_id,qtl_lines[3]]]
			gene_dict_caviar[gene_name]=[[gcta_id,qtl_lines[-1]]]
			gene_dict[gene_name]=1
		elif gene_name in gene_dict:
		#	gcta_id=chr+"_"+pos+"_"+allel1+"_"+allel2+"_b38"
			gcta_id=snp_id
			gene_dict_gcta[gene_name].append(gcta_id)
		#	gene_dict_caviar[gene_name].append([gcta_id,qtl_lines[3]])
			gene_dict_caviar[gene_name].append([gcta_id,qtl_lines[-1]])
			gene_dict[gene_name]=gene_dict[gene_name]+1

	d = dict((k, tuple(v)) for k, v in gene_dict_gcta.iteritems())
	caviar = dict((k, tuple(v)) for k, v in gene_dict_caviar.iteritems())
	print(caviar)
	for current_gene in d:
		gene_directory="/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project/2022-10-05-heritability/input/split_gene/"+tissue_name+"/"+current_gene
		if not os.path.exists(gene_directory):
			os.makedirs(gene_directory)

		output_file_caviar=gene_directory+"/"+tissue_name+"_3aQTL_all_ttest.txt"
		output_file_caviar_file=open(output_file_caviar, "w")

		output_file_gcta=gene_directory+"/"+tissue_name+"_3aQTL_all.snp"
		output_file_gcta_file=open(output_file_gcta, "w")
		for qtl_id in caviar[current_gene]:
			output_file_caviar_file.writelines(qtl_id[0] + "\t" + qtl_id[1] + "\n")

		for qtl_id in d[current_gene]:
			output_file_gcta_file.writelines(qtl_id+"\n")


if __name__ == '__main__':
	tissue_name = sys.argv[1]
	split_qtl(tissue_name)
