#!/usr/bin/bash

#==========#
# Get args #
#==========#

GTFfile=$1 # GTF (not gzipped)
leafcutterOutDir=$2 # Path to Leafcutter2 output directory
leafcutterOutPrefix=$3 # Prefix of Leafcutter2 output files
sampleLookupFile=$4
selectedSampleList=$5
leafcutterScriptDir=$6 # Path to leafcutter `script` foldr
exonFile=$7 # This gzipped TSV file contains annotation information for exons, with each row corresponding to a single exon entry. The file includes the following columns:chr: Chromosome name; start: 0-based genomic start position of the exon; end: 0-based exclusive genomic end position of the exon; strand: Strand of the exon on the chromosome, indicated by "+"; gene_id: Unique identifier for the gene associated with the exon (e.g., Ensembl gene ID: ENSG00000123456); gene_name: Common/gene symbol name of the associated gene (e.g., TP53, BRCA1); transcript_id: Unique identifier for the transcript that the exon belongs to (e.g., Ensembl transcript ID: ENST00000345678), linking the exon to its parent transcript.
geneinfoFile=$8 # This TSV file contains annotation information for genes, with each row corresponding to a single gene. The file includes the following columns:chr: gene_name: Common/gene symbol name of the associated gene (e.g., TP53, BRCA1);gene_type: biotype of gene; gene_id: Unique identifier for the gene associated with the exon (e.g., Ensembl gene ID: ENSG00000123456); Chromosome name; start: 1-based genomic start position of the exon; end: 1-based inclusive genomic end position of the exon; strand: Strand of the exon on the chromosome, indicated by "+".
outDir=$9 


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts


phenotypeScript=$scriptDir/leafcutter2.cluster_prepare_tensorqtl.py



#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outDir ]]; then
	mkdir -p $outDir
fi



#====================================================================#
# Map splicing clusters to genes and make phenotype for sQTL calling #
#====================================================================#

python $phenotypeScript \
    $wkpath/tissue_group/alllist \
    $exonFile \
    $GTFfile \
    ${leafcutterOutPrefix} \
    $sampleLookupFile \
    --leafcutter_dir $leafcutterScriptDir \
    --geneinfo $geneinfoFile \
    -o $outpath
fi


