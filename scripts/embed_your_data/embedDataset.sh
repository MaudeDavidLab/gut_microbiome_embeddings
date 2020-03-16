# You will be starting with a repseqs.fasta file

###############################################################
####### Step 1: Modify all pathnames to suit your needs  ######
###############################################################

data_dir="/nfs3/PHARM/David_Lab/christine/quality_vectors_git/data/IBD_huttenhower/" #This should contain your repseqs, seqtab
out_dir="/nfs3/PHARM/David_Lab/christine/quality_vectors_git/data/IBD_huttenhower/embed/" #Where do you want your files to go?
r_software_dir="/nfs3/PHARM/David_Lab/christine/software/R-3.6.0/bin"
blast_software_dir="../../software/ncbi-blast-2.9.0+/bin/"
blast_db="blastdb"
seqtab_file_name="seqtab_final.txt" #Labeled by full ASVs, not ids. Columns are asvs
dim="100"
qual_vec_file_name="/nfs3/PHARM/David_Lab/christine/quality_vectors_git/embed/embed_.07_100id""$dim""dim.txt" #Labled by ids, not full ASVs

###############################################################################
####   Step 2: Use blast to align repseqs to embedding database sequences  ####
###############################################################################
echo "Run blast"
"$blast_software_dir/blastn" -db "$blast_db/embedding_db_.07" -query "$data_dir/repseqs.fasta" -out "$out_dir/blast_hits.tsv"  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"

##############################################################################################
#######  Step 3: Filter blast hits so that we only have the best one per sequence  ###########
##############################################################################################
echo "Filter hits"
cat "$out_dir/blast_hits.tsv" | sort -k1,1 -k5,5g -k6,6nr | sort -u -k1,1 --merge > "$out_dir/best_hits.tsv"

#####################################################
#######  Step 3: Get embedded data  #################
#####################################################
#Arguments:
#asv_file_name
#best_hits_file_name
#embedded_trasmat_file_name
#echo "Embed Data"
#"$r_software_dir/Rscript" embed_asv_table.R "$data_dir/$seqtab_file_name" "$data_dir/repseqs_ml.fasta" "$out_dir/best_hits.tsv" "$qual_vec_file_name" "$out_dir"

#mv "$out_dir/embedded.rds" "$out_dir/embedded_""$dim"".rds"
#mv "$out_dir/embedded.txt" "$out_dir/embedded_""$dim"".txt"

