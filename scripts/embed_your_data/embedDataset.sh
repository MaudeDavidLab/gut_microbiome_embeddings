# You will be starting with a repseqs.fasta file

###############################################################
####### Step 1: Modify all pathnames to suit your needs  ######
###############################################################

data_dir="../../data/halfvarson" #This should contain your repseqs, seqtab
out_dir="../../data/halfvarson/embed"
r_software_dir="/nfs3/PHARM/David_Lab/christine/software/R-3.6.0/bin"
blast_software_dir="../../software/ncbi-blast-2.9.0+/bin/"
blast_db="../../data/blastdb"
seqtab_file_name="seqtab_t.txt"
qual_vec_file_name="../../data/embed/embed_.07_100dim.txt" #Labled by ids, not full ASVs


###############################################################################
####   Step 2: Use blast to align repseqs to embedding database sequences  ####
###############################################################################
echo "Run blast"
"$blast_software_dir/blastn" -db "$blast_db/embedding_db_.07" -query "$data_dir/repseqs.fasta" -out "$out_dir/blast_hits.tsv"  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore"

##############################################################################################
#######  Step 3: Filter blast hits so that we only have the best one per sequence  ###########
##############################################################################################
echo "Filter hits"
cat "$out_dir/blast_hits.tsv" | sort -k1,1 -k5,5g -k6,6nr | sort -u -k1,1 --merge > "$out_dir/best_hits.tsv"

#####################################################
#######  Step 3: Get embedded data  #################
#####################################################
#Arguments:
#asv_file_name - taxa are rows
#best_hits_file_name
#embedded_trasmat_file_name
echo "Embed Data"
"$r_software_dir/Rscript" embed_asv_table.R "$data_dir/$seqtab_file_name" "$out_dir/best_hits.tsv" "$qual_vec_file_name" "$out_dir"

