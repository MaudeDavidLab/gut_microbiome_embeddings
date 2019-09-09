# 1. Download fastqs from EBI:
  - Select columns "Run accession" and "Sample Title" to create err-to-qid.txt
  - Grab ftp links for each fastq from the appropriate column in manifest PRJEB11419
  - paste "ftp://" in front of each link
  - wget -i ftp_fastqs.txt
  
# 2. Process using Dada2:
  - dada2_filtering.R
  - dada2_inferVariants.R
  - dada2_chimeras.R

# 3. Filter taxa for prevalence
  - filterTaxaPrevalence.R
 
####Switch to Python, where it's much easier to work with larger matrices

# 4. Make Glove input
  - make_glove_input.py
  - Only use samples designated for training. Save the sample names of the test set. 
  
# 5. Run Glove algorithm
  - runGlove.sh
  
# 6. Label glove output for clarity later on:
  - label_qualvec_transform_mat.R
  
# 7. Process the tables together
  - Match ASV sample names with mapping sample names (err-to-qid.txt)
  - Drop ASVs in ASV table that are not present in transformation matrix
  - Include only samples that have dietary information
  - Include only samples that have IBD diagnosis
  - Split into train and test sets. Use test sample names saved in step 4. 
  - save as python objects otu_train/test_.07.obj and map_train/test_.07.obj
  
# 8. Robustness over hyperparameter test (Fig. 2A):
  - robustness_test.py 
  
# 9. Test model performance on test set (Fig. 2B):
  - predict_AGP_test.py (
  
# 10. Train and test model on independent dataset (Fig. 2C)
  - data from this paper by Halfvarson et. al: https://www.ncbi.nlm.nih.gov/pubmed/28191884
  - halfvarson.predict.py
  
# 11. Run piphillan
  - Get nearest 3-letter genome code for each ASV 

# 12. Get associated pathway information for each 3-letter genome code
   - KEGG_genomes.R
   
# 13. Run correlation test between embedding dimensions and metabolic pathways
  - metabolic_correlations.R
