# 1. Download fastqs from EBI for the American Gut Dataset:
  - Select columns "Run accession" and "Sample Title" to create err-to-qid.txt
  - Grab ftp links for each fastq from the appropriate column in manifest PRJEB11419
  - paste "ftp://" in front of each link
  - In a consol, run `wget -i ftp_fastqs.txt`
  
# 2. Process using Dada2:
### Because the number of samples to process was so large, they are processed in 3 steps, with objects saved
### 1. Change the path name in dada2_filtering.R to point to the folder that holds the downloaded fastqs
### 2. Run `Rscript dada2_filtering.R` . The script will create two folders, filtered_100 and filtered_150, based on the number of nucleotides per read. There was a large variability in the quality and length of reads.
### 3. Run `Rscript dada_inferVariants.R  --filtpath --errorFile --outfile `. The script will output and object named using the outfile passed in as an arguement.
### 4. Run `Rscript dada2_chimeras.R`. The script loads all of the seqtab objects from the previous step, merges them, and removes chimeric sequences. The script outputs the file seqtab_final.txt, which is the sample x ASV table. It will also output seqtab_filterxx.txt files where ASVs are filtered to include only those present in 0.01%, 0.03%, 0.05%, 0.07%, 0.09% of samples.

# 3. Make Glove input
### 1. Change data_dir to path where the seqtab text files are located, and the outfile appropriately.
### 2. Run `python make_glove_input.py`. The script will output a text file called glove_input_filterXX.txt (as dictated by your outfile variable). This is a text file in which each line contains a list of the full length ASVs that are present in a single samples. The script also outputs a file called test_samples that includes the sample ids of the samples NOT used to make embeddings. These will be used as the test set later on, for figure 3. 
  
# 4. Run Glove algorithm
### 1. Install GloVe, using the instructions here: https://nlp.stanford.edu/projects/glove/
### 2. Specify the folder where Glove is installed, and the filtration level (created in step 3).
### 2. Run `bash runGlove.sh`. The script will output glove_emb_filter"$filter" "_$vector_size" where filter is the filtration level selected and vector_size is the number of properties. This is the embedding transformation matrix (Fig. 1 C)
  
# 5. Label glove output for clarity later on:
### Sometimes full ASV labels will be required, sometimes sequence ids. This script produces equivalent embedding transformation matrices with ids instead of ASV labels.
### `Rscript label_qualvec_transform_mat.R`
  
# 6. Process the tables together
  ### Run `python match_otu_mapping.py`
  ### The script will:
  - Match ASV sample names with mapping sample names (err-to-qid.txt)
  - Drop ASVs in ASV table that are not present in transformation matrix
  - Include only samples that have dietary information
  - Include only samples that have IBD diagnosis
  - Split into train and test sets. Use test sample names saved in step 4. 
  - save as python objects otu_train/test_.07.obj and map_train/test_.07.obj
  
# 7. Robustness over hyperparameter test (Fig. 2):
### Open `robustness_test_ASV.ipynb` with jupyter notebook and run each cell. 
### The script will use cross validation to check the performance of data transformation (ASVs, embedding, PCA) over hyperparameters of a random forest model. It plots each hyperparamter combination as a dot in the boxplot of the appropriate data transformation.
  
# 8. Test model performance on test set (Fig. 3):
### 1. Change the data_dir and fig_dir directories in the script
### 2. Run `python predict_AGP_testset.py`
### The script trains a random forest model on each data transformation, and then tests that model performance on a held out test set (those samples in test_samples.txt from step 3). It outputs figures in the specified figure directory under the name "curves_AGP_test_*data_transformation*.pdf"
  
# 9. Train and test model on independent dataset (Fig. 4)
### 1. Download data from this paper by Halfvarson et. al: https://www.ncbi.nlm.nih.gov/pubmed/28191884
### 2. Change the data_dir, fig_dir, and id_thresh (% similarity of ASV in query to embedding set to be considered a match)
### 3. Run `python halfvarson_predict.py`
### The script trains random forest models using cross validation to predict IBD status. It includes incrementally more samples in the training set, to observe the effect of higher training set sizes on performance, and compare that effect across data transformations. It outputs intermediate csv and pdf files into the figure directory given at the beginning of the file, and outputs a plot showing the AUC ROC and AUC PR over training set size for each transformation.

# 10. Train model on AGP, and test on Halfvarson (Fig. 5B)
### 1. Change the data_dir, fig_dir, and id_thresh (% similarity of ASV in query to embedding set to be considered a match)
### 2. Open `train_ag_test_half.ipynb` in a jupyter notebook, and run all cells. 
Results including accuracy, precision, recall, f1, and f2 scores will be printed in the last two cells for an ASV-based model and an embedding-based model

# 11. Train model on AGP, and test on Schirmer( Fig. 5C)
### 1. Change the data_dir, fig_dir, and id_thresh (% similarity of ASV in query to embedding set to be considered a match)
### 2. Open `train_ag_test_huttenhower.ipynb` in a jupyter notebook, and run all cells. 
Results including accuracy, precision, recall, f1, and f2 scores will be printed in the last two cells for an ASV-based model and an embedding-based model
  
# 12. Run piphillan
### 1. Go to the following url: http://piphillin.secondgenome.com/
### 2. Fill in the form with the appropriate information. "OTU abundance table" is the seqtab_final.txt ASV table from step 3. "Representative Sequence File" is the repseqs.fasta file from step 3. Select KEGG database version Oct 2018, and a 97% similarity cutoff
### 3. Wait to receive an email with the results, and save into a data directory. We will use the file 200p_otu_genome_hit_table.txt which contained the nearest 3-letter genome code

# 13. Get associated pathway information for each 3-letter genome code
### Run `Rscript KEGG_genomes.R`. The script will create an ASV by pathway (kegg id) binary table called otu_pathway_table.txt that indicates which pathways are present in which ASVs according to KEGG and Piphillan.
   
# 14. Run correlation test between embedding dimensions and metabolic pathways
### 1. Change data_dir and fig_dir in the beginning of the script.
### Run `Rscript metabolic_correlations.R`.
### The script will create correlation matrices between properties and metabolic pathways, and plot these. It will also calculate all the significantly correlation pathways per property using a permutation test, and save these in pathways/property_pathway_dict_allsig.txt. Lastly, it will visualize each property as a combination of broad metabolic categories, and save this in properties_visualized.pdf
