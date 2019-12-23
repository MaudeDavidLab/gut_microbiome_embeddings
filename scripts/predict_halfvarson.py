import pandas as pd
import numpy as np
import helper_functions as hf
import pickle
import matplotlib.pyplot as plt
import importlib
import math
import copy
import re
import sys
from sklearn import preprocessing
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
import random
importlib.reload(hf)
import warnings
from sklearn.model_selection import train_test_split
warnings.filterwarnings('ignore')

from collections import Counter
import multiprocessing as mp

###################################
##### UPDATE VARIABLES  ###########
###################################
print("Starting program")
data_dir = "../data/"
fig_dir = "../data/halfvarson/figs/"
id_thresh = 100

qual_vecs, embed_ids, embed_seqs = hf.getQualVecs(data_dir)
#seed = sys.argv[1]
#print("Seed: " + str(seed))

otu = pd.read_csv(data_dir + "halfvarson/seqtab.txt", sep = "\t")
otu.shape

mapping = pd.read_csv(data_dir + "halfvarson/mapping.txt", sep = "\t")
mapping = mapping.set_index("sample_name")

#Keep only diagnoses we're interested in right now
mapping['sample_name'] = mapping.index.values
keep = [i in ["CD", "UC", "HC"] for i in mapping["diagnosis_full"].values]
mapping = mapping.loc[keep, :]
mapping = mapping.loc[[i in otu.index.values for i in mapping.index.values], :]
otu_use = otu.loc[mapping.index, :]

mapping["sample_id"] = mapping.index.values



best_hits = pd.read_csv(data_dir + "halfvarson/embed_2/best_hits.tsv", header = None, sep = "\t")
best_hits.columns = ["query_id", "hit_id", "query_seq", "hit_seq", "evalue", "bitscore","length", "perc_id"]
best_hits = best_hits.set_index('query_id')
keep = [i < 1E-29 for i in best_hits['evalue'] ]
best_hits = best_hits.loc[keep, :]
keep = [i >= id_thresh for i in best_hits['perc_id']]
best_hits = best_hits.loc[keep, :]
best_hits = best_hits.loc[[i in otu_use.columns.values for i in best_hits['query_seq']], :]

def embed(otu_use):
    print("embedding dataset")
    global best_hits
    global qual_vecs
    #Embed asv table
    #Get only those ASVs that have a match in the embedding set
    otu_use = otu_use.loc[:, best_hits['query_seq'].values]
    print(best_hits.shape)
    
    #Assign id of best match and reorder columns
    otu_use.columns = best_hits['hit_id'].values
    
    #Put transformation matrix in order to be dotted with the ASV table
    qual_vecs_tmp = qual_vecs.loc[otu_use.columns.values, :]
    
    print("OTU shape: " + str(otu_use.shape))
    print("Qual vec shape: " + str(qual_vecs_tmp.shape))
    embedded = pd.DataFrame(np.dot(hf.asinh(otu_use), qual_vecs_tmp))
    embedded.columns = qual_vecs_tmp.columns.values
    
    return(embedded)


def getDiagnosisInfo(otu_train, otu_test, map_train, map_test, pos_class, neg_class):
    map_train = map_train.loc[[i in neg_class or i in pos_class for i in map_train.diagnosis_full], :]
    map_test = map_test.loc[[i in neg_class or i in pos_class for i in map_test.diagnosis_full], :]
    otu_train = otu_train.loc[map_train.index.values, :]
    otu_test = otu_test.loc[map_test.index.values, :]
    y_train = [i in pos_class for i in map_train.diagnosis_full.values]
    y_test = [i in pos_class for i in map_test.diagnosis_full.values]
    return(otu_train, otu_test, y_train, y_test)



def trainHyperParameters(X_train, y_train):

    aucs = []
    precisions = []
    f1s = []
    
    depths = [2, 3, 5, 7]
    n_estimators = [50, 70, 90, 110, 130]
    weights = [1, 3, 5, 10, 15]
    
    #depths = [2, 3]
    #n_estimators = [50, 70]
    #weights = [1, 3]
    
    df = np.zeros((len(depths) * len(n_estimators) * len(weights), 6))
    i = 0
    for max_depth in depths:
        for n_est in n_estimators:
            for weight in weights:
                #print(max_depth, n_est, weight, end = "\t")
                auc_crossVal, auc_prec_crossVal, f1_crossVal, _ = hf.crossValPrediction(X_train, y_train,
                                                                                        max_depth = max_depth,
                                                                                        n_estimators = n_est,
                                                                                        weight = weight ,
                                                                                        folds = 3)
                df[i, :] = [np.mean(auc_crossVal), np.mean(auc_prec_crossVal), np.mean(f1_crossVal), 
                                   max_depth, n_est, weight]
                i += 1
   
    return(df)




def getParamDfs(X_train, y_train):

	asin_params = trainHyperParameters(hf.asinh(X_train), y_train)
	embed_params = trainHyperParameters(embed(X_train), y_train)

	pca = PCA(n_components= np.min([qual_vecs.shape[1], X_train.shape[0]]))
	pca.fit(hf.asinh(X_train))
	otu_pca = pd.DataFrame(pca.transform(hf.asinh(X_train)))
	pca_params = trainHyperParameters(otu_pca, y_train)

	return(asin_params, embed_params, pca_params)



def getBestParams(params):
    #f1 metric
    tmp = params[np.argmax(params[:, 2]), :]
    print(tmp)
    return([int(i) for i in tmp[3:6]])


importlib.reload(hf)

def drawFigures(X_train, X_test, y_train, y_test, asin_params, embed_params, pca_params):
    plot = False
    f = plt.figure(figsize=(15,5))
    _, roc_auc_asin, _, _, _, pr_auc_asin, _, _, _ = hf.predictIBD(X_train = hf.asinh(X_train),
                                                                                X_test = hf.asinh(X_test),
                                                                                y_train = y_train,
                                                                                y_test = y_test,
                                                                                max_depth = asin_params[0],
                                                                                n_estimators = asin_params[1],
                                                                                weight = asin_params[2],
                                                                                plot = plot, plot_pr = plot)
    print("Asinh")
    #print(otu_use.shape)



    f = plt.figure(figsize=(15,5))
    _, roc_auc_embed, _, _, _, pr_auc_embed, _, _, _ = hf.predictIBD(X_train = embed(X_train),
                                                                     X_test = embed(X_test),
                                                                     y_train = y_train,
                                                                     y_test = y_test,
                                                                     max_depth = embed_params[0],
                                                                     n_estimators = embed_params[1],
                                                                     weight = embed_params[2],
                                                                     plot = plot, plot_pr = plot)
    print("Embed")




    f = plt.figure(figsize=(15,5))
    from sklearn.decomposition import PCA
    pca = PCA(n_components= np.min([qual_vecs.shape[1], X_train.shape[0]]))
    pca.fit(hf.asinh(X_train))
    otu_pca = pd.DataFrame(pca.transform(hf.asinh(X_train)))
    otu_pca_test = pd.DataFrame(pca.transform(hf.asinh(X_test)))
    _, roc_auc_pca, _, _, _, pr_auc_pca, _, _, _ = hf.predictIBD(X_train = otu_pca,
                                                               X_test = otu_pca_test, 
                                                               y_train = y_train,
                                                               y_test = y_test,
                                                               max_depth = pca_params[0],
                                                               n_estimators = pca_params[1],
                                                               weight = pca_params[2],
                                                                plot = plot, plot_pr = plot)
    print("PCA")
    #print(otu_pca.shape)

    return(roc_auc_asin, pr_auc_asin, roc_auc_embed, pr_auc_embed, roc_auc_pca, pr_auc_pca)
from matplotlib import pyplot as plt





def runTest(seed, otu_use, mapping):
	try: 
		print("SEED: " + str(seed))
		train_sizes = [10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250]
		#train_sizes = [10, 30]
		n = len(train_sizes)
		df = pd.DataFrame({'train_size' : np.zeros(n), 'roc_asin': np.zeros(n), 'pr_asin':np.zeros(n),
		                                               'roc_embed': np.zeros(n), 'pr_embed':np.zeros(n),
		                                                'roc_pca': np.zeros(n), 'pr_pca':np.zeros(n)})
		i = 0

		for i in range(n):
		    print(i)
		    np.random.seed(seed)
		    train_patient_nums = np.random.choice(np.unique(mapping.patientnumber), train_sizes[i])
		    test_patient_nums = np.unique(mapping.patientnumber[[(not i in train_patient_nums) for i in  mapping.patientnumber]])

		    mapping = mapping.loc[[i in otu_use.index.values for i in mapping.index.values], :]
		    map_train = mapping.loc[[i in train_patient_nums for i in mapping.patientnumber], :]
		    map_test = mapping.loc[[i in test_patient_nums for i in mapping.patientnumber], :]
		    otu_train = otu_use.loc[map_train.index.values, :]
		    otu_test = otu_use.loc[map_test.index.values, :]
		    otu_test.shape[0] / (otu_train.shape[0] + otu_test.shape[0])

		    df.train_size[i] = otu_train.shape[0]
		    
		    otu_train_hc_ibd, otu_test_hc_ibd, y_train_hc_ibd, y_test_hc_ibd = getDiagnosisInfo(otu_train, otu_test, 
		                                                                                    map_train, map_test, 
		                                                                                    pos_class = ["CD" , "UC"], neg_class = "HC")


		    asin_params_hc_ibd , embed_params_hc_ibd, pca_params_hc_ibd = getParamDfs(otu_train_hc_ibd, y_train_hc_ibd)
		    asin_params_hc_ibd_best = getBestParams(asin_params_hc_ibd)
		    embed_params_hc_ibd_best = getBestParams(embed_params_hc_ibd)
		    pca_params_hc_ibd_best = getBestParams(pca_params_hc_ibd)


		    roc_auc_asin, pr_auc_asin, roc_auc_embed, pr_auc_embed, roc_auc_pca, pr_auc_pca = drawFigures(otu_train_hc_ibd, otu_test_hc_ibd, y_train_hc_ibd, y_test_hc_ibd,
		                                                                                    asin_params_hc_ibd_best, embed_params_hc_ibd_best, pca_params_hc_ibd_best)
		    
		    df.roc_asin[i] = roc_auc_asin
		    df.pr_asin[i] = pr_auc_asin
		    df.roc_embed[i] = roc_auc_embed
		    df.pr_embed[i] = pr_auc_embed
		    df.roc_pca[i] = roc_auc_pca
		    df.pr_pca[i] = pr_auc_pca
		    
		    print("Train size: " + str(otu_train.shape[0]))
		    print("ROC asin: " + str(roc_auc_asin))
		    print("PR asin: " + str(pr_auc_asin))
		    print("ROC embed: " + str(roc_auc_embed))
		    print("PR embed: " + str(pr_auc_embed))
		    print("ROC pca: " + str(roc_auc_pca))
		    print("PR pca: " + str(pr_auc_pca))


		df.to_csv(fig_dir + "/halfvarson_df_fullasv_" + str(id_thresh) + "id" + str(seed) + ".csv")

		f = plotLineGraph(df)
		f.savefig(fig_dir + '/halfvarson_performance_trainsize_fullasv_' + str(id_thresh) + 'id' + str(seed) + '.pdf')
	except Exception as e:
		return e



#Step 1: Init multiprocessing.Pool()
pool = mp.Pool(mp.cpu_count())

# Step 2: `pool.apply`
results = [pool.apply(runTest, args=(seed, otu_use, mapping)) for seed in range(0, 100)]

pool.close()


###########################################
####### Read in files to plot results #####
###########################################

def plotLineGraph(df, df_std):
    f = plt.figure(figsize=(15,5))
    flatui = ["#9b59b6", "#3498db", "#e74c3c", "#2ecc71"]
    plt.subplot(1,2,1)
    plt.plot(df.train_size, df.roc_asin, color = flatui[0], marker = 'o')
    plt.plot(df.train_size, df.roc_embed, color = flatui[1], marker = 'o')
    plt.plot(df.train_size, df.roc_pca, color = flatui[2], marker = 'o')
    #plt.errorbar(df.train_size, df.roc_asin, yerr = df_std.roc_asin, color = flatui[0])
    #plt.errorbar(df.train_size, df.roc_embed, yerr = df_std.roc_embed,  color = flatui[1])
    #plt.errorbar(df.train_size, df.roc_pca, yerr = df_std.roc_pca, color = flatui[2])
    
    plt.legend(labels = ['Asin', 'Embed', 'PCA'], loc = 'lower right')

    plt.subplot(1, 2, 2)
    plt.plot(df.train_size, df.pr_asin, color = flatui[0], marker = 'o')
    plt.plot(df.train_size, df.pr_embed, color = flatui[1], marker = 'o')
    plt.plot(df.train_size, df.pr_pca, color = flatui[2], marker = 'o')
    
    plt.legend(labels = ['Asin', 'Embed', 'PCA'], loc = 'lower right')
    return(f)

import glob
import pandas as pd
files = glob.glob(fig_dir + '/' + str(id_thresh) + "id/*.csv")

dfs = []
for f in files:
    df = pd.read_csv(f, index_col = 0)
    dfs.append(df)
print(len(dfs))
df_avg = pd.concat(dfs).groupby(level=0).mean()
df_std = pd.concat(dfs).groupby(level=0).std()
df_std
df_avg = pd.concat(dfs).groupby(level=0).mean()
print(df_avg)
f = plotLineGraph(df_avg, df_std)
f.save_fig(fig_dir + '/halfvarson_performance.png')