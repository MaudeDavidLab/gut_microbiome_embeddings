import pandas as pd
import numpy as np
import helper_functions as hf
import pickle
import matplotlib.pyplot as plt
from sklearn import preprocessing
import importlib
importlib.reload(hf)
import math
import copy
import random




def getQualVecObj(file, otu):
    qual_vecs = pd.read_csv(file, sep = " ", index_col = 0, header=None, dtype = {0:str})
    qual_vecs = qual_vecs.loc[otu.columns.values]
    return(qual_vecs)

def getQualVecs(filt, otu_train):
    filt = str(filt)[1:]
    qual_vecs_50 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_newfilter" + filt + "_50.txt", otu_train)
    qual_vecs_100 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_newfilter" + filt + "_100.txt", otu_train)
    qual_vecs_250 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_newfilter" + filt + "_250.txt", otu_train)
    qual_vecs_500 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_newfilter" + filt + "_500.txt", otu_train)
    qual_vecs_750 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_newfilter" + filt + "_750.txt", otu_train)
    return(qual_vecs_50, qual_vecs_100, qual_vecs_250, qual_vecs_500, qual_vecs_750)



def getQualVecs2(filt, otu_train):
    qual_vecs_50 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_new07perc_feces" + "_50.txt", otu_train)
    qual_vecs_100 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_new07perc_feces" + "_100.txt", otu_train)
    qual_vecs_250 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_new07perc_feces" + "_250.txt", otu_train)
    qual_vecs_500 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_new07perc_feces" + "_500.txt", otu_train)
    qual_vecs_750 = getQualVecObj("../data/AG_new/feces/glove_emb_AG_new07perc_feces" + "_750.txt", otu_train)
    return(qual_vecs_50, qual_vecs_100, qual_vecs_250, qual_vecs_500, qual_vecs_750)


def trainHyperParameters(X_train_list, y_train_list, X_val_list, y_val_list):
    #Need to change this to cross-validation
    depths = [2, 3, 5, 10, 15]
    n_estimators = [50, 65, 80, 95, 110, 125, 140, 155, 170, 200]
    weights = [1,5,20,30]

    #depths = [2, 3]
    #n_estimators = [50, 65]
    #weights = [1,5]


    aucs = np.zeros((len(depths) * len(n_estimators) * len(weights), 6))
    aucs_train = np.zeros((len(depths) * len(n_estimators) * len(weights), 6))
    i = 0
    for depth in depths:
        for trees in n_estimators:
            for weight in weights:
                auc_crossVal = []
                auc_train_crossVal = []
                precision_crossVal = []
                f1_crossVal = []
                for X_train, y_train, X_val, y_val in zip(X_train_list, y_train_list, X_val_list, y_val_list):
                    auc, auc_train, fpr, tpr, precision, f1, _ = hf.predictIBD(X_train, y_train, X_val, y_val, "Embedding weighted by averaging taxa",
                                             max_depth = depth, n_estimators = trees, weight = weight, plot = False)
                    auc_crossVal.append(auc)
                    auc_train_crossVal.append(auc_train)
                    precision_crossVal.append(precision)
                    f1_crossVal.append(f1)
                    
                    
                avg_auc_crossVal = np.mean(auc_crossVal)  
                avg_auc_train_crossVal = np.mean(auc_train_crossVal)
                avg_precision_crossVal = np.mean(precision_crossVal)
                avg_f1_crossVal = np.mean(f1_crossVal)
                
                aucs[i, :] = [avg_auc_crossVal, depth, trees, weight, avg_precision_crossVal, avg_f1_crossVal]
                aucs_train[i, :] = [avg_auc_train_crossVal, depth, trees, weight, avg_precision_crossVal, avg_f1_crossVal]

                print(depth, trees, weight, avg_auc_train_crossVal, avg_auc_crossVal, avg_precision_crossVal, avg_f1_crossVal)
                i = i + 1
                
                 
    return(aucs, aucs_train)  




def getEmbedPlotData(otu_train, otu_test, map_train, map_test, filt, target = "IBD", folds = 10, save_dir = ""):
    print(filt)
    if save_dir == "":
        print("Provide save directory")
    #Get appropriate quality vector objects from file
    if filt is None:
        qual_vecs_50, qual_vecs_100, qual_vecs_250, qual_vecs_500, qual_vecs_750 = getQualVecs2(filt, otu_train)
        
    else:
        qual_vecs_50, qual_vecs_100, qual_vecs_250, qual_vecs_500, qual_vecs_750 = getQualVecs(filt, otu_train)
       
    
    #Train hyperparameters on embeddings at each dimension
    dims = [50, 100, 250, 500, 750]
    i = 0
    embed_aucs = []
    embed_aucs_train = []
    for qual_vecs in [qual_vecs_50, qual_vecs_100, qual_vecs_250, qual_vecs_500, qual_vecs_750]:
        print("Embed " + str(dims[i]))
        X_train_list, X_val_list, X_test, y_train_list, y_val_list, y_test = hf.getCrossValMlInput(otu_train, otu_test, map_train, map_test,target = target, embed = True, qual_vecs = qual_vecs, folds = folds)
        
        embed_aucs_tmp, embed_aucs_train_tmp = trainHyperParameters(X_train_list, y_train_list, X_val_list, y_val_list)
        embed_aucs.append(embed_aucs_tmp)
        embed_aucs_train.append(embed_aucs_train_tmp)
        i = i + 1
        
        
    #Save object
    ##f = open(save_dir + "embed_aucs.obj", "wb")
    #pickle.dump(embed_aucs, f)
    #f.close()
    return(embed_aucs)
    
def getPcaPlotData(otu_train, otu_test, map_train, map_test, filt, target = "IBD", folds = 10, save_dir =""):
    dims = [50, 100, 250, 500, 750]
    if save_dir == "":
        print("Provide save directory")
    #Train hyperparameters on pca reduced at each dimension
    pca_aucs = []
    pca_aucs_train = []
    for dim in dims:
        print("PCA " + str(dim))
        X_train_list, X_val_list, X_test, y_train_list, y_val_list, y_test = hf.getCrossValMlInput(otu_train, otu_test, map_train, map_test, 
                                                                target = target, pca_reduced = True, numComponents = dim, folds = folds)
        pca_aucs_tmp, pca_aucs_train_tmp = trainHyperParameters(X_train_list, y_train_list, X_val_list, y_val_list)
        pca_aucs.append(pca_aucs_tmp)
        pca_aucs_train.append(pca_aucs_train_tmp)
 
    #Save object
    #f = open(save_dir + "pca_aucs.obj", "wb")
    #pickle.dump(pca_aucs, f)
    #f.close()
    return(pca_aucs)

    
def getAsinPlotData(otu_train, otu_test, map_train, map_test, filt, target = "IBD", folds = 10, save_dir = ""):
    print("Asin")
    if save_dir == "":
        print("Provide save directory")
    X_train_list, X_val_list, X_test_list, y_train_list, y_val_list, y_test = hf.getCrossValMlInput(otu_train, otu_test, map_train, map_test, 
                                                                target = target, asinNormalized = True, folds = folds)
    asin_aucs, asin_aucs_train = trainHyperParameters(X_train_list, y_train_list, X_val_list, y_val_list)
    
    #Save object
    #f = open(save_dir + "asin_aucs.obj", "wb")
    #pickle.dump(asin_aucs, f)
    #f.close()
    return(asin_aucs)

def getPlotData(asin_aucs, embed_aucs, pca_aucs):
    aucs_plot = pd.concat([pd.DataFrame(asin_aucs), 
                           pd.DataFrame(embed_aucs[0]), pd.DataFrame(pca_aucs[0]),
                           pd.DataFrame(embed_aucs[1]), pd.DataFrame(pca_aucs[1]),
                           pd.DataFrame(embed_aucs[2]), pd.DataFrame(pca_aucs[2]),
                           pd.DataFrame(embed_aucs[3]), pd.DataFrame(pca_aucs[3]),
                           pd.DataFrame(embed_aucs[4]), pd.DataFrame(pca_aucs[4])])
    

    aucs_plot.columns = ["aucs", "depth", "numTrees", "weight", "precision", "f1"]
    
    aucs_plot["method"] =  ["asin"] * asin_aucs.shape[0] \
    + ["embed_50"] * embed_aucs[0].shape[0] + ["pca_50"] * pca_aucs[0].shape[0] \
    + ["embed_100"] * embed_aucs[1].shape[0] + ["pca_100"] * pca_aucs[1].shape[0] \
    + ["embed_250"] * embed_aucs[2].shape[0] + ["pca_250"] * pca_aucs[2].shape[0] \
    + ["embed_500"] * embed_aucs[3].shape[0] + ["pca_500"] * pca_aucs[3].shape[0] \
    + ["embed_750"] * embed_aucs[4].shape[0] + ["pca_750"] * pca_aucs[4].shape[0] \


    aucs_plot["method_color"] =  ["asin"] * asin_aucs.shape[0] \
    + ["embed"] * embed_aucs[0].shape[0] + ["pca"] * pca_aucs[0].shape[0] \
    + ["embed"] * embed_aucs[1].shape[0] + ["pca"] * pca_aucs[1].shape[0] \
    + ["embed"] * embed_aucs[2].shape[0] + ["pca"] * pca_aucs[2].shape[0] \
    + ["embed"] * embed_aucs[3].shape[0] + ["pca"] * pca_aucs[3].shape[0] \
    + ["embed"] * embed_aucs[4].shape[0] + ["pca"] * pca_aucs[4].shape[0] \
    
   


    data = aucs_plot
    cat_dtype = pd.api.types.CategoricalDtype(categories=["asin", "pca_50", "embed_50", "pca_100", "embed_100",\
                                                         "pca_250", "embed_250", "pca_500", "embed_500",\
                                                         "pca_750", "embed_750"], ordered=True)
    
    
    data.method = data.method.astype(cat_dtype)



    return(aucs_plot)



############################
############################
############################
import seaborn as sns
from matplotlib import rc

def plotRobustness(aucs_plot,y, pvals, numTaxa = 0):
    import seaborn as sns
    fig = plt.figure(figsize=(25,10))
    ax = fig.add_subplot(111)

    ax2 = ax.twinx()
    if y == "aucs":
        ax.set(ylim=(0.40, 0.8))
        ax2.set(ylim=(0.40, 0.8))
    if y == "precision":
        ax.set(ylim=(0, 0.3))
        ax2.set(ylim=(0, 0.3))
    
    swarm = sns.swarmplot(x = "method", y = y, 
                        data = aucs_plot,
                        hue = "weight",
                        alpha = 0.75, zorder = 1,
                        ax = ax2, palette = sns.color_palette("Set2"))
    
    ax2.set(xticklabels = [])
    flatui = ["#9b59b6", "#3498db", "#e74c3c", "#2ecc71"]
    box = sns.boxplot(x = "method", y = y, data = aucs_plot, hue = "method_color", ax = ax, palette = flatui)
    ax.set(xticklabels = [])
    
    handles, _ = ax2.get_legend_handles_labels()
    chartBox = ax2.get_position()
    ax2.legend(loc='upper center', bbox_to_anchor=(1.075, 0.8), shadow=True, ncol=1, title = "Weight on Positive Class")

    handles, _ = ax.get_legend_handles_labels()
    chartBox = ax.get_position()
    ax.legend(loc='upper center', bbox_to_anchor=(1.075, 0.6), shadow=True, ncol=1, title = "Dim. Red. Method")

    ax.set_ylabel("AUC", fontsize = 30)
    ax2.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticklabels([])
    ax2.set_xticklabels([])
    ax.set_xticks([])
    ax2.set_xticks([])
    #ax.margins(y = 0.1)
    
    
    # Colored x labels
    vert_offset = 0.09
    
    
    x1  = 0
    z, h, col = aucs_plot['aucs'].min() - vert_offset, 0.001, 'k'
    ax.text(x1, z+h, "Asin", ha='center', va='bottom', color= flatui[0], fontsize = 28)
    
    for i in range(5):
        x1  = 2*i + 1
        z, h, col = aucs_plot['aucs'].min() - vert_offset, 0.001, 'k'
        ax.text(x1, z+h, "PCA", ha='center', va='bottom', color= flatui[2], fontsize = 28)

        x1  = 2*i + 2
        z, h, col = aucs_plot['aucs'].min() - vert_offset, 0.001, 'k'
        ax.text(x1, z+h, "Embed", ha='center', va='bottom', color= flatui[1], fontsize = 28)
        
        
        
        #x1  = 2*i + 1.5
        #z, h, col = aucs_plot['aucs'].min() - p_offset, 0.001, 'k'
        #ax.text(x1, z+h, "p = " + str(round(pvals[i], 3)), ha='center', va='bottom', color= "black" , fontsize = 28)


    
    
    #ax.set_xticklabels(["Asin", "PCA", "Embed", "PCA", "Embed", "PCA", "Embed","PCA", "Embed","PCA", "Embed"],
    #                   fontsize = 30, 
    #                   color = "black") 
                       #flatui[0] + flatui[[1,2,1,2,1,2,1,2,1,2]]
    # statistical annotation
    vert_offset = 0.05
    p_offset = -0.03
    p_fontsize = 24
    x1 = 0
    z, h, col = aucs_plot[y].max() + vert_offset, 0.001, 'k'
    ax.text(x1, z + h, numTaxa, ha='center', va='bottom', color=col, fontsize = 28)
    
    
    x1, x2 = 1, 2   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    z, h, col = aucs_plot[y].max() + vert_offset, 0.001, 'k'
    ax.plot([x1, x1, x2, x2], [z, z+h, z+h, z], lw=1.5, c=col)
    ax.text((x1+x2)*.5, z+h, "50 dim", ha='center', va='bottom', color=col, fontsize = 28)
    ax.text((x1+x2)*.5, z + p_offset,  "p = " + '%.3e' % pvals[0], ha='center', va='bottom', color=col, fontsize = p_fontsize)

    x1, x2 = 3,4   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    z, h, col = aucs_plot[y].max() + vert_offset, 0.001, 'k'
    ax.plot([x1, x1, x2, x2], [z, z+h, z+h, z], lw=1.5, c=col)
    ax.text((x1+x2)*.5, z+h, "100 dim", ha='center', va='bottom', color=col, fontsize = 28)
    ax.text((x1+x2)*.5, z + p_offset,  "p = " + '%.3e' % pvals[1], ha='center', va='bottom', color=col, fontsize = p_fontsize)

    x1, x2 = 5,6   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    z, h, col = aucs_plot[y].max() + vert_offset, 0.001, 'k'
    ax.plot([x1, x1, x2, x2], [z, z+h, z+h, z], lw=1.5, c=col)
    ax.text((x1+x2)*.5, z+h, "250 dim", ha='center', va='bottom', color=col, fontsize = 28)
    ax.text((x1+x2)*.5, z + p_offset,  "p = " + '%.3e' % pvals[2], ha='center', va='bottom', color=col, fontsize = p_fontsize)

    x1, x2 = 7,8  # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    z, h, col = aucs_plot[y].max() + vert_offset, 0.001, 'k'
    ax.plot([x1, x1, x2, x2], [z, z+h, z+h, z], lw=1.5, c=col)
    ax.text((x1+x2)*.5, z+h, "500 dim", ha='center', va='bottom', color=col, fontsize = 28)
    ax.text((x1+x2)*.5, z + p_offset,  "p = " + '%.3e' % pvals[3], ha='center', va='bottom', color=col, fontsize = p_fontsize)
        
    x1, x2 = 9,10   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    z, h, col = aucs_plot[y].max() +vert_offset, 0.001, 'k'
    ax.plot([x1, x1, x2, x2], [z, z+h, z+h, z], lw=1.5, c=col)
    ax.text((x1+x2)*.5, z+h, "750 dim", ha='center', va='bottom', color=col, fontsize = 28)
    ax.text((x1+x2)*.5, z + p_offset,  "p = " + '%.3e' % pvals[4], ha='center', va='bottom', color=col, fontsize = p_fontsize)
    
    return(fig)
    #handles_box, _ = ax.get_legend_handles_labels()
    #chartBox = ax.get_position()
    #ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    #ax.legend(loc='upper center', bbox_to_anchor=(1.2, 0.8), shadow=True, ncol=1, title = "Dim. Red. Method")