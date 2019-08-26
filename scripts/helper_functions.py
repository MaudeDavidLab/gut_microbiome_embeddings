import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import pickle
import random

from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import confusion_matrix
import math
import copy
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.metrics import precision_recall_curve


from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import average_precision_score
from inspect import signature
from sklearn.metrics import f1_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import fbeta_score


def asinh(otu):
    return(np.arcsinh(otu))

def log_normalize(otu):
    return(np.log(otu + 1))

def is_number(s):
    try:
        if np.isnan(float(s)):
            return False
        else: 
            return True
    except ValueError:
        return False
    
def getQualVecs(data_dir):
    qual_vec_file = data_dir + "embed/embed_.07_100dim.txt"
    qual_vecs = pd.read_csv(qual_vec_file, sep = " ", index_col = 0, dtype = {0:str})
    qual_repseqs = pd.read_csv(data_dir + "embed/seqs_.07_embed.fasta", sep = "\t", header = None)
    
    import re
    ids = qual_repseqs.iloc[range(0, qual_repseqs.shape[0], 2), 0]
    ids = [re.sub(">", "", i) for i in ids.values]

    seqs = qual_repseqs.iloc[range(1, qual_repseqs.shape[0], 2), 0]
    seqs = [str(i) for i in seqs.values]

    #Drop <unk> character
    ids = ids[0: len(ids)-1]
    seqs = seqs[0: len(seqs)-1]
    qual_vecs = qual_vecs.iloc[0: len(seqs), :]

    print(len(ids))
    print(len(seqs))
    print(qual_vecs.shape)
    return(qual_vecs, ids, seqs)

def matchOtuQual(otu, embed_ids, embed_seqs):
    otu_reorder = otu.loc[:, embed_seqs]
    if np.sum([i==j for i,j in zip(otu_reorder.columns.values, embed_seqs)]) == len(embed_seqs):
        print("all good")
    else:
        print("There's a problem, stop")
    otu_reorder.columns = embed_ids
    return(otu_reorder)

def match_otu_map(otu, mapping):
    map_clean = mapping.loc[otu.index.values,:] #Keep samples if present in otu
    otu = otu.loc[map_clean.index.values, :]
    map_clean = map_clean.reindex(sorted(map_clean.index.values))
    otu = otu.reindex(sorted(otu.index.values), axis = 0)
    
    keep = map_clean.HMP_SITE == "FECAL"
    map_clean = map_clean.loc[keep,:]
    otu = otu.loc[keep,:]
    return(otu, map_clean)

def filterForMetadata(otu_clean, map_clean, number_criteria, cat_criteria):
    keep = [True] * map_clean.shape[0]
    print("Samples originally: " + str(sum(keep)))
    for criteria in cat_criteria:
        #Keep sample if it has desired metadata available
        keep_tmp = [( (i != "Unknown") and (i != "Unspecified") and (i!="other" ) and (i != "unspecified") and (isinstance(i, str)) )  for i in map_clean[criteria]] 
   
        keep = [(i and j) for (i,j) in zip(keep, keep_tmp)]
    print("Samples after categorical filter: " + str(sum(keep)))

    for criteria in number_criteria:
        keep_tmp = [is_number(i) for i in map_clean[criteria]]
        keep = [i and j for (i,j) in zip(keep, keep_tmp)] 
    print("Samples after numerical filter: " + str(sum(keep)))

        
    otu_keep = otu_clean.loc[keep, :]
    map_keep = map_clean.loc[keep, cat_criteria + number_criteria]   
   

    print("Filter for desired metadata present")
    print("Samples: " + str(otu_keep.shape[0]) + "  Taxa: " + str(otu_keep.shape[1]))
    return(otu_keep, map_keep)


def splitTrainTest(otu_keep, map_keep, test_samples):
    test_samples = test_samples[[i in otu_keep.index.values for i in test_samples]] #Only include ids that we haven't dropped in previous steps
    otu_train = otu_keep.loc[[not(i in np.array(test_samples)) for i in otu_keep.index.values], :]
    otu_test = otu_keep.loc[test_samples, :]

    map_train = map_keep.loc[[not(i in np.array(test_samples)) for i in map_keep.index.values], :]
    map_test = map_keep.loc[test_samples, :]

    print("OTU TRAIN: " + str(otu_train.shape))
    print("MAP TRAIN: " + str(map_train.shape))
    print("OTU TEST: " + str(otu_test.shape))
    print("MAP TEST: " + str(map_test.shape))
    return(otu_train, otu_test, map_train, map_test)

def getDataAG(otu_file, mapping_file, qual_vec_file, test_samples_file, number_criteria, cat_criteria):
    
    otu = pd.read_csv(otu_file, sep = "\t", index_col= 0)
    otu.index = otu.index.map(str)
    otu.head()
    otu = otu.T
    print("Original data dimensions")
    print("Samples: " + str(otu.shape[0]) + "  Taxa: " + str(otu.shape[1]))
    
    #Format quality vector matrices
    qual_vecs = pd.read_csv(qual_vec_file, sep = " ", index_col = 0, header=None, dtype = {0:str})
    otu_clean, qual_vecs_clean = match_otu_qual(otu, qual_vecs)
    

    mapping = pd.read_csv(mapping_file, sep = "\t", index_col=0)
    otu_clean, map_clean = match_otu_map(otu_clean, mapping)
    otu_clean, map_clean = filterForMetadata(otu_clean, map_clean, number_criteria, cat_criteria)

    #Make train/test set
    f = open(test_samples_file, "rb")
    test_samples = pickle.load(f)
    f.close

    otu_train, otu_test, map_train, map_test = splitTrainTest(otu_clean, map_clean, test_samples)

    return(otu_train, otu_test, qual_vecs_clean, map_train, map_test)


def normalize(otu):
    #Normalize
    sample_sums = otu.sum(axis=1)
    otu_norm = otu.div(sample_sums, axis=0)
    return(otu_norm)


def embed(otu, qual_vecs):
    qual_vecs_use = qual_vecs.loc[list(otu.columns.values)]
    df = pd.DataFrame(np.dot(otu, qual_vecs_use), index = otu.index.values)
    return(df)
    
    
def plotLineOfBestFit(xi, y, title = ""):
    #plt.figure(figsize=(15,5))
    print(np.max(y))

    #y = y.values
    slope, intercept, r_value, p_value, std_err = stats.linregress(xi,y)
    line = slope*xi+intercept
    plt.plot(xi,y,'o', color = "#AC76C8")
    plt.plot(xi, line, color = "orange")
    
    #perfect line
    #plt.plot(
    
    #plt.xlabel("Predicted Value")
    #plt.ylabel("True Value")
    #plt.title(title)
    plt.ylim((np.min(y),np.max(y)))
    plt.xlim((np.min(y), np.max(y)))
    print("Slope: " + str(slope))
    print("R value: " + str(r_value))
    
    
    
def makeNumeric(var, dictionary, map_train, map_test, map_train_correct, map_test_correct):
    map_train_correct[var] = [dictionary[i] for i in map_train[var]]
    map_test_correct[var] = [dictionary[i] for i in map_test[var]]
    return(map_train_correct, map_test_correct)

def makeMappingNumeric(map_train, map_test, number_criteria, cat_criteria):
    map_train_correct = copy.deepcopy(map_train)
    map_test_correct = copy.deepcopy(map_test)
    
    if len(number_criteria)>0:
        for num_var in number_criteria:
            print(num_var)
            map_train_correct[num_var] = [float(i) for i in map_train[num_var]]
            map_test_correct[num_var] = [float(i) for i in map_test[num_var]]
  
    
    #Deal with variables that have to do frequency
    freq_dict = {"Never": 0, "Rarely (a few times/month)":1, "Occasionally (1-2 times/week)": 2,
             "Regularly (3-5 times/week)":3, "Daily":4, "Rarely (less than once/week)" : 1}
    
    freq_vars = ["EXERCISE_FREQUENCY", "ONE_LITER_OF_WATER_A_DAY_FREQUENCY", "SEAFOOD_FREQUENCY", "PROBIOTIC_FREQUENCY", 
             "OLIVE_OIL", "FRUIT_FREQUENCY", "SUGAR_SWEETENED_DRINK_FREQUENCY", "MILK_CHEESE_FREQUENCY", "RED_MEAT_FREQUENCY",
            "MEAT_EGGS_FREQUENCY", "VEGETABLE_FREQUENCY"]
    keep = [(i in cat_criteria) for i in freq_vars]
    freq_vars = np.array(freq_vars)[keep]
    for var in freq_vars:
        print(var)
        map_train_correct, map_test_correct = makeNumeric(var, freq_dict, map_train, map_test, map_train_correct, map_test_correct)

    #Deal with sleep
    if "SLEEP_DURATION" in cat_criteria:
        print("SLEEP_DURATION")
        sleep_dict = {"Less than 5 hours":1, "5-6 hours":2, "6-7 hours":3, "7-8 hours":4, "8 or more hours":5 }
        map_train_correct, map_test_correct = makeNumeric("SLEEP_DURATION", sleep_dict, map_train, map_test, map_train_correct, map_test_correct)

    if "SEX" in cat_criteria:
        print("SEX")
        sex_dict = {"male": 0, "female":1}
        map_train_correct, map_test_correct = makeNumeric("SEX", sex_dict, map_train, map_test, map_train_correct, map_test_correct)

    if "IBD" in cat_criteria:
        print("IBD")
        ibd_dict = {'I do not have this condition':0, 'Self-diagnosed':1, 
                    'Diagnosed by a medical professional (doctor, physician assistant)':1,
                    'Diagnosed by an alternative medicine practitioner':1}
        map_train_correct, map_test_correct = makeNumeric("IBD", ibd_dict, map_train, map_test, map_train_correct, map_test_correct)
        
    if "IBS" in cat_criteria:
        print("IBS")
        ibd_dict = {'I do not have this condition':0, 'Self-diagnosed':1, 
                    'Diagnosed by a medical professional (doctor, physician assistant)':1,
                    'Diagnosed by an alternative medicine practitioner':1}
        map_train_correct, map_test_correct = makeNumeric("IBS", ibd_dict, map_train, map_test, map_train_correct, map_test_correct)
    
    if "GLUTEN" in cat_criteria:
        print("GLUTEN")
        gluten_dict = {"No":0, 'I do not eat gluten because it makes me feel bad':1, 'I was diagnosed with celiac disease':2,
                       'I was diagnosed with gluten allergy (anti-gluten IgG), but not celiac disease':2}
        map_train_correct, map_test_correct = makeNumeric("GLUTEN", gluten_dict, map_train, map_test, map_train_correct, map_test_correct)
    return(map_train_correct, map_test_correct)


def predictIBD(X_train, y_train, X_test, y_test, graph_title = "", max_depth = 12, n_estimators = 140, plot = False, plot_pr = False, weight = 20, feat_imp = False, flipped = False):
    weights = {0:1, 1:weight}
    m = RandomForestClassifier(max_depth= max_depth, random_state=0, n_estimators= n_estimators, class_weight = weights)
    m.fit(X_train, y_train)
    probs = m.predict_proba(X_test)
    probs_train = m.predict_proba(X_train)
    roc_auc, fpr, tpr, precision, f1, f2 = computeMLstats(m, data = X_test, y = y_test, plot = plot, plot_pr = plot_pr, graph_title = graph_title, flipped = flipped)
    

    feat_imp_sort = getFeatureImportance(m, data = X_train, y = y_train)
    
    return(m, roc_auc, None, fpr, tpr, precision, f1, f2, feat_imp_sort)

 
def getFeatureImportance(m, data, y):
    feat_imp = m.feature_importances_
    feat_imp_labeled = zip(data.columns.values, feat_imp)
    feat_imp_sort = sorted(feat_imp_labeled, key = lambda t: t[1], reverse = True)
    return(feat_imp_sort)

def computeMLstats(m, data, y, plot = False, plot_pr = False, graph_title = None, flipped = False):
    probs = m.predict_proba(data)
    
    #Flip for opposite class imbalance
    if flipped:
        y = [1 - i for i in y]
        probs = 1 - probs
    
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y, probs[:, 1])
    roc_auc = auc(fpr, tpr)

    
    #Compute precision-recall
    precision, recall, _ = precision_recall_curve(y, probs[:,1])

    #avg_pr = average_precision_score(precision, recall)
    average_precision = average_precision_score(y, probs[:,1])
    
    f1 = f1_score(y, np.argmax(probs, axis = 1))
    f2 = fbeta_score(y, np.argmax(probs, axis = 1), beta = 2)
    
    if plot:
        plt.subplot(1, 2, 1)
        plt.plot(fpr, tpr, lw=2, alpha=0.3, label='AUC ROC = %0.2f' %  roc_auc)
        #'AUC PR = %0.2f' % pr_avg_pr
        
        plt.legend(loc="lower right")
        x = np.linspace(0, 1, 10)
        plt.plot(x, x)
        plt.title(graph_title)
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        
    if plot_pr:
        plt.subplot(1,2,2)
        # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
        step_kwargs = ({'step': 'post'}
                       if 'step' in signature(plt.fill_between).parameters
                       else {})
        plt.step(recall, precision, color='b', alpha=0.2,
                 where='post', label='AUC PR = %0.2f' %  average_precision)
        plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)
        plt.legend(loc="lower right")
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
    
    return(roc_auc, fpr, tpr, average_precision, f1, f2)

from scipy import interp
from itertools import cycle
def predict_multiclass(X_train, y_train, X_test, y_test, graphTitle = "", max_depth = 12, n_estimators = 140, plot = True, weight = 20):
    weights = {0:1, 1:weight, 2:weight}
    
    #y manipulating
    y_train = label_binarize(y_train, classes=[0, 1, 2])
    y_test = label_binarize(y_test, classes=[0, 1, 2])
    
    m = OneVsRestClassifier(RandomForestClassifier(max_depth= max_depth, random_state=0, n_estimators= n_estimators))
    y_score = m.fit(X_train, y_train).predict_proba(X_test)
    
    
    probs = m.predict_proba(X_test)
    probs_train = m.predict_proba(X_train)

    # Compute ROC curve and area the curve
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    n_classes = 3
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    plt.figure()
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
    for i, color in zip(range(n_classes), colors):
        plt.plot(fpr[i], tpr[i], color=color,
                 label='ROC curve of class {0} (area = {1:0.2f})'
                 ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Some extension of Receiver operating characteristic to multi-class')
    plt.legend(loc="lower right")
    plt.show()

    return()

def getPCAReduced(X_train, X_val, X_test, components = 500):
    pca = PCA(n_components= components)
    pca.fit(X_train)
    X_train_pca = pca.transform(X_train)
    X_val_pca = pca.transform(X_val)
    X_test_pca = pca.transform(X_test)
    return(X_train_pca, X_val_pca, X_test_pca)

def plotPCA(table, otu_raw, components):
    pca = PCA(n_components= components)
    pca = pca.fit(table)
    table_pca = pca.transform(table)
    table_pca = table_pca / np.max(table_pca)
    df = pd.DataFrame(table_pca, index = table.index.values)
    sample_sums_table = otu_raw.sum(axis = 1)
    plt.scatter(df.iloc[:,0], df.iloc[:,1], c = sample_sums_table, cmap='viridis')
    plt.colorbar()
    plt.xlabel(pca.explained_variance_ratio_[0])
    plt.ylabel(pca.explained_variance_ratio_[1])
    
    
 #Create input for ML alg including otu abundances and metadata
def combineData(microbe_data, mapping_data, names = []):
    micro_norm = preprocessing.scale(microbe_data)
    map_norm = preprocessing.scale(mapping_data)
    data = pd.concat([pd.DataFrame(micro_norm), pd.DataFrame(map_norm)], axis = 1)
    if not names == []:
        data.columns = np.concatenate((names, [i for i in mapping_data.columns.values]))
    return(data)

def setTarget(mapping, target = ""):
    y = [float(i) for i in mapping[target]]
    mapping = mapping.drop(target, axis = 1)
    return(mapping, y)

def getMlInput(otu_train, otu_test, map_train, map_test, target, 
               embed = False, pca_reduced = False, asinNormalized = False, percNormalized = False, pathwayEmbed = False,
               qual_vecs = None, numComponents = 250, names = []):
    
    #split training set again to get some validation data for training hyperparameters
    otu_train_train = otu_train.sample(frac = 0.9, random_state = 10)
    otu_val = otu_train.drop(otu_train_train.index.values)
    map_train_train = map_train.loc[otu_train_train.index.values]
    map_val = map_train.drop(otu_train_train.index.values)
    
    map_train_train, y_train = setTarget(map_train_train, target = target)
    map_val, y_val = setTarget(map_val, target = target)
    map_test, y_test = setTarget(map_test, target = target)
    
    if embed:
        X_train = combineData(embed_average(otu_train_train, qual_vecs), map_train_train, names = qual_vecs.columns.values)
        X_val = combineData(embed_average(otu_val, qual_vecs), map_val, names = qual_vecs.columns.values)
        X_test = combineData(embed_average(otu_test, qual_vecs), map_test, names = qual_vecs.columns.values)
    elif pca_reduced:
        pca_train, pca_val, pca_test = getPCAReduced(otu_train_train, otu_val, otu_test, components = numComponents)
        X_train = combineData(pca_train, map_train_train, names = names)
        X_val = combineData(pca_val, map_val, names = names)
        X_test = combineData(pca_test, map_test, names = names)
    elif asinNormalized:
        X_train = combineData(asinh(otu_train_train), map_train_train, names = names)
        X_val = combineData(asinh(otu_val), map_val, names = names)
        X_test = combineData(asinh(otu_test), map_test, names = names)
    elif percNormalized: 
        X_train = combineData(otu_train_train.div(otu_train_train.sum(axis=1), axis=0), map_train_train, naming = naming)
        X_val = combineData(otu_val.div(otu_val.sum(axis=1), axis=0), map_val, names = names)
        X_test = combineData(otu_test.div(otu_test.sum(axis=1), axis=0), map_test, names = names)
    elif pathwayEmbed:
        X_train = combineData(embed_average(otu_train_train, pathway_table), map_train_train, names = names)
        X_val = combineData(embed_average(otu_val, pathway_table), map_val, names = names)
        X_test = combineData(embed_average(otu_test, pathway_table), map_test, names = names)
    
    return(X_train, X_val, X_test, y_train, y_val, y_test)  



def getCrossValMlInput(otu_train, otu_test, map_train, map_test, target, 
               embed = False, pca_reduced = False, asinNormalized = False, percNormalized = False, pathwayEmbed = False,
               qual_vecs = None, numComponents = 250, naming = "topics", folds = 10):
    
    map_test, y_test = setTarget(map_test, target = target)
    
    #split training set again to get some validation data for training hyperparameters
    random.seed(1)
    kf = KFold(n_splits = folds)
    kf.get_n_splits(otu_train)
    X_train_list = []
    X_val_list = []

    y_train_list = []
    y_val_list = []

    i = 0
    for train_index, val_index in kf.split(otu_train):
        otu_train_train = otu_train.iloc[train_index, :]
        otu_val = otu_train.iloc[val_index, :]
        map_train_train = map_train.iloc[train_index, :]
        map_val = map_train.iloc[val_index, :]
    
        map_train_train, y_train = setTarget(map_train_train, target = target)
        map_val, y_val = setTarget(map_val, target = target)
        

        if embed:
            X_train = combineData(embed_average(otu_train_train, qual_vecs), map_train_train, naming = naming)
            X_val = combineData(embed_average(otu_val, qual_vecs), map_val, naming = naming)
            if i == 0:
                X_test = combineData(embed_average(otu_test, qual_vecs), map_test, naming = naming)
        elif pca_reduced:
            pca_train, pca_val, pca_test = getPCAReduced(otu_train_train, otu_val, otu_test, components = numComponents)
            X_train = combineData(pca_train, map_train_train, naming = naming)
            X_val = combineData(pca_val, map_val, naming = naming)
            if i == 0:
                X_test = combineData(pca_test, map_test, naming = naming)
        elif asinNormalized:
            X_train = combineData(asinh(otu_train_train), map_train_train, naming = naming)
            X_val = combineData(asinh(otu_val), map_val, naming = naming)
            if i == 0:
                X_test = combineData(asinh(otu_test), map_test, naming = naming)
        elif percNormalized: 
            X_train = combineData(otu_train_train.div(otu_train_train.sum(axis=1), axis=0), map_train_train, naming = naming)
            X_val = combineData(otu_val.div(otu_val.sum(axis=1), axis=0), map_val, naming = naming)
            if i == 0:
                X_test = combineData(otu_test.div(otu_test.sum(axis=1), axis=0), map_test, naming = naming)
        elif pathwayEmbed:
            X_train = combineData(embed_average(otu_train_train, pathway_table), map_train_train, naming = naming)
            X_val = combineData(embed_average(otu_val, pathway_table), map_val, naming = naming)
            if i == 0:
                X_test = combineData(embed_average(otu_test, pathway_table), map_test, naming = naming)
        X_train_list.append(X_train)
        X_val_list.append(X_val)
       
        y_train_list.append(y_train)
        y_val_list.append( y_val)
        i = i + 1
      
    return(X_train_list, X_val_list, X_test, y_train_list, y_val_list, y_test) 


def trainHyperParameters(X_train, y_train, X_val, y_val):
    depths = [2, 3, 5, 10]
    n_estimators = [50, 65, 80, 95, 110, 125, 140, 155]
    aucs = np.zeros((len(depths), len(n_estimators)))
    aucs_train = np.zeros((len(depths), len(n_estimators)))
    for depth in depths:
        for trees in n_estimators:
            auc, auc_train = predictIBD(X_train, y_train, X_val, y_val, "Embedding weighted by averaging taxa",
                                     max_depth = depth, n_estimators = trees, plot = False)
            aucs[depths.index(depth), n_estimators.index(trees)] = np.mean(auc)
            aucs_train[depths.index(depth), n_estimators.index(trees)] = np.mean(auc_train)
            print(depth, trees, np.mean(auc_train), np.mean(auc))
            
    plt.figure(figsize=(15,5))
    for i in range(aucs.shape[0]):
        plt.subplot(2,3, i + 1)
        plt.plot(n_estimators, aucs[i,:])
        plt.plot(n_estimators, aucs_train[i, :])
        plt.title("Depth: " + str(depths[i]))
        plt.xticks(n_estimators)

def embed_average(otu, qual_vecs):
    if(np.sum([i == j for i,j in zip(otu.columns.values, qual_vecs.index.values)]) == otu.shape[1]):
        print("all good")
    else:
        print("There's a problem")
    df = pd.DataFrame(np.dot(asinh(otu), qual_vecs), index = otu.index.values)
    return(df)

from sklearn.model_selection import StratifiedShuffleSplit
def crossValPrediction(otu_use, y, max_depth = 10, n_estimators = 65, weight = 5, plot = False, plot_pr = False, folds = 5):
    kf = StratifiedShuffleSplit(n_splits = folds)
    kf.get_n_splits(otu_use, y)
    
    auc_crossVal = []
    auc_prec_crossVal = []
    f1_crossVal = []
    feat_imp_crossVal = []
    i = 0
    for train_index, val_index in kf.split(otu_use, y):
        otu_train = otu_use.iloc[train_index, :]
        otu_val = otu_use.iloc[val_index, :]
        y_train = np.array(y)[train_index]
        y_val = np.array(y)[val_index]
        
        plt.subplot(1, 2, 1)
        m, auc, auc_train, fpr, tpr, prec, f1, f2, feat_imp = predictIBD(otu_train, y_train, otu_val, y_val,
                  max_depth = max_depth, n_estimators = n_estimators, weight = weight, plot = plot, plot_pr = plot_pr, feat_imp = True)
        auc_crossVal.append(auc)
        auc_prec_crossVal.append(prec)
        f1_crossVal.append(f1)
        feat_imp_crossVal.append(feat_imp)
        
        i = i + 1
    return(auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_crossVal)

#########################################
##### Feature Importance Functions  #####
#########################################
def getFeatImpDf(feat_imp):
    keys = [i[0] for i in feat_imp[0]]
    feat_imp_total = {key: 0 for key in keys}
    feat_imp_total

    for i in range(len(feat_imp)):
        feat_imp_sort = sorted(feat_imp[i], key = lambda t: t[0])
        keys = [j[0] for j in feat_imp_sort]
        values = [j[1] for j in feat_imp_sort]
        for key,value in zip(keys, values):
            feat_imp_total[key] += value
    cum_importance_sort = sorted(feat_imp_total.items(), key=lambda x: x[1], reverse = True)

    feat_imp_df = pd.DataFrame({'properties': [i[0] for i in cum_importance_sort],
                 'cum_importance': [i[1] for i in cum_importance_sort]})

    feat_imp_df = feat_imp_df.set_index("properties")
    feat_imp_df.iloc[0:5, :]
    return(feat_imp_df)


def associate_one_level(feature, label_left, label_right):
    neg_associations[feature][label_left] += 1
    pos_associations[feature][label_right] += 1
    
def left_side_propogate(feature, label_left, label_right):
    if label_left == label_right:
        neg_associations[features[0]][label_left] += 1 #Doesn't matter which label, they're equal
    else:
        associate_one_level(feature, label_left, label_right)

def right_side_propogate(feature, label_left, label_right):
    if label_left == label_right:
        pos_associations[features[0]][label_left] += 1
    else:
        associate_one_level(feature, label_left, label_right)
