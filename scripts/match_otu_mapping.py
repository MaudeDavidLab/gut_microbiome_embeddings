import pandas as pd
import numpy as np
import helper_functions as hf
import pickle
import matplotlib.pyplot as plt
from sklearn import preprocessing
import importlib
import math
import copy

data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new"
otu_file = data_dir + "/seqtab_final_filter.07.txt"
otu = pd.read_csv(otu_file, sep = "\t", index_col= 0)
print("Samples: " + str(otu.shape[0]) + "  Taxa: " + str(otu.shape[1]))

mapping_file = data_dir + "/AG_mapping.txt"
mapping = pd.read_csv(mapping_file, sep = "\t", index_col=0)

err_qid = pd.read_csv(data_dir + "/err-to-qid.txt", sep = "\t", index_col = 0)

convert_sample_ids = err_qid.loc[otu.index.values, :]
otu = otu.set_index(convert_sample_ids.sample_title)

otu_clean, qual_vecs_clean = hf.match_otu_qual(otu, qual_vecs)
otu_clean, map_clean = hf.match_otu_map(otu_clean, mapping)

number_criteria = []
cat_criteria = ["IBD", "EXERCISE_FREQUENCY", "SEX", "ONE_LITER_OF_WATER_A_DAY_FREQUENCY", 
        "SEAFOOD_FREQUENCY", "PROBIOTIC_FREQUENCY", "OLIVE_OIL", "FRUIT_FREQUENCY", 
         "SLEEP_DURATION", "SUGAR_SWEETENED_DRINK_FREQUENCY", "MILK_CHEESE_FREQUENCY",
         "RED_MEAT_FREQUENCY","MEAT_EGGS_FREQUENCY", "VEGETABLE_FREQUENCY", "BODY_SITE"]

otu_clean, map_clean = hf.filterForMetadata(otu_clean, map_clean, number_criteria, cat_criteria)



#Make train/test set
test_samples_file = data_dir + "/test_samples.txt"
with open(test_samples_file) as f:
    test_samples = f.read().split()

err_qid = pd.read_csv(data_dir + "/err-to-qid.txt", sep = "\t", index_col = 0)
test_samples = err_qid.loc[test_samples, "sample_title"]
test_samples = test_samples[ test_samples == test_samples] #delete Nan values
test_samples = test_samples[[test_samples[i] in otu_clean.index.values for i in range(len(test_samples))]]


otu_train, otu_test, map_train, map_test = hf.splitTrainTest(otu_clean, map_clean, test_samples)

map_train = map_train.drop('BODY_SITE', axis = 1)
map_test = map_test.drop('BODY_SITE', axis = 1)
map_train, map_test = hf.makeMappingNumeric(map_train, map_test, number_criteria, cat_criteria)

filt = ".07"
f = open(data_dir + "/otu_train_" + str(filt) + ".obj", "wb")
pickle.dump(otu_train, f)
f.close()

f = open(data_dir + "/otu_test_" + str(filt) + ".obj", "wb")
pickle.dump(otu_test, f)
f.close()

f = open(data_dir + "/map_train_" + str(filt) + ".obj", "wb")
pickle.dump(map_train, f)
f.close()

f = open(data_dir + "/map_test_" + str(filt) +  ".obj", "wb")
pickle.dump(map_test, f)
f.close()
