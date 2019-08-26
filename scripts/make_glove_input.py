import csv
import numpy as np
import random

data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new/filter_.05/"

f = open(data_dir + "seqtab_filter5e-04.txt", 'r')
outfile = open(data_dir + 'glove_input_filter5e-04.txt', mode = 'w')
test_samples_file = open(data_dir + '../test_samples.txt', 'w')
print("seqtab_filter5e-04")


writer = csv.writer(outfile, delimiter = "\t", quoting = csv.QUOTE_NONE, escapechar = '')

taxa_names = f.readline()
taxa_names = taxa_names.strip().split("\t")
i = 0
j = 0
test_samples = []
random.seed(0)
for line in f:
    vec = line.split("\t")
    sample_id = vec[0]
    if random.random() > 0.15:
        if i == 3261 or i == 3260:
            print([i for i in vec[1:]])
        present = [float(i) > 0 for i in vec[1:]]
        writer.writerow(np.array(taxa_names)[present])
        print(i, end = '\t')
        i = i + 1
    else:
        test_samples.append(sample_id)
    j += 1
	
test_samples_file.write("\t".join(test_samples))
f.close()
outfile.close()
