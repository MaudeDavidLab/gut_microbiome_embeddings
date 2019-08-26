glove_dir="/nfs3/PHARM/David_Lab/christine/software/GloVe-1.2/build"
filter=test

echo "filter"
echo "$filter"

echo "vocab count"
$glove_dir/vocab_count -verbose 2 -max-vocab 50000 -min-count 10 < glove_input_07perc_feces.txt > vocab_filter"$filter".txt

echo "Coooccurence count"
$glove_dir/cooccur -verbose 2 -symmetric 1 -window-size 50000 -vocab-file vocab_filter"$filter".txt -memory 20 -overflow-file overflow_filter"$filter" < glove_input_07perc_feces.txt > cooccur_filter"$filter".bin

echo "\n Cooccurent shuffle"
$glove_dir/shuffle -verbose 2 -memory 20 < cooccur_filter"$filter".bin  > cooccur_filter"$filter".shuf.bin

echo "Run Glove 50"
$glove_dir/glove -input-file cooccur_filter"$filter".shuf.bin -vocab-file vocab_filter"$filter".txt -save-file glove_emb_filter"$filter""_$vector_size" -gradsq-file gradsq_filter"$filter" -verbose 2 -vector-size 50 -threads 10 -iter 100

echo "Run Glove 100"
$glove_dir/glove -input-file cooccur_filter"$filter".shuf.bin -vocab-file vocab_filter"$filter".txt -save-file glove_emb_filter"$filter""_$vector_size" -gradsq-file gradsq_filter"$filter" -verbose 2 -vector-size 100 -threads 10 -iter 100

echo "Run Glove 250"
$glove_dir/glove -input-file cooccur_filter"$filter".shuf.bin -vocab-file vocab_filter"$filter".txt -save-file glove_emb_filter"$filter""_$vector_size" -gradsq-file gradsq_filter"$filter" -verbose 2 -vector-size 250 -threads 10 -iter 100

echo "Run Glove 500"
$glove_dir/glove -input-file cooccur_filter"$filter".shuf.bin -vocab-file vocab_filter"$filter".txt -save-file glove_emb_filter"$filter""_$vector_size" -gradsq-file gradsq_filter"$filter" -verbose 2 -vector-size 500 -threads 10 -iter 100

echo "Run Glove 750"
$glove_dir/glove -input-file cooccur_filter"$filter".shuf.bin -vocab-file vocab_filter"$filter".txt -save-file glove_emb_filter"$filter""_$vector_size" -gradsq-file gradsq_filter"$filter" -verbose 2 -vector-size 750 -threads 10 -iter 100

