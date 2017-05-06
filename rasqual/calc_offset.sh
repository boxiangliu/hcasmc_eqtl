# Calculate RASQUAL offset:
Rscript rasqual/reorder_gcc_by_Y.R ../processed_data/rasqual/expression/Y.txt ../processed_data/rasqual/expression/gcc.exon.txt ../processed_data/rasqual/expression/Y.tidy.txt ../processed_data/rasqual/expression/gcc.exon.tidy.txt
cut -f2-2 ../processed_data/rasqual/gcc.exon.tidy.txt > ../processed_data/rasqual/gcc.exon.tidy.cut.txt
cp /srv/persistent/bliu2/tools/rasqual/R/{makeOffset.R,gcCor.R} rasqual
R --vanilla --quiet --args ../processed_data/rasqual/expression/Y.tidy.txt ../processed_data/rasqual/expression/gcc.exon.tidy.cut.txt ../processed_data/rasqual/expression/K.txt < rasqual/makeOffset.R
