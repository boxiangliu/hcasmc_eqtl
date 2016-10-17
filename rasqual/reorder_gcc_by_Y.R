library(data.table)
args=commandArgs(T)
Y_file=args[1]
gcc_file=args[2]
Y_out=args[3]
gcc_out=args[4]

# Y_file='../processed_data/rasqual/Y.txt'
# gcc_file='../processed_data/rasqual/gcc.exon.txt'
# Y_out='../processed_data/rasqual/Y.tidy.txt'
# gcc_out='../processed_data/rasqual/gcc.exon.tidy.txt'


Y=fread(Y_file)
gcc=fread(gcc_file)



# remove any element in Y that is not in gcc: 
Y=Y[Y$V1%in%gcc$V1,]


# remove any element in Y that has zero count:
Y=Y[which(rowSums(Y[,2:ncol(Y),with=F])!=0),]


# reorder gcc according to Y:
gcc=gcc[match(Y$V1,gcc$V1),]


# sanity check: 
stopifnot(gcc$V1==Y$V1)


# output: 
write.table(Y,Y_out,quote=F,row.names=F,col.names=F,sep='\t')
write.table(gcc,gcc_out,quote=F,row.names=F,col.names=F,sep='\t')