# Select background variants
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15
library(foreach)
library(doMC)
registerDoMC(10)


read_snpsnap=function(){
	snpsnap=fread('/srv/persistent/bliu2/shared/SNPsnap/kb1000_collection.dist_to_intron.tab',select=c(1,2,3,4,5,7,25:30))
	snpsnap[,chr:=str_split_fixed(snpID,':',2)[,1]]
	snpsnap=snpsnap[chr%in%c(1:22),]
	return(snpsnap)
}

select_background_variants=function(x,snpsnap,bg_size=500){
	# INPUT: 
	# x: a list of snp IDs (chr:pos)

	background_ls=foreach(i=1:length(x))%dopar%{

		snpid=x[i]
		print(sprintf('INFO - %s',snpid))
		if (!snpid%in%snpsnap$snpID){
			return(NA)
		}
		tmp=snpsnap[snpID==snpid,list(freq_bin,dist_nearest_gene_snpsnap_protein_coding,friends_ld07,gene_count)]
		freq=tmp$freq_bin
		dist=tmp$dist_nearest_gene_snpsnap_protein_coding
		ld07=tmp$friends_ld07
		count=tmp$gene_count

		step=0
		tmp=data.frame()
		while (nrow(tmp)<=bg_size){
			step=step+1
			ub=1+0.1*step
			lb=1-0.1*step
			# tmp=snpsnap[(freq_bin<=freq+step)&(freq_bin>=freq-step)&(gene_count<=ub*count)&(gene_count>=lb*count)&(dist_nearest_gene_snpsnap_protein_coding<=ub*dist)&(dist_nearest_gene_snpsnap_protein_coding>=lb*dist)&(friends_ld07<=ub*ld07)&(friends_ld07>=lb*ld07),]
			tmp=snpsnap[(freq_bin<=freq+step)&(freq_bin>=freq-step)&(dist_nearest_intron<=ub*dist)&(dist_nearest_intron>=lb*dist)&(friends_ld07<=ub*ld07)&(friends_ld07>=lb*ld07),]
			print(sprintf('INFO - %s background variants selected',nrow(tmp)))
		}
		print(sprintf('INFO - tolerance: %s',step))
		set.seed(42)
		tmp=tmp[!snpID%in%x,]
		# message('BUG - ', nrow(tmp))
		if (nrow(tmp)==bg_size){
			background=tmp
		} else {
			background=tryCatch(tmp[sample(1:nrow(tmp),bg_size),],error=function(e){tmp})
		}
		return(background)
	}
	names(background_ls)=x
	return(background_ls)
}


create_index_set=function(background_ls,snpsnap){
	
	index_set=foreach(foreground_variant=names(background_ls),.combine='rbind')%dopar%{
		message('INFO - ',foreground_variant)
		if (!foreground_variant%in%snpsnap$snpID){
			message('Not in SNPsnap. Skipped.')
			return(data.table())
		}
		tmp=rbind(snpsnap[snpID==foreground_variant,],background_ls[[foreground_variant]])
		tmp$foreground_variant=foreground_variant
		return(tmp)
	}
	index_set[,chr:=str_split_fixed(snpID,':',2)[,1]]
	index_set[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]
	return(index_set)
}



