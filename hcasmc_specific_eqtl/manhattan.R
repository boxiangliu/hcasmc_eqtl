# make manhattan plot of eQTLs, using color to indicate hcasmc specific eqtls

#### library:
library(lattice)
library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)


#### function:
#' make manhattan plot: 
manhattan.plot<-function(chr, pos, pvalue, 
        sig.level=NA, annotate=NULL, ann.default=list(),
        should.thin=T, thin.pos.places=2, thin.logp.places=2, 
        xlab="Chromosome", ylab=expression(-log[10](p-value)),
        col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
 
        if (length(chr)==0) stop("chromosome vector is empty")
        if (length(pos)==0) stop("position vector is empty")
        if (length(pvalue)==0) stop("pvalue vector is empty")
 
        #make sure we have an ordered factor
        if(!is.ordered(chr)) {
                chr <- ordered(chr)
        } else {
                chr <- chr[,drop=T]
        }
 
        #make sure positions are in kbp
        if (any(pos>1e6)) pos<-pos/1e6;
 
        #calculate absolute genomic position
        #from relative chromosomal positions
        posmin <- tapply(pos,chr, min);
        posmax <- tapply(pos,chr, max);
        posshift <- head(c(0,cumsum(posmax)),-1);
        names(posshift) <- levels(chr)
        genpos <- pos + posshift[chr];
        getGenPos<-function(cchr, cpos) {
                p<-posshift[as.character(cchr)]+cpos
                return(p)
        }
 
        #parse annotations
        grp <- NULL
        ann.settings <- list()
        label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                col=NULL, fontface=NULL, fontsize=NULL, show=F)
        parse.label<-function(rawval, groupname) {
                r<-list(text=groupname)
                if(is.logical(rawval)) {
                        if(!rawval) {r$show <- F}
                } else if (is.character(rawval) || is.expression(rawval)) {
                        if(nchar(rawval)>=1) {
                                r$text <- rawval
                        }
                } else if (is.list(rawval)) {
                        r <- modifyList(r, rawval)
                }
                return(r)
        }
 
        if(!is.null(annotate)) {
                if (is.list(annotate)) {
                        grp <- annotate[[1]]
                } else {
                        grp <- annotate
                } 
                if (!is.factor(grp)) {
                        grp <- factor(grp)
                }
        } else {
                grp <- factor(rep(1, times=length(pvalue)))
        }
 
        ann.settings<-vector("list", length(levels(grp)))
        ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
 
        if (length(ann.settings)>1) { 
                lcols<-trellis.par.get("superpose.symbol")$col 
                lfills<-trellis.par.get("superpose.symbol")$fill
                for(i in 2:length(levels(grp))) {
                        ann.settings[[i]]<-list(pch=pch, 
                                col=lcols[(i-2) %% length(lcols) +1 ], 
                                fill=lfills[(i-2) %% length(lfills) +1 ], 
                                cex=cex, label=label.default);
                        ann.settings[[i]]$label$show <- T
                }
                names(ann.settings)<-levels(grp)
        }
        for(i in 1:length(ann.settings)) {
                if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
                ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                        parse.label(ann.settings[[i]]$label, levels(grp)[i]))
        }
        if(is.list(annotate) && length(annotate)>1) {
                user.cols <- 2:length(annotate)
                ann.cols <- c()
                if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
                        ann.cols<-match(names(annotate)[-1], names(ann.settings))
                } else {
                        ann.cols<-user.cols-1
                }
                for(i in seq_along(user.cols)) {
                        if(!is.null(annotate[[user.cols[i]]]$label)) {
                                annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                        levels(grp)[ann.cols[i]])
                        }
                        ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                annotate[[user.cols[i]]])
                }
        }
        rm(annotate)
 
        #reduce number of points plotted
        if(should.thin) {
                thinned <- unique(data.frame(
                        logp=round(-log10(pvalue),thin.logp.places), 
                        pos=round(genpos,thin.pos.places), 
                        chr=chr,
                        grp=grp)
                )
                logp <- thinned$logp
                genpos <- thinned$pos
                chr <- thinned$chr
                grp <- thinned$grp
                rm(thinned)
        } else {
                logp <- -log10(pvalue)
        }
        rm(pos, pvalue)
        gc()
 
        #custom axis to print chromosome names
        axis.chr <- function(side,...) {
                if(side=="bottom") {
                        panel.axis(side=side, outside=T,
                                at=((posmax+posmin)/2+posshift),
                                labels=levels(chr), 
                                ticks=F, rot=0,
                                check.overlap=F
                        )
                } else if (side=="top" || side=="right") {
                        panel.axis(side=side, draw.labels=F, ticks=F);
                }
                else {
                        axis.default(side=side,...);
                }
         }
 
        #make sure the y-lim covers the range (plus a bit more to look nice)
        prepanel.chr<-function(x,y,...) { 
                A<-list();
                maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
                A$ylim=c(0,maxy);
                A;
        }
 
        xyplot(logp~genpos, chr=chr, groups=grp,
                axis=axis.chr, ann.settings=ann.settings, 
                prepanel=prepanel.chr, scales=list(axs="i"),
                panel=function(x, y, ..., getgenpos) {
                        if(!is.na(sig.level)) {
                                #add significance line (if requested)
                                panel.abline(h=-log10(sig.level), lty=2);
                        }
                        panel.superpose(x, y, ..., getgenpos=getgenpos);
                        if(!is.null(panel.extra)) {
                                panel.extra(x,y, getgenpos, ...)
                        }
                },
                panel.groups = function(x,y,..., subscripts, group.number) {
                        A<-list(...)
                        #allow for different annotation settings
                        gs <- ann.settings[[group.number]]
                        A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
                        A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
                        A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
                        A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
                        A$x <- x
                        A$y <- y
                        do.call("panel.xyplot", A)
                        #draw labels (if requested)
                        if(gs$label$show) {
                                gt<-gs$label
                                names(gt)[which(names(gt)=="text")]<-"labels"
                                gt$show<-NULL
                                if(is.character(gt$x) | is.character(gt$y)) {
                                        peak = which.max(y)
                                        center = mean(range(x))
                                        if (is.character(gt$x)) {
                                                if(gt$x=="peak") {gt$x<-x[peak]}
                                                if(gt$x=="center") {gt$x<-center}
                                        }
                                        if (is.character(gt$y)) {
                                                if(gt$y=="peak") {gt$y<-y[peak]}
                                        }
                                }
                                if(is.list(gt$x)) {
                                        gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
                                }
                                do.call("panel.text", gt)
                        }
                },
                xlab=xlab, ylab=ylab, 
                panel.extra=panel.extra, getgenpos=getGenPos, ...
        );
}


#' parse the `geno` column and add `chr` and `pos` columns
format=function(x) {
	tmp=str_split_fixed(x$geno,'_',n=4)
	tmp[,1]=str_replace(tmp[,1],'chr','')
	x$chr=as.integer(tmp[,1])
	x$pos=as.integer(tmp[,2])
	return(x)
}

#### main 
# read eQTL data: 
eqtl=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.sid_parsed.txt')
setnames(eqtl,c('pheno','chr','pos','ref','alt','dist','pval','beta','varbeta'))


# format eQTL data:
eqtl=eqtl%>%mutate(geno=paste(chr,pos,ref,alt,'b37',sep='_'))
eqtl=eqtl%>%select(chr,pos,pval,pheno,geno)
eqtl[,id:=paste0(pheno,'_',geno)]

# read hcasmc specific eQTL ids 
hcasmc_specific_eqtl=fread('../processed_data/hcasmc_specific_eqtl/hcasmc_specific_eqtl.autosomes.txt',header=F)
setnames(hcasmc_specific_eqtl,c('pheno','geno'))


# format hcasmc specific eQTL ids:   
hcasmc_specific_eqtl=format(hcasmc_specific_eqtl)
hcasmc_specific_eqtl[,id:=paste0(pheno,"_",geno)]


# annotate variants based on whether it passes genome wide threshold, and whether it has a HCASMC specific eQTL, :
ann=rep(1, length(eqtl$pval))
# ann[eqtl$geno%in%hcasmc_specific_eqtl$geno]=2
ann[eqtl$id%in%hcasmc_specific_eqtl$id]=2
ann=factor(ann, levels=1:2, labels=c(""," "))
# table(ann) 
       #          HCASMC_specific 
       # 90286677            5362
# ann[eqtl$pheno%in%eqtl[eqtl$pval<5e-8,pheno]]=3
# ann=factor(ann, levels=1:3, labels=c("","HCASMC_specific",'eQTL'))


# make manhattan plot:
fig_file='../figures/hcasmc_specific_eqtl/manhattan.gene.2.png'
png(fig_file, width=950, height=500)
print(manhattan.plot(chr=eqtl$chr,pos=eqtl$pos,pvalue=eqtl$pval,annotate=ann))
dev.off()


# find example for HCASMC-specific eQTLs:  
hcasmc_specific_eqtl2=eqtl[which(ann==" "),]
hcasmc_specific_eqtl2[pval==min(hcasmc_specific_eqtl2$pval),]

hcasmc_specific_eqtl2[pval==min(hcasmc_specific_eqtl2[chr==6,pval]),]
# ENSG00000137331.11_6_30507577_A_AC_b37 IER3
# This gene functions in the protection of cells from Fas- or tumor necrosis factor type alpha-induced apoptosis. Partially degraded and unspliced transcripts are found after virus infection in vitro, but these transcripts are not found in vivo and do not generate a valid protein.

hcasmc_specific_eqtl2[pval==min(hcasmc_specific_eqtl2[chr==16,pval]),]
# ENSG00000169688.10_16_56666848_G_A_b37 MT1B
# MT1B (Metallothionein 1B) is a Protein Coding gene. Among its related pathways are Metallothioneins bind metals and Metabolism.

hcasmc_specific_eqtl2[pval==min(hcasmc_specific_eqtl2[chr==17,pval]),]
# ENSG00000120071.8_17_43654468_C_T_b37 KANSL1
# This gene encodes a nuclear protein that is a subunit of two protein complexes involved with histone acetylation, the MLL1 complex and the NSL1 complex. The corresponding protein in Drosophila interacts with K(lysine) acetyltransferase 8, which is also a subunit of both the MLL1 and NSL1 complexes.


# ACTA2
hcasmc_specific_eqtl2[str_detect(pheno,'ENSG00000107796'),]
#    chr      pos        pval             pheno                geno
# 1:  10 90666952 8.61969e-06 ENSG00000107796.8 10_90666952_C_T_b37
#                                       id
# 1: ENSG00000107796.8_10_90666952_C_T_b37
# The protein encoded by this gene belongs to the actin family of proteins, which are highly conserved proteins that play a role in cell motility, structure and integrity. Alpha, beta and gamma actin isoforms have been identified, with alpha actins being a major constituent of the contractile apparatus, while beta and gamma actins are involved in the regulation of cell motility. This actin is an alpha actin that is found in skeletal muscle. Defects in this gene cause aortic aneurysm familial thoracic type 6. Multiple alternatively spliced variants, encoding the same protein, have been identified. [provided by RefSeq, Nov 2008]
