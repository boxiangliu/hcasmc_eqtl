#!/usr/bin/env Rscript
# bosh liu
# 2016/05/11
# durga
# make plots of ibd
library(scales)
input_file='../processed_data/160511_calc_ibd_plink/indep_pairwise_50_5_0.2.genome'
input=fread(input_file,header=T)


# plot PI_HAT:
pdf('../figures/160511_analyze_ibd/pi_hat.pdf')
with(input,plot(PI_HAT,main='PI_HAT',ylab='p(IBD=1)+2*P(IBD=2)', breaks= c(0,1), col=ifelse(FID1%in%c('1848','1858','24635') | FID2%in%c('1848','1858','24635'),'red','black')))
text(input$PI_HAT,labels=ifelse(input$PI_HAT>0.5,paste(input$FID1,input$FID2,sep="_"),""),srt=90,adj=1.1)
abline(h=0.125,col='red',lty=2)
axis(2,at=0.125,labels=c('0.125'))
legend('left',legend=c('1848 or 1858 or 24635 in comparison', 'otherwise'),col=c('red','black'),pch=1)
dev.off()


# plot p(Z=1):
pdf('../figures/160511_analyze_ibd/Z1.pdf')
plot(input$Z1,main='Z1',ylab='p(IBD=1)')
dev.off()


# 24156:
pdf('../figures/160511_analyze_ibd/24156.pdf')
with(input,plot(PI_HAT,main='PI_HAT',ylab='p(IBD=1)+2*P(IBD=2)', breaks= c(0,1), col=ifelse(FID1%in%c('24156') | FID2%in%c('24156'),'red',alpha('green',0.3))))
legend('left',legend=c('24156 in comparison', 'otherwise'),col=c('red','green'),pch=1)
dev.off()