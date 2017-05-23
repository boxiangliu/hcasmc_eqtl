library(cowplot)
library(MASS)

# data from http://www.pubhealth.org.cn/cn/magazine_data/zhyf2001/0102pdf/010208.pdf. retrieved 09/30/2016
# Trend of changes in mortality of cardiovascular diseases in some areas of Beijing during 1984 to 1998 WU Guixian et al, Beijing Institute of Heart, Lung and Blood Vessel Disease
beijing=c(7.4,8.1,8.6,7.9,9.4,9.1,11.0,11.9,9.0,10.4,7.7,9.6,7.9,10.2,12.9)
# men=c(8.2,8.4,9.9,8.9,10.5,9.8,11.5,11.9,9.6,10.2,8.5,17.1,14.9,18.0,14.8)
# women=c(5.9,7.8,6.7,7.9,8.1,10.5,9.1,8.1,10.6,7.6,11.4,11.6,17.4,9.6)
to_plot=data.frame(year=seq(1984,1998),pct=beijing)
fit=rlm(pct~year,to_plot)
p1=ggplot(to_plot,aes(x=year,y=pct))+geom_point()+theme_bw()+stat_smooth(method='rlm')+xlab('Year')+ylab('Percentage')+annotate(geom='text',x=1988,y=13,label=paste0('Annual change: ', prettyNum(coefficients(fit)['year'],digits=2),'%'))
save_plot('../figures/prevalence_of_CAD/beijing_total_1984to1998.pdf',p1)


# data from https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6040a1.htm. retrieved 09/30/2016
# Prevalence of Coronary Heart Disease --- United States, 2006--2010, CDC, 
us=c(6.7,6.2,6.3,5.8,6.0)
to_plot=data.frame(year=seq(2006,2010),pct=us)
fit=rlm(pct~year,to_plot)
p2=ggplot(to_plot,aes(x=year,y=pct))+geom_point()+theme_bw()+stat_smooth(method='rlm')+xlab('Year')+ylab('Percentage')+annotate(geom='text',x=2008,y=7,label=paste0('Annual change: ', prettyNum(coefficients(fit)['year'],digits=2),'%'))
save_plot('../figures/prevalence_of_CAD/us_total_2006to2010.pdf',p2)
