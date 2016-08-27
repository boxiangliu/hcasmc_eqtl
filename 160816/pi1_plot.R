args=commandArgs(T)
# in_file='../processed_data/160816/replication/pi1.txt'
# out_fig='../figures/160816/pi1.pdf'
in_file=args[1]
out_fig=args[2]
pi1=fread(args[1])%>%rename(tissue=V1,pi1=V2)
p=ggplot(pi1,aes(x=reorder(tissue,-pi1),y=pi1))+geom_bar(stat='identity')+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('Tissue')
save_plot(out_fig,p,base_width=6,base_height=6)
