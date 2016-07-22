# libraries:
library(gap)

# function to perform equivalence test for multiple groups: 
equivalence=function(response,group,alpha,eps){
	x=data.table(group=group,response=response)
	N=nrow(x)
	k=length(unique(x$group))
	n_bar=N/k
	x[,n:=.N,by='group']
	x[,mu:=sum(response)/n,by='group']
	df=N-k
	S_within=x[,sum((response-mu)^2)/df]
	x_summ=x%>%dplyr::select(group,n,mu)%>%unique()
	mu_dot=x_summ[,sum(n*mu)/sum(n)]
	temp=x%>%group_by(group)%>%summarize(s=(n/n_bar)*(mu-mu_dot)^2)%>%unique()%>%mutate(S=sum(s))
	S_between=unique(temp$S)
	psi2=S_between/S_within
	c=qf(alpha,df1=k-1,df2=N-k,ncp=n_bar*eps^2)
	cr=(k-1)/n_bar*qf(alpha,df1=k-1,df2=N-k,ncp=n_bar*eps^2)
	res=list()
	res$n_bar=n_bar
	res$k=k
	res$N=N
	res$alpha=alpha
	res$eps=eps
	res$psi2=psi2
	res$cr=cr
	return(res)
}