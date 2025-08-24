library(ggplot2)

p0=0
r0=0
beta=1/4
alpha=1/4
time_increment=0.1
p_r_matrix=matrix(0, nrow=100, ncol=3)

p_r_matrix[1,1]=p0
p_r_matrix[1,2]=r0
p_r_matrix[1,3]=0
x=2

for(time in seq(time_increment*2, 10, by=time_increment)){
	p0=p_r_matrix[x-1,1]
	r0=p_r_matrix[x-1,2]

	pn=p0 + time_increment*(alpha*r0 - beta*p0)
	rn=r0 + time_increment*(1/(1+p0)-r0)
	
	p_r_matrix[x,1]=pn
	p_r_matrix[x,2]=rn
	p_r_matrix[x,3]=time
	x=x+1
}
print(p_r_matrix[100,1:3])

p_r_data=as.data.frame (p_r_matrix)
p_r_plot=ggplot(p_r_data)+
geom_line(aes(x=V3, y=V1, color='protein'))+
geom_line(aes(x=V3, y=V2, color='mRNA'))
head(p_r_data)
p_r_plot
