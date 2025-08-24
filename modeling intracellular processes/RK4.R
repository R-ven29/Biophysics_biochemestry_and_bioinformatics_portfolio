library(ggplot2)

Ec_p <- function(R, P){
	inc=(alpha*R-beta*P)
	return(inc)
}

Ec_r <- function(R, P){
	inc=(1/(1+P)-R)
	return(inc)
}

p0=0
r0=0
beta=1
alpha=1
time_increment=0.1
p_r_matrix=matrix(0, nrow=100, ncol=3)

p_r_matrix[1,1]=p0
p_r_matrix[1,2]=r0
p_r_matrix[1,3]=0
x=2

for(time in seq(time_increment*2, 10, by=time_increment)){
	p0=p_r_matrix[x-1,1]
	r0=p_r_matrix[x-1,2]

	k1_r=time_increment*Ec_r(R=r0, P=p0)
	k1_p=time_increment*Ec_p(R=r0, P=p0)

	k2_r=time_increment*Ec_r(R=r0+k1_r/2, P=p0+k1_p/2)
	k2_p=time_increment*Ec_p(R=r0+k1_r/2, P=p0+k1_p/2)
	
	k3_r=time_increment*Ec_r(R=r0+k2_r/2, P=p0+k2_p/2)
	k3_p=time_increment*Ec_p(R=r0+k2_r/2, P=p0+k2_p/2)

	k4_r=time_increment*Ec_r(R=r0+k3_r, P=p0+k3_p)
	k4_p=time_increment*Ec_p(R=r0+k3_r, P=p0+k3_p)

	pn=p0+(k1_p+2*k2_p+2*k3_p+k4_p)/6
	rn=r0+(k1_r+2*k2_r+2*k3_r+k4_r)/6

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
