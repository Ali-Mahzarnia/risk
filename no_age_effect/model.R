library(ggplot2)
load("connectivity_plain.rda")
load("response.rda")


v2m = function(v){
  x=length(v)
  n=(1+sqrt(1+8*x))/2
  mat <- matrix(0, nrow = n, ncol =  n)
  mat[upper.tri(mat, diag = FALSE)] <- v
  mat = mat + t(mat)
  return(mat)
}


m2v = function(m){
  x=dim(m)[1]
  n=(x-1)*x/2
  v <- vector("double", n)
  v = m[upper.tri(m, diag = FALSE)] 
  return(v)
}


#for example
# v<-c(1,2,3,4,5,6,7,8,9,10)
# m = matrix( seq(1,9),3,3)
# v= m2v(m)
# m= v2m(v)
len = length( m2v(connectivity[,,1]))


image = matrix(NA, nrow = dim(connectivity)[3] , len)
for (i in 1:dim(connectivity)[3]) {
  image[i,] =m2v(connectivity[,,i])
  
}


k =10 # degree of polynomial of t in f(t)
time = as.data.frame(matrix( NA , dim(response)[1] , k))
for (kk in 1:k) {
  time[,kk] = response$age^kk
  colnames(time)[kk] = paste0("t",kk)
}

##### regressing out f(t) from connectome
image_raw = image
for (j in 1:dim(image)[2]) {
  lm =lm(image[,j]~as.matrix(time) )
  image[,j] = lm$residuals
  }



ad = response$risk_for_ad
ad [ad==1 | ad ==0] = -1
ad [ad ==2 | ad ==3 ] = 1 
data = as.data.frame(cbind(image,time, ad))

library(glmnet)
cv = cv.glmnet(x= as.matrix(cbind(image,time)), y=data$ad, alpha=1, family="binomial")
plot(cv, xvar="lambda")
fit = glmnet(x= as.matrix(cbind(image,time)), y=data$ad, alpha=1, family="binomial", lambda = cv$lambda.min)
fit$beta

# now prediction of X* pick someone normal as X*, Y* at age 35 what happens to them at age 80 , star 
star=71
goal_age = 80
timestar = as.data.frame(matrix( NA , 1 , k))
for (kk in 1:k) {
  timestar[,kk] = goal_age^kk
  colnames(timestar)[kk] = paste0("t",kk)
}
datastar = cbind (data[star,1:len ], timestar)

pstar = predict(fit, as.matrix(datastar), type="response")
data$t1[star ] #their current age
pstar #probability of developing AD


## function of time for each subject
ages = seq(min(response$age), 100)
# ages = seq(60, 100)

# star=1

# plot=ggplot()
data_pstar = vector(mode = "list" , dim(connectivity)[3])
for (star in 1:dim(connectivity)[3]) {
# for (star in 1:2) {

  pstars = matrix(NA, length(ages) )
  
for (i in 1:length(ages)) {
timestar = as.data.frame(matrix( NA , 1 , k))
for (kk in 1:k) {
  timestar[,kk] = ages[i]^kk
  colnames(timestar)[kk] = paste0("t",kk)
}
datastar = cbind (data[star,1:len ], timestar)
pstar = predict(fit, as.matrix(datastar), type="response")
pstars[i] = pstar #probability of developing AD
}

data_pstar[[star]] = as.data.frame(cbind(ages,pstars))
# data_pstar[[star]] = as.data.frame(cbind(ages,pstars, rep(response$sex[star],length(ages)), rep(response$genotype[star],length(ages)), rep(response$Weight[star],length(ages)), rep(response$risk_for_ad[star],length(ages)), rep(response$age[star],length(ages))  ))
# colnames(data_pstar[[star]] ) = c("ages", "Probability" ,  "Sex", "Genotype" , "Weight" , "Risk", "Current_Age")

cat("done", star,"\n")
#p=plot(ages, pstars , xlab="age", ylab="Probability of developing MCI or AD", main=paste0("Subject ", response$Subject[star]) )

#data$t1[star ] #their current age
# ggsave(paste0('/Users/ali/Desktop/may23/risk/code/figures/',response$Subject[star],".png" ) , plot = p, device = "png")

}
# ggsave(paste0('/Users/ali/Desktop/may23/risk/code/figures/.png" ) , plot = p, device = "png")
# data_pstar = as.data.frame(data_pstar)
library(reshape2)

melt = melt(data_pstar, id.vars=c("ages") )

ggplot(melt,                            
       aes(x = ages,
           y = value,
           col = L1, group=L1)) +
  geom_line()

ggsave(paste0('/Users/ali/Desktop/may23/risk/no_age_effect/code/figures/all.png' ) , plot = last_plot(), device = "png")

melt$sex = response$sex[melt$L1]
ggplot(melt,                            
       aes(x = ages,
           y = value,
           col = sex, group=L1, alpha=0.7)) +
  geom_line()
ggsave(paste0('/Users/ali/Desktop/may23/risk/no_age_effect/code/figures/sex.png' ) , plot = last_plot(), device = "png")

melt$current_age = response$age[melt$L1]
ggplot(melt,                            
       aes(x = ages,
           y = value,
           col = current_age, group=L1, alpha= 1)) + 
   scale_colour_gradient(low = "red", high = "blue")+
  geom_line()
ggsave(paste0('/Users/ali/Desktop/may23/risk/no_age_effect/code/figures/current_age.png' ) , plot = last_plot(), device = "png")


melt$weight = response$Weight[melt$L1]
ggplot(melt,                            
       aes(x = ages,
           y = value,
           col = weight, group=L1, alpha=0.7)) +
  scale_colour_gradient(low = "red", high = "blue")+
  geom_line()
ggsave(paste0('/Users/ali/Desktop/may23/risk/no_age_effect/code/figures/weight.png' ) , plot = last_plot(), device = "png")


melt$genotype = response$genotype[melt$L1]
melt$genotype[melt$genotype == "APOE23"] = "APOE33"
melt$genotype[melt$genotype == "APOE34"] = "APOE44"

ggplot(melt,                            
       aes(x = ages,
           y = value,
           col = genotype, group=L1, alpha=1)) +
  geom_line()
ggsave(paste0('/Users/ali/Desktop/may23/risk/no_age_effect/code/figures/genotype.png' ) , plot = last_plot(), device = "png")


melt$risk = response$risk_for_ad[melt$L1]
ggplot(melt,                            
       aes(x = ages,
           y = value,
           col = as.factor(risk), group=L1, alpha=1)) +
  geom_line()
ggsave(paste0('/Users/ali/Desktop/may23/risk/no_age_effect/code/figures/risk.png' ) , plot = last_plot(), device = "png")



## edges 
edges = fit$beta[1:len]
connectivitvals = v2m(edges)

library('igraph');
connectivitvalsones=connectivitvals
connectivitvalsones[lower.tri(connectivitvalsones,diag = F)] =0
t=which(connectivitvalsones!=0, arr.ind=TRUE)
t <- cbind(t, connectivitvals[which(connectivitvals!=0,arr.ind=TRUE)]) 
t.graph=graph.data.frame(t,directed=F)
E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'red','blue') 
#t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
minC <- rep(-Inf, vcount(t.graph))
maxC <- rep(Inf, vcount(t.graph))
minC[1] <- maxC[1] <- 0
l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
                    miny=minC, maxy=maxC)      

pathnames='/Users/ali/Desktop/may23/risk/no_age_effect/code/anatomyInfo_whiston_new.csv'
datanmes=read.csv(pathnames, header = TRUE, sep = ",", quote = "")
datanmes$ROI

#noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab

#datanmes=datanmes[-noreadcsf]

datanmess=datanmes$ROI # remove csf
#datanmess=datanmes$ROI



par(mfrow=c(1,1))

#set.vertex.attribute(t.graph, "name", value=datanmes$ROI   )


png("nets.png", units="cm", width=20, height=20, res=300)  

plot(t.graph, layout=l, 
     rescale=T,
     asp=0,
     edge.arrow.size=0.1, 
     vertex.label.cex=0.8, 
     vertex.label.family="Helvetica",
     vertex.label.font=4,
     #vertex.label=t.names,
     vertex.shape="circle", 
     vertex.size=7, 
     vertex.color="white",
     vertex.label.color="black", 
     #edge.color=E(t.graph)$color, ##do not need this since E(t.graph)$color is already defined.
     edge.width=as.integer(cut(abs(E(t.graph)$V3), breaks = 5)))
dev.off()
# connectivitvals=connectivitvals+t(connectivitvals) #symetric


getwd()
filename=paste(getwd(), "/", "valandpos.mat", sep = "")
#writeMat(filename, nonzeroposition = nonzeroposition, connectivitvals = connectivitvals , oddzeroposition=indexofzeros)




subnets=groups(components(t.graph))
subnetsresults=vector(mode = "list", length = length(subnets))
colsumabs=colSums(abs(connectivitvals))
colsum=colSums(connectivitvals)

leftright=datanmes$Bigpart




####################3
for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,9,length(temp) )
  net[1,]=as.numeric(temp)
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  #indofleftright=tt>=164
  #net[5,][indofleftright]="Right"
  #net[5,][!indofleftright]="Left"
  
  
  net[2,]=datanmess[temp]
  net[5,]=leftright[temp]
  net[1,]=paste(net[1,],net[5,])
  net[3,]= as.numeric( colsum[temp]   )
  net[4,]= as.numeric( colsumabs[temp]   )
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  for (j in 1:length( net[8,])) {
    tempindex=which(datanmes$ROI %in% net[2,j]  )
    if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
  }
  net[9,]=sum(abs(as.numeric(net[3,])))/length(net[3,])
  
  subnetsresults[[i]]=net 
}

#install.packages("xlsx")
library(xlsx)


for (i in 1:length(subnetsresults)){
  net=t(subnetsresults[[i]])
  write.xlsx2(net, "nets.xlsx", sheetName =  paste0(i), append=TRUE )
}




net_new=matrix(NA, length(subnetsresults),5)


for (j in 1:dim(net_new)[1]) {
  temps=subnetsresults[[j]]
  net_new[j,1]=j
  net_new[j,2]= paste(temps[8,], collapse = ", ")
  net_new[j,3] = paste(paste(temps[5,],temps[2,]), collapse = ", ")
  net_new[j,4] = paste(temps[7,1])
  net_new[j,5] = paste(temps[6,1])
}
colnames(net_new)=c("Sub-Network", "Region Number", "Region Name", "Sub-Network Weight", "Sub_Network Average of reduced sum")


write.xlsx2(net_new, "net_new.xlsx" )








# install.packages("vioplot")
library("vioplot")


for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,8,length(temp) )
  net[2,]=datanmess[temp]
  net[1,]=as.numeric(temp)
  net[3,]= as.numeric( colsumabs[temp]   )
  net[4,]= as.numeric( colsum[temp]   )
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  indofleftright=tt>=164
  net[5,][indofleftright]="Right"
  net[5,][!indofleftright]="Left"
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  for (j in 1:length( net[8,])) {
    tempindex=which(datanmes$ROI %in% net[2,j]  )
    if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
  }
  subnetsresults[[i]]=net
}





for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,8,length(temp) )
  net[2,]=datanmess[temp]
  net[1,]=as.numeric(temp)
  net[3,]= as.numeric( colsumabs[temp]   )
  net[4,]= as.numeric( colsum[temp]   )
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  indofleftright=tt>=164
  net[5,][indofleftright]="Right"
  net[5,][!indofleftright]="Left"
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  for (j in 1:length( net[8,])) {
    tempindex=which(datanmes$ROI %in% net[2,j]  )
    if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
  }
  subnetsresults[[i]]=net 
}



#for (i in 1:length(subnetsresults)) {
#  net=subnetsresults[i]
#  print(net[[1]][1:2,])
#}


for (i in 1:length(subnetsresults)) {
  net=subnetsresults[i]
  cat( i,'th sub-net: the summation of all edges in this sub-net is' ,sum(as.numeric(net[[1]][4,])), 'and summation of absolut values of all edges in this subnet is', sum(as.numeric(net[[1]][3,])),'\n')
  cat(  'the fsirst row is the Region #, second row is the name of Region, the third row is the sum of absulote values of the edges of each region, and the last row is the sum of edges of each region \n')
  print(net)
  cat( '\n \n \n')
}


capture.output(subnetsresults, file = "subnet.txt")


write.csv(subnetsresults, row.names = T)




#install.packages("xlsx")
library(xlsx)

# 
# for (i in 1:length(subnetsresults)){
#   net=subnetsresults[[i]]
#   write.xlsx2(net, "ssir.xlsx", sheetName =  paste0(i), append=TRUE )
# }



