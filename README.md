
# Exploring how genomic architecture may present as balancing selection
## Evan Durland
### 1/24/2022
### questions? contact: durlandevan@gmail.com

## Introduction
This is a walkthrough of a simulation method to explore how contrasting genomic architectures in oysters may create changes in pooled minor allele frequency (MAF) over larval development.  

Take a scenario when you have a single genetic locus with two alleles (C/G) and it is linked to two genes with developmental importance to oysters:

....geneA....(**C**).....geneB....

....geneA....(**G**).....geneB....

Now assume that geneA becomes mutated (**x**), with detrimental effects early in development (between time 1 and 2) and it origininally showed up with the 'C' allele of the marker.  The population now has three genetic variants:

....geneA....(**C**).....geneB....

....geneA**x**....(**C**).....geneB....

....geneA....(**G**).....geneB....

in a pool, changes in the frequency of the 'C' allele reflect the death of individuals where 'C' is associated with 'geneA**x**'. Now assume a second mutation shows up for geneB which kills larvae late in development and is associated with the 'G' allele.  The population has four variants:

....geneA....(**C**).....geneB....

....geneA**x**....(**C**).....geneB....

....geneA....(**G**).....geneB....

....geneA....(**G**).....geneB**x**....

In this scenario, you can have a relative decrease in allele '**C**' early on, due to some of those mutations and then it would actually rise in frequency again as '**G**' allele individuals with the 'geneB**x**' variants die off.  

The question is, in a 5 x 19 factorial cross how frequent would these genes, alleles and alternate linkages have to be in order to produce these weird 'flip-flop' patterns of allele frequency change.  Lets simulate things to find out.  First, here's some functions for simulating parents based on starting allele frequencies, creating crosses among them and calculating how the frequency of alleles changes from time 1 to 3

```{r}
library(tidyverse)
# simulate parental genotypes:
gtyp_parents <- function(af_major, af_hidden, n_male, n_fem, print_parents = "N"){ # here we assign 'normal' associations to 'af_major' and 'hidden' mutations to 'af_hidden'
  n_par <- (n_male + n_fem)
  n_hap <-  n_par * 2 #number of haplotypes
  n_A <- ceiling(n_hap * af_major[1]) #number of 'A' alleles
  n_B <- ceiling(n_hap * af_major[2]) #number of 'B' alleles
  n_Ax <- ceiling(n_A * af_hidden[1]) #number of 'Ax' alleles
  n_Bx <- ceiling(n_B * af_hidden[2]) #number of 'Ax' alleles
  #print(paste("A: ",n_A,"B: ", n_B,"Ax: ",n_Ax, "Bx: ", n_Bx))
  allele_pool <- c(rep("A",n_A),rep("B",n_B),rep("Ax",n_Ax),rep("Bx",n_Bx))
  gtps <- c()
  for(par in 1:n_par){
    idx <- sample(1:length(allele_pool),2,replace = FALSE) #sample random alleles by position
    gtp <- paste(allele_pool[idx[1]],allele_pool[idx[2]],sep="_") #assign to a specific genotype
    gtps <- append(gtps, gtp) #append to list
    allele_pool <- allele_pool[-idx] #remove these alleles from the pool
  }
  #print(gtps)
  males <<- gtps[1:n_male] #assign first genotypes to males
  females <<- gtps[n_male+1:n_fem] #assign rest to females
  if(print_parents =="Y"){
    print("males:")
    print(males)
    print("females: ")
    print(females)
  }
  #return(males, females)
}

# simulate the crosses in a factorial design:
cross_sim<-function(dads,moms){
  genFreq<-c()
  crosses<-c()
  for(dad in 1:length(dads)){
    sire<-unlist(strsplit(dads[dad],split="_"))
    for(mom in 1:length(moms)){
      dam<-unlist(strsplit(moms[mom],split="_"))
      c1<-paste(sire[1],dam[1],sep="_")
      c2<-paste(sire[1],dam[2],sep="_")
      c3<-paste(sire[2],dam[1],sep="_")
      c4<-paste(sire[2],dam[2],sep="_")
      #add crosses to the list:
      crosses<-append(crosses,c1)
      crosses<-append(crosses,c2)
      crosses<-append(crosses,c3)
      crosses<-append(crosses,c4)
      #add all parental alleles to the list:
      genFreq<-append(genFreq,sire[1])
      genFreq<-append(genFreq,sire[2])
      genFreq<-append(genFreq,dam[1])
      genFreq<-append(genFreq,dam[2])
    }
  }
  #now, get summary stats of the pool:
  #realized frequencies of each allele:
  f_A <- length(genFreq[grep("A$",genFreq)])/length(genFreq)
  f_B <- length(genFreq[grep("B$",genFreq)])/length(genFreq)
  f_Ax <- length(genFreq[grep("Ax",genFreq)])/length(genFreq)
  f_Bx <- length(genFreq[grep("Bx",genFreq)])/length(genFreq)
  est_fA <- f_A + f_Ax
  
  #each genotype proportion:
  crosses <- gsub("B_A$","A_B",crosses)
  crosses <- gsub("Ax_A$","A_Ax",crosses)
  crosses <- gsub("B_Ax$","Ax_B",crosses)
  crosses <- gsub("Bx_A$","A_Bx",crosses)
  crosses <- gsub("Bx_Ax$","Ax_Bx",crosses)
  crosses <- gsub("Bx_B$","B_Bx",crosses)
  #crosses <<- crosses
  p_AA <- length(crosses[crosses=="A_A"])/length(crosses)
  p_BB <- length(crosses[crosses=="B_B"])/length(crosses)
  p_AxAx <- length(crosses[crosses=="Ax_Ax"])/length(crosses)
  p_BxBx <- length(crosses[crosses=="Bx_Bx"])/length(crosses)
  
  p_AB <- length(crosses[crosses=="A_B"])/length(crosses)
  p_AAx <- length(crosses[crosses=="A_Ax"])/length(crosses)
  p_ABx <- length(crosses[crosses=="A_Bx"])/length(crosses)
  p_AxB <- length(crosses[crosses=="Ax_B"])/length(crosses)
  p_BBx <- length(crosses[crosses=="B_Bx"])/length(crosses)
  p_AxBx <- length(crosses[crosses=="Ax_Bx"])/length(crosses)
  
  return(c("f_A_est" = est_fA,
           "f_A" = f_A,
           "f_B" = f_B,
           "f_Ax" = f_Ax,
           "f_Bx" = f_Bx,
           "p_AA"= p_AA,
           "p_BB"= p_BB,
           "p_AB"= p_AB,
           "p_AxAx" = p_AxAx,
           "p_BxBx" = p_BxBx,
           "p_AAx"= p_AAx,
           "p_ABx"= p_ABx,
           "p_AxB"= p_AxB,
           "p_BBx"= p_BBx,
           "p_AxBx"= p_AxBx))
}

# simulate the change in frequency given fitness (w) of each hidden allele:
sim_chg <- function(prev_gen, w_Ax, w_Bx){
  #extract genotype frequencies
  p_AA <- prev_gen[6]
  p_BB <- prev_gen[7]
  p_AB <- prev_gen[8]
  p_AxAx <- prev_gen[9]
  p_BxBx <- prev_gen[10]
  p_AAx <- prev_gen[11]
  p_ABx <- prev_gen[12]
  p_AxB <- prev_gen[13]
  p_BBx <- prev_gen[14]
  p_AxBx <- prev_gen[15]
  #calculate average fitness
  avg_fit <- sum(p_AA,p_BB,p_AB) + 
    sum(c(p_AxAx,p_AAx,p_AxB)* w_Ax) +
    sum(c(p_BxBx,p_ABx,p_BBx)* w_Bx) +
    (p_AxBx * w_Ax * w_Bx)
  #calculate change in genotype frequencies:
  n_AA <- round(p_AA/avg_fit, 4)
  n_AB <- round(p_AB/avg_fit, 4)
  n_BB <- round(p_BB/avg_fit, 4)
  n_AxAx <- round(p_AxAx * w_Ax / avg_fit, 4)
  n_BxBx <- round(p_BxBx * w_Bx / avg_fit, 4)
  n_AAx <- round(p_AAx * w_Ax / avg_fit, 4)
  n_ABx <- round(p_ABx * w_Bx / avg_fit, 4)
  n_AxB <- round(p_AxB * w_Ax / avg_fit, 4)
  n_BBx <- round(p_BBx * w_Bx / avg_fit, 4)
  n_AxBx <- round(p_AxBx * w_Ax * w_Bx / avg_fit, 4)
  #calculate new allele frequencies:
  f_A <- round(sum(n_AA,c(n_AB,n_AAx,n_ABx)*0.5),4)
  f_B <- round(sum(n_BB,c(n_AB,n_AxB,n_BBx)*0.5),4)
  f_Ax <- round(sum(n_AxAx,c(n_AAx,n_AxB,n_AxBx)*0.5),4)
  f_Bx <- round(sum(n_BxBx,c(n_ABx,n_BBx,n_AxBx)*0.5),4)
  est_fA <- f_A + f_Ax
  #print out
  return(c("f_A_est" = est_fA,
           "f_A" = f_A,
           "f_B" = f_B,
           "f_Ax" = f_Ax,
           "f_Bx" = f_Bx,
           "p_AA"= n_AA,
           "p_BB"= n_BB,
           "p_AB"= n_AB,
           "p_AxAx" = n_AxAx,
           "p_BxBx" = n_BxBx,
           "p_AAx"= n_AAx,
           "p_ABx"= n_ABx,
           "p_AxB"= n_AxB,
           "p_BBx"= n_BBx,
           "p_AxBx"= n_AxBx))
  
}

maf_cat<- function(x){
  if(x <= 0.1){
    maf_cat="0<MAF<10%"
  }
  if(0.1< x & x <=0.2){
    maf_cat = "10<MAF<20%"
  }
  if(0.2< x & x <=0.3){
    maf_cat = "20<MAF<30%"
  }
  if(0.3< x & x <=0.4){
    maf_cat = "30<MAF<40%"
  }
  if(0.4< x & x <=0.5){
    maf_cat = "40<MAF<50%"
  }
  if(x > 0.5 ){
    maf_cat = "50%<MAF"
  }
  return(maf_cat)
}
```
Next, we put these together for our scenario:
```{r}
#put it together:
scenarios <-list()
maf<-seq(from=0.05, to=0.5,by=0.02)

maf_df <- data.frame("A"=maf,
                     "B"=1-maf)
rm(maf)
for(maf in 1:nrow(maf_df)){
  af_major <- maf_df[maf,]
  for(hid in seq(1:5)/10){
    af_hidden <- c(hid,hid)
    for(a in 1:10){
      gtyp_parents(af_major, af_hidden, 5,20, print_parents = "N")
      Gen0 <- cross_sim(males,females)
      w_Ax = 0.3
      w_Bx = 1
      Gen0<-rbind(Gen0,sim_chg(Gen0,w_Ax,w_Bx))
      w_Ax = 1
      w_Bx = 0.3
      Gen0<-rbind.data.frame(Gen0,sim_chg(Gen0[2,],w_Ax,w_Bx))
      out_df <- Gen0[,1:5]
      row.names(out_df)<-c(1,2,3)
      out_df$MAF <- af_major[[1]]
      out_df$hidden <- hid
      out_df$mafcat <- maf_cat(out_df$f_A_est[1])
      scenarios[[paste(maf,hid,a)]]= out_df
      #scenarios[[paste(maf,a)]]=as.data.frame(Gen0[,1:5])
    }
  }
}
```
Finally, we will bind the list object we created into a dataframe and plot our simulations:
```{r}
df <-bind_rows(scenarios,.id = "run")
df$time <- rep(c(1,2,3),nrow(df)/3)
unique(df$mafcat)
df$mafcat <- factor(df$mafcat,levels=unique(df$mafcat))
head(df)
ggplot(df, aes(time,f_A_est))+
  geom_point(alpha=0.05)+
  geom_line(aes(group=as.factor(run)),alpha=0.05)+
  stat_smooth(aes(color=as.factor(hidden*100)),se=FALSE,size=1.5)+
  facet_wrap(~mafcat)+
  scale_x_continuous("time",breaks=c(1,2,3),labels = c(1,2,3))+
  #scale_color_discrete("Alternate vQTL frequency")+
  labs(y="minor allele frequency (estimated)",
         title=paste("f(a) trajectories for variable starting minor allele frequencies",
                       "and 10-50% disjointed vQTL abundance","fitness of vQTL=0.3",sep= "\n"))+
  scale_color_manual(paste("ALternate vQTL" , "frequency",sep= "\n"),values=c("green3","red","blue","purple","black"))+
  theme_bw()
```
![](https://github.com/E-Durland/balancing-oysters/blob/main/example.png)
So we can see that even at high abundance (MAF ~50%), and with a large proportion of alternate linkages (~50%) the overall effect is relatively little!
