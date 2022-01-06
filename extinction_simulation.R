
source("scripts/R/functional_redundancy.R")
library(tidyverse)
library(doParallel)



filename='data/simulations/extinction_simulations.csv'
n.express<-c(1:100)#possible n for trait possessing taxa
n.nonexpress<-(0:100) #possible n for taxa not possessing a trait
q<-0.5 #diversity order
a.range<-seq(0,0.4,by=0.01) # a values range for defining the evenness of a lognormal community
inc<-1#resolution of trait level
trait.min<-1 #minimum trait level if a trait possessing taxa (normalized later)
trait.max<-1000 #maximum trait level if a trait possessing taxa (normalized later)
it<-1000 #number of artificial communities
it2<-100 #number of extinction simulations
S0<-1 #maximum value for abundance distribution (normalized later)
loss.thres<-0.99 #proportion of trait needed to be removed for extinction simulation to end
ncore=3#number of coresd for parallelization

cl<-makeCluster(2) #create cluster object
registerDoParallel(cores = ncore) #assign number of cores

tmp.trait.sample<-seq(trait.min,trait.max,by=inc) #vector containing possible trait levels

#parallelize extinction simulations
community.trait<-foreach (i=1:it, .combine = rbind, .packages=c('e1071','tidyverse')) %dopar% {
  
  n.e<-sample(n.express,1) #randomly choose number of trait possesing taxa
  n.n<-sample(n.nonexpress,1) #randomly choose taxa not possessing trait
  n.tot=n.e+n.n #calculate community richness
  
  #create data.frame representing temporary community
  d_tmp_community<-data.frame(relative.abundance=rep(NaN,n.tot),
                              trait=rep(NaN,n.tot),
                              # scaled.trait10=rep(NaN,n),
                              scaled.trait=rep(NaN,n.tot),
                              member=1:n.tot)
  
  R<-c(0:(n.tot-1)) #octaves
  a<-sample(a.range,1) #choose a parameter for generating lognormal community
  SR<-S0*exp(-a^2*R^2) #generate abundances of lognormal community
  SR<-SR/sum(SR) #convert to proportions
  d_tmp_community$relative.abundance<-SR #assign relative abundances to taxa in community data.frame
  
  trait.rand<-sample(tmp.trait.sample,n.e,replace = TRUE) #randomly choose, with replacement, trait levels

  trait.rand<-sample(c(trait.rand,rep(0,n.n))) #randomize the order of trait levels; random trait vector is concatenated with a 0's vector of length n.n
  d_tmp_community$trait<-trait.rand #assign trait levels
  d_tmp_community$scaled.trait<-d_tmp_community$relative.abundance*d_tmp_community$trait #caluclate tau
  
  
  tmp_trait_scaled<-d_tmp_community %>%
    dplyr::select(-trait,-relative.abundance) %>%
    mutate(scaled.trait=scaled.trait/sum(scaled.trait),
              member=member) %>%
    spread(key="member",value="scaled.trait") #convert artifical community into long format, used by royalty_functional_redundancy
 

  
  d.trait.loss<-data.frame(member.loss=rep(NaN,it2)) #a dataframe for determining the proportion of a community that goes extinct

  #perform extinction simulations,it2 replicates 
  for (k in 1:it2){
    trait.loss<-0 #intial trait loss from community aggregated parameter
    tmp_trait_scaled.express<-tmp_trait_scaled #create dummy data.frame which updates after each iteration
    member.loss<-0
   while(trait.loss<loss.thres){ #keep looping until trait level removed exceeds threshold
      indx<-sample(1:length(tmp_trait_scaled.express),1) #randomly select community member
      trait.loss<-trait.loss+tmp_trait_scaled.express[indx] #calculate total trait level removed from community aggregated parameter
      tmp_trait_scaled.express<-tmp_trait_scaled.express[-indx] #remove member from community
      member.loss<-member.loss+1 #monitor number of extinctions
    }
    d.trait.loss$member.loss[k]<-member.loss #assign total members loss during replicate to data.frame 
  }


  tmp_community<-d_tmp_community %>%
    dplyr::select(-trait,-scaled.trait) %>%
    spread(key="member",value="relative.abundance") #relative abundances, long format, for ricotta method
  
  tmp_trait<-d_tmp_community %>%
    dplyr::select(-relative.abundance,-scaled.trait) %>%
    mutate(trait=trait/max(trait)) %>%
    spread(key="member",value="trait") #trait levels, long format, for ricotta method
  
 
  Q.rao<-ricotta_fr(tmp_community,t(tmp_trait)) #calculate functional redundancy using rao quadratic diversity, ricotta method
 
  fr.fraction<-length(tmp_trait_scaled[tmp_trait_scaled>0])/n.tot #calculate fraction of community with trait
  fr.hester<-hester_functional_redundancy(tmp_trait_scaled) #calculate median overlap in pathways with trait
  fr.mouillot<-n.tot/2 #divide richness by number of functional entities--absent/present results in two FE
  fr.ricotta<-Q.rao$R  
  fr.royalty<-royalty_functional_redundancy(tmp_trait_scaled,q)[1,1] #caluclate functional redundancy using method here
  
  loss.frac=mean(d.trait.loss$member.loss)/n.tot #calculate proportion of community that went extinct
  tmp_community_out<-data.frame(fr=c(fr.royalty,fr.ricotta,fr.mouillot,fr.hester,fr.fraction),
                                method=c("Here; q=0.5","Ricotta","Mouillot","Hester","Here; q=0"),
                                loss.frac=rep(loss.frac,5),
                                express.level=rep(n.e/n.tot,5)) #create data.frame used for appending replicate simulations together, express level is proportion of community with a trait
  
}

stopCluster(cl)#stop cluster object

ggplot()+
  geom_point(data=community.trait,aes(x=fr,y=loss.frac,color=express.level))+
  geom_smooth(data=community.trait,aes(x=fr,y=loss.frac),method="gam")+
  scale_color_distiller(palette = "RdYlBu")+
  facet_wrap(.~method,scales="free") +
  ylim(0.5,1)+
  ylab('Trait Stability')+
  xlab('Functional Redundancy')
#write model outputs
write.csv(community.trait,file =filename,quote = FALSE,row.names = FALSE)
