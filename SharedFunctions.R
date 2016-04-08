## Shared functions between R markdown files



SimulateChromatid <- function(num.markers=10, recomb.prob=0.3) {
  paternal<-sample(c(0,1),num.markers,rep=TRUE,prob=c(0.5,0.5))
  maternal<-sample(c(0,1),num.markers,rep=TRUE,prob=c(0.5,0.5))
  recombination.event<-sample(c(TRUE,FALSE),num.markers,prob=c(recomb.prob,1-recomb.prob),rep=TRUE)
  ## creating empty vectors  
  child<-c();
  parent.path<-c();
  current.genome<-sample(c('paternal','maternal'),1,prob=c(0.5,0.5))
  ## main loop to generate child chromatid genotype
  for (base.index in 1:length(paternal)) {
  
    if (recombination.event[base.index]){ #switching current parental genome
      if(current.genome=='paternal'){
        current.genome<-'maternal'
      } else {
          current.genome<-'paternal'
      }      
    }   
   
    if(current.genome=='paternal'){
       child[base.index]<-paternal[base.index]
    } else {
       child[base.index]<-maternal[base.index]	
    }
    parent.path[base.index]<-current.genome
  }
  return(list(maternal=maternal,paternal=paternal,parent.path=parent.path,child=child))
}


GetEmissionProb<-function(mom.maternal,mom.paternal,dad.maternal,dad.paternal){
  row.1<-c(as.numeric(mom.maternal+dad.maternal==0),
           as.numeric(mom.maternal+dad.maternal==1),
           as.numeric(mom.maternal+dad.maternal==2))
  row.2<-c(as.numeric(mom.maternal+dad.paternal==0),
           as.numeric(mom.maternal+dad.paternal==1),
           as.numeric(mom.maternal+dad.paternal==2)) 
  row.3<-c(as.numeric(mom.paternal+dad.maternal==0),
           as.numeric(mom.paternal+dad.maternal==1),
           as.numeric(mom.paternal+dad.maternal==2))
  row.4<-c(as.numeric(mom.paternal+dad.paternal==0),
           as.numeric(mom.paternal+dad.paternal==1),
           as.numeric(mom.paternal+dad.paternal==2))      
  
  emission.matrix<-rbind(row.1,row.2,row.3,row.4)
  colnames(emission.matrix)<-c('0','1','2')
  rownames(emission.matrix)<-c('mom.maternal/dad.maternal','mom.maternal/dad.paternal',
                               'mom.paternal/dad.maternal','mom.paternal/dad.paternal')
  
  return(emission.matrix)
}


GetTransitionProb<-function(r){
  # r is the recombination rate
  transition.matrix<-matrix(rep(r/3,16),ncol=4);
  transition.matrix[1,1]=1-r;
  transition.matrix[2,2]=1-r;
  transition.matrix[3,3]=1-r;
  transition.matrix[4,4]=1-r;
  rownames(transition.matrix)<-c('mom.maternal/dad.maternal','mom.maternal/dad.paternal',
                                 'mom.paternal/dad.maternal','mom.paternal/dad.paternal')
  colnames(transition.matrix)<-rownames(transition.matrix)      
  return(transition.matrix)
}




Viterbi2 <- function(obs,start.p, trans.p,mom.chromatid,dad.chromatid) {
  #initializing
  v <- matrix(NA, nr=length(obs), nc=dim(trans.p)[1])
  emit.p<-GetEmissionProb(mom.chromatid$maternal[1],
                          mom.chromatid$paternal[1],
                          dad.chromatid$maternal[1],
                          dad.chromatid$paternal[1])
  for(i in 1:4){
    v[1,i]=start.p[i]*emit.p[i,obs[1]+1]
  }
  
  for(i in 2:length(obs)) {# from observation 2 to t
    emit.p<-GetEmissionProb(mom.chromatid$maternal[i],
                            mom.chromatid$paternal[i],
                            dad.chromatid$maternal[i],
                            dad.chromatid$paternal[i])
    for (l in 1:dim(trans.p)[1]) { #from state 1 to N (4 states)
      v[i,l] <- emit.p[l,obs[i]+1] * max(v[(i-1),] * trans.p[l,])
    }
  }   
   return(v) 
}

Viterbi3 <- function(obs,start.p, trans.p,mom.chromatid,dad.chromatid) {
  #initializing
  v <- matrix(NA, nr=length(obs), nc=dim(trans.p)[1])
  emit.p<-GetEmissionProb(mom.chromatid$maternal[1],
                          mom.chromatid$paternal[1],
                          dad.chromatid$maternal[1],
                          dad.chromatid$paternal[1])
  emit.p<-emit.p+10^(-12)                        
  for(i in 1:4){
    v[1,i]=log(start.p[i],2)+log(emit.p[i,obs[1]+1],2)
  }
  
  for(i in 2:length(obs)) {# from observation 2 to t
     emit.p<-GetEmissionProb(mom.chromatid$maternal[i],
                             mom.chromatid$paternal[i],
                             dad.chromatid$maternal[i],
                             dad.chromatid$paternal[i])
     for (l in 1:dim(trans.p)[1]) { #from state 1 to N (4 states)
       v[i,l] <- log(emit.p[l,obs[i]+1]+10^(-12)) + max(v[(i-1),] + log( trans.p[l,],2))
     }
  }   
  
  return(v) 
}

ViterbiSimulator <- function(num.markers,num.simulations,start.p,recombination.rate){ 
  sim.results<-c()
  trans.p<-GetTransitionProb(recombination.rate)
  for (i in 1:num.simulations){
    child1.chromatid.mom<-SimulateChromatid(num.markers,recombination.rate) #100 markers
    child1.chromatid.dad<-SimulateChromatid(num.markers,recombination.rate)
    obs<-child1.chromatid.mom$child+child1.chromatid.dad$child
    v.2<-Viterbi3 (obs,start.p, trans.p,child1.chromatid.mom,child1.chromatid.dad)
    vitrowmax.2 <- apply(v.2, 1, function(x) which.max(x))
    simulated<-paste("mom.",child1.chromatid.mom$parent.path,"/","dad.",child1.chromatid.dad$parent.path,sep="")
    estimated<-rownames(trans.p)[vitrowmax.2]
    sim.results[i]<-sum(simulated==estimated)
  }
  return(sim.results)
}

ViterbiSimulator2 <- function(num.markers,num.simulations,start.p,recombination.rate){ 
  sim.results<-c()
  sim.recomb.events<-c()
  est.recomb.events<-c()
  trans.p<-GetTransitionProb(recombination.rate)
  for (i in 1:num.simulations){
    child1.chromatid.mom<-SimulateChromatid(num.markers,recombination.rate) #100 markers
    child1.chromatid.dad<-SimulateChromatid(num.markers,recombination.rate)
    obs<-child1.chromatid.mom$child+child1.chromatid.dad$child
    v.2<-Viterbi3 (obs,start.p, trans.p,child1.chromatid.mom,child1.chromatid.dad)
    vitrowmax.2 <- apply(v.2, 1, function(x) which.max(x))
    simulated<-paste("mom.",child1.chromatid.mom$parent.path,"/","dad.",child1.chromatid.dad$parent.path,sep="")
    estimated<-rownames(trans.p)[vitrowmax.2]
    sim.results[i]<-sum(simulated==estimated)
    sim.recomb.events[i]<-CountRecombEvents(simulated);
    est.recomb.events[i]<-CountRecombEvents(estimated);
  }
  return(list(sim.results=sim.results,sim.recomb.events=sim.recomb.events,est.recomb.events=est.recomb.events))
}  
 

CountRecombEvents<-function(parent.path){
  current<-parent.path[1];
  counter<-0;
  for(val in parent.path){
    if(val!=current){
      counter<-counter+1
      current<-val
    }
  }
  return (counter)
} 
