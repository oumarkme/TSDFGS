####################################################
#                             \-- r-score (RPT)    # 
#              \-- target ----\-- pev-score (PPT)  # 
#              \              \-- CD-score (CPT)   #
#   \-- pop ---\                                   #
#   \          \              \-- r-score(RPN)     #
#   \          \-- untarget --\-- pev-score(PPN)   #
#   \                         \-- CD-score(CPN)    #
#   \                                              #
#   \                         \-- r-score (RNT)    #
#   \          \-- target ----\-- pev-score (PNT)  #
#   \          \              \-- CD-score (CNT)   #
#   \-- npop --\                                   #
#              \              \-- r-score (RNN)    #
#              \-- untarget --\-- pev-score(PNN)   #
#                             \-- CD-score (CNN)   #
####################################################

#' Optimal training set determination
#' 
#' This function is designed for determining optimal training set.
#' 
#' @author Jen-Hsiang Ou
#' 
#' @param geno A numeric matrix of principal components (rows: individuals; columns: PCs).
#' @param cand An integer vector of which rows of individuals are candidates of the training set in the geno matrix.
#' @param n.train The size of the target training set. This could be determined with the help of the ssdfgp function provided in this package.
#' @param subpop A character vector of sub-population's group name. The algorithm will ignore the population structure if it remains NULL.
#' @param test An integer vector of which rows of individuals are in the test set in the geno matrix. The algorithm will use an un-target method if it remains NULL.
#' @param method Choices are rScore, PEV and CD. rScore will be used by default.
#' @param min.iter Minimum iteration of all methods can be appointed. One should always check if the algorithm is converged or not. A minimum iteration will set by considering the candidate and test set size if it remains NULL.
#' 
#' @return This function will return 3 information including OPTtrain (a vector of chosen optimal training set), TOPscore (highest scores of before iteration), and ITERscore (criteria scores of each iteration).
#' 
#' @export
#' @examples
#' data(geno)
#' \dontrun{optTrain(geno, cand = 1:404, n.train = 100)}
#' 
optTrain = function(geno, cand, n.train, subpop=NULL, test=NULL, method="rScore", min.iter=NULL)
{  
  if(!method%in%c("rScore","PEV","CD")){stop("Method not found. Please choose one from (rScore, PEV, CD)")}
  n=n.train; N=nrow(geno); Nc=length(cand)
  geno = as.matrix(geno)
  
  if(is.null(subpop)){
    if(is.null(test)){
      if(method=="PEV"){
        ## npop untarget pev
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*33)}
        sol = sample(cand,n)
        score = pev_score(geno[sol,], geno[cand[!cand%in%sol],])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          new.sol[sample(n,1)] = sample(cand[!cand%in%sol],1)
          new.score=pev_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          if(new.score < score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score, new.score)
            top.score = c(top.score, new.score)
          }else{
            iter.score=c(iter.score, new.score)
            top.score=c(top.score, score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }else if(method=="CD"){
        ## npop untarget CD
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*65)}
        sol = sample(cand,n)
        score = cd_score(geno[sol,], geno[cand[!cand%in%sol],])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          new.sol[sample(n,1)] = sample(cand[!cand%in%sol],1)
          new.score=cd_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score, new.score)
            top.score = c(top.score, new.score)
          }else{
            iter.score=c(iter.score, new.score)
            top.score=c(top.score, score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }else{
        ## npop untarget r-score
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*80)}
        sol = sample(cand,n)
        score = r_score(geno[sol,], geno[cand[!cand%in%sol],])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          new.sol[sample(n,1)] = sample(cand[!cand%in%sol],1)
          new.score=r_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score, new.score)
            top.score = c(top.score, new.score)
          }else{
            iter.score=c(iter.score, new.score)
            top.score=c(top.score, score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }
    }else{
      if(method=="PEV"){
        ## npop target pev
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*65)}
        sol=matrix(NA, n, 30); for(i in 1:30){sol[,i] = sample(cand, n)}
        score=rep(NA,30); for(i in 1:30){score[i]=pev_score(geno[sol[,i],], geno[test,])}
        iter.score=score; top.score=min(score)
        stop=0; iter=1
        while(stop==0){
          elite = which(rank(score)<4)
          del = sample(seq(30)[-elite], 3, prob=(score[-elite]/sum(score[-elite])))
          for(i in 1:length(del)){
            par = sample(seq(30)[-del],2)
            sol[,del[i]] = sample(unique(c(sol[,par[1]],sol[,par[2]])),n)
          }
          for(i in 1:30){
            new.sol=sol[,i]
            if(i %in% del){
              sol[sample(n,ceiling(n*.03)),i]=sample(cand[!cand%in%sol[,i]],ceiling(n*.03))
              score[i] = pev_score(geno[sol[,i],], geno[test,])
            }else{
              new.sol[sample(n,ceiling(n*.03))]=sample(cand[!cand%in%new.sol],ceiling(n*.03))
              new.score = pev_score(geno[new.sol,], geno[test,])
              if(new.score < score[i]){sol[,i] = new.sol; score[i]=new.score}
            }
          }
          iter.score=c(iter.score,mean(score)); top.score=c(top.score, min(score))
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
        sol = sol[,which(score==min(score))[1]]
      }else if(method=="CD"){
        ## npop tagget cd
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*40)}
        sol = sample(cand, n)
        score = cd_score(geno[sol,], geno[test,])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          new.sol[sample(n,1)] = sample(cand[!cand%in%new.sol],1)
          new.score=cd_score(geno[new.sol,], geno[test,])
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score, new.score)
            top.score=c(top.score,new.score)
          }else{
            iter.score=c(iter.score, new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-10){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }else{
        ## npop target r-score
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*40)}
        sol=matrix(NA, n, 30); for(i in 1:30){sol[,i] = sample(cand, n)}
        score=rep(NA,30); for(i in 1:30){score[i]=r_score(geno[sol[,i],], geno[test,])}
        iter.score=score; top.score=max(score)
        stop=0; iter=1
        while(stop==0){
          elite = which(rank(score)>27)
          del = sample(seq(30)[-elite], 3, prob=((1/score[-elite])/sum(1/score[-elite])))
          for(i in 1:length(del)){
            par = sample(seq(30)[-del],2)
            sol[,del[i]] = sample(unique(c(sol[,par[1]],sol[,par[2]])),n)
          }
          for(i in 1:30){
            new.sol=sol[,i]
            if(i %in% del){
              sol[sample(n,ceiling(n*.03)),i]=sample(cand[!cand%in%sol[,i]],ceiling(n*.03))
              score[i] = r_score(geno[sol[,i],], geno[test,])
            }else{
              new.sol[sample(n,ceiling(n*.03))]=sample(cand[!cand%in%new.sol],ceiling(n*.03))
              new.score = r_score(geno[new.sol,], geno[test,])
              if(new.score > score[i]){sol[,i] = new.sol; score[i]=new.score}
            }
          }
          iter.score=c(iter.score,mean(score)); top.score=c(top.score, max(score))
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
        sol = sol[,which(score==max(score))[1]]
      }
    }
  }else{
    subpop=as.character(subpop)
    if(is.null(test)){
      if(method=="PEV"){
        ## pop untarget pev
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
        if(length(cand)!=length(subpop)){stop("Input data is not correct.")}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*(as.numeric(table(subpop))/sum(as.numeric(table(subpop)))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{stop=1}
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop==pops[i]],pop.ratio[i]))
        }
        score=pev_score(geno[sol,], geno[cand[!cand%in%sol],])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1)
          new.sol[sample(which(subpop[sol]==pops[p]),1)]=sample(which(subpop==pops[p])[!which(subpop==pops[p])%in%sol],1)
          new.score=pev_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          if(new.score < score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,new.score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
        
      }else if(method=="CD"){
        ## pop untarget cd
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
        if(length(cand)!=length(subpop)){stop("Input data is not correct.")}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*(as.numeric(table(subpop))/sum(as.numeric(table(subpop)))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{stop=1}
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop==pops[i]],pop.ratio[i]))
        }
        score=cd_score(geno[sol,], geno[cand[!cand%in%sol],])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1)
          new.sol[sample(which(subpop[sol]==pops[p]),1)]=sample(which(subpop==pops[p])[!which(subpop==pops[p])%in%sol],1)
          new.score=cd_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,new.score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-10){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }else{
        ## pop untarget r-score
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*70)}
        if(length(cand)!=length(subpop)){stop("Input data is not correct.")}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*(as.numeric(table(subpop))/sum(as.numeric(table(subpop)))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{stop=1}
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop==pops[i]],pop.ratio[i]))
        }
        score=r_score(geno[sol,], geno[cand[!cand%in%sol],])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1)
          new.sol[sample(which(subpop[sol]==pops[p]),1)]=sample(which(subpop==pops[p])[!which(subpop==pops[p])%in%sol],1)
          new.score=r_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,new.score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }
    }else{
      if(method=="PEV"){
        
        ## pop target pev
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*as.numeric(table(subpop))/sum(as.numeric(table(subpop))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{
            stop=1
          }
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
        }
        score=pev_score(geno[sol,], geno[test,])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1,prob=pop.ratio)
          new.sol[as.numeric(sample(as.character(which(subpop[new.sol]==pops[p])),1))]=cand[sample(which(subpop[cand]==pops[p]),1)]
          new.score=pev_score(geno[new.sol,],geno[test,])
          if(new.score<score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,score)
            top.score=c(top.score,score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }else if(method=="CD"){
        ## pop target cd
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*60)}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*as.numeric(table(subpop))/sum(as.numeric(table(subpop))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{
            stop=1
          }
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
        }
        score=cd_score(geno[sol,], geno[test,])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1,prob=pop.ratio)
          new.sol[as.numeric(sample(as.character(which(subpop[new.sol]==pops[p])),1))]=cand[sample(which(subpop[cand]==pops[p]),1)]
          new.score=cd_score(geno[new.sol,],geno[test,])
          if(new.score>score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,score)
            top.score=c(top.score,score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }else{
        ## pop target r-score
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*60)}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*as.numeric(table(subpop))/sum(as.numeric(table(subpop))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{
            stop=1
          }
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
        }
        score=r_score(geno[sol,], geno[test,])
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1,prob=pop.ratio)
          new.sol[as.numeric(sample(as.character(which(subpop[new.sol]==pops[p])),1))]=cand[sample(which(subpop[cand]==pops[p]),1)]
          new.score=r_score(geno[new.sol,],geno[test,])
          if(new.score>score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,score)
            top.score=c(top.score,score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
      }
    }
  }
  ret = list(
    OPTtrain=as.numeric(sol),
    TOPscore=as.numeric(top.score[-1]),
    ITERscore=as.numeric(iter.score[-1])
  )
  return(ret)
}



