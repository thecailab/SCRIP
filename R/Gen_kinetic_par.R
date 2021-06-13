
#############################################################
# Master Equation Related Functions
#############################################################

#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from 
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @param evf_res the EVFs generated for cells
#' @param is_in_module a vector of length ngenes. 0 means the gene is not in any gene module. In case of non-zero values, genes with the same value in this vector are in the same module.
#' @return a list of 3 matrices, each of dimension ngenes * nevf
GeneEffects <- function(ngenes,nevf,randseed,prob,geffect_mean,geffect_sd, evf_res, is_in_module){
  set.seed(randseed)
  
  gene_effects <- lapply(c('kon','koff','s'),function(param){
    effect <- lapply(c(1:ngenes),function(i){
      nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
      nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=geffect_mean,sd=geffect_sd)
      return(nonzero)
    })
    return(do.call(rbind,effect))
  })
  
  
  mod_strength <- 0.8
  if (sum(is_in_module) > 0){ # some genes need to be assigned as module genes; their gene_effects for s will be updated accordingly
    for (pop in unique(is_in_module[is_in_module>0])){ # for every population, find the 
      #which(evf_res[[2]]$pop == pop)
      for (iparam in c(1,3)){
        nonzero_pos <- order(colMeans(evf_res[[1]][[iparam]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[iparam]]),
                             decreasing = T)[1:ceiling(nevf*prob)]
        geffects_centric <- rnorm(length(nonzero_pos), mean = 0.5, sd=0.5)
        for (igene in which(is_in_module==pop)){
          for (ipos in 1:length(nonzero_pos)){
            gene_effects[[iparam]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
          }
          gene_effects[[iparam]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
        }
      }
      
      nonzero_pos <- order(colMeans(evf_res[[1]][[2]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[2]]),
                           decreasing = T)[1:ceiling(nevf*prob)]
      geffects_centric <- rnorm(length(nonzero_pos), mean = -0.5, sd=0.5)
      for (igene in which(is_in_module==pop)){
        for (ipos in 1:length(nonzero_pos)){
          gene_effects[[2]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
        }
        gene_effects[[2]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
      }
    }
  }
  return(gene_effects)
}

#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample,den_fun){
  probs <- den_fun$y/sum(den_fun$y)
  bw <- den_fun$x[2]-den_fun$x[1]
  bin_id <- sample(size=nsample,x=c(1:length(probs)),prob=probs,replace=T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)),function(j){
    runif(n=counts[j],min=(den_fun$x[sampled_bins[j]]-0.5*bw),max=(den_fun$x[sampled_bins[j]]+0.5*bw))
  })
  samples <- do.call(c,samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from 
#' @param bimod the bimodality constant
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params <- function(gene_effects,evf,match_param_den,bimod,scale_s){
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    # X=matrix(data=c(1:10),ncol=2) 
    # this line is to check that the row and columns did not flip
    temp <- alply(X, 1, function(Y){Y})
    values <- do.call(c,temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample=max(ranks),den_fun=match_param_den[[i]]))
    temp3 <- matrix(data=sorted[ranks],ncol=length(X[1,]),byrow=T)
    return(temp3)
  })
  
  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x}))*scale_s
  
  return(scaled_params)
}


#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param param_realdata the fitted parameter distribution to sample from 
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params2 <- function(gene_effects,evf,bimod,ranges){
  params <- lapply(gene_effects,function(X){evf %*% t(X)})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    temp <- apply(X,2,function(x){1/(1+exp(-x))})
    temp2 <- temp*(ranges[[i]][2]-ranges[[i]][1])+ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
}



#' Generating EVFs for cells sampled along the trajectory of cell development
#' @param phyla tree for cell developement
#' @param ncells number of cells
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param impulse if the impluse model should be used instead of Brownian motion
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param tip The leaf that the path with impulse lead to
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @param plotting Whether to plot the trajectory or not
#' @param plotname The string to be used in the output file name
#' @param seed the random seed 
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
ContinuousEVF <- function(phyla,ncells,n_nd_evf,n_de_evf,impulse=F,evf_center=1,vary='s',
                          Sigma,plotting=T,plotname='cont_evf.pdf',seed){
  set.seed(seed)
  edges <- cbind(phyla$edge,phyla$edge.length)
  edges <- cbind(c(1:length(edges[,1])),edges)
  connections <- table(c(edges[,2],edges[,3]))
  root <- as.numeric(names(connections)[connections==2])
  tips <- as.numeric(names(connections)[connections==1])
  internal <- as.numeric(names(connections)[connections==3])
  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)    
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)    
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }
  neutral <- SampleSubtree(root,0,evf_center,edges,ncells)
  param_names <- c("kon", "koff", "s")
  evfs <- lapply(c(1:3),function(parami){
    nd_evf <- lapply(c(1:N_ND_evfs[parami]),function(ievf){
      rnorm(ncells,evf_center,Sigma)
    })
    nd_evf <- do.call(cbind,nd_evf)
    if(N_DE_evfs[parami]!=0){
      #if there is more than 1 de_evfs for the parameter we are looking at
      if(impulse==T){
        pdf(file = plotname,width=15,height=5)
        tip <- rep(tips,ceiling(N_DE_evfs[parami]/length(tips)))
        de_evf <- lapply(c(1:N_DE_evfs[parami]),function(evf_i){
          impulse <-ImpulseEVFpertip(phyla, edges,root,tips,internal, neutral, tip[evf_i],Sigma,evf_center=evf_center)
          if(plotting==T){PlotRoot2Leave(impulse,tips,edges,root,internal)}
          re_order <- match(
            apply(neutral[,c(1:3)],1,function(X){paste0(X,collapse='_')}),
            apply(impulse[,c(1:3)],1,function(X){paste0(X,collapse='_')}))            
          return(impulse[re_order,])
        })
        dev.off()
      }else{
        de_evf <- lapply(c(1:N_DE_evfs[parami]),function(evf_i){
          SampleSubtree(root,0,evf_center,edges,ncells,neutral=neutral)    
        })
        # pdf(file = plotname,width=15,height=5)
        # if(plotting==T){PlotRoot2Leave(cbind(neutral,tips,edges,root,internal)}
        # dev.off()
      }
      de_evf <- lapply(de_evf,function(X){X[,4]})
      de_evf <- do.call(cbind,de_evf)
      de_evf <- de_evf[c(1:ncells),] # the de_evf for more cells are generated. Here we need to truncate the last ones.
      evfs <- cbind(nd_evf,de_evf)
      colnames(evfs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(de_evf[1,])),c(1:length(de_evf[1,])),sep='_'))
    }else{
      evfs <- nd_evf
      colnames(evfs)<-paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_')
    }
    return(evfs)
  })
  meta <- data.frame(pop=apply(neutral[,c(1:2)],1,function(X){paste0(X,collapse='_')}),depth=neutral[,3])
  return(list(evfs,meta[c(1:ncells),]))
  # note previously the number of sampled evfs and meta isn't necessarily ncells? 
}

#' Generating EVFs for cells sampled from tip populations from a tree
#' @param phyla tree for cell developement
#' @param ncells_total number of cells from all populations
#' @param min_popsize size of the rarest population
#' @param i_minpop to specify which population has the smallest size
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @param n_nd_evf number of non-Diff EVFs
#' @param n_de_evf number of Diff EVFs
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param evf_center the value used to generated evf means. Suggested value is 1
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the population each cell comes from (pop)
DiscreteEVF <- function(phyla, ncells_total, min_popsize, i_minpop, Sigma, n_nd_evf, n_de_evf, 
                        vary, evf_center, seed){
  set.seed(seed)
  npop <- length(phyla$tip.label)
  # set the number of cells in each population: first give each population min_popsize cells
  # then randomly distribute the rest of cells to all populations except the smallest one
  ncells_pop <- rep(min_popsize, npop)
  if (ncells_total < min_popsize*npop) {
    stop("The size of the smallest population is too big for the total number of cells")}
  larger_pops <- setdiff(1:npop, i_minpop)
  ncells_pop[larger_pops] <- floor((ncells_total-min_popsize)/length(larger_pops))
  leftover <- ncells_total-sum(ncells_pop)
  if (leftover > 0){
    temp <- sample(larger_pops, leftover, replace = F); ncells_pop[temp] <- ncells_pop[temp] + 1
  }
  
  vcv_evf_mean <- vcv.phylo(phyla,cor=T)
  param_names <- c("kon", "koff", "s")
  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)    
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)    
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }
  
  if (sum(N_DE_evfs) < 5) {warning("The number of DE evfs is less than 5; in the case of a small number of DE evfs, the structure of generated data 
	                       may not closely follow the input tree. One can either increase nevf or percent_DEevf to avoid this warning.")}
  
  evfs <- lapply(1:3, function(iparam){
    if (N_ND_evfs[iparam] > 0) {
      pop_evf_nonDE <- lapply(c(1:npop),function(ipop){
        evf <- sapply(c(1:(N_ND_evfs[iparam])),function(ievf){rnorm(ncells_pop[ipop],evf_center,Sigma)})
        return(evf)
      })
      pop_evf_nonDE <- do.call(rbind,pop_evf_nonDE)
      colnames(pop_evf_nonDE) <- rep('nonDE',N_ND_evfs[iparam])
    } else {pop_evf_nonDE <- NULL}
    if (N_DE_evfs[iparam] > 0){
      pop_evf_mean_DE <- mvrnorm(N_DE_evfs[iparam],rep(evf_center,npop),vcv_evf_mean)
      pop_evf_DE <- lapply(c(1:npop),function(ipop){
        evf <- sapply(c(1:N_DE_evfs[iparam]),function(ievf){rnorm(ncells_pop[ipop],pop_evf_mean_DE[ievf,ipop],Sigma)})
        return(evf)
      })
      pop_evf_DE <- do.call(rbind,pop_evf_DE)
      colnames(pop_evf_DE) <- rep('DE',N_DE_evfs[iparam])
    } else {pop_evf_DE <- NULL}
    
    evfs_per_param <- cbind(pop_evf_nonDE,pop_evf_DE)
    colnames(evfs_per_param) <- sprintf("%s_%s_evf%d", param_names[iparam],colnames(evfs_per_param), 
                                        1:(N_ND_evfs[iparam]+N_DE_evfs[iparam]))
    return(evfs_per_param)
  })
  meta <- data.frame(pop=do.call(c,lapply(c(1:npop),function(i){rep(i,ncells_pop[i])})))
  return(list(evfs,meta))
}

#' Generate both evf and gene effect and simulate true transcript counts
#' @param ncells_total number of cells
#' @param min_popsize the number of cells in the smallest population
#' @param i_minpop specifies which population has the smallest size
#' @param ngenes number of genes
#' @param evf_center the value which evf mean is generated from
#' @param nevf number of evfs
#' @param evf_type string that is one of the following: 'one.population','discrete','continuous'
#' @param n_de_evf number of differential evfs between populations
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param impulse use the impulse function when generating continuous population or not. Default is F. 
#' @param Sigma parameter of the std of evf values within the same population
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param param_realdata pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.
#' @param gene_effect_prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero gene effect sizes are sampled from 
#' @param gene_effect_sd the standard deviation of the normal distribution where the non-zero gene effect sizes are sampled from 
#' @param bimod the amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param gene_module_prop proportion of genes which are in co-expressed gene module
#' @param prop_hge the proportion of very highly expressed genes
#' @param mean_hge the parameter to amplify the gene-expression levels of the very highly expressed genes
#' @param randseed random seed
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
#' @export
SimulateTrueCounts <- function(ncells_total,min_popsize,i_minpop=1,ngenes, 
                               evf_center=1,evf_type="one.population",nevf=10,
                               phyla, randseed, n_de_evf=0,vary='s',Sigma=0.4,
                               geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                               bimod=0,param_realdata="zeisel.imputed",scale_s=1,impulse=F,
                               gene_module_prop=0,prop_hge=0.015, mean_hge=5){
  set.seed(randseed)
  n_nd_evf=nevf-n_de_evf
  seed <- sample(c(1:1e5),size=2)
  param_names <- c("kon", "koff", "s")
  if(evf_type=='one.population'){
    evf_mean=rep(evf_center,nevf); evf_sd=rep(Sigma,nevf)
    evfs <- lapply(1:3, function(iparam){
      evfs_per_param <- lapply(c(1:ncells_total),function(celli){
        evf <- sapply(c(1:nevf),function(ievf){rnorm(1,evf_mean[ievf],evf_sd[ievf])})
        return(evf)})
      evfs_per_param <- do.call(rbind, evfs_per_param)
      colnames(evfs_per_param) <- sprintf("%s_evf%d", param_names[iparam], 1:nevf)
      return(evfs_per_param)})
    evf_res <- list(evfs=evfs, meta=data.frame(pop=rep(1, ncells_total)))
  } else if(evf_type=='discrete'){
    evf_res <- DiscreteEVF(phyla,ncells_total,min_popsize,i_minpop=i_minpop,Sigma,
                           n_nd_evf, n_de_evf, vary=vary, evf_center=evf_center, seed=seed[1])
  }else if(evf_type=='continuous'){
    evf_res <- ContinuousEVF(phyla,ncells_total,n_nd_evf=nevf-n_de_evf,n_de_evf=n_de_evf,
                             evf_center=evf_center,vary=vary,impulse=impulse,
                             Sigma,plotting=T,seed=seed[1])    
  }
  
  if (gene_module_prop > 0 & evf_type == "continuous"){
    stop("the gene modules are developed only for evf_type to be discrete")
  }
  is_in_module <- numeric(ngenes) # assign modules automatically and return module idx for genes. 
  if (gene_module_prop > 0){
    total_module_genes <- ceiling(ngenes*gene_module_prop)
    min_module_size <- 10
    if (total_module_genes < min_module_size) {
      warning("The proportion of module genes is too low. Will not assign any module genes.")
    } else { # assign module genes to populations. The module genes for population i is expected to be highly expressed in i but not always.
      pops <- unique(evf_res[[2]]$pop); npop <- length(pops)
      nblocks <- floor(total_module_genes/min_module_size)
      nmod_genes_perpop <- numeric(npop)
      nmod_genes_perpop[1:npop] <- floor(nblocks/npop)*min_module_size
      extra_blocks <- nblocks - npop*floor(nblocks/npop)
      if ( extra_blocks > 0) {
        nmod_genes_perpop[1:extra_blocks] <- nmod_genes_perpop[1:extra_blocks] + min_module_size
      }
      extra_genes <- total_module_genes-sum(nmod_genes_perpop)
      if (extra_genes > 0){
        nmod_genes_perpop[which.min(nmod_genes_perpop)] <- nmod_genes_perpop[which.min(nmod_genes_perpop)] + extra_genes
      }
      # now assign values to is_in_module according to how many module genes are assigned to each population
      temp <- sample(ngenes); temp2 <- c(0, cumsum(nmod_genes_perpop))
      for (pop in 1:npop){
        is_in_module[temp[(temp2[pop]+1): temp2[pop+1]]] <- pop
      }
    }
  }
  
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd, evf_res=evf_res, is_in_module=is_in_module)
  
  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      data(param_realdata.zeisel.imputed)
    } else {stop("wrong input for parameter param_realdata")}
    
    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params[,3]=log(base=10,match_params[,3])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod,scale_s=scale_s)
  }else{
    params <- Get_params2(gene_effects,evf_res[[1]],bimod,ranges)
  }
  
  return(list(gene_effects=gene_effects, kinetic_params=params))
}

