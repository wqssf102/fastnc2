##此代码copy于micrntwk包，作者：吴求生，邮箱：565715597@qq.com
##完整的代码在：https://github.com/wqssf102/micrntwk
#####
calcor <- function(codt=codt,r_th=0.6,p_th=0.05,type=c("cor","mine"),n.cores=4,group=NULL){
  type <- match.arg(type, several.ok = TRUE)
  type <- unique(type)
  ######spearman
  calcora <- function(codt=codt,r_th=0.6,p_th=0.05){
    correst <-  Hmisc::rcorr(as.matrix(codt), type = 'spearman')
    ###R value
    cor_R  <-  correst$r %>% as.data.frame()
    diag(cor_R) <- NA
    cor_R[upper.tri(cor_R)] <- NA
    nrr <- as.data.frame(cbind(rownames(cor_R),cor_R))
    names(nrr)[1] <- "spnamer"
    longr <- reshape2::melt(nrr)
    nlongr <- na.omit(longr)###remove NA
    
    ###P value
    cor_p0 <-  correst$P
    myp <- reshape2::melt(cor_p0)
    #P adjust：c("holm", "hochberg", "hommel",
    ##"bonferroni", "BH", "BY","fdr", "none")
    cor_p1 <- p.adjust(cor_p0, method ="fdr") %>% as.data.frame()
    colnames(cor_p1) <- "pzhi"
    ####
    newdat <- cbind(hang=myp[,1],lie=myp[,2],cor_p1) %>% as.data.frame()
    newdat$pzhi <- round(newdat$pzhi,digits = 4)
    newdat0 <- reshape2::dcast(newdat,hang~lie,value.var = 'pzhi')
    rownames(newdat0) <- newdat0[,1]
    php <- newdat0[,-1] %>% as.data.frame()
    php[upper.tri(php)] <- NA
    npp <- as.data.frame(cbind(rownames(php),php))
    names(npp)[1] <- "spnamep"
    longp <- reshape2::melt(npp)
    nlongp <- na.omit(longp)###
    names(nlongp)[3] <- "p_value"
    ##R P merge
    res_dt <- as.data.frame(cbind(nlongr,nlongp$p_value))
    names(res_dt) <- c("Source","Target","R","P")
    #########
    res_dt$R <- round(res_dt$R,digits=2)
    res_dt$P <- round(res_dt$P,digits=2)
    res_dt <- res_dt[-(which((res_dt[,4]>p_th)|abs(res_dt[,3])<r_th)),]
    res_dt_gg <- res_dt[,c("Source","Target","R")]
    res_dt_gg$edcol <- ifelse(res_dt_gg$R>0,1,-1)
    gg <- graph.data.frame(res_dt_gg,directed = FALSE)
    E(gg)$weight <- abs(E(gg)$R)
    # gg <- simplify(gg,remove.multiple = F)
    resgg <- list(gg,res_dt_gg)
    names(resgg) <- c("gg","corgg")
    set.seed(123456)
    resgg$modularity <- modularity(cluster_fast_greedy(resgg$gg))
    return(resgg)
  }
  #####
  
  ####MIC
  calcorb <- function(codt=codt,r_th=0.6,n.cores=8){
    message("The mine method has no P value")
    codt <- codt[,colSums(codt)!=0]
    res <- minerva::mine(codt,n.cores=n.cores,alpha = 0.6,var.thr=1e-20)
    mc <- res$MIC
    mc <- reshape2::melt(mc)
    mc <- mc[(mc$Var1)!=(mc$Var2),]
    mc <- mc[mc$value>r_th,]
    gg <- graph.data.frame(mc[,c(1:2)],directed = FALSE)
    E(gg)$weight <- mc$value
    gg <- simplify(gg)
    resgg <- list(gg,mc)
    names(resgg) <- c("gg","corgg")
    set.seed(123456)
    resgg$modularity <- modularity(cluster_fast_greedy(resgg$gg))
    return(resgg)
  }
  ##################################
  ########
  if(is.null(group)){
    if(type=="cor"){
      mycormt <- calcora(codt=codt,r_th=r_th,p_th=p_th)
    }else{
      mycormt <- calcorb(codt=codt,r_th=r_th,n.cores=n.cores)
    }
    rescor <- mycormt
  }else{
    if(type=="cor"){
      codt <- split(codt,group$group)
      rescor <- lapply(codt,calcora,r_th=r_th,p_th=p_th)
    }else{
      codt <- split(codt,group$group)
      rescor <- lapply(codt,calcorb,r_th=r_th,n.cores=n.cores)
    }
  }
  rescor$group <- group
  return(rescor)
}
