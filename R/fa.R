

var_fa <- function(model){
  ASM <- MrBean:::fa.asreml( model , trunc.char = NULL)
  L.star = ASM$gammas[[1]]$`rotated loads`
  psi = ASM$gammas[[1]]$`specific var`
  VarTot = sum(diag(L.star %*% t(L.star))) / sum(diag(L.star %*% t(L.star) + diag(psi) ))

  paf.site <- ASM$gammas[[1]]$`site %vaf`

  VarGenEnv <- diag(L.star %*% t(L.star) + diag(psi) )
  TotGenVar <- sum(VarGenEnv)

  VarFA1 <- sum(VarGenEnv*paf.site[,1])/100
  VarFA2 <- sum(VarGenEnv*paf.site[,2])/100

  PerVarFA1 <- round(VarFA1/TotGenVar*100,1)
  PerVarFA2 <- round(VarFA2/TotGenVar*100,1)

  return(perc = c(PerVarFA1,PerVarFA2))
}


"circleFun" <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

require(ggrepel)

biplot_fa2 <- function(model,
                       predictions,
                       one = -1,
                       second = -1,
                       fscore = 2,
                       fmult = 10,
                       alpha = 0.3,
                       alpha_site = 0.5,
                       alpha_ind = 0.5,
                       subtitle=NULL,
                       gen = NULL,
                       size_ind_biplot = 3){

  ASM <- MrBean:::fa.asreml( model , trunc.char = NULL)
  L.star = ASM$gammas[[1]]$`rotated loads`
  L.star[,1] <- L.star[,1]*one
  L.star[,2] <- L.star[,2]*second
  psi = ASM$gammas[[1]]$`specific var`
  Gvar <- ASM$gammas[[1]]$Gmat
  Cmat <- ASM$gammas[[1]]$Cmat
  Env.means <- predictions
  names(Env.means)[1:2] <- c("site", "BLUE")
  faComp <- data.frame(site = rownames(L.star), fa1 = L.star[,1], fa2 = L.star[,2], psi = psi, Vg = diag(Gvar), BLUE = Env.means$BLUE)

  percentg <- var_fa(model)

  # Without Standardize Loadings
  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  loadings = ggplot(faComp, aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) +
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha)

  # Standardize Loadings
  faCompR <- faComp
  faCompR[,2:3] <- diag(1/sqrt(diag(Gvar))) %*% L.star
  d <- data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=faCompR[,2], vy=faCompR[,3])
  circle <- circleFun(c(0,0),2,npoints = 100)
  loading_C <- ggplot(faCompR, aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) +
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha=alpha) +
    geom_path(data=circle, aes(x,y))


  # Biplot
  fa12_scores = ASM$blups[[1]]$scores
  names(fa12_scores)[2:3] <- c("comp", "Genotype")
  fa12_scores %<>% select(-blup) %>%  spread(. ,"comp", value = "blupr")
  names(fa12_scores) = c("Genotype", "fa1", "fa2")
  fa12_scores$fa1 <- fa12_scores$fa1*one
  fa12_scores$fa2 <- fa12_scores$fa2*second
  fa12_scores$Score <-  ifelse(sqrt(fa12_scores$fa1^2+fa12_scores$fa2^2)>fscore,1,0)

  if(!is.null(gen)){

    fa12_scores[fa12_scores$Genotype %in% gen, "Score" ] <- 1

  }

  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  biplot = ggplot(faComp, aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1,  alpha = alpha_site) +
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha) +
    geom_label_repel(data = subset(fa12_scores, Score==1),
                     aes(label = Genotype, x = fmult*fa1 , y= fmult*fa2),
                     colour = "red",segment.colour = "red" , size=size_ind_biplot, alpha = alpha_ind)


  return(list(g1=loadings, g2= loading_C, g3 = biplot))
}
