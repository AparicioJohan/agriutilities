#' ran
#'
#' @param var factor to be taken as random
#'
#' @return string
#' @export
#'
#' @examples
#' # ran("rep")
ran <- function(var){
  effect <- paste0("(",1,"|",var,")")
  return(effect)
}

#' varG
#'
#' @param model lmer model
#' @param comp String: random component to be extracted
#'
#' @return
#' @export
#'
#' @examples
VarG <- function(model, comp){
  v <- as.data.frame(VarCorr(model))
  v <- v[v$grp==comp,"vcov"]
  return(v)
}

#' VarE
#'
#' @param model lmer model
#'
#' @return residual variance
#' @export
#'
#' @examples
#' # VarE(model)
VarE <- function(model){
  v <- as.data.frame(VarCorr(model))
  v <- v[v$grp=="Residual","vcov"]
  return(v)
}

#' Cullis heritability
#'
#' @param model lmer model
#' @param gen  String genotype
#'
#' @return value
#' @export
#'
#' @examples
#' # in progress
h.cullis <- function(model, gen){
  aveped <- mean(attr(ranef(model,drop=T)[[gen]],"postVar"))
  vc.g <- as.data.frame(VarCorr(model))
  vc.g <- vc.g[vc.g$grp==gen,"vcov"]
  ifelse(vc.g==0, 0 , round(1-aveped/vc.g,3) )
}

#' varG.pvalue
#'
#' @param model lmer model
#' @param gen string genotype
#'
#' @return pvalue significance of random component
#' @export
#'
#' @examples
#'  # in progress
varG.pvalue <- function(model, gen){
  table <- try(suppressWarnings(broom.mixed::tidy(lmerTest::ranova(model))), silent = T)
  if(length(class(table))==1){
    return(NA)
  } else{
    term <- grepl(gen, x = table$term)
    as.numeric(table[term,"p.value"])
  }
}

lme4_res <- function(model, return=F){
  res <-  residuals(model, scaled=TRUE)
  data <- model@frame
  data$res <- res
  data$Classify <- NA
  data$Classify[which(abs(data$res)>=3)] <- "Outlier"
  data$Classify[which(abs(data$res)<3)  ] <- "Normal"
  ix = ifelse(length(which( abs(res)>3 ))>=1, length(which( abs(res)>3 )) , 0  )
  if (return) {
    return(data)
  } else{
    return(ix)
  }
}

mult_summary <- function(models, gen = "Name", y = "response"){
  exp <- names(models)
  gv <- unlist(lapply(models, VarG, gen))
  ev <- unlist(lapply(models, VarE))
  he <- unlist(lapply(models, h.cullis, gen ))
  out <- unlist(lapply(models, lme4_res ))
  summ <- data.frame(Experiment=exp, y = y ,varG = gv, varE = ev, h2 = he, outliers=out , row.names = NULL)
  return(summ)
}

mult_lme4 <- function(data, equation, var_sub ){
  models <- list()
  data[,var_sub] <- as.factor(data[,var_sub])
  for (exp in levels(data[,var_sub])) {
    tmp_dt <- dplyr::filter(data,.data[[var_sub]]%in%exp)
    model <-  try(lmerTest::lmer(equation,data=tmp_dt, na.action = na.omit), silent = T)
    if (class(model)=="try-error") {
      models[[exp]] <- NULL
    } else{
      models[[exp]] <- model
    }
  }
  return(models)
}


lme4_BLUPs <- function(model, genotype){
  BLUPS <- ranef(model)[[genotype]]
  BLUPS <- data.frame(as.factor(row.names(BLUPS)),BLUPS[,1])
  colnames(BLUPS) <- c("Genotype","Effect")
  BLUPS <- dplyr::arrange(BLUPS,desc(Effect))
  BLUPS <- data.frame(BLUPS[,1],round(BLUPS[,2],2))
  names(BLUPS) <- c("Line","BLUPs")
  d <- broom.mixed::augment(ranef(model))
  d <- d[d$grp==genotype,c("level","std.error")]
  d <- data.frame(level=d[,1],std.error=round(d[,2],2))
  BLUPS <- merge(BLUPS,d,by.x="Line",by.y="level")
  BLUPS
}


lme4_plotly <- function(blups){
  BLUPS <- blups
  BLUPS$Lu <- BLUPS[,2]-1.645*BLUPS[,3]
  BLUPS$Ls <- BLUPS[,2]+1.645*BLUPS[,3]
  v <- as.character(BLUPS[order(BLUPS[,2],decreasing = TRUE),1])
  names(BLUPS)[2] <- "predicted.value"
  p <- ggplot(BLUPS,aes(x=Line , y=predicted.value))+
    geom_point(size = 1) +
    geom_errorbar(aes(ymax = Ls, ymin = Lu))+
    theme_bw() +
    geom_hline(yintercept = mean(BLUPS[,2]), linetype=2 ,color="red")+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    ylab(names(BLUPS)[2])+scale_x_discrete(limits=v)
  plotly::ggplotly(p)
}


lme4_ggplot <- function(blups, title = NULL , subtitle = NULL){
  BLUPS <- blups
  BLUPS$Lu <- BLUPS[,2]-1.645*BLUPS[,3]
  BLUPS$Ls <- BLUPS[,2]+1.645*BLUPS[,3]
  v <- as.character(BLUPS[order(BLUPS[,2],decreasing = TRUE),1])
  names(BLUPS)[2] <- "predicted.value"
  p <- ggplot(BLUPS,aes(x=Line , y=predicted.value))+
    geom_point(size = 1) +
    geom_errorbar(aes(ymax = Ls, ymin = Lu))+
    theme_bw() +
    geom_hline(yintercept = mean(BLUPS[,2]), linetype=2 ,color="red")+
    theme_ipsum(base_size = 10) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 70, hjust = 1),
          axis.ticks.x=element_blank())+
    labs(x="", y=names(BLUPS)[2],
         title=title,
         subtitle=subtitle) +
    scale_x_discrete(limits=v)
  return(p)
}



res_data_lme4 <- function(Model){
  if(class(Model)=="lm"){
    Data <- Model$model
    VarE <- sigma(Model)^2
  } else {
    Data <- Model@frame
    VarE <- VarE(Model)
  }
  Data$Index <- 1:length(residuals(Model))
  Data$Residuals <- residuals(Model)
  u <- +3*sqrt(VarE)
  l <- -3*sqrt(VarE)
  Data$Classify <- NA
  Data$Classify[which(abs(Data$Residuals)>=u)] <- "Outlier"
  Data$Classify[which(abs(Data$Residuals)<u)  ] <- "Normal"
  Data$l <- l
  Data$u <- u
  Data$fit <-  fitted.values(Model)
  return(Data)
}

res_qqplot <- function(data_out, title = NULL){
  q <- dplyr::filter(data_out,!is.na(Classify)) %>%
    ggpubr::ggqqplot(x="Residuals",
                     fill="Classify",
                     ggtheme=theme_ipsum(),
                     ylab = "Sample Quantile",
                     xlab = "Theoretical Quantile", title =title )
  return(q)
}

