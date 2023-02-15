#! /usr/bin/env Rscript
# (c) Konstantin Riege

suppressMessages({
  library(survival)
  library(survminer)
  library(stringr)
})
args = commandArgs(TRUE);

odir = args[4];
cat(paste0("about to generate survival curves in ",odir,"\n"))

clinical = read.table(args[1], header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="")
max_time = as.double(args[2])
scores = read.table(args[3], header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="")
clinical = merge(clinical, scores, by="barcode")

covariates = c()
if(length(args)>4){
  covariates = args[5:length(args)]
}

if (nrow(clinical)==0){
  cat("empty clinical file. make sure clinical and scores patients do match.\n", file=stderr())
  quit(status=1)
}

clinical = clinical[! (is.na(clinical$days_to_death) & is.na(clinical$days_to_last_follow_up)) , ]
clinical = cbind(clinical,t(apply(clinical,1,function(l) {
  time_orig = max(as.numeric(l["days_to_death"]),as.numeric(l["days_to_last_follow_up"]),na.rm=T)
  event_orig = as.numeric(l["vital_status"]=="Dead")
  time_revised = min(time_orig,max_time)
  if (time_revised<time_orig) event_revised=0 else event_revised=event_orig
  return(c(time_orig=time_orig, event_orig=event_orig, time_revised=time_revised, event_revised=event_revised))
})))

survplot = function(path){
  groups = levels(clinical$group)

  # check factors - important for single tumor plots
  coeffs=c()
  for (c in covariates){
    v=clinical[,c]
    if (length(levels(as.factor(v))) > 1){
      coeffs = c(coeffs,c)
    }
  }
  f = as.formula(paste("surv~",paste(c("group",coeffs),collapse=" + ")))

  # cox fit
  cox = coxph(f,clinical)
      # new = data.frame(group = groups)
      # fitx = survfit(cox,newdata=new,data=clinical)
      # tsvx=summary(fitx)
      # tsvx=data.frame(cbind(tsvx$surv,tsvx$time))
      # colnames(tsvx)=c(groups,"time")
  cox=summary(cox)
  # OR KM fit
  fit = survfit(surv~group,clinical)
  tsv=summary(fit)
  tsv=data.frame(surv=tsv$surv,time=tsv$time,group=unlist(str_split(tsv$strata, "="))[c(F,T)])

  f = as.formula(paste("surv~",paste(c("score",coeffs),collapse=" + ")))
  coxv=summary(coxph(f,clinical))
  #print(coxv)

  #p1=paste("Cox Likelihood ratio test (variable)\n   global p-value = ", round(coxv$logtest[[3]],5) ,sep="")
  p1=paste("Cox Likelihood ratio test (variable)\n   global p-value = ", coxv$logtest[[3]] ,sep="")

  sdf=survdiff(surv~group,data=clinical,rho = 0)
  #print(sdf)
  #p2=paste("Log-rank test (groups)\n   global p-value = ", round(1-pchisq(sdf$chisq, length(sdf$n) - 1),5) ,sep="")
  p2=paste("Log-rank test (groups)\n   global p-value = ", 1-pchisq(sdf$chisq, length(sdf$n) - 1) ,sep="")

  #print(cox)
  #p3=paste("Cox Likelihood ratio test (groups)\n   global p-value = ", round(cox$logtest[[3]],5) ,sep="")
  p3=paste("Cox Likelihood ratio test (groups)\n   global p-value = ", cox$logtest[[3]] ,sep="")

  coxdf = as.data.frame(cox$coefficients)
  pv = coxdf[grep("^group",rownames(coxdf)),ncol(coxdf)]
  #p4=paste("   ",groups[1] , " vs. " , groups[2:length(groups)] , " p-value = " , round(pv,5) , sep="", collapse="\n")
  p4=paste("   ",groups[1] , " vs. " , groups[2:length(groups)] , " p-value = " , pv , sep="", collapse="\n")

  labels <- paste(groups, " (", unname(xtabs(~group, data=clinical)) , ")", sep="")

  g = ggsurvplot(fit=fit, risk.table = F, pval = paste(p1,p2,p3,p4,sep="\n"), conf.int = F, xlim = NULL,
                 title = "Kaplan-Meier plot", xlab = "Time since diagnosis (days)",
                 legend.title = "", legend.labs = labels, legend="top", palette = "jco", pval.size=4)
  #g$plot
  suppressMessages(ggsave(g$plot, filename = paste(path,"KM.pdf",sep="."), width=16, height=10))
  write.table(tsv, row.names=F, file=paste(path,"KM.tsv",sep="."), quote=F, sep="\t")

      # g = ggsurvplot(fit=fitx, risk.table = F, pval = paste(p1,p2,p3,p4,sep="\n"), conf.int = F, xlim = NULL,
      #                title = "Cox PH regression plot", xlab = "Time since diagnosis (days)",
      #                legend.title = "", legend.labs = labels, legend="top", palette = "jco", pval.size=4)
      # #g$plot
      # suppressMessages(ggsave(g$plot, filename = paste(path,"COX.pdf",sep="."), width=16, height=10))
      # write.table(tsvx, row.names=F, file=paste(path,"COX.tsv",sep="."), quote=F, sep="\t")
}

for (set in colnames(scores)[2:ncol(scores)]){
  clinical$score = clinical[,colnames(clinical) == set]
  clinical = clinical[order(clinical$score),]

  clinical$group = cut(seq(1,nrow(clinical)), breaks=2, label=c("low","high"), include.lowest=T)
  surv <- Surv(clinical$time_revised, clinical$event_revised)
  survplot(file.path(odir,paste0(set,".2.censored")))
  surv <- Surv(clinical$time_orig, clinical$event_orig)
  survplot(file.path(odir,paste0(set,".2.full")))

  if(! is.null(clinical$TP53_status) &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]$group )) == 2 &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]$group )) == 2 ){

    tmp = clinical
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".2.full.p53_wt")))

    clinical = tmp
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".2.full.p53_nv")))
    clinical = tmp
  }

  clinical$group = cut(seq(1,nrow(clinical)), breaks=3, label=c("low","medium","high"), include.lowest=T)
  surv <- Surv(clinical$time_revised, clinical$event_revised)
  survplot(file.path(odir,paste0(set,".3.censored")))
  surv <- Surv(clinical$time_orig, clinical$event_orig)
  survplot(file.path(odir,paste0(set,".3.full")))

  if(! is.null(clinical$TP53_status) &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]$group )) == 3 &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]$group )) == 3 ){

    tmp = clinical
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".3.full.p53_wt")))

    clinical = tmp
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".3.full.p53_nv")))
    clinical = tmp
  }

  clinical$group = cut(seq(1,nrow(clinical)), breaks=4, label=c("low","low-med","med-high","high"), include.lowest=T)
  surv <- Surv(clinical$time_revised, clinical$event_revised)
  survplot(file.path(odir,paste0(set,".4.censored")))
  surv <- Surv(clinical$time_orig, clinical$event_orig)
  survplot(file.path(odir,paste0(set,".4.full")))

  if(! is.null(clinical$TP53_status) &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]$group )) == 4 &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]$group )) == 4 ){

    tmp = clinical
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".4.full.p53_wt")))

    clinical = tmp
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".4.full.p53_nv")))
    clinical = tmp
  }

  clinical$group = cut(seq(1,nrow(clinical)), breaks=5, label=c("low","low-med","medium","med-high","high"), include.lowest=T)
  surv <- Surv(clinical$time_revised, clinical$event_revised)
  survplot(file.path(odir,paste0(set,".5.censored")))
  surv <- Surv(clinical$time_orig, clinical$event_orig)
  survplot(file.path(odir,paste0(set,".5.full")))

  if(! is.null(clinical$TP53_status) &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]$group )) == 5 &&
    length(unique( clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]$group )) == 5 ){

    tmp = clinical
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="wildtype",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".5.full.p53_wt")))

    clinical = tmp
    clinical=clinical[! is.na(clinical$TP53_status) & clinical$TP53_status=="mutated",]
    surv <- Surv(clinical$time_orig, clinical$event_orig)
    survplot(file.path(odir,paste0(set,".5.full.p53_nv")))
    clinical = tmp
  }
}

