# Power calculation for genotype/phenotype association work in pulmonary arterial hypertensio
# Tim Dawes (2018)

install.pckages("lme4")
library(lme4)


# Calculate the variances for the variables used in the LMM
    l1<- list.files(pattern="_1.txt")
    l2<- list.files(pattern="_2.txt")
    diff<- matrix(0, nrow=18028, ncol=400)
    l1.1D.ED<- list.files(pattern="Scan1_ED.txt")
    l1.1D.ES<- list.files(pattern="Scan1_ES.txt")
    l2.1D.ED<- list.files(pattern="Scan2_ED.txt")
    l2.1D.ES<- list.files(pattern="Scan2_ES.txt")

# Set up some matrices to store the output
    data1D.1.volumes<- matrix(0,nrow=20,ncol=10)
    data1D.2.volumes<- matrix(0,nrow=20,ncol=10)
    colnames(data1D.1.volumes)<- c("LVEDV","LVMED","RVEDV","RVMED", "LVESV","LVMES","RVESV","RVMES","RVSV","RVEF")
    rownames(data1D.1.volumes)<- substr(l1.1D.ED,1,9)
    colnames(data1D.2.volumes)<- c("LVEDV","LVMED","RVEDV","RVMED", "LVESV","LVMES","RVESV","RVMES","RVSV","RVEF")
    rownames(data1D.2.volumes)<- substr(l1.1D.ED,1,9)

# Read in the data
    options(warn=-1)
    for (i in 1:20)
    {
      data1D.1.volumes[i,c(1:3,5:7)]<- unlist(c(read.table(l1.1D.ED[i]),read.table(l1.1D.ES[i])))
      data1D.2.volumes[i,c(1:3,5:7)]<- unlist(c(read.table(l2.1D.ED[i]),read.table(l2.1D.ES[i])))
      cat("\n",l1.1D.ED[i],"...",sep="")
    }
    options(warn=0)


# Calculate the RVEF / RVEDV data
    data1D.1.volumes<- data1D.1.volumes / 1000
    data1D.2.volumes<- data1D.2.volumes / 1000
    rvsv.1<- data1D.1.volumes[,3] - data1D.1.volumes[,7]
    rvef.1<- rvsv.1 / data1D.1.volumes[,3] * 100
    
    rvsv.2<- data1D.2.volumes[,3] - data1D.2.volumes[,7]
    rvef.2<- rvsv.2 / data1D.2.volumes[,3] * 100
    
    data1D.1.volumes[,9:10]<- c(rvsv.1,rvef.1)
    data1D.2.volumes[,9:10]<- c(rvsv.2,rvef.2)
    diffs<- data1D.1.volumes - data1D.2.volumes
    rownames(diffs)<- substr(l1.1D.ED,1,9)

    diffs.absolute<- abs(diffs)

# Read in the 3D data
    motion1<- motion2<- matrix(0, 20*18028,20, dimnames=list(paste(rep(paste("ID",1:20,sep=""),each=no.points), rep(paste("P",1:18028,sep=""),times=20),sep=""), paste("Time",1:20,sep="")))

# Read in function for 1st scan
    Repro.3D.1<- list.files(pattern="rv_wholecyclefunction_Scan1_")
    for (i in 1:20)
    {
      cat(i)
      start<- (i-1)*18028+1
      end<- start+18027
      motion1[start:end,]<- data.matrix(read.table(Repro.3D.1[i], colClasses = c(rep("NULL", 60),rep("numeric", 20)), header = TRUE))
      rownames(motion1[start:end,])<- rep(substr(Repro.3D.1[i],3,nchar(Repro.3D.1[i])),18028) 
    }

# Read in function for 2nd scan
    Repro.3D.2<- list.files(pattern="rv_wholecyclefunction_Scan2_")
    for (i in 1:20)
    {
      cat(i)
      start<- (i-1)*18028+1
      end<- start+18027
      motion2[start:end,]<- data.matrix(read.table(Repro.3D.2[i], colClasses = c(rep("NULL", 60),rep("numeric", 20)), header = TRUE))
      rownames(motion2[start:end,])<- rep(substr(Repro.3D.2[i],3,nchar(Repro.3D.2[i])),18028)  
    }

# Train partial least squares model based on 1st scan data
    train1<- data.frame(RVEF = rep(sapply(1:20, function(x) {rvef.1[match(substr(Repro.3D.1,3,11)[x], substr(l1.1D.ED,1,9))]}), each=18028), X = motion1)
    pls.model1<- mvr(RVEF ~ ., data=train1, validation="CV")
    pls.model1$scores[,1]

# Train partial least squares model based on 2nd scan data
    train2<- data.frame(RVEF = rep(sapply(1:20, function(x) {rvef.2[match(substr(Repro.3D.2,3,11)[x], substr(l2.1D.ED,1,9))]}), each=18028), X = motion2)
    pls.model2<- mvr(RVEF ~ ., data=train2, validation="CV")
    mean(abs(pls.model1$scores[,1] - pls.model2$scores[,1]))

# Variance in the principal components from the mvr function used on the reproducibility cohorts (Verror)
  var(pls.model1$scores[,1] - pls.model2$scores[,1]) 



diff<- matrix(0, nrow=18028, ncol=20)

for (i in 1:20)
  {
  cat(i)
  firstscan.start<- (i-1)*18028+1
  firstscan.end<- firstscan.start+18027
  
  secondscan.start<- firstscan.start + (18028*20)
  secondscan.end<- secondscan.start + 18027
  
  first.scan<- pc.motion.mat[firstscan.start:firstscan.end,]
  second.scan<- pc.motion.mat[secondscan.start:secondscan.end,]
  
  diff.vector<- c(first.scan[,1] - second.scan[,1])
  diff[,i]<- diff.vector
}
  
  
  nSubjects.at.start<- sum(na.omit(as.numeric(substr(c(d[,8]),1,4)))<2017)
  nScans.at.start<- sum(na.omit(as.numeric(sapply(8:17, function(x) {substr(c(d[,x]),1,4)}))<2017))
  nSubjects.now<- length(na.omit(unlist(c(d[,8:17]))))
  
  # One SD of RVEF = 3.5-4.0%
  # So, bBMPR2 = 1 implies reduction of RVEF of 3.5-4.0% per year



# Set up simulation function
  
    sim1 <- function(bSEX, bAGE, b0, bBMPR2, Vsubj, Verror, nSubjects.at.start, nScans.at.start, i) {
      nSubjects<- nSubjects.at.start + i
      nScans<- nScans.at.start + i
      
      Subject <- sample(1:nSubjects, nScans, replace=TRUE)
      AGE <- rnorm(nSubjects,58.9,16.7)[Subject]
      SEX<- c("M","F")[round(runif(nSubjects,0,1),0)+1][Subject]
      TimeFromDiagDate <- runif(nSubjects,0,5)[Subject]
      
      BMPR2<- rep(0,nSubjects)
      BMPR2[sample(1:nSubjects,0.24*nSubjects)]<- 1
      BMPR2[sample(1:nSubjects,0.05*nSubjects)]<- 2
      BMPR2<- BMPR2[Subject]
      BMPR2[1:i]<- 2
      
      
      # random effects per subject
      S.re <- sample(rnorm(nSubjects, 0, sqrt(Vsubj)), nScans, replace=TRUE)
      
      # epsilons
      eps <- sample(rnorm(nSubjects, 0, sqrt(Verror)), nScans, replace=TRUE)
      
      # put it all together
      PC1 <- b0 + bSEX*(SEX=='M') + bAGE*AGE + bBMPR2*sd(pls.model1$scores[,1])*TimeFromDiagDate*BMPR2 + S.re + eps # 
      # put into a data frame
      mydata <- data.frame(Subject = paste('s',Subject, sep=''), AGE = AGE, SEX = SEX, BMPR2 = BMPR2, TimeFromDiagDate = TimeFromDiagDate, PC1 = PC1)
      
      # analyze looking at interaction term with LR test
      return(mydata)
    }
 

# Do simulation!
    nReps<- 500
    s<- c(1,seq(10,130,10))
    out<-rep(0, nReps)
    cols<- colorRampPalette(c("white","red"))(6)[-1]
    # BMPR2 carriers have 9% lower RVEF at diagnosis (van der Bruggen, 2016, Circ); assuming 2 year disease duration at Dx --> 4.5% higher drop per year = ~1 sd/year
    # Annual drop in RVEF: Std-beta = -0.38 per 5 years
    # Effect.size is the number of sds for BMPR2, 1sd = 11.5% change in global 3D function

    effect.size<- seq(0.01,0.05,length.out=5)
    power<- matrix(0, nrow=length(effect.size),ncol=length(s), dimnames=list(effect.size,s))
    

    
for (j in 1:length(effect.size)) {
  for (i in 1:length(s)) {cat(s[i],",",sep="")
    for (k in 1:nReps)
      {
        mydata <- sim1(bSEX=-0.576, bAGE=-0.076, b0=0, bBMPR2=effect.size[j]*sd(rvef.1), Vsubj=.1, Verror=14.7, nSubjects.at.start, nScans.at.start, s[i])
        
        fit1<- tryCatch(lme(PC1 ~ BMPR2 + TimeFromDiagDate + AGE + SEX, mydata, random = ~ 1 | Subject, correlation=corCAR1(), control = lmeControl(opt='optim', msTol=(10^-6)), method="ML"),
                        error = function(e) {TRUE})
        fit2<- tryCatch(lme(PC1 ~ BMPR2 + TimeFromDiagDate + AGE + SEX + BMPR2*TimeFromDiagDate, mydata, random = ~ 1 | Subject, correlation=corCAR1(), control = lmeControl(opt='optim', msTol=(10^-6)), method="ML"),
                        error = function(e) {TRUE})
        
        if (is.logical(fit1)==TRUE || is.logical(fit2)==TRUE) {
          fit1<- lmer(PC1 ~ BMPR2 + TimeFromDiagDate + BMPR2*TimeFromDiagDate + AGE + SEX + (1|Subject), mydata, REML=F)
          fit2<- lmer(PC1 ~ BMPR2 + TimeFromDiagDate + AGE + SEX + (1|Subject), mydata, REML=F)
          out[k]<- anova(fit2,fit1)[2,8]} else {out[k]<- anova(fit1,fit2)[2,9]}
        
      }
    
          # mean SD = 39.2, median = 24, bAGE based on drop of 1% per year in RVEF, van de Veerdonk, JACC, 2011
      power[j,i]<- mean(out<0.05)
     
      plot(1:ncol(power), 100*power[1,], pch=19, col=cols[1], ylim=c(0,100), type='p', xlab="Number of additional subjects", ylab="Power (%)", xaxt='n')
        lines(lowess(1:ncol(power),100*power[1,], f=1), col=cols[1],lwd=4)
        axis(side=1, at=1:length(s), labels=s)
        if (j>1) {for (k in 2:j) {points(1:ncol(power), 100*power[k,], pch=19, col=cols[k], type='p');
          lines(lowess(1:ncol(power),100*power[k,], f=1), col=cols[k],lwd=4, lend=2) }}
        
      # Add lines to graph
          abline(h=80, lty=2, col="lightblue", lwd=4)
          abline(h=90, lty=2, col="darkblue", lwd=4)
          #v<- s*c(0,43,55,95,126,126)/130, function(x) {which.min(abs(x-s))})
          v<- 1+(ncol(power)-1)*c(0,43,55,95,126,126)/129
          for (x in v) {segments(x,.4,x,99.6, lty=1, col="darkblue", lwd=4, lend=2)}
          # 43 scans done 1/2/17 - 7/8/18, with 12 more booked as 7/8/18 --> date and paid for (?)
          # Power if 43 (done) + 12 (booked+paid) + 40 (collabs) + 31 (AMS grant) = 125 scans done (incl. 3% drop out in HH future)
          
          for (x in 1:(length(v)-1)) {
            text(v[x]+0.2, 97, labels=x, cex=2);
            rect(v[x],0,v[x+1],100, col="darkblue", border=NA, density=(x*10))}
          text(tail(v,1)+0.2, 97, labels=length(v)-1, cex=2);
          for (x in 1:5) {text(0.75,c(8,15,25,40,62)[x],labels=paste(x,"%",sep=""))}
          
          
          
          
  }
}

write.table(power, file="Results.txt", col.names=T, row.names=T)


        