ipr.densplotUnimod <- function(df,hyper,hypo,main="Unimodal Methylation by Sample",ymax=10){
  plot(2,1,type="b",xlim=c(0,1),ylim=c(0,ymax),
       xlab="Beta-val",ylab="Density",main=main)
  plotcol <- colours()[sample(nrow(df),nrow(df))]
  for(samplei in 1:nrow(df)){
    lines(density(ipr.dmpModmethy(betas=df[samplei,],hyper=hyper,hypo=hypo),na.rm=TRUE),col=plotcol[samplei])
  }
}
