#' Title of the function
#'
#' Description of what Est0 does
#' @export
Est0 <- function(data) {
  # 함수 내용
  print("Est0 function is working!")
}

Est0<-function(train0){
#  use_package("survival")
  ###### Right censored data for cure rate model ############################

  cens=train0$cens
  TT=train0$TT;
  train0$ys0<-ys0<-ifelse(cens==0,0,1)
  delta=ys0

  n<-nrow(train0)
  print (n)

  fit1<-glm(ys0~z1+z2+z3, data=train0, family=binomial)
  y0<-fit1$fitted

  fit0<-survival::survfit(survival::Surv(TT,cens)~1, data=train0)
  ss=ipred::getsurv(fit0,TT)
  tm<-fit0$time; m=length(tm)
  delta=cens

  if(m==0)  {weights<-rep(0,n) ; break}

  ss0=rep(0,n)

  m2=m-1
  for(jj in 1:m2){
    for(i in 1:n) {
      if(tm[jj]<=TT[i]&TT[i]<=tm[jj+1]) ss0[i]=ss[jj]
    }}

  ys=rep(0,n)
  for(i in 1:n) {
    ys[i]=ifelse(cens[i]==1,1,y0[i]*ss0[i]/(y0[i]*ss0[i]+1-y0[i]))
  }

  ys<-ifelse(ys>=1,1,ys)

  if(m==1) break

  ###############################################
  err1=0
  ss1=rep(0,m)
  ss02=rep(0,n)

  fitt<-glm(ys~z1+z2+z3,data=train0, family=binomial)
  y0<-fitt$fitted


  for (kk in 1:50){
    for(i in 1:n) {
      ys[i]=ifelse(cens[i]==1,1,y0[i]*ss0[i]/(y0[i]*ss0[i]+1-y0[i]))
    }

    ys<-ifelse(ys>1,1,ys)

    tot=0; tot1=1
    ss22=rep(0,m)

    for(jj in 1:m){
      dd<-sum(tm[jj]==TT&cens==1)
      rr=sum(ifelse(tm[jj]<=TT,1,0)*ys)
      lam=dd/rr
      tot=tot+lam
      clam=tot
      ss1[jj]=exp(-clam)
      qqq=(1-dd/rr)
      tot1=tot1*qqq
      ss1[jj]=tot1
      #print (c(ss1[jj],ss22[jj]))
    }


    m2=m-1
    for(jj in 1:m2){
      for(i in 1:n) {
        if(tm[jj]<=TT[i]&TT[i]<=tm[jj+1]) ss02[i]=ss1[jj]
      }}
    err1=sum(abs(ss0-ss02))
    ss0=ss02
    if(kk>1&err1<0.005) break
  }

  return(list("surv"=ss1,"tm"=tm,"ys"=ys))

}
