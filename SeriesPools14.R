## function for modeling a series of m pools, including carbon (mass) and 14C ##

SeriesPools14<- function 
(t,      
 m.pools, 
 ki,	
 Tij, 
 C0,	
 F0_Delta14C,  
 In,    
 xi=1,  
 inputFc,
 lambda=-0.0001209681, 
 lag=0, 
 solver=deSolve.lsoda.wrapper,  
 pass=FALSE  
)	
{ 
  t_start=min(t)
  t_end=max(t)
  if(length(ki)!=m.pools) stop("ki must be of length = m.pools")
  if(length(C0)!=m.pools) stop("the vector with initial conditions must be of length = m.pools")
  if(length(In)==m.pools){
    inputFluxes=BoundInFluxes(
      function(t){matrix(nrow=m.pools,ncol=1,In)},
      t_start,
      t_end
    )
  }
  if(class(In)=="data.frame"){
    x=In[,1]  
    y=In[,2]  
    inputFlux=splinefun(x,y)
    inputFluxes=BoundInFluxes(
      function(t){matrix(nrow=m.pools,ncol=1,c(inputFlux(t),rep(0,m.pools-1)))},
      min(x),
      max(x)
    )
  }
  A=-1*abs(diag(ki))
  a=abs(ki[-length(ki)])*Tij
  ij=matrix(c((2:m.pools),(1:(m.pools-1))),ncol=2)
  A[ij]=a
  if(length(xi)==1) fX=function(t){xi}
  if(class(xi)=="data.frame"){
    X=xi[,1]
    Y=xi[,2]
    fX=splinefun(X,Y)
  }
  Af=BoundLinDecompOp(
    function(t){fX(t)*A},
    t_start,
    t_end
  )
  Fc=BoundFc(inputFc,lag=lag,format="Delta14C")
  Mod=GeneralModel_14(t=t,A=Af,ivList=C0,initialValF=ConstFc(F0_Delta14C,"Delta14C"),inputFluxes=inputFluxes,inputFc=Fc,di=lambda,pass=pass)
  return(Mod)
}
