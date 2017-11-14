.fN_ratio = function(lambda,mu,gamma,ts,te,rho,sigma) {
  if(ts<te) stop("Error, starting time cannot be smaller than ending time")
  if(ts-te<1e-10) return(1)
  const = .constants(lambda,mu,gamma,sigma)
  A1 = lambda*(1-rho) + const[2]
  A2 = lambda*(1-rho) + const[3]
  res = exp(const[1]*(te-ts))*((A2*exp(-const[1]*te)-A1)/(A2*exp(-const[1]*ts)-A1))^2
  res
}

.p0 = function(lambda,mu,gamma,t,rho,sigma,init_time=0, init_value = NULL) {
  if(init_time<1e-10) init_value=1-rho
  const = .constants(lambda,mu,gamma,sigma)
  A1 = lambda*init_value + const[2]
  A2 = lambda*init_value + const[3]
  res = (A1*const[3]-A2*const[2]*exp(-const[1]*(t-init_time)))/(lambda*(A2*exp(-const[1]*(t-init_time))-A1))
  res
}

#Formula from Stadler 2010
#Assumptions: gamma = sigma = 0 in unresolved parts, invalid otherwise!
.pn = function(lambda,mu,nspecies,tcut,rho=1) {
  d = lambda - mu
  res = rho*d^2*exp(-d*tcut)*(rho*lambda*(1-exp(-d*tcut)))^(nspecies-1)
  res = res/(rho*lambda+((1-rho)*lambda-mu)*exp(-d*tcut))^(nspecies+1)
  res
}

.fN_ratio_decay = function(lambda0,z,mu,gamma,ts,te,t0,rho,sigma,ival=NULL) {
  if(ts<te) stop("Error, starting time cannot be smaller than ending time")
  if(ts-te<1e-10) return(1)
  
  lbda = .average_rate(lambda0,z,ts,te,t0)
  const = .constants(lbda,mu,gamma,sigma)
  A1 = lbda*ival + const[2]
  A2 = lbda*ival + const[3]
  res = exp(const[1]*(te-ts))*((A2-A1)/(A2*exp(-const[1]*(ts-te))-A1))^2
  res
}

.calc_p_decay = function(lambda0,z,mu,gamma,grid,t0,rho,sigma,ival) {
  values = c()
  
  #steps in the grid need to be recorded
  values[length(grid)] = ival
  if(length(grid)==2) return(values)
  for(i in (length(grid)-1):2) {
    l = .average_rate(lambda0,z,grid[i],grid[i+1],t0)
    values[i] = .p0(l,mu,gamma,grid[i],rho,sigma,init_time = grid[i+1],init_value = values[i+1])
  }
  values
}

.initial_calc_p_decay = function(lambda0,z,mu,gamma,grid,t0,rho,sigma,stepsize) {
  ord = order(grid)
  
  values = c()
  ival = 1-rho
  grid = c(0,grid[ord])
  for(g in 2:length(grid)) {
    te = grid[g-1]
    ts = grid[g]
    if(ts > 1e-10 && (ts-te)>1e-10) {
      nsteps = floor((ts-te)/stepsize)
      ts_step = ts - nsteps*stepsize
      te_step = te
      for(i in 1:(nsteps+1)) {
        l = .average_rate(lambda0,z,ts_step,te_step,t0)
        ival = .p0(l,mu,gamma,ts_step,rho,sigma,init_time = te_step,init_value = ival)
        te_step = ts_step
        ts_step = ts_step + stepsize
      }
      #if(nsteps>0 && ts_step != ts) browser()
    }
    values[g-1] = ival
  }
  values[order(ord)]
}

.average_rate = function(val0, rate, ts, te, t0) {
  if(rate==0) return(val0)
  if(rate*(ts-te)<1e-10) return(val0*exp(rate*(ts-t0)))
  av = (val0/(rate*(ts-te)))*(exp(rate*(ts-t0))-exp(rate*(te-t0)))
  av
}

.constants = function(lambda,mu,gamma,sigma=0) {
  c = sqrt((gamma+mu+lambda)^2-4*mu*lambda*(1-sigma))
  x1 = -(gamma+mu+lambda+c)/2
  x2 = -(gamma+mu+lambda-c)/2
  const = c(c,x1,x2)
  const
}

.get_desc = function(tree) {
  desc = list()
  length(desc) = tree$Nnode
  ntips = length(tree$tip.label)
  mat2 = tree[[1]]
  apply(mat2, 1, function(v) {
    a=v[1]-ntips
    desc[[a]] <<- c(v[2],desc[[a]])
  })
  
  desc
}