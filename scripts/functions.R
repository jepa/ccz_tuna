
# This script contains the subset of functions needed to run the analysis
# All functions were created and taken from the code provided by Bell et al., 2021
# Original files can be found in ... 

# Functions are presented in order of usage for clarity

# Functions needed for the data_extraction() function

# From read-sea-var.R script
read.var.dym <- function(file.in,t0.user=NULL,tfin.user,region=c(NA,NA,NA,NA),dt=30,apply.mask=FALSE){
  
  #1-reading
  message("Reading file ",file.in,"...")	    
  con<-file(file.in,"rb") 
  file.type<-readChar(con,4)
  grid.id<-readBin(con,integer(),size=4)
  minval<-readBin(con,numeric(),size=4)
  maxval<-readBin(con,numeric(),size=4)
  nlon<-readBin(con,integer(),size=4)
  nlat<-readBin(con,integer(),size=4)
  nlevel<-readBin(con,integer(),size=4)
  t0.file<-readBin(con,numeric(),size=4)
  tfin.file<-readBin(con,numeric(),size=4)
  
  xlon<-array(0,c(nlat,nlon))
  ylat<-array(0,c(nlat,nlon))
  for(i in 1:nlat){
    xlon[i,]<-readBin(con,numeric(),n=nlon,size=4)
  }
  for(i in 1:nlat){
    ylat[i,]<-readBin(con,numeric(),n=nlon,size=4)
  }
  tvect<-readBin(con,numeric(),n=nlevel,size=4) 
  
  mask<-array(0,c(nlat,nlon))
  for(i in 1:nlat){
    mask[i,]<-readBin(con,integer(),n=nlon,size=4)
  }
  
  #2-time vector
  bytestoskip<-0
  dates<-tvect
  if (!is.null(t0.user)){ # extract sub-time vector
    if ((length(t0.user)==2 | length(tfin.user)==2) & dt==30){
      t0.user<-c(t0.user[1:2],15)
      tfin.user<-c(tfin.user[1:2],15)
    }
    if ((length(t0.user)==2 | length(tfin.user)==2) & dt!=30){
      message("Warning: the startdate and enddate do not contain day, will use first of month!")
      t0.user<-c(t0.user[1:2],1)
      tfin.user<-c(tfin.user[1:2],1)
      print(t0.user)
      print(tfin.user)
    }
    t0.user.date<-as.Date(paste(t0.user,collapse="-"))
    tfin.user.date<-as.Date(paste(tfin.user,collapse="-"))
    if (dt==30)
      dates<-gen.monthly.dates(year.month.sea(t0.file),year.month.sea(tfin.file))
    
    if (dt!=30)
      dates<-get.date.sea(t0.file)+seq(0,(nlevel-1)*dt,dt)
    
    ind<-which(dates>=t0.user.date & dates<=tfin.user.date)
    
    if (length(ind)==0 | any(is.na(ind))){
      message("Problem with dates!")
      print(ind)
      return()
    }
    if (length(ind)>0 & all(!is.na(ind))){
      dates<-dates[ind]		    
      message("Extracting data from ",dates[1]," to ",dates[length(dates)])		    
      tvect<-tvect[ind]	
      bytestoskip<-(ind[1]-1)*nlon*nlat*4
      nlevel<-length(tvect)
    }
  }
  
  if (bytestoskip>0){
    message("Skipping ",bytestoskip/(nlon*nlat*4)," matrices...")		
    pos<-seek(con,bytestoskip,"current")
  }
  
  data<-array(0,c(nlevel,nlat,nlon))
  for(ti in 1:nlevel){
    for(i in 1:nlat){
      data[ti,i,]<-readBin(con,numeric(),n=nlon,size=4)
    }
  }
  #convert invalid values in DYM to NA to avoid 
  #errors while treating these values as valid ones
  data<-ifelse(data==-999,NA,data)
  
  close(con)
  
  #	do.warning<-function(ind,limit)
  
  if (!any(is.na(region))){#extract sub-region
    x1<-region[1]; x2<-region[2]		
    y1<-region[3]; y2<-region[4]		
    dx<-xlon[1,2]-xlon[1,1]
    dy<-ylat[1,1]-ylat[2,1]
    i1<-xtoi.dym(x1,xlon[1,1],dx)
    i2<-xtoi.dym(x2,xlon[1,1],dx)
    if (i1<1) i1<-1; if (i1>nlon) i1<-nlon
    if (i2<1) i2<-1; if (i2>nlon) i2<-nlon
    j2<-ytoj.dym(y1,ylat[1,1],dy)
    j1<-ytoj.dym(y2,ylat[1,1],dy)
    if (j1<1) j1<-1; if (j2>nlat) j2<-nlat
    message("Extracting data from ",xlon[1,i1]," to ",xlon[1,i2]," and from ",ylat[j2,1]," to ",ylat[j1,1])
    mask<-mask[j1:j2,i1:i2]
    data<-data[,j1:j2,i1:i2]
    xlon<-xlon[j1:j2,i1:i2]
    ylat<-ylat[j1:j2,i1:i2]; 
  }
  nlat<-nrow(ylat)
  #3-apply mask, flip and transpose
  landmask.na<-ifelse(mask==0,NA,1)
  if (length(tvect)>1){
    if (apply.mask){		
      #couldn't find how make 3d array from landmask.na to multiply data on. 
      for (n in 1:length(dates)) data[n,,]<-data[n,,]*landmask.na 
    }
    data<-data[,nlat:1,]
    if (length(tvect)>1) data<-apply(data,3:2,t) 
    if (length(tvect)==1) data<-t(data[nlat:1,]) 
  }
  if (length(tvect)==1){
    message("HERE")
    data<-data*landmask.na	 
    data<-data[nlat:1,]
    data<-t(data) 
  }
  
  return(list(x=xlon[1,],y=rev(ylat[,1]),t=dates,var=data,landmask=mask))
} #end  read.var.dym


# Function needed to run read.var.dym()
# From the Utilities script
gen.monthly.dates <- function(t0,tfin){
  
  years<-seq(t0[1],tfin[1],1)
  for (y in years){
    dyear<-as.Date(paste(y,1:12,15,sep="-"))
    is.out<-(dyear<as.Date(paste(t0[1],t0[2],15,sep="-"))|
               dyear>as.Date(paste(tfin[1],tfin[2],15,sep="-")))
    if (any(is.out))
      dyear<-dyear[!is.out]
    if (y==years[1]) dates<-dyear
    if (y!=years[1]) dates<-c(dates,dyear)
  }
  return(dates)
}

# Function needed to run read.var.dym()
# From the read_varDYM.R script
year.month.sea <- function(ndat){
  year  <- trunc(ndat)
  days <- trunc((ndat - year)*365);
  date<-as.Date(paste(year,1,1,sep="-"))+days-1
  month<-as.integer(format(date,"%m"))
  return(c(year,month))
}


# From read_varDYM script
xtoi.dym <- function(x,xmin,dx) {
  return(round((x-xmin)/dx,digits=0)+1)
}
# From read_varDYM script
ytoj.dym <- function(y,ymin,dy){
  return(round((ymin-y)/dy,digits=0)+1)
}


# From the Utilities.R script
sum1d <- function(vars){
  
  nx<-dim(vars)[2]
  ny<-dim(vars)[3]
  SUM<-array(NA,c(nx,ny))
  
  for (i in 1:nx){
    for (j in 1:ny){
      SUM[i,j]<-sum(vars[,i,j],na.rm=TRUE)
    }
  }
  
  SUM<-ifelse(SUM==0,NA,SUM)
  return(SUM)
}