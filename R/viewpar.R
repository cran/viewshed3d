#' Compute visibility in 3D TLS scene with parallel computation
#'
#' @description The \strong{view.par} function is similar to the \code{\link{view}} with paralel computation capabilities.
#'
#' @param data a data.frame containing the xyz coordinates of a TLS point cloud
#' @param param a data.frame containing the directions that will be used to explore the point cloud
#' @param scene.center a vector containing the xyz coordinates of the scene center. Default is 0,0,0.
#' @param shape  the shape of the scene radius reshaping (see \code{\link{reshape.scene}} function for more details). Default = "3d".
#' @param scene.radius the cut-off distance relative to the \emph{scene.center}
#' @param plot.result logical. If \emph{TRUE} the \% of visibility vs. distance is ploted
#' @param ncores the number of cores to use for parallel computation. Default is N detected cores-1
#'
#' @return a data.frame containing the distance, the number of points found at this distance and the \% of remaining visibility and the distance
#'
#' @note some existing parameters are already provided with the viewshed3d package with parameters for an angular resolution ranging from 1 to 10 : \link{param_1},
#' \link{param_2}, \link{param_3}, \link{param_4}, \link{param_5}, \link{param_6}, \link{param_7}, \link{param_8}, \link{param_9}, \link{param_10}.
#'
#'
#' @export
#'
#' @import VoxR rgl doParallel foreach parallel iterators graphics
#'
#' @examples
#' \donttest{
#' library(viewshed3d)
#' data(TLSrecons)
#' data(param_10)
#'
#' center=c(mean(TLSrecons[,1]),mean(TLSrecons[,2]),mean(TLSrecons[,3]))
#' # ncores is set to 1 for the example
#' view.data=view.par(TLSrecons,param = param_10,scene.center = center, scene.radius = 2 , ncores = 1)
#'
#' head(view.data)
#' ]

view.par=function(data,param,scene.center,shape,scene.radius,plot.result,ncores){

  requireNamespace('VoxR')
  requireNamespace('rgl')
  requireNamespace('doParallel')
  requireNamespace('foreach')
  requireNamespace('parallel')
  requireNamespace('iterators')
  requireNamespace('graphics')

  if(missing(param)){stop("Please provide a data frame containing the directions")}
  if(missing(scene.center)){scene.center=c(0,0,0) ; print("Defaiult scene center is 0,0,0")}
  if(missing(scene.radius)){scene.radius=5 ; print("Default scene radius = 5")}
  if(missing(plot.result)){plot.result=T}
  if(missing(ncores)){ncores=detectCores()-1 ; print(paste("Parallel computing launched with",ncores,"cores"))}

  angular.res=param[1,1]
  print(paste("Angular resolution =",angular.res))
  param=param[-1,]
  if(missing(scene.center)==F){
    data[,1]=data[,1]-scene.center[1] ; data[,2]=data[,2]-scene.center[2] ; data[,3]=data[,3]-scene.center[3]
  }
  nc=ncol(data)
  if(missing(scene.radius)==F){
    if(missing(shape)){shape="3d" ; print("Shape = 3d")}
    if(shape=="2d"){
      data[,nc]=point.distance(data.frame(data[,1:2],0),point=c(0,0,0))
      data=subset(data,data[,nc]<=scene.radius)
    }
    if(shape=="3d"){
      data[,nc]=point.distance(data,point=c(0,0,0))
      data=subset(data,data[,nc]<=scene.radius)
    }
  }
  row.names(data)=1:nrow(data)
  near.dist=c()
  print("Computing nearest voxels in directions")

  cl <- parallel::makeCluster(ncores, outfile = "")
  registerDoParallel(cl)

  near.dist<-foreach(i = icount(nrow(param)), .combine=c, .packages = c("VoxR","rgl")) %dopar%{
    rotdat=as.matrix(data[,1:3])
    rotdat=rotate3d(rotdat, pi*param[i,3]/180, 0, 0, 1)
    rotdat=subset(rotdat,rotdat[,2]<0)
    rotdat=rotate3d(rotdat, pi*param[i,2]/180, 1, 0, 0)
    inangX=as.data.frame(rotdat)
    inangX[,4]=axis.angle(inangX,"Z")
    inangX=subset(inangX,inangX[,4]<=angular.res/2)
    if(nrow(inangX)>0){
      inangX[,4]=point.distance(inangX,c(0,0,0))
      nearest=subset(inangX,inangX[,4]==min(inangX[,4]))[1,]
      return(nearest[,4])#################################################
    }
  }
  print("Computing visibility")
  view=data.frame(round(near.dist,2),1)
  view=aggregate(view[,2],by=list(view[,1]),FUN=sum)
  seen=1
  view[,3]=NA
  for(i in 1:nrow(view)){
    see=seen-(view[i,2]/nrow(param))
    view[i,3]=see
    seen=see
  }
  if(plot.result==T){
    plot(view[,3]*100~view[,1],type="l",lwd=3,xlab="Distance from center (m)",ylab="Percent visibility (%)",cex.lab=1.5,ylim=c(0,100))
  }
  names(view)=c("Distance from center","number of voxels","% remaining visibility")
  print("Done")
  stopCluster(cl)
  return(view)
}
