#' Reconstruct the missing portions of a TLS scene's ground
#'
#' @description The \strong{reconstruct.ground.par} function is similar to the \code{\link{reconstruct.ground}} function with parallel computation capabilities.
##'
#'
#' @param data a data.frame containing the xyz and point class of ground points segmented with the \code{\link{class.ground}} and optionally with the \code{\link{ground.filter}} functions
#' @param grid.res the resolution of a grid used to locate the new points to add. Default = 0.05
#' @param scene.radius optional. Define scene radius relative to the 0,0,0 coordinate. Should be used after the scene was reshaped with the \code{\link{reshape.scene} function}
#' @param scene.center a vector containing the xyz coordinates of the scene center for reshaping
#' @param shape the shape of the scene radius reshaping (see \code{\link{reshape.scene}} function for more details)
#' @param method the method to find the neighboring points.
#' If = "distance" all the points located within a distance \emph{d} are used to estimate \emph{z}.
#' If = "k-nearest" only the \emph{k} nearest neighbors are used.
#' @param d set the distance of neighboring points if \emph{method} = "distance". Default = 1.
#' @param k set the number of neighboring points if \emph{method} = "k-nearest". Default = 5.
#' @param reconstruct.all logical. If \emph{TRUE} the entire ground surface is reconstructed. If FALSE only the ground
#' portions with missing points are reconstructed
#' @param ncores the number of cores to use for parallel computation. Default = N detected cores-1
#'
#'
#' @return a data.frame containing the xyz coordinates plus the class of the reconstructed points
#' @export
#'
#' @import VoxR doParallel foreach parallel iterators
#'
#' @examples
#'\donttest{
#' library(viewshed3d)
#' data(TLSclass)
#' ground=subset(TLSclass,TLSclass[,4]=='ground')
#' # increase ncores for parallel computation
#' reconstructed=reconstruct.ground.par(ground,grid.res = 0.2, ncores=2)'
#' library(rgl)
#' open3d()
#' plot3d(reconstructed,add=TRUE)
#' }


reconstruct.ground.par <- function(data,grid.res,scene.radius,scene.center,shape,method,d,k,reconstruct.all,ncores){
  requireNamespace('VoxR')
  requireNamespace('doParallel')
  requireNamespace('foreach')
  requireNamespace('parallel')
  requireNamespace('iterators')

  if(missing(method)){method = "k-nearest" ; k=5 ; print("Reconstruct ground with k-nearest neighboors and k=5")}
  if(missing(grid.res)){grid.res=0.05; print("Default voxel resolution = 0.05")}
  if(missing(reconstruct.all)){reconstruct.all=T}
  if(missing(ncores)){ncores=detectCores()-1 ; print(paste("Parallel computing launched with",ncores,"cores"))}
  if(method=="distance" & missing(d)){d=1 ; k=NA ; print("Default distance (d) is set to 1")}
  if(method=="k-nearest" & missing(k)){k=5 ; d=NA ; print("Default number of neighboors (k) is set to 5")}
  if(method=="k-nearest" & k<=1){stop("k should be greater than 1")}

  data=data.frame(round(data[,1]/grid.res)*grid.res,round(data[,2]/grid.res)*grid.res,data[,3])
  data=aggregate(data[,3],by=list(data[,1],data[,2]),FUN=mean)
  nc=ncol(data)+1
  if(missing(scene.radius)==F){
    if(missing(scene.center)){scene.center=c(0,0,0) ; print("Scene.center was set to 0,0,0")}
    if(missing(shape)){shape="3d" ; print("Shape = 3d")}
    if(shape=="2d"){
      data[,nc]=point.distance(data.frame(data[,1:2],0),point=c(scene.center[1],scene.center[2],0))
      data=subset(data,data[,nc]<=scene.radius)
      if(nrow(data)==0){stop("No remaining ground points, scene.radius is to small or scene.center is wrong")}
    }
    if(shape=="3d"){
      data[,nc]=point.distance(data,point=scene.center)
      data=subset(data,data[,nc]<=scene.radius)
      if(nrow(data)==0){stop("No remaining ground points, scene.radius is to small or scene.center is wrong")}
    }
  }
  full=expand.grid(seq(min(data[,1]),max(data[,1]),grid.res),seq(min(data[,2]),max(data[,2]),grid.res))
  full[,3]=mean(data[,3])
  if(reconstruct.all==F){
    names(full)=names(data)[1:3]
    full=rbind(data[,1:3],full[,1:3])
    full[,4]=1
    full=aggregate(full[,4],by=list(full[,1],full[,2]),FUN=sum)
    full=subset(full,full[,3]==1)
  }
  full[,3]=mean(data[,3])
  if(missing(scene.radius)==F){
    if(shape=="2d"){
      full[,nc]=point.distance(data.frame(full[,1:2],0),point=c(scene.center[1],scene.center[2],0))
      full=subset(full,full[,nc]<=scene.radius)
      if(nrow(full)==0){stop("No remaining ground points, scene.radius is to small or scene.center is wrong")}
    }
    if(shape=="3d"){
      full[,nc]=point.distance(full,point=scene.center)
      full=subset(full,full[,nc]<=scene.radius)
      if(nrow(full)==0){stop("No remaining ground points, scene.radius is to small or scene.center is wrong")}
    }
  }
  cl <- parallel::makeCluster(ncores, outfile = "")
  registerDoParallel(cl)

  print("reconstructing the ground")
  i=NULL
  if(method=="distance"){
    returns<-foreach(i = icount(nrow(full)), .combine=rbind, .packages = c("VoxR")) %dopar% {
      x=full[i,1] ; y=full[i,2]
      data[,4]=point.distance(data.frame(data[,1],data[,2],0),c(x,y,0))
      near=subset(data,data[,4]<d)
      return(data.frame(sum(near[,3]*(near[,4]/sum(near[,4])))))
    }
    full[,3]=returns
  }
  if(method=="k-nearest"){
    returns<-foreach(i = icount(nrow(full)), .combine=rbind, .packages = c("VoxR")) %dopar% {
      x=full[i,1] ; y=full[i,2]
      data[,4]=point.distance(data.frame(data[,1],data[,2],0),c(x,y,0))
      near=data[order(data[,4])[1:k],]
      return(data.frame(sum(near[,3]*(near[,4]/sum(near[,4])))))
    }
    full[,3]=returns
  }
  if(reconstruct.all==F){
    names(full)=names(data)[1:3]
    full=rbind(full[,1:3],data[1:3])
  }
  full=data.frame(full[,1:3],"g") ; names(full)=c("x","y","z","class")
  print("Done")
  stopCluster(cl)
  return(full)
}
