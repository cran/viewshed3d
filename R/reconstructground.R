#' Reconstruct the missing portions of a TLS scene's ground
#'
#' @description The \strong{reconstruct.ground} function is used to reconstruct the missing portions of the ground after segmentation of the TLS scene
#' with the \code{\link{class.ground}} function. It creates a grid of pixels of user defined resolution and attributes to each pixel an elevation
#' value (i.e., a \emph{z} value) calculated as the average value (weighted by distance) of the neighboring ground points. Two methods for neighboring points
#' selection are provided: k nearest neighbors or absolute cartesian distance.
#'
#' @seealso \code{\link{reconstruct.ground.par}} for more effecient computation.
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
#'
#'
#' @return a data.frame containing the xyz coordinates plus the class of the reconstructed points
#' @export
#' @import VoxR stats
#' @examples
#'
#' library(viewshed3d)
#' data(TLSclass)
#' ground=subset(TLSclass,TLSclass[,4]=='ground')
#' reconstructed=reconstruct.ground(ground,grid.res = 0.2)
#'
#' library(rgl)
#' open3d()
#' plot3d(reconstructed,add=TRUE)

#'
reconstruct.ground <- function(data,grid.res,scene.radius,scene.center,shape,method,d,k,reconstruct.all){

  requireNamespace('VoxR')
  requireNamespace('stats')

  if(missing(method)){method = "k-nearest" ; k=5 ; print("Reconstruct ground with k-nearest neighboors and k=5")}
  if(missing(grid.res)){grid.res=0.05; print("Default voxel resolution = 0.05")}
  if(missing(reconstruct.all)){reconstruct.all=T}
  if(method=="distance" & missing(d)){d=1 ; k=NA ; print("Default distance (d) is set to 1")}
  if(method=="k-nearest" & missing(k)){k=5 ; d=NA ; print("Default number of neighboors (k) is set to 5")}
  if(method=="k-nearest" & k<=1){stop("k should be greater than 1")}

  print("Ground reconstruction init.")
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
  pb <- txtProgressBar(min=1,max=nrow(full),style=3)
  print("reconstructing the ground")
  if(method=="distance"){
    for(i in 1:nrow(full)){
      x=full[i,1] ; y=full[i,2]
      data[,4]=point.distance(data.frame(data[,1],data[,2],0),c(x,y,0))
      near=subset(data,data[,4]<d)
      full[i,3]=sum(near[,3]*(near[,4]/sum(near[,4])))
      setTxtProgressBar(pb, i)
    }
  }
  if(method=="k-nearest"){
    for(i in 1:nrow(full)){
      x=full[i,1] ; y=full[i,2]
      data[,4]=point.distance(data.frame(data[,1],data[,2],0),c(x,y,0))
      near=data[order(data[,4])[1:k],]
      full[i,3]=sum(near[,3]*(near[,4]/sum(near[,4])))
      setTxtProgressBar(pb, i)
    }
  }
  if(reconstruct.all==F){
    names(full)=names(data)[1:3]###############
    full=rbind(full[,1:3],data[1:3])
  }
  full=data.frame(full[,1:3],"g") ; names(full)=c("x","y","z","class")
  print("Done")
  return(full)
}
