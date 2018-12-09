#' Compute visibility in 3D TLS scene
#'
#' @description The \strong{view} function explores a TLS point cloud in all directions of the 3D space and records the nearest point in each direction.
#' A single direction is thus assumed to be a sligthline that ends as soon as an object is encountered (see package description for more details).
#' The view function requires a parameter data.frame produced with the \code{\link{view.param}} function, and the angle of a single slightline is thus defined at this step.
#'
#' @seealso \code{\link{view.par}} for more effecient computation.
#'
#' @param data a data.frame containing the xyz coordinates of a TLS point cloud
#' @param param a data.frame containing the directions that will be used to explore the point cloud
#' @param scene.center a vector containing the xyz coordinates of the scene center. Defauls is 0,0,0.
#' @param shape  the shape of the scene radius reshaping (see \code{\link{reshape.scene}} function for more details). Default = "3d"
#' @param scene.radius the cut-off distance relative to the \emph{scene.center}
#' @param plot3d logical. If \emph{TRUE} a 3D view of the point cloud exploration is plotted
#' @param plot.result logical. If \emph{TRUE} the \% of visibility vs. distance is ploted
#'
#' @return a data.frame containing the distance, the number of point found at this distance and the \% of remaining visibility.
#'
#' @note some existing parameters are already provided with the viewshed3d package with parameters for an angular resolution ranging from 1 to 10 : \link{param_1},
#' \link{param_2}, \link{param_3}, \link{param_4}, \link{param_5}, \link{param_6}, \link{param_7}, \link{param_8}, \link{param_9}, \link{param_10}.
#'
#' @export
#'
#' @import VoxR rgl graphics
#'
#' @examples
#' \donttest{
#' library(viewshed3d)
#' data(TLSrecons)
#' data(param_10)
#'
#' center=c(mean(TLSrecons[,1]),mean(TLSrecons[,2]),mean(TLSrecons[,3]))
#'
#' view.data=view(TLSrecons,param = param_10,scene.center = center, scene.radius = 2, plot3d = TRUE)
#'
#' head(view.data)
#' }

view=function(data,param,scene.center,scene.radius,shape,plot3d,plot.result){
  requireNamespace('VoxR')
  requireNamespace('rgl')
  requireNamespace('graphics')

  if(missing(param)){stop("Please provide a data frame containing the directions")}
  if(missing(scene.center)){scene.center=c(0,0,0) ; print("Defaiult scene center is 0,0,0")}
  if(missing(scene.radius)){scene.radius=5 ; print("Default scene radius = 5")}
  if(missing(plot3d)){plot3d=F}
  if(missing(plot.result)){plot.result=T}

  angular.res=param[1,1]
  print(paste("Angular resolution =",angular.res))
  param=param[-1,]

  if(plot3d==T){open3d()}

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
  if(plot3d==T){plot3d(data,size=0.7,add=T)}
  print("Computing nearest voxels in directions")
  pb <- txtProgressBar(min=1,max=nrow(param),style=3)
  for(i in 1:nrow(param)){
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
      near.dist=c(near.dist,nearest[,4])
      if(plot3d==T){
        plot3d(data[as.numeric(row.names(nearest)),],add=T,size=6,col="red")
      }
    }
    setTxtProgressBar(pb, i)
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
  return(view)
}
