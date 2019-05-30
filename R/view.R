#' Compute visibility in 3D TLS scene
#'
#' @description The \strong{view} function explores a TLS point cloud in all directions of the 3D space and records the nearest point in each direction.
#' A single direction is thus assumed to be a sligthline that ends as soon as an object is encountered (see package description for more details).
#' The view function requires a parameter data.frame produced with the \code{\link{view.param}} function, and the angle of a single slightline is thus defined at this step.
#' The visibility is computed from a user defiend location.
#'
#' @param data a data.frame containing the xyz coordinates of a TLS point cloud
#' @param param a data.frame containing the directions that will be used to explore the point cloud
#' @param location a vector containing the xyz coordinates of the point representing the animal location. Default is 0,0,0.
#' @param plot3d logical. If \emph{TRUE} a 3D view of the point cloud exploration is plotted
#' @param plot.result logical. If \emph{TRUE} the \% of visibility vs. distance is ploted
#'
#' @return a list containing the remaining visibility as function of the distance from the scene center ($visibility) and the 3D point cloud of
#' the portion of data seen from the scene center ($points).
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
#' view.data=view(TLSrecons,param = param_1,location = center, plot3d = TRUE)
#' }

view=function(data,param,location,plot3d,plot.result){
  requireNamespace('VoxR')
  requireNamespace('rgl')
  requireNamespace('graphics')
  requireNamespace('data.table')

  if(missing(plot3d)){plot3d=F}
  if(missing(plot.result)){plot.result=T}
  if(missing(location)){location = c(0,0,0)}

  dat=data.frame(data[,1]-location[1],data[,2]-location[2],data[,3]-location[3])

  nr=nrow(param)
  a=param[-1,]
  res.ang=param[1,1]
  param=data.frame(a[,2],1)
  param=aggregate(param[,2],by=list(param[,1]),FUN=sum)

  dat=dat[,1:3]
  dat=setDT(dat)
  names(dat)=c("X","Y","Z")
  dat[,V4:=axis.angle(dat[,1:3],"Z")]
  dat[,V4:=round(dat[,4]/res.ang)*res.ang]
  dat[,V5:=0]
  dat[,V6:=0]
  ang=unique(dat[,4])
  pb <- txtProgressBar(min = 1, max = nrow(ang), style = 3)
  for(j in 1:nrow(ang)){
    setTxtProgressBar(pb, j)
    inang=dat[V4==as.numeric(ang[j])]
    ndiv=param[param[,1]==as.numeric(ang[j]),2]
    res=360/ndiv
    inang[,V5:=axis.angle(inang[,1:3],"Y",projected=TRUE,plan="xy")]
    inang[,V5:=round(inang[,5]/(360/ndiv))*(360/ndiv)]
    inang[,V6:=point.distance(inang,c(0,0,0))]
    near=inang[,min(V6),by=V5]
    points=inang[V6%in%near[,V1] & V5 %in% near[,V5],]
    if(j==1){
      out=data.frame(points)
    }else{
      out=rbind(out,data.frame(points))
    }
  }
  view=data.table(out[,6],1)
  view[,V1:=round(V1,digits=2)]
  view=view[,sum(V2),by=V1]
  view=view[order(V1)]
  names(view)=c("distance_from_center","visibility")
  view[,visibility:=(1-cumsum(visibility)/nr)*100]
  view=rbind(list(0,100),view)
  if(plot3d==T){
    open3d()
    plot3d(dat,add=T,size=0.8)
    spheres3d(c(0,0,0),add=T,col="green",radius=0.2)
    plot3d(out,add=T,size=2,col="red")
  }
  if(plot.result==T){
    plot(view,type="l",lwd=3,xlab="Distance from animal location", ylab="% remaining visibility",ylim=c(0,100))
  }
  out=data.frame(out[,1]+location[1],out[,2]+location[2],out[,3]+location[3])
  names(out)=c("X","Y","Z")
  ret=list(visibility=data.frame(view),points=data.frame(out[,1:3]))
  return(ret)
}
