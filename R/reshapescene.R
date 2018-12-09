#' Reshape a scene by applying a cut off distance and recentering a TLS scene
#'
#'@description The \strong{reshape.scene} function enables reshaping a TLS scene at any time in the data preparation process. It enables application of a
#'cut off distance and recentering a TLS scene by providing a user defined scene center, or by allowing users to interactively select a region of the point cloud to be used as the scene center.
#'
#' @param data a data.frame containing the xyz coordinates of a TLS scene
#' @param point (optional) a vector containing the xyz coordinates of the scene center. If no point is provided the user had to interactively select a region of the point cloud
#' @param radius the scene radius. Default = 5
#' @param shape If = "2d" the cut off distance is applied on the x and y axes of the scene and the scene edge is a cylinder.
#' If = "3d" the cut off distance is applied on the three dimensions of the scene and the scene edge is a sphere. Default = "2d"
#' @param plot3d logical. If \emph{TRUE} a 3d plot of the result is provided.
#'
#' @return a dataframe containing the reshaped TLS scene
#' @export
#'
#' @import tcltk2 tcltk rgl VoxR
#'
#' @examples
#'
#' library(viewshed3d)
#' data(TLSclass)
#'
#' # when a point is provided, the scene is automatically reshaped
#' center=c(mean(TLSclass[,1]),mean(TLSclass[,2]),mean(TLSclass[,3]))
#' reshaped=reshape.scene(TLSclass,point = center,radius = 2, plot3d = TRUE)
#'
#' # if no point is provided, the user has to select a region (draw a squared region with the
#' # right click) of the point cloud that will be used as the scene center
#' if(interactice()){
#' reshaped=reshape.scene(TLSclass,radius = 2, plot3d = TRUE)
#, }


reshape.scene=function(data,point,radius,shape,plot3d){

  requireNamespace('VoxR')
  requireNamespace('rgl')
  requireNamespace('tcltk2')
  requireNamespace('tcltk')

  if(missing(plot3d)){plot3d=F}
  if(missing(shape)){shape="2d"}
  if(missing(radius)){radius=5}
  if(missing(point)){
    if(ncol(data)==3){
      open3d()
      plot3d(data,add=T)
    }
    if(ncol(data)>=4){
      if(unique(data[,4]=="ground"|data[,4]=="vegetation")){
        ground=subset(data,data[,4]=="ground")
        veg=subset(data,data[,4]=="vegetation")
        open3d()
        plot3d(veg[,1:3],add=T,col="darkgreen")
        plot3d(ground[,1:3],col="brown",add=T)
      }
    }
    run=T
    while(run==T){
      print("Select a region of the point cloud (right cick) that will be used as center of the reshaped scene")
      xyz=selectpoints3d(value = T, closest = F, multiple = F,button="right")
      robert=plot3d(xyz,col="red",size=10,add=T)
      xyz=c(mean(xyz[,1]),mean(xyz[,2]),mean(xyz[,3]))
      rep=tkmessageBox(message = "Is it the right location ?",
                       icon="question", type = "yesno", default = "yes")
      if(as.character(rep)=="yes"){
        rgl.close()
        run=F
      }else{rgl.pop(id=robert)}
    }
  }else{
    xyz=point
  }
  nc=ncol(data)+1
  if(shape=="2d"){
    data[,nc]=point.distance(data.frame(data[,1:2],0),point=c(xyz[1:2],0))
    data=subset(data,data[,nc]<=radius)
  }
  if(shape=="3d"){
    data[,nc]=point.distance(data,point=xyz)
    data=subset(data,data[,nc]<=radius)
  }
  data[,1]=data[,1]-xyz[1] ; data[,2]=data[,2]-xyz[2] ; data[,3]=data[,3]-xyz[3]
  if(plot3d==T){
    open3d()
    plot3d(data[,1:3],add=T)
    plot3d(data.frame(0,0,0),size=10,col="red",add=T)
  }
  data=data[,1:(nc-1)]
  print(paste("Selected point coordinates = ",xyz[1],xyz[2],xyz[3]))
  return(data)
}
