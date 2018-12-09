#' Additionnal filters after first ground and vegetation segmentation
#'
#' @description The \strong{ground.filter} function provides the ability to filter an output from the \code{\link{class.ground}} function.
#' The first filter removes from the ground point cloud every point with a \emph{z} value above \emph{pix.percentile} threshold within a single pixel of \emph{grid.res} resolution.
#' The second filter is a global percentile filter that removes from the ground point cloud the points with a within-pixel \emph{z} value above \emph{global.percentile} threshold.
#' @param data a data.frame containing the xyz and point class of ground points segmented with the \code{\link{class.ground}} function
#' @param grid.res resolution of a grid's pixel. Default = 1
#' @param pix.percentile the percentile filter to use in each pixel. Default = 50
#' @param global.percentile the percentile filter to use at the plot scale. Default = 95
#'
#' @return a data.frame containing the xyz coordinates plus the point class: "ground" or "reclass to veg." if the point was classed as vegetation.
#'
#' @export
#'
#' @examples
#' library(viewshed3d)
#' data(TLSclass)
#'
#' ground=subset(TLSclass,TLSclass[,4]=='ground')
#' resegmented=ground.filter(ground)
#'
#' ground = subset(resegmented,resegmented[,4]=='ground')
#' vegetation = subset(resegmented,resegmented[,4]=='reclass to veg.')
#'
#' library(rgl)
#' open3d()
#' plot3d(ground,col="brown",add=TRUE)
#' plot3d(vegetation,col="darkgreen",add=TRUE)


ground.filter <- function(data,grid.res,pix.percentile,global.percentile){

  if(missing(pix.percentile)){pix.percentile=50 ; print("Filtering ground points with median filter in pixels")}
  if(missing(global.percentile)){global.percentile=95 ; print("Filtering ground points within 95th percentile")}
  if(missing(grid.res)){grid.res=1 ; print("Default grid resolution = 1")}

  row.names(data)=1:nrow(data)

  data[,5]=round(data[,1]/grid.res)*grid.res
  data[,6]=round(data[,2]/grid.res)*grid.res
  data[,5]=paste(data[,5],data[,6])
  pix=unique(data[,5])
  pb <- txtProgressBar(min=1,max=length(pix),style=3)
  for(i in 1:length(pix)){
    inpix=data[data[,5]==pix[i],] ; inpix=na.omit(inpix)
    inq=inpix[inpix[,3]>=quantile(inpix[,3],seq(0,1,0.01),na.rm=T)[pix.percentile],]
    data[as.numeric(row.names(inq)),4]="reclass to veg."
    data[as.numeric(row.names(inpix)),6]=inpix[,3]-min(inpix[,3])
    setTxtProgressBar(pb, i)
  }
  inground=subset(data,data[,4]=="ground")
  recal=inground[inground[,6]>=quantile(inground[,6],seq(0,1,0.01),na.rm=T)[global.percentile],]    #subset(data,data[,4]=="r")
  data[as.numeric(row.names(recal)),4]="reclass to veg."
  data=data[,1:4]
  return(data)
}
