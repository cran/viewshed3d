#' Ground and vegetation segmentation in TLS scene
#'
#' @description The \strong{class.ground} function combines three filters to segment the soil and the vegetation in a TLS scene of a forest.
#' In each pixel of a grid of \emph{grid.res} resolution, the points with a \emph{z} value between \emph{min(z)} and \emph{min(z)} + \emph{ground.thickness} are first classified as ground.
#' Then, in each pixel, the ground points with \emph{z} above a \emph{pix.percentile} percentile threshold are reclassified as vegetation. During this process, the \emph{z} values of each point relative
#' to the lower \emph{z} value of all points within the pixel are recorded. This enables application of a last percentile filter apllied at the plot scale, such that all points with a \emph{z} value higher than \emph{global.percentile} are reclassified as vegetation.
#' @param data a data.frame containing xyz coordinates
#' @param grid.res resolution of the grid's pixel. Default = 1.
#' @param ground.thickness the thickness of the soil layer used for primary segmentation. Default = 0.2
#' @param pix.percentile the percentile filter to use in each pixel. Default = 100 (i.e., no filter)
#' @param global.percentile the percentile filter to use at the plot scale. Default = 100 (i.e., no filter)
#' @param plot3d logical, if \emph{TRUE} the result of the segmentation is plotted with the soil in brown and the vegetation in green
#'
#'@details
#' The two percentile based filters were successfully tested in:
#' Muir, J., Goodwin, N., Armston, J., Phinn, S., & Scarth, P. (2017). An accuracy assessment of derived digital elevation models from terrestrial laser scanning in a sub-tropical forested environment. Remote Sensing, 9(8), 843.
#'
#' @return a data frame containing xyz coordinates plus the point class: "ground" or "vegetation"
#' @export
#'
#' @import VoxR rgl utils stats
#'
#' @examples
#' library(viewshed3d)
#' data(TLSscene)
#' segmented.scene=class.ground(TLSscene, plot3d = TRUE)

class.ground <- function(data,grid.res,ground.thickness,pix.percentile,global.percentile,plot3d){

  requireNamespace('VoxR')
  requireNamespace('utils')
  requireNamespace('stats')

  if(missing(pix.percentile)){pix.percentile=100 ; print("Ground segmentation without pixel percentile filter")}
  if(missing(global.percentile)){global.percentile=100 ; print("Ground segmentation without global percentile filter")}
  if(missing(grid.res)){grid.res=1 ; print("Default grid resolution = 1")}
  if(missing(ground.thickness)){ground.thickness=0.2 ; print("Default ground thickness = 0.2")}
  if(missing(plot3d)){plot3d=F}

  row.names(data)=1:nrow(data)

  vox_dat=data.frame(data)

  vox_dat[,4]=round(vox_dat[,1]/grid.res)*grid.res
  vox_dat[,5]=round(vox_dat[,2]/grid.res)*grid.res
  vox_dat[,4]=paste(vox_dat[,4],vox_dat[,5])
  vox_dat=vox_dat[,1:4]
  pix=unique(vox_dat[,4])
  vox_dat[,5]="vegetation"
  row.names(vox_dat)=1:nrow(vox_dat)
  print("Segmenting the soil in grid pixels")
  pb <- txtProgressBar(min=1,max=length(pix),style=3)
  for(i in 1:length(pix)){
    inpix=subset(vox_dat,vox_dat[,4]==pix[i])
    ground=subset(inpix,inpix[,3]<(min(inpix[,3])+ground.thickness))
    ground=ground[ground[,3]<=quantile(ground[,3],seq(0,1,0.01),na.rm=T)[pix.percentile+1],]
    vox_dat[as.numeric(row.names(ground)),5]="ground"
    vox_dat[as.numeric(row.names(inpix)),6]=inpix[,3]-min(inpix[,3])
    setTxtProgressBar(pb, i)
  }
  inground=subset(vox_dat,vox_dat[,5]=="ground")
  recal=inground[inground[,6]>=quantile(inground[,6],seq(0,1,0.01),na.rm=T)[global.percentile+1],]
  vox_dat[as.numeric(row.names(recal)),5]="vegetation"
  vox_dat=vox_dat[,c(1,2,3,5)] ; names(vox_dat)=c("x","y","z","class")
  print("Done")
  if(plot3d==T){
    requireNamespace('rgl')
    vegetation=subset(vox_dat,vox_dat[,4]=="vegetation")
    ground=subset(vox_dat,vox_dat[,4]=="ground")
    open3d()
    plot3d(vegetation[,1:3],add=T,col="darkgreen",size=2)
    plot3d(ground[,1:3],col="brown",add=T,size=3)
  }
  return(vox_dat)
}
