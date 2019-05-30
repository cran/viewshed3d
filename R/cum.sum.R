#' Compute cumulated viewsheds in a 3D TLS scene from multiple location
#'
#' @param data a data.frame containing the xyz coordinates of a TLS point cloud
#' @param locations a data.frame with each row containing the xyz coordinates of one point representing the animal location at a single point in time.
#' @param param a data.frame containing the directions that will be used to explore the point cloud
#' @param vox_res the voxel resolution
#'
#' @return a data frame containing the x,y,z coordinates of a voxel cloud plus the number of locations from which each voxel could be seen
#' @export
#'
#' @import data.table VoxR
#'
#' @examples
#' \donttest{
#' library(viewshed3d)
#' data(TLSrecons)
#' data(param_2)
#'
#' locations=data.frame(c(12,12),c(58,60),c(0,0))
#'
#' cum=cum.view(TLSrecons,param = param_2,location=locations,vox_res = 0.2)
#'
#' library(viridis)
#' library(rgl)
#' col=cividis(max(cum[,4])+1)
#'
#' plot3d(cum,col=col[cum[,4]+1],add=TRUE,size=8)
#' points3d(locations,col="red",size=5)
#'}
cum.view=function(data,locations,param,vox_res){
  requireNamespace("data.table")
  requireNamespace("VoxR")

  if(missing(locations)){stop("Please provide the locations")}
  if(missing(vox_res)){vox_res=0.5}

  vox_dat=round(data[,1:3]/vox_res)*vox_res
  vox_dat[,4]=0

  a=param[-1,]
  res.ang=param[1,1]
  param=data.frame(a[,2],1)
  param=aggregate(param[,2],by=list(param[,1]),FUN=sum)
  print(paste("Computing cumulative viewsheds from",nrow(locations),"locations"))
  if(nrow(locations)>1){pb <- txtProgressBar(min = 1, max = nrow(locations), style = 3)}
  for(i in 1:nrow(locations)){
    if(nrow(locations)>1){setTxtProgressBar(pb, i)}
    dat=data.frame(data[,1]-locations[i,1],data[,2]-locations[i,2],data[,3]-locations[i,3])
    dat=setDT(dat)
    names(dat)=c("X","Y","Z")
    dat[,V4:=axis.angle(dat[,1:3],"Z")]
    dat[,V4:=round(dat[,4]/res.ang)*res.ang]
    dat[,V5:=0]
    dat[,V6:=0]
    ang=unique(dat[,4])
    for(j in 1:nrow(ang)){
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
    out[,1]=out[,1]+locations[i,1] ; out[,2]=out[,2]+locations[i,2] ; out[,3]=out[,3]+locations[i,3]
    vox_out=unique(data.frame(round(out[,1]/vox_res)*vox_res,round(out[,2]/vox_res)*vox_res,round(out[,3]/vox_res)*vox_res))
    vox_out[,4]=1
    names(vox_dat)=names(vox_out)
    vox_dat=rbind(vox_out,vox_dat)
  }
  print("Ranging cumulative viewshed data")
  output=aggregate(vox_dat[,4],by=list(vox_dat[,1],vox_dat[,2],vox_dat[,3]),FUN=sum)
  names(output) = c("X","Y","Z","Nshed")
  return(data.frame(output))
}
