#' Computing parameters for the \code{\link{view}} function
#'
#' @description The \strong{view.param} function computes the parameters for the \code{\link{view}} function. Each line of the output is a direction corresponding
#' to a single slightline. The calibration is done by optimizing the filling of the surface of a sphere with circles, each direction being recorded as the 3D
#' polar coordinates of the circle center.
#' @param angular.res the angular resolution of  slightline that will be used in the \code{\link{view}} funtion. Default = 10.
#' @param method in some cases, when no space is available to insert a new circle, a hole might remain after the sphere surface filling.
#' Two options are available to correct this issue. If \emph{method = 'fill'} the gap is filled with a portion of circle and a new diection
#' is added but its angle with the neighboring direction is lower than \emph{angular.res}. If \emph{method = 'regularize} (the default value)
#' the directions are corrected so that they all have a similar angle between them. If \emph{method = 'n'} no action is taken.
#'
#' @return a data.frame containing the parameters to use as input in the \code{\link{view}} function.
#'
#' @note some existing parameters are already provided with the viewshed3d package with parameters for an angular resolution ranging from 1 to 10 : \link{param_1},
#' \link{param_2}, \link{param_3}, \link{param_4}, \link{param_5}, \link{param_6}, \link{param_7}, \link{param_8}, \link{param_9}, \link{param_10}.
#'
#' @export
#'
#' @import VoxR rgl
#'
#' @examples
#' \donttest{
#' library(viewshed3d)
#'
#' #- anglura.res = 40 is way to big. Just good for the example
#' param=view.param(angular.res = 40)
#'
#' head(param)
#' }
#'
view.param=function(angular.res,method){
  requireNamespace('VoxR')
  requireNamespace('rgl')

  if(missing(angular.res)){angular.res=10}
  if(missing(method)){method='regularize'}

  print("Creating the calibration sphere")
  a=seq(-4,4,0.02)
  sphere=expand.grid(a,a,a)
  sphere[,4]=point.distance(sphere,c(0,0,0))
  sphere=subset(sphere,sphere[,4]<=4)
  sphere=subset(sphere,sphere[,4]>3.98)
  data=sphere

  print("Start directions calculation")
  data[,4]=round(axis.angle(data,"Z")/angular.res)*angular.res
  data[,5]=0
  curID=1
  curang=0
  ang=sort(unique(data[,4]))
  row.names(data)=1:nrow(data)
  ang_coord=data.frame(NA,NA,NA) ; names(ang_coord)=c("ID","Xrot","Zrot") ; ang_coord=ang_coord[-1,]
  anginc=angular.res/10
  for(i in 1:length(ang)){
    inangZ=subset(data,data[,4]==ang[i])
    if(i == 1 | i==length(ang)){
      data[as.numeric(row.names(inangZ)),5]=curID
      ang_coord[nrow(ang_coord)+1,]=data.frame(curID,ang[i],0)
      curID=curID+1
    }else{
      run=T
      while(run==T){
        rotdat=as.matrix(inangZ[,1:3])
        rotdat=rotate3d(rotdat, pi*curang/180, 0, 0, 1)
        rotdat=subset(rotdat,rotdat[,2]<0)
        rotdat=rotate3d(rotdat, pi*ang[i]/180, 1, 0, 0)
        inangX=as.data.frame(rotdat)
        inangX[,4]=axis.angle(inangX,"Z")
        inangX=subset(inangX,inangX[,4]<=angular.res/2)
        inangX=data[as.numeric(row.names(inangX)),]
        test=max(inangX[,5])

        if(test==0){
          data[as.numeric(row.names(inangX)),5]=curID
          ang_coord[nrow(ang_coord)+1,]=data.frame(curID,ang[i],curang)
          if(curID%%50==0){print(paste(curID,"directions recorded."))}
          curID=curID+1
        }else{
          null=subset(inangX,inangX[,5]==0)
          nonull=subset(inangX,inangX[,5]!=0)
          if(min(nonull[,5])!=max(nonull[,5])){
            if(method=="fill"){
              data[as.numeric(row.names(null)),5]=curID
              ang_coord[nrow(ang_coord)+1,]=data.frame(curID,ang[i],curang)
            }
            run=F
            curID=curID+1
          }
        }
        curang=curang+anginc
      }
    }
  }
  if(method=="regularize"){
    print('correcting directions')
    row.names(ang_coord)=1:nrow(ang_coord)
    ang=unique(ang_coord[,2])
    for(i in 1:length(ang)){
      inang=subset(ang_coord,ang_coord[,2]==ang[i])
      if(nrow(inang)>1){
        a=c(0) ; a[2:nrow(inang)]=360/(nrow(inang))
        a=cumsum(a)
        ang_coord[as.numeric(row.names(inang)),3]=a+min(inang[,3])
      }
    }
  }
  head=data.frame(angular.res,NA,NA) ; names(head)=c("ID","Xrot","Zrot")
  ang_coord=rbind(head,ang_coord)
  print("Done")
  return(ang_coord)
}
