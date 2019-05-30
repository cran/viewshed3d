\name{viewshed3d-package}
\alias{viewshed3d-package}
\alias{viewshed3d}
\docType{package}
\title{
The \strong{viewshed3d} package to compute visibility in 3D scenes of ecosystems acquired with Terrestrial Laser Scanner.
}
\description{
For many animals, the ability to visually assess the environment and detect approaching predators is an important part of anti-predator strategies.  Because this can occur across spatial scales, estimation of the viewshed can help to quantify visibility as a continuous variable around animal locations and facilitate studies of habitat selection and predator-prey interactions. \cr
The viewshed3d package provides tools to compute viewshed in 3D scenes of an ecosystem acquired with a Terrestrial Laser Scanner (TLS). The main function of the package, \code{\link{view}}, is designed to sample the point cloud in every direction of the 3D space from a single user-defined location and to record the distance to the nearest point in each direction. Each direction is thus considered as a slightline - of a user defined angle - that is assumed to end when an object - of undefined nature and size - is encountered. The \code{\link{cum.view}} function computes the overlap between viewsheds calculated from different locations. The other functions of the package are designed to properly prepare the data before the visibility calculation. First, the ground has to be reconstructed as occlusion might induce holes that would result in infinite slightlines to continue below the ground surface. To do so, the \code{\link{class.ground}} and \code{\link{ground.filter}} functions can be used to segment a TLS scene into ground and vegetation points. The \code{\link{reconstruct.ground}} and r\code{\link{reconstruct.ground.par}} functions can then be used to fill the holes by reconstructing the ground point cloud. Additionaly, the \code{\link{reshape.scene}} function can be used at any time during the process to reshape the point cloud and define its central location (that will be used as the starting point for slightlines) in order to speed up the subsequent processes if the original point cloud is very large. Finally, the \code{\link{view.param}} function is provided to build the parameter data required to run the \code{\link{view}} function. Considering the relatively long processing time required to run the \emph{view.param} function, some predefined parameter data are provided for different slightline angles.
}
\author{
Bastien Lecigne and Jan Eitel

Maintainer : Bastien Lecigne <lecignebastien@gmail.com>
}

\keyword{package}
\keyword{3d point cloud}
\keyword{T-LiDAR}
\keyword{Terrestrial Laser Scanner}
\keyword{Viewshed}
\keyword{ground vegetation segmentation}

