% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct_ground.R
\name{reconstruct_ground}
\alias{reconstruct_ground}
\title{Optimal ground reconstruction for visibility computation}
\usage{
reconstruct_ground(
  data,
  ground_res,
  position,
  angular_res,
  method,
  full_raster,
  ...
)
}
\arguments{
\item{data}{LAS class object containing a 3d point cloud + a Classification
field that classes points as ground and non-ground, as provided by the
\code{\link[lidR]{classify_ground}} function from the
\code{\link[lidR]{lidR-package}}.}

\item{ground_res}{numeric. The grid resolution to reconstruct the ground on
the entire 3D scene. Default = 0.05. NOTE: a if needed, second grid may be added
with smaller (internally computed) resolution.}

\item{position}{vector of length 3 containing the xyz coordinates of the
animal position when computing the visibility with the
\code{\link{visibility}} function. Default = c(0,0,0).}

\item{angular_res}{numeric. The angular resolution of sightlines when
computing the visibility with the \code{\link{visibility}} function.
Default = 1.}

\item{method}{which algorithm to use for spatial interpolation. Can be
"knnidw", "tin" or "kriging". See documentation from the
\code{\link[lidR]{lidR-package}} for \code{\link[lidR]{knnidw}},
\code{\link[lidR]{tin}} and \code{\link[lidR]{kriging}}.}

\item{full_raster}{should the entire raster be interpolated for the ground
portion around the animal position? Parameter passed to the
\code{\link[lidR]{grid_terrain}} function available in the
\code{\link[lidR]{lidR-package}}.}

\item{...}{other arguments to pass to the spatial interpolation algorithm.
See documentation from \code{\link[lidR]{knnidw}}, \code{\link[lidR]{tin}}
and \code{\link[lidR]{kriging}}}
}
\value{
A LAS class object containing the 3D point cloud coordinates
with the ground reconstructed to be passed directly to the
\code{\link{visibility}} function. Note: the Classification field is
preserved.
}
\description{
Reconstructs the ground surface with a grid resolution defined
by the user and adds a second grid around the animal position with an
optimal resolution so that no sightline can pass through the ground when
computing visibility with the \code{\link{visibility}} function.
}
\examples{
\donttest{
#- import the tree_line_plot dataset
file <- system.file("extdata", "tree_line_plot.laz", package="viewshed3d")
tls <- lidR::readLAS(file,select="xyz")

#- class ground and vegetation points
class <- lidR::classify_ground(tls, lidR::csf(rigidness = 1L,
                                        class_threshold = 0.2,
                                        sloop_smooth = FALSE))

#- reconstruct the ground. Here the ground is reconstructed with the user
#- defined resolution only.
recons <- viewshed3d::reconstruct_ground(data=class,position = c(0,0,3),
                                         ground_res = 0.05,
                                         angular_res = 2,
                                         method="knnidw")

lidR::plot(recons,color="Classification",
           colorPalette = c("darkgreen","chocolate4"))

#- when the position is closer to the ground, the user defined resolution is
#- not sufficient and a second grid is added with the optimal resolution so
#- that no sightline can pass trough the ground when computing visibility.
#- In this example, full_raster = TRUE was used as a portion of the ground
#- near the animal location is not reconstructed because of a data
#- gap around a TLS scan position when using full_raster = FALSE.
recons <- viewshed3d::reconstruct_ground(data=class,position = c(0,0,1),
                                         ground_res = 0.05,
                                         angular_res = 2,
                                         method="knnidw",
                                         full_raster = TRUE)

lidR::plot(recons,color="Classification",
           colorPalette = c("darkgreen","chocolate4"))
}
}
