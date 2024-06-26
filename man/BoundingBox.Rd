% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundingBox.R
\name{BoundingBox}
\alias{BoundingBox}
\title{Create a bounding box for a study region.}
\usage{
BoundingBox(shapefile, d = 2, rotate = FALSE, verbose = FALSE)
}
\arguments{
\item{shapefile}{Spatial feature that defines the boundary of the area to define a bounding
box over.}

\item{d}{Dimension of the new Master Sample, at this stage we only work with d=2.}

\item{rotate}{Boolean of whether or not to randomly rotate the bounding box. This parameter
is not supported at this time.}

\item{verbose}{Print the rotation and random seed when it is generated.}
}
\value{
bounding box for a study area.
}
\description{
Randomly generate a seed from 10,000 possible values in right now 2 dimensions.
Note that in van Dam-Bates et al. (2018) we required that the random seed falls into main
object shape, such as one of the islands in New Zealand, or within marine environment for
BC west coast. However, with a random rotation, we are able to ignore that detail. If this
function is used without a random rotation, we recommend running it until
the first master sample point does indeed fall within the largest scale of the master sample use.
}
\details{
This function was first written by Paul van Dam-Bates for the
package BASMasterSample and later ported to this package, spbal.
}
\examples{
# Create a bounding box for the Gates, North Carolina study area -------------
# Use the North Carolina shapefile supplied in the sf R package.
shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
shp_gates <- shp_file[shp_file$NAME == "Gates",]
# Vertically aligned master sample bounding box.
bb <- spbal::BoundingBox(shapefile = shp_gates)
bb

}
