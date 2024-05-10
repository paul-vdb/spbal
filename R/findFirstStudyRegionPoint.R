# findFirstStudyRegionPoint.R

#' @name findFirstStudyRegionPoint
#'
#' @title Get a randomly chosen Halton point from within the study area and the associated seeds.
#'
#' @description This function repeatedly calls function spbal::getBASSample
#' to generate the Halton frame sample. This function selects the first point at random from those
#' points in the study area. This point and the seeds used to generate the sample are returned to
#' the caller.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) of the study area(s).
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly using function \code{generateUVector}.
#' @param bb Bounding box which defines the Master Sample. A bounding box must be
#' supplied.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing three variables:
#'
#' \itemize{
#' \item \code{seeds} The u1 and u2 seeds used to generate the first point.
#' \item \code{k} The index of the first point in the initial sample.
#' }
#'
# 1. Set J1 = 4 and J2 = 3.
# 2. Generate B = 2^J1 x 3^J2 points from a random-start Halton sequence H
#    with a random seed (u1, u2).
# 3. Find points from H in the study area. Call this set S. If S is empty,
#    increment J1 and J2 and go to step 2.
# 4. Randomly choose a point from S. Let xk be this point where k is the 'site index'
#    (I think that's what we call it).
# 5. Set the seed to (u1 + k - 1, u2 + k - 1).
# 6. Re-number ID by subtracting k (re-generate sample using seeds for 5 - first point must also be ID=1)

# For example, let (u1, u2) = (1, 5) and S = {x2, x6, x7}.
# If x6 is randomly chosen, then the new seed is (1 + 6 - 1, 5 + 6 - 1) = (6, 10)
# (the sixth point in H).

# The only difference is that the random-start Halton sequence must be length B.
#' @keywords internal
findFirstStudyRegionPoint <- function(shapefile, bb, seeds, verbose = FALSE){
  # must not be called without seeds! (also checked in getBASSample).
  if(base::is.null(seeds)){
    msg <- "spbal(findFirstStudyRegionPoint) The seeds parameter must not be NULL."
    msgs <- base::sprintf(msg)
    base::stop(msgs)
  }

  # Initialise variables.
  J <- base::c(4, 3)
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)

  # default number of sample points to find.
  n <- (bases[1]^J[1]) * (bases[2]^J[2])

  pts_in_intersection <- 0
  call.getBASSample.cnt <- 0

  while(pts_in_intersection < 1){
    # shapefile, bb, n, seeds
    call.getBASSample.cnt <- call.getBASSample.cnt + 1
    result <- getBASSample(shapefile = shapefile, bb = bb, n = n, seeds = seeds)
    diff_ <- result$sample
    seeds <- result$seed

    # find number of points within our study area/bb intersection.
    pts_in_intersection <- base::length(diff_$SiteID)
    n <- n * 2
  }

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Needed %s call(s) to locate first study area point."
    msgs <- base::sprintf(msg, call.getBASSample.cnt)
    base::message(msgs)
  }

  # select a point in the study area at random.
  base::set.seed(seeds[1] + seeds[2])
  k <- base::sample(pts_in_intersection, 1)
  k <- diff_$SiteID[k]

  # select our random first point. # return SiteID = 1, know what k is.
  #first.pt <- diff_ #[1,]

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Random point selected: %s."
    msgs <- base::sprintf(msg, k)
    base::message(msgs)
  }

  result <- base::list(seeds = seeds,
                       k     = k)
  return(result)
}

## Is this very efficient? Could definitely do something a bit smarter.
## Keep for now to have a basic simple random sample of a polygon.
SRSPoly <- function(n = 1, shapefile, bb, verbose = FALSE){
  bb.bounds <- sf::st_bbox(bb)
  n_found <- 0
  ndraw <- n + 10
  while(n_found < n){
    xy <- cbind(runif(ndraw, bb.bounds["xmin"], bb.bounds["xmax"]), runif(ndraw, bb.bounds["ymin"], bb.bounds["ymax"]))
    pts.coord <- sf::st_as_sf(base::data.frame(SiteID = 1:ndraw, xy), coords = c(2, 3))
    
    sf::st_crs(pts.coord) <- sf::st_crs(bb)
    # find the intersection. Generates the same as sf::st_intersection(pts.coord, shapefile)
    if(n_found == 0) {
      pts.intersect <- pts.coord[shapefile,]
    }else{ 
      pts.intersect <- rbind( pts.intersect, pts.coord[shapefile,] )
    }
    n_found <- nrow(pts.intersect)
    ndraw <- ndraw*2
  }
  pts.intersect$SiteID <- 1:n_found
  return(pts.intersect[1:n,])
}

SRSHaltonGrid <- function(n = 1, shapefile, bb, J = c(3,2), verbose = FALSE){

  bases <- c(2,3)
  Bxy <- bases^J

  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]
  
  ## Choose centre points of a Halton Grid to sample: 
  ## This ensures that the centroid of the random seed falls into the polygon.
  grdx <- ((1:Bxy[1])/Bxy[1] - 0.5/Bxy[1])*scale.bas[1] + shift.bas[1]
  grdy <- ((1:Bxy[2])/Bxy[2] - 0.5/Bxy[2])*scale.bas[2] + shift.bas[2]

  n_found <- 0
  ndraw <- n + 5
  call.intersect.cnt <- 0

  while(n_found < n & call.intersect.cnt < 50){
    xy <- cbind(sample(grdx, ndraw, replace = TRUE), sample(grdy, ndraw, replace = TRUE))
    pts.coord <- sf::st_as_sf(base::data.frame(SiteID = 1:ndraw, xy), coords = c(2, 3))
    
    sf::st_crs(pts.coord) <- sf::st_crs(bb)
    # find the intersection. Generates the same as sf::st_intersection(pts.coord, shapefile)
    if(n_found == 0) {
      pts.intersect <- pts.coord[shapefile,]
    }else{ 
      pts.intersect <- rbind( pts.intersect, pts.coord[shapefile,] )
    }
    call.intersect.cnt <- call.intersect.cnt+1
    n_found <- nrow(pts.intersect)
    ndraw <- ndraw*2

    if(verbose){
      msg <- "spbal(SRSHaltonGrid) after SRSHaltonGrid n_samples = %s. num_samples = %s"
      msgs <- base::sprintf(msg, call.intersect.cnt, n_found)
      base::message(msgs)
    }
    
  }
  ## Brute force if it didn't work the first time before giving up.
  if(n_found == 0){
    if(verbose){
      msg <- "spbal(SRSHaltonGrid) Using brute force to find a initial point."
      base::message(msg)
    }
  
    grd <- expand.grid(x = grdx, y = grdy)
    xy <- cbind(sample(grdx, ndraw, replace = TRUE), sample(grdy, ndraw, replace = TRUE))
    pts.coord <- sf::st_as_sf(base::data.frame(SiteID = 1:ndraw, xy), coords = c(2, 3))
    sf::st_crs(pts.coord) <- sf::st_crs(bb)
    pts.intersect <- pts.coord[shapefile,]
    if(nrow(pts.intersect) == 0) stop("Warning: Not able to generate any sample sites in this polygon.")
  }
  if(verbose){
    msg <- "spbal(SRSHaltonGrid) Needed %s call(s) to obtain %s samples."
    msgs <- base::sprintf(msg, call.intersect.cnt, n_found)
    base::message(msgs)
  }    
  pts.intersect$SiteID <- 1:n_found
  return(pts.intersect[1:n,])
}

## This should work well for both HaltonFrames and BAS.
setHaltonSeed <- function(shapefile, bb, J = NULL, verbose = FALSE){
  
  ## Check to see if we want to make sure the grid centre point falls in object.
  haltonFrame <- TRUE
  
  # Initialise variables.
  if( is.null(J) ){
    J <- base::c(7, 5)
    haltonFrame <- FALSE ## BAS if J in NULL
  }
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)
  Bxy <- bases^J
  
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  ## Get a single random sample from the polygon.
  ## We ensure this is a good seed for the Halton Grid by sampling it explicitly.
  pts.unif <- SRSHaltonGrid(n=1, shapefile = shapefile, bb = bb, J = J, verbose) %>% 
    st_coordinates()
  upts <- (pts.unif - shift.bas)/scale.bas

  pts_in_intersection <- 0
  call.getBASSample.cnt <- 0
  seeds0 <- c(0,0)

  ## If the first random point isn't in the polygon, then subdivide again by increasing J.
  while(pts_in_intersection < 1){
    Axy <- cbind(floor((upts[,1] + 2*.Machine$double.eps)*bases[1]^J[1]), 
                 floor((upts[,2] + 2*.Machine$double.eps)*bases[2]^J[2]))
    ix <- round(cppBASpts(n = Bxy[1], seeds = 0, bases = bases[1], FALSE)$pts*Bxy[1],0)  ## Make sure it's an integer
    iy <- round(cppBASpts(n = Bxy[2], seeds = 0, bases = bases[2], FALSE)$pts*Bxy[2],0)
    seeds0[1] <- which(ix == Axy[1]) - 1  ## Starts at 0.
    seeds0[2] <- which(iy == Axy[2]) - 1
    
    U <- sample(1000, 2, replace = TRUE)
    seeds <- seeds0 + U*Bxy
    
    call.getBASSample.cnt <- call.getBASSample.cnt + 1
    result <- spbal:::getBASSample(shapefile = shapefile, bb = bb, n = n, seeds = seeds)
    diff_ <- result$sample
    seeds <- result$seed
    
    pts_in_intersection <- base::length(diff_$SiteID)
    
    ## Keep boxes reasonably similar in size.
    if(bases[1]^J[1] <= bases[2]^J[2]){
      J[1] <- J[1] + 1
    }else{
      J[2] <- J[2] + 1
    }
    Bxy <- bases^J
   }

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Needed %s call(s) to locate first study area point."
    msgs <- base::sprintf(msg, call.getBASSample.cnt)
    base::message(msgs)
  }

  result <- base::list(seeds = seeds,
                       k     = 1) ## Keep k for backward compatability.
  return(result)
}