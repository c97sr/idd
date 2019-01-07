# Copyright Steven Riley (sr@stevenriley.net)
#
# This file is part of the library idd.
#
# idd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This work is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

#' Proposes parameter updates for an MCMC algorithm
#'
#' Uses either linear or log-linear random walks to propose 
#' multivariate jumps in parameter space.
#'
#' @param pt a table of parameter names, values and min and maxes. 
#' See Details.
#' @param fm a vector of strings of the parameters within the table 
#' to be fitted.
#' The default value is to assume all the parameters are to be fitted
#'
#' @details The table of parameter names must have the correct columns for
#' this to work... 
#'
#' @return A vector of proposed parameter values. The vector is the 
#' same length as ptab so it copies in the values of the parameters 
#' that are not being updated.
mh_update <- function(pt,fmask=1:(dim(pt)[1])) {

	bps <- pt[,"val"]
	rtn <- bps
	nofit <- length(rtn)

	for (i in 1:nofit) {

		# Set up and transform to unit scale
		rv <- runif(1)
		rv <- (rv-0.5)*pt[i,"step"]
		x <- bps[fmask[i]]
		x <- SR_to_unit(x,min=pt[i,"min"],max=pt[i,"max"],logflag=pt[i,"log"])
		x <- x + rv

		# Cyclical boundary conditions
		if (x < 0) x <- 1 + x
		if (x > 1) x <- x - 1

		# Test for errors and return to originl scales
		if (x < 0 || x > 1) stop("problem here")
		rtn[fmask[i]] <- SR_from_unit(
		  x,min=pt[i,"min"],max=pt[i,"max"],logflag=pt[i,"log"]
		)

	}

	rtn

}

