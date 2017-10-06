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

#' Solves the basic individual-based serial interval model
#'
#' Simulates a sequence of individual infection events.
#' Infections are drawn from simple offspring and serial
#' interval distributions.
#'
#' @param Tg generation time
#' @param R0 the basic reproductive number
#' @param Nmax the maximum number if infections in a cluster
#'
#' @return Returns a list of times of infection and indices.
#' The position in the list of infection times corresponds to
#' the index given in the second argument.
#'
#' XXXX This function not fully implemented yet XXXX
#'
#' @examples ind.tau()
ind_tau <- function(
	Tg=5,
	R0=2.6,
	Nmax=30
) {

  # Define the outputs to be returned
  rtn_times <- vector(mode="numeric",length=Nmax)
  rtn_infectors <- vector(mode="numeric",length=Nmax)

  # Initiate the
  rtn_times[1] <- 0
  rtn_infectors[1] <- 9999

  completed_infector <- 0
  current_infector <- 1

  # Main loop to go through every infector
  while (completed_infector <= Nmax) {

    # Not yet completed
    completed_infector <- completed_infector + 1
  }

  # Return list of infectors and times of infection
  list(t=rtn_times, i = rtn_infectors)

}
