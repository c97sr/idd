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

#' Solves the classic relationship between R0 and cumulative attack rate
#'
#' @param R0 basic reproductive number
#'
#' @example attack_rate(1.8)
#' @example attack_rate(4.0)
#' @example attack_rate(3.0)
attack_rate <- function(R0=1.8) {
  f <- function(a){1-exp(-R0*a)-a}
  tmp <- uniroot(f,c(0.000001,0.9999999),tol=0.001)
  tmp$root
}
