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

# Translates a value from a unit scale to either a log or a linear scale
# x is the number on the unit scale
# min is the minimum value of the non-unit scale
# max is the maximum value of the non-unit scale
# logbase is the base of the log scale if used
# log is a boolean with TRUE for a log scale and FALSE for a linear scale
# Returns the value on the non-unit scale
unit_scale_from <- function(x,min=1,max=100,logbase=10,logflag=FALSE) {
  if (logflag) {
    rtn <- min*logbase^(x*(log(max,logbase)-log(min,logbase)))
  } else {
    rtn <- min + (max-min)*x
  }
  rtn
}
