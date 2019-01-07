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
# GNU General Public License for more details

# Translates a number from a log or linear scale to a unit scale
# x is the number on the non-unit scale
# min is the assumed minimum value of the non-unit scale
# max is the assumed maximum value of the non-unit scale
# logbase is the base of the log scale if used
# log is a boolean for whether the scale is log or linear
# Returns the value on the unit scale
unit_scale_to <- function(y,min=1,max=100,logbase=10,logflag=FALSE) {
  if (logflag) {
    rtn <- (log(y,logbase)-log(min,logbase))/(log(max,logbase)-log(min,logbase))
  } else {
    rtn <- (y-min)/(max-min) 
  }
  rtn
}