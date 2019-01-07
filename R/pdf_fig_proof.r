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

#' Opens a pdf device with sensible options
#'
#' Starts a new pdf device with the width and height of the plotting
#' area specified in real distance units and with a specified margin
#' measured from the plotting area to the edge of the pdf bounding box.
#' This should make it easier to compose multi part figure. The device
#' that gets opened with this still needs to be closed.
#'
#' @param findex is a name index to make it easier to generate many plots
#' using the same basic file stem.
#' @param file the filename for the pdf
#' @param pw plotting area width in cm
#' @param ph plotting area height in cm
#' @param textscale scale factor for axis labels
#'
#' @example pdf_fig_proof()
#' 
#' @export 
pdf_fig_proof <- function(
  findex=1,
  file=paste("./pdf_sub_",findex,".pdf",sep=""),
  pw=7,
  ph=7,
  textscale=0.6,
  xpd=NA) {
    plotwidth <- pw/cm(1)
    plotheight <- ph/cm(1)
    margin <- 5/cm(1)
    pdfwidth <- plotwidth+2*margin
    pdfheight <- plotheight+2*margin
    posplot = c(
      margin/pdfwidth,
      (margin+plotwidth)/pdfwidth,
      margin/pdfheight,
      (margin+plotheight)/pdfheight
    )
    pdf(file=file,height=pdfheight,width=pdfwidth,useDingbats=FALSE)
    par(
      cex=textscale,
      fig=posplot,
      mai=c(0,0,0,0),
      xpd=xpd)
}
