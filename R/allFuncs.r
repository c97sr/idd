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
circle_plot <- function(wdir='~/png') {

  require("shape")
  frames = 110
  rands <- data.frame(matrix(nrow=0,ncol=3))
  for ( i in 1:frames) {
    # creating a name for each plot file with leading zeros
    if (i < 10) {name = paste(wdir,'/000',i,'plot.png',sep='')}
    if (i < 100 && i >= 10) {name = paste(wdir,'/00',i,'plot.png', sep='')}
    if (i >= 100) {name = paste(wdir,'/0', i,'plot.png', sep='')}

    png(name)
    plot( x = NULL,xlim=c(log10(0.05), log10(1)), ylim=c(1,3), xlab =expression('p'['c']), ylab =expression("R"[0]),
          xaxt='n')
    dr <- 0.05
    shape::filledcircle(r1 = 0.25, mid = c( log10(0.0753), 1.35), dr = dr,
                 col = colorRampPalette(c(rgb(0,1,1,1), rgb(0,1,1,0)), alpha = TRUE)(10))
    text(x= log10(0.0753),y= 1.35, "H1N1 2009")
    shape::filledcircle(r1 = 0.45, mid = c( log10(0.6), 1.8), dr = dr,
                 col = colorRampPalette(c(rgb(0,1,0,1), rgb(0,1,0,0)), alpha = TRUE)(10))
    text(x= log10(.5),y= 1.8, "Ebola 2014")
    shape::filledcircle(r1=0.5,mid = c( log10(0.9),2.5), dr = dr,
                 col = colorRampPalette(c(rgb(1,1,0,1), rgb(1,1,0,0)), alpha = TRUE)(10))
    text(x= log10(0.8),y=2.6,"SARS 2002")
    shape::filledcircle(r1=0.4,mid = c( log10(0.2),1.1), dr = dr,
                 col = colorRampPalette(c(rgb(0.4,0,1,1), rgb(0.4,0,1,0)), alpha = TRUE)(10))
    text(x= log10(0.2),y=1.1,"MERS 2012")
    axis(side = 1, at = c( log10(0.01),log10(0.05),log10(0.1),log10(0.5), log10(1)),labels=c(0.01,0.05,0.1,0.5,1))

    if (i < 2)
    {

      rands= rbind(rands,data.frame(x=runif(1,min=log10(0.05),max=log10(0.1)),y=runif(1,1.5,2),start=(i-1)))
      #points(x=rands[,1], y=rands[,2], pch = 20, col = rgb(1,0,1,0.5),"cex" = 0.5)
    }

    else if (i < 5 ) {
      rands= rbind(rands,data.frame(x=runif(1,min=log10(0.1),max=log10(0.2)),y=runif(1,1.5,2),start=(i-1)))
      #points(x=rands[,1], y=rands[,2], pch = 20, col = rgb(1,0,1,0.5),"cex" = 0.5)
    }
    else if (i < 100 ){
      if ( i %% 4 == 1 )
        rands= rbind(rands,data.frame(x=runif(1,min=log10(0.2),max=log10(0.5)),y=runif(1,1.5,2),start=(i-1)))

      #points(x=rands[,1], y=rands[,2], pch = 20, col = rgb(1,0,1,0.5),"cex" = 0.5)
    }
    for ( j in seq_len(nrow(rands)) ) {
      #mid <- c(rands[j,1],rands[j,2])
      r <- i - rands[j,3]
      if (is.element(r, c(1:10))) {
        rands <- runif(1,0.8,1.1)
        shape::filledcircle(r1=((1)/(r)), mid = c(rands[j,1],rands[j,2]), dr = dr,
                     col = colorRampPalette(c(rgb(1,0,1,0.5), rgb(1,0,1,0)), alpha = TRUE)(10))
      }
      else
        points(x=rands[j,1], y=rands[j,2], pch = 20, col = "red","cex" = 0.9)

    }

    if( i > 105) {
      rows <- nrow(rands)
      shape::filledcircle(r1=(0.25),mid = c((sum(rands[,1])/rows),(sum(rands[,2])/rows)), dr = dr,
                   col = colorRampPalette(c(rgb(1,0,1,(i-104)/6), rgb(1,0,1,0)), alpha = TRUE)(10))
      text(x= log10(0.3),y=1.825,"New Pathogen")

    }
    dev.off()
  }


}

# Some script lines to be used for development
if (FALSE) {

  # This file can probably go now
  
  # Reinstall and reload the library
  install( ".",force=TRUE)
  library("idd")
}
