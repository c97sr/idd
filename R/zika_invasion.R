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

#' Illustration of possible zika invasion dynamics
#'
#' Simulates a stochastic compartmental model to reproduce an
#' illustration of how two very similar communities might have
#' very different initial epidemics.
#'
#' @param its number of iterations for the second population
#' @param file name of file to write pdf or "screen"
#' @param N size of population
#'
#' @examples 
#' zika_invasion(its=100)
#' zika_invasion(its=100,N=1000000)
#' @export
zika_invasion <- function(
    its=10,
    ## file="screen",
    N=1000000,
    seed=19283719,
    highlights=c(1,2,3),
    highcolors=c("red","magenta","cyan"),
    maxts=360*2,
    R0axes=0:4,
    incaxes=seq(0,3.5,0.5)) {
  set.seed(seed)
  nosens <- length(highlights)
  senstab <- array(dim=c(its,maxts+1))
  transval <- vector(mode="numeric",length=maxts+1)
  De <- 11
  Tg <- 16
  Di <- Tg - De
  r1 <- comp_seir(
    De=De, Tg=Tg, A=0.9, I0 = 10, trickle=0/7,
    noTimeSteps=maxts,deterministic=FALSE,R0=2,N=N
  )
  ## if (file != "screen") pdf_fig_proof(file=file,pw=9)
  plot(r1$beta*Di,type="n",ylim=c(min(R0axes),max(R0axes)),axes=FALSE,xlab="",ylab="")
  axis(4,at=R0axes,col=grey.colors(10)[5],las=1)
  polygon(
      c(0:maxts,maxts:0),c(r1$beta*Di,rep(0,maxts+1)),
      col=grey.colors(10)[5],border=NA)
  par(new=TRUE)
  plot(r1$inf_inc/N*100,type="l",col="black",
       xlab="Day", ylab="Incidence",
       lwd=2,axes=FALSE,ylim=c(min(incaxes),max(incaxes)))
  axis(1)
  axis(2,las=1,at=incaxes)
  for (i in 1:its) {
    r2 <- comp_seir(
      De=De, Tg=Tg, A=0.9, I0 = 0, trickle=3/7,
      noTimeSteps=maxts,deterministic=FALSE,
      trickleStart=70,R0=2,N=N)
      points(r2$inf_inc/N*100,type="l",col=grey.colors(10)[9])
      senstab[i,] <- r2$inf_inc/N*100
  }
  for (i in 1:length(highlights)) {
      points(senstab[highlights[i],],type="l",col=highcolors[i],lwd=2)
  }
  ## if (file != "screen") dev.off()
}
