# Copyright Steven Riley (sr@stevenriley.net)

# This file is part of the library idsource.

# idsource is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This work is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this idsource.  If not, see <http://www.gnu.org/licenses/>. 

# You may need something simimilar to setwd("/Users/sriley/Dropbox/talks/20120309_MSc_Adv_Spatial/")
# A function to load up a file and prepare objects to be used for simulation
adm.load.sheep.prep.sim <- function(
		filename="msc_farm_data.csv",
		max_kernel_dist=15000,
		simsteps=0:100) {
	
	# Load up the file
	dat <- read.csv(filename)
	nofarms <- dim(dat)[1]
	
	# Add dummy rows for the different elements of the simulation
	dat <- cbind(seqid=1:nofarms,dat,g_i=rep(-1,nofarms))
	
	dat
	
}

adm.setup.aux.mat <- function(df) {
	nf 	<- dim(df)[1]
	rtn <- matrix(nrow=nf,ncol=nf)
	rtn[,] <- 0
	for (i in 1:(nf-1)) {
		for (j in (i+1):nf) {
			x1 		<- df[i,"easting"]
			y1 		<- df[i,"northing"]
			x2 		<- df[j,"easting"]
			y2 		<- df[j,"northing"]
			dist <- sqrt((x1-x2)^2+(y1-y2)^2)
			rtn[i,j] <- dist
			rtn[j,i] <- dist
		}
	}
	rtn
}

# Function to see farms randomly near to the coast
# Returns the value of seeds as UIDs
adm.seed.sheep <- function(tab,nseed=1,sdsq=c(0,10000000,0,1000000)) {
	subtab <- tab[	tab$easting > sdsq[1] &
					tab$easting < sdsq[2] &
					tab$northing > sdsq[3] &
					tab$northing < sdsq[4],]
	nopotseeds <- dim(subtab)[1]
	if (nseed > nopotseeds) stop("adm.seed.sheep: not enough potential seeds")
	x <- c()
	while (length(x) < nseed) {
		randno <- ceiling(runif(1)*nopotseeds)
		randindex <- subtab$seqid[randno]
		if (!(randindex %in% x)) x <- c(x,randindex)
	}
	x
}

adm.plot.farm.and.seeds <- function(tab,vecseeds) {
	subtab <- tab[tab$seqid %in% vecseeds,]
	plot(tab$easting,tab$northing)
	points(subtab$easting,subtab$northing,col="red",cex=1,pch=19)
}

# Assumes the existance of farmdf in an enclosing environment
# Essentially a void function so returns a 0 for successful execution
adm.apply.seed <- function(s,df,newgen=0) {
	nf <- dim(df)[1]
	rtn <- rep(-1,nf)
	rtn[s] <- newgen
	rtn
}

# plot((1:100)*1000,sapply((1:100)*1000,adm.offset.kernel),type="l",log="x")
adm.offset.kernel <- function(d,p=adm.sheep.params()) {
	rtn <- 0
	if (d < p["cutoff"]) rtn <- 1/(1+(d/p["offset"])^p["power"])
	rtn
}

# This function assumes the existance of a farmdf object in the calling scope
adm.apply.gen.model <- function(p,df,distmatrix=NULL) {
	nf <- dim(df)[1]
	epsilon <- 1e-10
	rtn <- df[,"g_i"]
	lastgen <- max(rtn)
	gen <- lastgen + 1
	for (i in 1:nf) {
		if (abs(df[i,"g_i"] - (gen-1)) < epsilon) {
			for (j in 1:nf) {
				if ((j != i) & (df[j,"g_i"] < 0)) {
					if (is.null(distmatrix)) {
						x1 		<- farmdf[i,"easting"]
						y1 		<- farmdf[i,"northing"]
						x2 		<- farmdf[j,"easting"]
						y2 		<- farmdf[j,"northing"]
						dist <- sqrt((x1-x2)^2+(y1-y2)^2)
					}
					else {
						dist <- distmatrix[i,j]
					}
					if (dist < p["cutoff"]) {
						if (p["beta"]*adm.offset.kernel(dist) > runif(1)) {
							rtn[j] <- gen
							# cat(i,j,"\n")
							# flush.console()
						}
					}
				}
			}
		}
	}	
	rtn
}


adm.plot.gen <- function(datf,gens=c(0,1,2,3),cols=c("red","blue","green","orange")) {
	plot(datf[,"easting"],datf[,"northing"])
	for (i in 1:length(gens)) points(datf[datf$g_i==gens[i],"easting"],datf[datf$g_i==gens[i],"northing"],pch=19,col=cols[i])
}

# Function to setup a parameter object for the smv simulations
adm.sheep.params <- function() {
	rtn <- c(	offset=1000,
				power=2,
				cutoff=10000,
				beta=0.05)
	rtn
}
