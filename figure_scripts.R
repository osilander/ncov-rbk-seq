setwd("~/Documents/Manuscripts/Current/coronavirus")
library(RColorBrewer)
d1_files <- dir("./data/all_amplicon_depth", pattern="txt", full.names=T)
d1_file.names <- dir("./data/all_amplicon_depth", pattern="txt", full.names=F)
d2_files <- dir("./data/patient_depth", pattern="txt", full.names=T)

### Figure 1 compares the amplicons of
### different length
fig1 <- 0
### Figure 2 is compares the amplicon coverage as you increase sampling
fig2 <- 0
### Figure 3 compares the coverage for five 
### patient samples
fig3 <- 1
### Figure 6 plots fraction of Tp SNPs 
fig4 <- 0
### Figure 6 plots fraction of genome with at least 10X or 100X coverage at
### different read numbers 
fig5 <- 0
### Figure 5 plots number of N's at
### different read numbers 
fig6 <- 0
### Figure 6 plots coverage as template is decreased
fig7 <- 0

##########################
#### start fig 1 plotting
##########################
if(fig1) {
	draw.tiles <- 0
	plot.col <- brewer.pal(11, name="Spectral")
	plot.col <- plot.col[c(1,2,4,5,7,8,10,11)]
	## order depends on how dir is read. We order the coverage plots
	### as: 400bp contigs, 1200bp contigs, 1500 bp, 2000bp, with the S1 and S2
	#o <- c(7,8,1,2,3,4,5,6)
	o <- c(1,5,2,6,3,7,4,8)
	plot.amps <- c("400 bp","400 bp","1200 bp","1200 bp","1500 bp","1500 bp","2000 bp","2000 bp")
	plot.nums <- rep(c(20.3,31.2),4)
	
	### amplicon lists
	amplicon.list <- list()
	## p1_400
	amplicon.list[[1]] <- c(30,410,642,1028,1242,1651,1875,2269,2505,2904,3144,3531,3771,4164,4294,4696,4939,5321,5563,5957,6167,6550,6718,7117,7305,7694,7943,8341,8595,8983,9204,9585,9784,10171,10362,10763,10999,11394,11555,11949,12110,12490,12710,13096,13319,13699,13918,14299,14545,14926,15171,15560,15827,16209,16416,16833,17065,17452,17674,18062,18253,18672,18896,19297,19548,19939,20172,20572,20786,21169,21357,21743,21961,22346,22516,22903,23122,23522,23789,24169,24391,24789,24978,25369,25601,25994,26197,26590,26835,27227,27446,27854,28081,28464,28677,29063,29288,29693)
	#p2_400
	amplicon.list[[2]] <- c(320,726,943,1337,1573,1964,2181,2592,2826,3210,3460,3853,4054,4450,4636,5017,5230,5644,5867,6272,6466,6873,7035,7415,7626,8019,8249,8661,8888,9271,9477,9858,10076,10459,10666,11074,11306,11693,11863,12256,12417,12802,13005,13400,13599,13984,14207,14601,14865,15246,15481,15886,16118,16510,16748,17152,17381,17761,17966,18348,18596,18979,19204,19616,19844,20255,20472,20890,21075,21455,21658,22038,22262,22650,22797,23214,23443,23847,24078,24467,24696,25076,25279,25673,25902,26315,26520,26913,27141,27533,27784,28172,28394,28779,28985,29378,29486,29866)
	#p1_1200
	amplicon.list[[3]] <- c(30,1205,2153,3257,4167,5359,6283,7401,8253,9400,10343,11469,12450,13621,14540,15735,16624,17754,18596,19678,20553,21642,22511,23631,24633,25790,26744,27894,28677,29790)
	#p2_1200
	amplicon.list[[4]] <- c(1100,2266,3144,4262,5257,6380,7298,8385,9303,10451,11372,12560,13509,14641,15608,16720,17622,18706,19574,20698,21532,22612,23518,24736,25690,26857,27784,29007)
	#p1_1500
	amplicon.list[[5]] <- c(30,1476,2630,4147,5423,6903,8072,9432,10556,12034,13362,14786,16108,17458,18596,19981,21142,22610,23892,25324,26626,28007,28331,29790)
	#p2_1500
	amplicon.list[[6]] <- c(1372,2737,4046,5548,6747,8184,9317,10688,11922,13460,14686,16208,17350,18706,19877,21241,22511,24002,25213,26728,27906,29378)
	#p1_2000 
	amplicon.list[[7]] <- c(30,2079,3771,5585,7298,9123,10885,12719,14476,16401,18167,19981,21694,23631,25359,27164,27872,29790)
	#p2_2000
	amplicon.list[[8]] <- c(1955,3878,5472,7401,9010,11013,12619,14575,16290,18275,19877,21797,23518,25491,27050,29045)
	
	d1_files <- d1_files[o]
	d1_file.names <- d1_file.names[o]
	all.depths <- list()
	all.scaled.depths <- list()
	pdf(file="./figures/Fig_1.pdf", height=11, width=8)
	par(mfrow=c(10,1))
	par(mar=c(0,5,0,2))
	par(las=1)
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(1,2),ylim=c(1,2),xlab="",ylab="")
	lab.pos <- c(1.95, 1.85, 1.75, 1.9, 1.90, 1.7, 1.9, 1.75)*1e3

	for(i in 1:length(d1_files)) {
		if(i %% 2) {
			par(mar=c(0,5,1,2))
		}
		else {
			par(mar=c(1,5,0,2))
			draw.tiles <- 1
		}
		d <- read.table(d1_files[i])
		d.runmed <- d
		d.runmed[,3] <- runmed(d[,3],11)
		mean.d <- mean(d.runmed[,3])
		scaled.d <- 1e3/mean.d
		d.runmed[,3] <- d.runmed[,3]*scaled.d

		all.scaled.depths[[d1_file.names[i]]] <- d.runmed[,3]

		plot.lim <- 2
		plot(d.runmed[,2], d.runmed[,3], ty="l", xlab="Position on genome", ylab="", xlim=c(50, 29850), xaxt="n", yaxt="n", lwd=0.5, bty="n",ylim=c(1,3e3))
		polygon(c(d.runmed[,2],rev(d.runmed[,2])),c(d.runmed[,3],rep(0,length(d.runmed[,3]))), col=plot.col[i], bty="n")
		axis(2,at=c(0,plot.lim*1e3),labels=paste(c(0,plot.lim),"K",sep=""), cex.axis=1.1)
		text(1200,lab.pos[i], labels=bquote(paste('C'['q'], .(plot.nums[i]), '  ', .(plot.amps[i]))), cex=1.1)
		if(i==5) {
			par(las=0)
			mtext("Coverage depth", side=2, line=3, cex=1, adj=-2.0 )
			par(las=1)
		}
		if(draw.tiles) {
			top.tile <- amplicon.list[[i-1]]
			bottom.tile <- amplicon.list[[i]]
			
			for(x in seq(1,length(top.tile),by=2)) {
				rect(top.tile[x],2600,top.tile[x+1],2700, col="light grey",lwd=0.5)
				rect(bottom.tile[x],2800,bottom.tile[x+1],2900, col="light grey",lwd=0.5)
				text((top.tile[x]+top.tile[x+1])/2, 2400,cex=0.7,labels=(x), font=2)
			}			
		}
		draw.tiles <- 0
	}
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(50, 29850),ylim=c(1,2),xlab="",ylab="")
	axis(3,at=seq(0,30e3,by=5e3),labels=paste(seq(0,30,by=5),"K",sep=""), line=-0.8, tcl=0.5, mgp=c(-3,-2,0), cex.axis=1.2)
	text(15e3,1.5,labels="Genomic position", cex=1.2)
	dev.off()
	#write.table(sample.stats,file="./figures/table_5.txt",sep="\t")
	cat("Fig 1\n")
}
##########################
#### end fig 1 plotting
##########################

##########################
#### start fig 2 plotting
##########################
if(fig2) {
	plot.col <- brewer.pal(11, name="Spectral")
	## we adjust a bit :)
	plot.col <- plot.col[c(1,2,4,5,7,8,10,11)]
	fract.cov <- matrix(nrow=length(d1_files), ncol=10)
	fract.cov <- as.data.frame(fract.cov)
	o <- c(7,8,1,2,3,4,5,6)
	d1_files <- dir("./data", pattern="_d1", full.names=T)
	d1_files <- d1_files[o]
	
	for (i in 1:length(d1_files)) {
		d <- read.table(file=d1_files[i])
		
		d <- d[which(d[,2]>150),]
		d <- d[which(d[,2]<29.5e3),]
		
		total.bp <- 29.35e3*mean(d[,3])
		for (j in 1:20) {
			d.sub <- d[,3]*j*1e6/total.bp
			fract.cov[i,j] <- 1- length(which(d.sub<30))/length(d.sub)
		}
		
	}
	pdf(file="./figures/Fig_2.pdf", height=4, width=5)
	par(las=1, mar=c(5,5,2,1))
	plot(-1,-1,xli=c(0,20),ylim=c(0.90,1), xlab="Mbp of sequence data", ylab="Fraction of bp with > 30X coverage")
	
	### this is seriously fiddly
	### first we change the labels
	i.lab <- gsub("./data/", "",d1_files)
	i.lab <- gsub("_bp_d1_", " bp ",i.lab)
	i.lab <- gsub(".txt", "",i.lab)
	i.lab <- gsub("S2", "",i.lab)
	i.lab <- gsub("S1", "",i.lab)
	### then we adjust the label placements
	i.lab[3] <- gsub(" \nCq20", "",i.lab[3])
	i.pos <- c(1,1,2,1,1,1,2,2)
	i.xpos <- c(5,5,2,5,5,5,5,5)
	i.offset <- rep(0.5,8)
	i.offset[3] <- 0.2
	i.offset[1] <- 0.4
	### then we make some adjustments on the white rectangles under the writing
	i.rect <- rep(0.009,8)
	i.rect[2] <- 0.014
	i.rect[6] <- 0.007
	### do the actual plotting
	for (i in 1:length(d1_files)) {
		points(1:20,fract.cov[i,],ty="o", bg=plot.col[i], pch=21, lwd=1.2, cex=1.1)
		if(i==2 | i==6) {
			rect(i.xpos[i]-1.5, fract.cov[i,i.xpos[i]]-i.rect[i], 5+0.5, fract.cov[i,5]-0.003,col="white", border="white")
			text(i.xpos[i], fract.cov[i,i.xpos[i]], labels=i.lab[i], cex=0.7, offset=i.offset[i], pos=i.pos[i])	
			text(i.xpos[i], fract.cov[i,i.xpos[i]]-0.006, labels=bquote(paste('C'['q']*'31')),cex=0.7, pos=i.pos[i], offset=0.4)
		}
		# again need to adjust labelling
		if(i==1 | i==3) {
			text(i.xpos[i], fract.cov[i,i.xpos[i]], labels=i.lab[i],cex=0.7, pos=i.pos[i], offset=i.offset[i])
			text(i.xpos[i], fract.cov[i,i.xpos[i]]-0.006, labels=bquote(paste('C'['q']*'20')),cex=0.7, pos=i.pos[i], offset=0.4)	
		}
	}
	dev.off()
	cat("Fig 2\n")

}
##########################
#### end fig 2 plotting
##########################



##########################
#### start fig 3 plotting
##########################
if(fig3) {
	plot.col <- brewer.pal(5, name="Spectral")
	
	amplicon.list <- list()
#p1_1200
	amplicon.list[[1]] <- c(30,1205,2153,3257,4167,5359,6283,7401,8253,9400,10343,11469,12450,13621,14540,15735,16624,17754,18596,19678,20553,21642,22511,23631,24633,25790,26744,27894,28677,29790)
#p2_1200
	amplicon.list[[2]] <- c(1100,2266,3144,4262,5257,6380,7298,8385,9303,10451,11372,12560,13509,14641,15608,16720,17622,18706,19574,20698,21532,22612,23518,24736,25690,26857,27784,29007)

	## order depends on how dir is read. We order the coverage plots by Cq
	### as: S1 S3 S5 S4 S2
	o <- c(1,3,5,4,2)
	plot.nums <- c(20.3,20.6,25.7,28.5,31.2)
	d2_files <- d2_files[o]
	pdf(file="./figures/Fig_3.pdf", height=6, width=8)
	par(mfrow=c(7,1))
	par(mar=c(0,5,0,2))
	par(las=1)
	plot(-1-100,xaxt="n",yaxt="n",bty="n",xlim=c(50, 29850),ylim=c(0,1000),xlab="",ylab="")


	lab.pos <- c(1.8e3, 800, 900, 750, 900)
	plot.lim <- c(3.5, 1.5, 1.5, 2, 1.5)
	for(i in 1:length(d2_files)) {
		d <- read.table(d2_files[i])
		d.runmed <- d
		d.runmed[,3] <- runmed(d[,3],11)
		
		plot(d.runmed[,2], d.runmed[,3], ty="l", xlab="Position on genome", ylab="", xlim=c(50, 29850), xaxt="n", yaxt="n", lwd=0.5, bty="n",ylim=c(1,plot.lim[i]*1e3))
		polygon(c(d.runmed[,2],rev(d.runmed[,2])),c(d.runmed[,3],rep(0,length(d.runmed[,3]))), col=plot.col[i], bty="n")
		axis(2,at=c(0,plot.lim[i]*1e3-500),labels=paste(c(0,plot.lim[i]),"K",sep=""), cex.axis=1.1)
		text(0,lab.pos[i],labels=bquote(paste('C'['q'], .(plot.nums[i]))), cex=1.1)
		
		# with sapce if you want
		#text(500,lab.pos[i],labels=bquote(paste('C'['q']*' ', .(plot.nums[i]))), cex=1.1)
		if(i==3) {
			par(las=0)
			mtext("Coverage depth", side=2, line=3.3, cex=0.8, adj=1.5 )
			par(las=1)
		}
		
		if(i==1) {
			### figure out tiles
		top.tile <- amplicon.list[[1]]
		bottom.tile <- amplicon.list[[2]]
		for(x in seq(1,length(top.tile),by=2)) {
			rect(top.tile[x],3000.,top.tile[x+1],3150, col="light grey",lwd=0.5)
			text((top.tile[x]+top.tile[x+1])/2, 2800,cex=0.7,labels=(x), font=2)
				
			rect(bottom.tile[x],3250,bottom.tile[x+1],3400, col="light grey",lwd=0.5)
			text((bottom.tile[x]+ bottom.tile[x+1])/2, 3050,cex=0.7,labels=(x+1), font=2)
				
	}
			
		}

	}
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(50, 29850),ylim=c(1,2),xlab="",ylab="")
	axis(3,at=seq(0,30e3,by=5e3),labels=paste(seq(0,30,by=5),"K",sep=""), line=-0.8, tcl=0.5, mgp=c(-3,-2,0), cex.axis=1.2)
	text(15e3,1.4,labels="Genomic position", cex=1.2)
	dev.off()
	cat("Fig 3\n")

}
##########################
#### end fig 3 plotting
##########################

##########################
#### start fig 4 plotting
##########################
if(fig4) {
	pdf(file="./figures/Fig_4.pdf", height=6, width=8)
	par(mfrow=c(2,2))
	min.depth <- c(10,30,50,100)
	if(0) {
		all.depth.files <- dir("./data", pattern="depth", full.names=T)
		barcodes <- paste("barcode0",1:5,sep="")
		#re-oreder by Cq
		barcodes <- barcodes[c(1,3,5,4,2)]
		bc.depth.files <- list()
		for(i in 1:length(barcodes)) {
			bc.depth.files[[i]] <- all.depth.files[grep(barcodes[i],all.depth.files)]
		}
		
		bc01 <- bc.depth.files[[1]]
		bc01.len <- gsub("./data/2020-05-05_2246_barcode01.","",bc01)
		bc01.len <- gsub(".sorted.depth.txt","",bc01.len)
		bc01.len <- as.numeric(bc01.len)
		o <- order(bc01.len)
		read.len <- bc01.len[o]
		### sort by length
		for(i in 1:length(barcodes)) {
			bc.depth.files[[i]] <- bc.depth.files[[i]][o]
		}
		
		bc.depths <- list()
		for(i in 1:length(barcodes)) {
			depth.mat <- matrix(nrow=length(read.len), ncol=(length(min.depth)+1))
			depth.mat[,1] <- read.len
			colnames(depth.mat) <- c("reads",paste("min_",min.depth,sep=""))
		
			for (m in 1:length(min.depth)) {
				for (f in 1:length(bc01)) {
					d <- read.table(file=bc.depth.files[[i]][f])
					d <- d[100:29500,]
					depth.mat[f,(m+1)] <- 1-length(which(d[,3] < min.depth[m]))/length(d[,3])
				}
			}
			bc.depths[[i]] <- depth.mat
		}
	}
	
	layout( matrix(c(1,1,3,3,1,2,3,4,5,5,7,7,5,6,7,8),byrow=T,nrow=4))
	par(mar=c(5,6,2,1))
	par(las=1)
	
	plot.col <- brewer.pal(length(barcodes), name="Spectral")
	y.lims <- c(0.96,0.55,0,0)
	legend.y <- c(0.995,0.95,0.9,0.85)
	for(m in 1:length(min.depth)) {
		par(mar=c(5,6,2,1))
	#### normal
		plot(-1,-1,xlim=c(-200,3e4),ylim=c(y.lims[m],1),ylab=paste("Fraction of basepairs\nwith > ",min.depth[m], "X coverage", sep=""),xlab="Number of subsampled reads",xaxt="n")
		plot.nums <- c(20.3,20.6,25.7,28.5,31.2)
		text.adj <- c(0.0002,-0.0005,-0.0015,0,0)
		for (i in 1:length(barcodes)) {
			points(bc.depths[[i]][,1], bc.depths[[i]][,m+1], pch=21, bg=plot.col[i], lwd=1.2, ty="o",cex=1.2)
			if(m==1) {
				text(500,bc.depths[[i]][1,m+1]+text.adj[i], labels=bquote(paste('C'['q'], .(plot.nums[i]))), cex=0.85, offset=0.23)
			}
		}
		axis(1,at=seq(0,3e4,by=5e3),labels=paste(seq(0,30,by=5),"K", sep=""))
		legend(20e3, legend.y[m], legend=plot.nums, bty="n",pch=21, pt.bg=plot.col, cex=0.9)

	### Zoom
	par(mar=c(7.5,3,0.,4))
	plot(-1,-1,xlim=c(0,3e4),ylim=c(0.995,1.0005),ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
	
	for (i in 1:length(min.depth)) {
		for (i in 1:length(barcodes)) {
			points(bc.depths[[i]][,1], bc.depths[[i]][,m+1], pch=21, bg=plot.col[i], lwd=1, ty="o",cex=1)
		}
	}
	axis(1,at=seq(0,3e4,by=1e4),labels=paste(seq(0,30,by=10),"K", sep=""),cex.axis=0.9,mgp=c(0.4,0.4,0))
	axis(2,at=c(0.995,1),labels=c("0.995","1"),cex.axis=0.9,mgp=c(0.6,0.6,0))
	#line=-0.8, tcl=0.5, mgp=c(-3,-2,0)
	}
	dev.off()
	cat("Fig 4\n")

}
##########################
#### end fig 4 plotting
##########################

##########################
#### start fig 5 plotting
##########################
if(fig5) {
	
	barcodes <- paste("barcode0",1:5,sep="")
	plot.col <- brewer.pal(length(barcodes), name="Spectral")
	n.count <- read.table(file="./data/ncount/all_counts.txt", header=T)
	pdf(file="./figures/Fig_5.pdf", height=3.5, width=5)
	
	layout( matrix(c(1,2,1,1),byrow=T,nrow=2))
	
	par(mar=c(5,6,1,1))
	par(las=1)
	plot(-1,-1e3,xlim=c(0,3e4),ylim=c(0,1.5e4),xaxt="n", yaxt="n",xlab="Number of subsampled reads", ylab="Number of ambiguous\nbases in genome")
	
	axis(1,at=seq(0,3e4,by=10e3),labels=paste(seq(0,30,by=10),"K", sep=""))
	axis(2,at=seq(0,15e3,by=5e3),labels=paste(seq(0,15,by=5),"K", sep=""))
	plot.nums <- c(20.3,20.6,25.7,28.5,31.2)
	text.adj <- c(0,-400,400,0,0)
	for(b in 1:length(barcodes)) {
		barcode.n <- subset(n.count, barcode==barcodes[b])
		o <- order(barcode.n$reads)
		barcode.n <- barcode.n[o,]
		points(barcode.n$reads, barcode.n$n.count, pch=21, bg=plot.col[b], lwd=1, ty="o",cex=1.1)
		text(800, barcode.n$n.count[1]+text.adj[b], labels=bquote(paste('C'['q'], .(plot.nums[b]))), cex=0.7)
	}
	
	
	### ZOOM
	par(mar=c(3,4,2.,3))
	plot(-1,-1e3,xlim=c(0,3e4),ylim=c(-2,1.5e2),xaxt="n",yaxt="n", ylab="",xlab="Reads",bty="n")
	
	axis(1,at=seq(0,3e4,by=10e3),labels=paste(seq(0,30,by=10),"K", sep=""), cex.axis=0.9,mgp=c(0.4,0.4,0))
	axis(2,at=seq(0,1.5e2,by=50),labels=seq(0,150,by=50), cex.axis=0.9,mgp=c(0.6,0.6,0))
	par(las=0)
	mtext("Ambiguous bases",2,line=1.9,cex=0.75)
	mtext("Reads",1,line=1.2,cex=0.75)
	for(b in 1:length(barcodes)) {
		barcode.n <- subset(n.count, barcode==barcodes[b])
		o <- order(barcode.n$reads)
		barcode.n <- barcode.n[o,]
		points(barcode.n$reads, barcode.n$n.count, pch=21, bg=plot.col[b], lwd=1, ty="o",cex=1)
		
	}
	dev.off()
	cat("Fig 5\n")

}
##########################
#### end fig 6 plotting
##########################

##########################
#### start fig 6 plotting
##########################
if(fig6) {
	pdf(file="./figures/Fig_6.pdf", height=3.5, width=4.5)
	d <- read.table("./data/all_collated_snps.txt")
	x <- d[,1]/20000*100
	### these are simply the pairwise comparisons
	colnames(d) <- c("1_2","1_4","1_5","2_1","2_4","2_5","4_1","4_2","4_5","5_1","5_2","5_4")
	par(las=1)
	par(mar=c(5,6,2,2))
	plot.col <- brewer.pal(11, name="Spectral")
	plot.col <- c(plot.col, plot.col[1])
	plot(x, d[,2], ylim=c(0,1), ty="o", cex=0.7, log="x", xlab="Percentage of contaminating reads", ylab="Fraction of true\npositive SNPs called", bg=plot.col[1], pch=21,xaxt="n")
	axis(1,at=c(0.2,0.5,1,2,5,10,20,50),labels=c("0.2","0.5","1","2","5","10","20","50"))
	for (i in 3:13) {
		points(x, d[,i], ty="o", cex=0.7, bg=plot.col[i-1], pch=21)
	}
	abline(v=6,lty=3,lwd=0.5)
	dev.off()
	cat("Fig 6\n")

}
##########################
#### end fig 6 plotting
##########################



##########################
#### start fig 7 plotting
### not used.
##########################
if(fig7) {
	plot.col <- colorRampPalette(c("light green", "light blue"))( 18 )
	plot.col <- sort(c(plot.col,plot.col))
	depth.files <- dir("~/Documents/Manuscripts/Current/coronavirus/data/depth",patter="depth", full.names=T)
	## order depends on how dir is read. We order the coverage plots by Cq
	### as: S1 S3 S5 S4 S2
	## Cq 24.3 is 4.2e4 copies into ARTIC (11ul)
	## 31.2 is almost exactly 120-fold lower, so 350 copies.
	## 10-fold dilution is 35
	## 100-fold dilution is 3.5
	## 1000-fold dilution is non-existent
	o <- c(1,3,5,4,2)
	plot.labels <- paste(rep(c("350","35","3.5","0.35","0.035","No"),3), "copies")
	pdf(file="./figures/Fig_7.pdf", height=6, width=14)
	
	par(mfrow=c(8,2))
	par(mar=c(0,5,0,2))
	par(las=1)
	
	#this is just empty space
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(1,2),ylim=c(1,2),xlab="",ylab="")
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(1,2),ylim=c(1,2),xlab="",ylab="")

	### plotting is in rows
	for(i in c(1,13,2,14,3,15,4,16,5,17,6,18)) {
		d <- read.table(depth.files[i])
		d.runmed <- d
		d.runmed[,3] <- runmed(d[,3],11)
		#mean.d <- mean(d.runmed[,3])
		#scaled.d <- 1e3/mean.d
		#d.runmed[,3] <- d.runmed[,3]*scaled.d
		#plot.lim <- 2
		plot(d.runmed[,2], d.runmed[,3], ty="l", xlab="Position on genome", ylab="", xlim=c(50, 29850), xaxt="n", lwd=0.5, bty="n",log="y",ylim=c(1,5000),yaxt="n")
		polygon(c(d.runmed[,2],rev(d.runmed[,2])),c(d.runmed[,3],rep(1,length(d.runmed[,3]))), col=plot.col[i], bty="n")
		axis(2,at=c(1,50,500),labels=c(1,50,500), cex.axis=1.1)
		text(2000,25,labels=plot.labels[i], cex=1.1)
	}
	
	### this is just empty space at the bottom of column 1
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(50, 29850),ylim=c(1,2),xlab="",ylab="")
	### but then we add a label and text
	axis(3,at=seq(0,30e3,by=5e3),labels=paste(seq(0,30,by=5),"K",sep=""), line=-0.8, tcl=0.5, mgp=c(-3,-2,0), cex.axis=1.2)
	text(15e3,1.4,labels="Genomic position", cex=1.2)
	
	### this is just empty space at the bottom of colun two
	plot(-1-1,xaxt="n",yaxt="n",bty="n",xlim=c(50, 29850),ylim=c(1,2),xlab="",ylab="")
	### but then we add a label and text
	axis(3,at=seq(0,30e3,by=5e3),labels=paste(seq(0,30,by=5),"K",sep=""), line=-0.8, tcl=0.5, mgp=c(-3,-2,0), cex.axis=1.2)
	text(15e3,1.4,labels="Genomic position", cex=1.2)
	dev.off()
	cat("Fig 7\n")

}
##########################
#### end fig 7 plotting
##########################