library(gplots)
library(rtracklayer)
epigenome=list()
chr=NA
start=NA
end=NA

make.color.panel=function(cols,ncols=256) {
	panel=NULL
	for (i in 2:length(cols)) {
		panel=c(panel,colorpanel(ncols,cols[i-1],cols[i]))
	}
	return (panel)
}

get.contact.values=function(values,breaks) {
	min.value=breaks[1]
	max.value=breaks[length(breaks)]
	values[values<min.value]=min.value
	values[values>max.value]=max.value
	return (values)
}

vals.to.cols=function(vals,breaks,ncols) {
	min=breaks[1]
	max=breaks[length(breaks)]
	n=length(breaks)-1
	cols=rep(-1,length(vals))
	for (i in 1:n) {
		ind = (breaks[i]<=vals) & (vals<=breaks[i+1])
		if(!any(ind)) {
			next
		}
		cols[ind]=(vals[ind]-breaks[i])/(breaks[i+1]-breaks[i])
		cols[ind]=(i-1)*ncols+cols[ind]*(ncols-1)+1
		cols[ind]=round(cols[ind])
	}
	return (cols)
}

get.contact.colors=function(values,breaks,colors,ncols,panel) {
	values=get.contact.values(values,breaks)
	if(length(breaks)!=length(colors)) {
		stop("number of breaks and colors must be equal\n")
	}
	col=vals.to.cols(values,breaks,ncols)
	if(!all(is.finite(col))) {
		stop("missing values\n")
	}
	panel[col]
}

parking=function(left,right) {
	y=rep(-1,length(right))
	lengths=right-left
	for (i in order(lengths,decreasing=TRUE)) {
		otherleft=left
		otherleft[i]=NA
		otherright=right
		otherright[i]=NA
		placed=FALSE
		y[i]=0
		while(placed==FALSE) {
			placed=sum((right[i]>otherleft[y==y[i]])&(left[i]<otherright[y==y[i]]),na.rm=TRUE)==0
			if(placed==FALSE) {
				y[i]=y[i]+1
			}
		}
	}
	return(y)
}

plot_genes=function(genome,chr,start,end) {
	gene=genome[genome$Chr==chr & genome$Start > start & genome$End < end,]
	if (dim(gene)[1]>0) {	
		y_plot=parking(gene$Start,gene$End)
		plot(c(start,end),c(1,-max(y_plot)-0.5),col="white",ylab="",xlab="",fg="white",col.axis="white",xaxs="i",yaxs="i")
		arrowHeads=pretty(start:end,n=50)
		for (i in 1:dim(gene)[1]) {
			x=c(gene$Start[i],arrowHeads[arrowHeads>gene$Start[i]&arrowHeads<gene$End[i]],gene$End[i])
			if(gene$Strand[i]=="-") {
				arrows(x[2:length(x)],-y_plot[i],x[1:length(x)-1],col="blue",length=0.08)
			}
			else {
				arrows(x[1:length(x)-1],-y_plot[i],x[2:length(x)],col="blue",length=0.08)
			}
			text(gene$Start[i],-y_plot[i]+0.4,adj=0,labels=gene$Name[i])
		}
	}
}

import_bw = function(dir,files) {
	for (i in 1:length(files)) {
		file=paste0(dir,files[i],".bw")
		cat(paste("Reading file:",file,"\n"))
		epigenome[[files[i]]]<<-import(file)
		cat("Finished reading file\n")
	}
}

plot_tracks = function(tracks,chr,min,max,plotlevels,plotcols) {
	levels=unique(plotlevels[order(plotlevels)])
	for (i in 1:length(levels)) {
		leveltracks=list()
		levelkey=list()
		plotlim=numeric(2)
		counter=1
		for (j in 1:length(tracks)) {
			if (plotlevels[j]==levels[i]) {
				tmp=epigenome[[tracks[j]]]
				tmp=tmp[as.character(seqnames(tmp))==paste0("chr",chr) & start(tmp)>=min & end(tmp) <=max,]
				plotlim[1]=min(c(plotlim[1],score(tmp)))
				plotlim[2]=max(c(plotlim[2],score(tmp)))
				leveltracks[[counter]]=tmp
				levelkey[[counter]]=tracks[j]
				counter=counter+1
			}
		}
		for (j in 1:length(leveltracks)) {
			#par(mar=c(0,2,0,2)+0.1)
			plot.new()
			plot.window(xlim=c(min,max),ylim=plotlim,xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
			segments(start(leveltracks[[j]]),0,end(leveltracks[[j]]),score(leveltracks[[j]]),col=plotcols[j])
			title(main=levelkey[[j]],col.main=plotcols[j],adj=0,line=-2)
		}
	}
}

plot_track=function(track,chr,start,end,name,col) {
	tmp=track[as.character(seqnames(track))==chr & start(track)>start & end(track)<end,]
	plot.new()
	plot.window(xlim=c(start,end),ylim=c(0,max(score(tmp))),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	segments(start(tmp),0,end(tmp),score(tmp),col=col)
	title(main=name,col.main=col,adj=0,line=-2)
}

#assumes tab is normalised matrix - chr,ID,coord_bin1,coord_bin2,from1,to1,from2,to2,[normvalue field]
#If arrowhead-called TADs included, assumes tad is data.frame - ID,chr,start,end
triangular=function(tab,id,field,eps=NA,png=NA,vmax,tad=NA,mode="hic",main) {
	breaks=NULL
	ncols=256
	colors=NULL
	if(mode=="hic") {
		breaks=c(0,(vmax/4),(vmax/2),vmax)
		colors=c("white","orange","red","black")
	} else if (mode=="diff") {
		breaks=c(-vmax,0,vmax)
		colors=c("blue","gray","red")
	} else {
		cat("Bad mode given\n")
		return(NULL)
	}
	panel=make.color.panel(colors,ncols)

	xx=NULL
	yy=NULL
	polcol=NULL
	
	tab=tab[tab$ID==id,]
	heads=names(tab)
	ind=which(heads==field)
	if(length(ind)!=1) {
		stop("Bad field name given - must be present and unique\n")
	}
	vals=tab[,ind]
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	xmin=min(xleft)
	xright=tab$to1
	xmax=max(xright)
	ybottom=tab$from2
	ytop=tab$to2

	for (i in 1:length(xleft)) {
		if(xleft[i]<ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}
	if(!is.na(eps)) {
		setEPS()
		postscript(eps)
	} else if(!is.na(png)) {
		png(png,width=1500,height=1000)
	}

	
	plot(xx,yy,xlim=c(xmin,xmax),type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F,main=main)
	polygon(x=xx,y=yy,col=polcol,border=NA)
	axis(1,at=pretty(xx),labels=as.character(pretty(xx)),line=0)

	if(!is.na(tad)) {
		tad=tad[tad$ID==id,]
		xlef=tad$start
		xrigh=tad$start
		ybott=tad$end
		yto=tad$end
		for (i in 1:length(xlef)) {
			polygon(x=c(xlef[i],(yto[i]+xlef[i])/2,yto[i]),y=c(0,(yto[i]-xlef[i])/2,0),col=NA,border="green")
		}
	}


	if(!is.na(eps) | !is.na(png)) {
		dev.off()
	}
}

bar=function(vmax,eps=NA,png=NA) {
	ncols=256
	colors=c("white","orange","red","black")
	panel=make.color.panel(colors,ncols)
	n=length(panel)
	if(!is.na(eps)) {
		setEPS()
		postscript(eps)
	} else if (!is.na(png)) {
		png(png,width=1500,height=50)
	}
	plot(rep(1,n),col=panel,type="h",xlim=c(0,n),ylim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	labs=signif(c(0,(vmax/4),(vmax/2),(vmax*0.75),vmax),digits=3)
	axis(1,at=c(0,(n/4),(n/2),(n*0.75),n),labels=as.character(labs),line=0)
}

#Plot insulation score domainogram, assumes ins data.frame with chr,ID,from,to,insulation scores from increasing window sizes
allins=function(ins,id,vmax) {
	ncols=256
	colors=c("blue","gray","red")
	panel=make.color.panel(colors,ncols)
	breaks=c(-vmax,0,vmax)
	tab=ins[ins$ID==id,]
	n=dim(tab)[2]
	mat=as.matrix(tab[,5:n])
	n=n-4
	plot(c(min(tab$from),max(tab$to)),c(0,n),type="n",xlab="Genomic coordinate",ylab="",yaxt="n",xaxs="i",yaxs="i")
	for (i in 1:n) {
		cols=NULL
		vals=mat[,i]
		for (j in 1:length(vals)) {
			if(is.na(vals[j])) {
				cols=c(cols,"white")
			} else {
				cols=c(cols,get.contact.colors(vals[j],breaks,colors,ncols,panel))
			}
		}
		rect(tab$from,i-1,tab$to,i,col=cols,border=NA)
	}
}

#Assumes values are always ninth column
pairwise=function(tab1,tab2,id,vmax,cond1,cond2) {
	breaks=c(0,(vmax/4),(vmax/2),vmax)
	ncols=256
	colors=c("white","orange","red","black")
	panel=make.color.panel(colors,ncols)

	layout(matrix(1:2,ncol=1),heights=c(100,100))
	par(mar=c(2,2,2,2))

	tab=tab1[tab1$ID==id,]
	vals=tab[,9]
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	xright=tab$to1
	ybottom=tab$from2
	ytop=tab$to2

	xx=NULL
	yy=NULL
	polcol=NULL

	for (i in 1:length(xleft)) {
		if(xleft[i]<ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}

	title=paste0(cond1," ",id," - chr",tab[1,"chr"])
	plot(xx,yy,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F,main=title)
	polygon(x=xx,y=yy,col=polcol,border=NA)
	axis(1,at=pretty(xx),labels=as.character(pretty(xx)),line=0)

	tab=tab2[tab2$ID==id,]
	vals=tab[,9]
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	xright=tab$to1
	ybottom=tab$from2
	ytop=tab$to2

	xx=NULL
	yy=NULL
	polcol=NULL

	for (i in 1:length(xleft)) {
		if(xleft[i]>ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}

	title=paste0(cond2," ",id," - chr",tab[1,"chr"])
	plot(xx,yy,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F,main=title)
	polygon(x=xx,y=yy,col=polcol,border=NA)
} 

plot.ins=function(ins,id,conds,cols) {
	ins=ins[ins$ID==id,]
	ins$mid=ins$from+((ins$to-ins$from)/2)
	ind=names(ins) %in% conds
	tot=unlist(ins[,ind])
	tot=tot[!is.na(tot)]
	yrange=range(tot)
	ind=which(ind)
	plot(ins$mid,ins[,ind[1]],type="l",col=cols[1],xlim=c(start,end),ylim=yrange,xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	for (i in 2:length(ind)) {
		lines(ins$mid,ins[,ind[i]],col=cols[i])
	}
	abline(h=0)
	axis(2,at=pretty(tot),labels=as.character(pretty(tot)))
}

tri=function(tab,id,vmax,main,x=FALSE,rev=FALSE,diff=FALSE) {
	breaks=NULL
	ncols=256
	colors=NULL
	if(diff==TRUE) {
		breaks=c(-vmax,0,vmax)
		colors=c("blue","gray","red")
	} else {
	breaks=c(0,(vmax/4),(vmax/2),vmax)
	colors=c("white","orange","red","black")
	}
	
	panel=make.color.panel(colors,ncols)

	xx=NULL
	yy=NULL
	polcol=NULL
	
	tab=tab[tab$ID==id,]
	chr <<- paste0("chr",tab$chr[1])
	
	if (diff==TRUE) {
		vals=tab$z
	} else {
		vals=tab[,9]
	}
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	start <<- min(xleft)
	xright=tab$to1
	end <<- max(xright)
	ybottom=tab$from2
	ytop=tab$to2
	
	if(rev==FALSE) {

	for (i in 1:length(xleft)) {
		if(xleft[i]<ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}
	par(mar=c(1,3,1,1))
	plot(xx,yy,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F,main=main)
	polygon(x=xx,y=yy,col=polcol,border=NA)
	if(x==T) {
		axis(1,at=pretty(xx),labels=as.character(pretty(xx)),line=0)
	}
	} else {
		for (i in 1:length(xleft)) {
		if(xleft[i]>ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}
	plot(xx,yy,type="n",xlim=c(start,end),xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F)
	polygon(x=xx,y=yy,col=polcol,border=NA)
	}
}


triple=function(tab1,tab2,tab3,id,cond1,cond2,cond3,vmax,genes) {
	layout(matrix(1:4,ncol=1),heights=c(50,50,50,10))
	tri(tab=tab1,id=id,vmax=vmax,main=cond1)
	tri(tab=tab2,id=id,vmax=vmax,main=cond2)
	tri(tab=tab3,id=id,vmax=vmax,main=cond3)
	par(mar=c(1,3,1,1))
	plot_genes(genes,chr,start,end)
}


pairwise.ins=function(tab1,tab2,id,vmax,ins1,ins2,genes,main) {
	breaks=c(0,(vmax/4),(vmax/2),vmax)
	ncols=256
	colors=c("white","orange","red","black")
	panel=make.color.panel(colors,ncols)

	layout(matrix(1:4,ncol=1),heights=c(50,10,10,50))
	par(mar=c(1,3,1,1))

	tab=tab1[tab1$ID==id,]
	vals=tab[,9]
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	xmin=min(xleft)
	xright=tab$to1
	xmax=max(xright)
	ybottom=tab$from2
	ytop=tab$to2

	xx=NULL
	yy=NULL
	polcol=NULL

	for (i in 1:length(xleft)) {
		if(xleft[i]<ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}

	plot(xx,yy,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F,main=main,xlim=c(xmin,xmax))
	polygon(x=xx,y=yy,col=polcol,border=NA)
	axis(1,at=pretty(xx),labels=as.character(pretty(xx)),line=0)

	par(mar=c(1,3,1,1))
	plot_genes(genome=genes,chr=paste0("chr",tab[1,"chr"]),start=xmin,end=xmax)

	ins1=ins1[ins1$ID==id,]
	ins2=ins2[ins2$ID==id,]
	ins1$mid=ins1$from+((ins1$to-ins1$from)/2)
	ins2$mid=ins2$from+((ins2$to-ins2$from)/2)
	tot=c(ins1$ins,ins2$ins)
	tot=tot[!is.na(tot)]
	
	par(mar=c(1,3,1,1))
	plot(range(ins$mid),range(tot),type="n",xlim=c(xmin,xmax),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
	lines(ins1$mid,ins1$ins)
	lines(ins2$mid,ins2$ins,lty=2)
	abline(h=0)
	axis(2,at=pretty(range(tot)),labels=as.character(pretty(range(tot))),line=0)

	tab=tab2[tab2$ID==id,]
	vals=tab[,9]
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	xright=tab$to1
	ybottom=tab$from2
	ytop=tab$to2

	xx=NULL
	yy=NULL
	polcol=NULL

	for (i in 1:length(xleft)) {
		if(xleft[i]>ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}
	par(mar=c(1,3,1,1))
	plot(xx,yy,type="n",xlim=c(xmin,xmax),xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F)
	polygon(x=xx,y=yy,col=polcol,border=NA)
} 

tab3ins=function(tab,id,vmax,ins1,ins2,ins3,cols=c("darkgreen","blue","red"),genes,main) {
	ins1=ins1[ins1$ID==id,]
	ins2=ins2[ins2$ID==id,]
	ins3=ins3[ins3$ID==id,]
	ins1$mid=ins1$from+((ins1$to-ins1$from)/2)
	ins2$mid=ins2$from+((ins2$to-ins2$from)/2)
	ins3$mid=ins3$from+((ins3$to-ins3$from)/2)
	insall=c(ins1$ins,ins2$ins,ins3$ins)
	insall=insall[!is.na(insall)]
	layout(matrix(1:3,ncol=1),heights=c(100,10,40))
	par(mar=c(1,3,1,1))
	tri(tab=tab,id=id,vmax=vmax,main=main,x=T)
	par(mar=c(1,3,1,1))
	plot_genes(genes,chr,start,end)
	par(mar=c(1,3,1,1))
	plot(ins1$mid,ins1$ins,type="l",col=cols[1],xlim=c(start,end),ylim=range(insall),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	lines(ins2$mid,ins2$ins,col=cols[2])
	lines(ins3$mid,ins3$ins,col=cols[3])
	abline(h=0)
	axis(2,at=pretty(range(insall)),labels=as.character(pretty(range(insall))),line=0)
}

inscall=function(tab,id,cell,vmax,main,genes,scores,threshold=NA) {
	layout(matrix(1:7,ncol=1),heights=c(50,rep(10,6)))
	par(mar=c(1,3,1,1))
	tri(tab,id,vmax,main,x=T)
	par(mar=c(1,3,1,1))
	plot_genes(genes,chr,start,end)
	sub=scores[[cell]][[id]]
	for (i in 1:length(sub)) {
		par(mar=c(1,3,1,1))
		plot(sub[[i]]$mid,sub[[i]]$score,type="h",xlim=c(start,end),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
		ypos=max(sub[[i]]$score)*0.8
		text(start+20000,ypos,names(sub)[i])
		if(!is.na(threshold)) {
			abline(h=threshold,lty=2)
		}
	}
}

tabz3ins=function(tab,id,vmax,ins1,ins2,ins3,cols=c("darkgreen","blue","red"),genes,main) {
	ins1=ins1[ins1$ID==id,]
	ins2=ins2[ins2$ID==id,]
	ins3=ins3[ins3$ID==id,]
	ins1$mid=ins1$from+((ins1$to-ins1$from)/2)
	ins2$mid=ins2$from+((ins2$to-ins2$from)/2)
	ins3$mid=ins3$from+((ins3$to-ins3$from)/2)
	insall=c(ins1$ins,ins2$ins,ins3$ins)
	insall=insall[!is.na(insall)]
	layout(matrix(1:3,ncol=1),heights=c(100,10,40))
	par(mar=c(1,3,1,1))

	ncols=256
	breaks=c(-vmax,0,vmax)
	colors=c("blue","gray","red")
	panel=make.color.panel(colors,ncols)
	xx=NULL
	yy=NULL
	polcol=NULL
	
	tab=tab[tab$ID==id,]
	ind=which(names(tab)=="z")
	if(length(ind)!=1) {
		stop("z field must be present and unique\n")
	}
	vals=tab[,ind]
	chr <<- paste0("chr",tab$chr[1])
	vals[is.na(vals)]=0
	vals=get.contact.values(vals,breaks)
	col=get.contact.colors(vals,breaks,colors,ncols,panel)
	xleft=tab$from1
	start <<- min(xleft)
	xright=tab$to1
	end <<- max(xright)
	ybottom=tab$from2
	ytop=tab$to2

	for (i in 1:length(xleft)) {
		if(xleft[i]<ybottom[i]) {next}
		xx=c(xx,NA,(xleft[i]+ybottom[i])/2, (xleft[i]+ytop[i])/2, (xright[i]+ytop[i])/2, (xright[i]+ybottom[i])/2)
		yy=c(yy,NA,(xleft[i]-ybottom[i])/2, (xleft[i]-ytop[i])/2, (xright[i]-ytop[i])/2, (xright[i]-ybottom[i])/2)
		polcol=c(polcol,col[i])
	}
	
	plot(xx,yy,xlim=c(start,end),type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",ylab="",xlab="",axes=F,main=main)
	polygon(x=xx,y=yy,col=polcol,border=NA)
	axis(1,at=pretty(xx),labels=as.character(pretty(xx)),line=0)

	par(mar=c(1,3,1,1))
	plot_genes(genes,chr,start,end)
	par(mar=c(1,3,1,1))
	plot(ins1$mid,ins1$ins,type="l",col=cols[1],xlim=c(start,end),ylim=range(insall),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	lines(ins2$mid,ins2$ins,col=cols[2])
	lines(ins3$mid,ins3$ins,col=cols[3])
	abline(h=0)
	axis(2,at=pretty(range(insall)),labels=as.character(pretty(range(insall))),line=0)
}


bordercall=function(tab,id,vmax,main,bor,cutoff=0.095,freq=2,ins,genes) {
	layout(matrix(1:3,ncol=1),heights=c(50,10,40))
	par(mar=c(1,3,1,1))
	tri(tab,id,vmax,main,x=T)
	bor=bor[bor$ID==id,]
	ind=grep("score",names(bor))
	mat=as.matrix(bor[,ind])
	x=rowSums(mat>cutoff,na.rm=T)
	abline(v=bor[x>freq,"mid"],lty=2)
	par(mar=c(1,3,1,1))
	plot_genes(genes,chr,start,end)
	par(mar=c(1,3,1,1))
	ins=ins[ins$ID==id,]
	ind=grep("w",names(ins))
	tot=unlist(c(ins[,ind]))
	tot=tot[!is.na(tot)]
	plot(ins$mid,ins[,ind[1]],type="l",xlim=c(start,end),ylim=range(tot),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	abline(h=0)
	for (i in 2:length(ind)) {
		lines(ins$mid,ins[,ind[i]],col=i)
	}
}


bordercalldiff=function(tab1,tab2,cond1,cond2,id,vmax=0.05,ins,bor,q,genes,inscols=c("blue","red")) {
	layout(matrix(1:4,ncol=1),heights=c(50,10,20,50))
	par(mar=c(1,3,1,1))
	tri(tab1,id=id,vmax=vmax,main=paste(cond1,id),x=T)
	bor=bor[bor$ID==id & bor$q<q,]
	bor$mid=bor$from+((bor$to-bor$from)/2)
	abline(v=bor$mid,lty=2)
	par(mar=c(1,3,1,1))
	plot_genes(genes,chr,start,end)
	par(mar=c(1,3,1,1))
	ins=ins[ins$ID==id,]
	ins$mid=ins$from+((ins$to-ins$from)/2)
	tot=c(ins[,cond1],ins[,cond2])
	tot=tot[!is.na(tot)]
	plot(c(start,end),range(tot),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
	lines(ins$mid,ins[,cond1],col=inscols[1])
	lines(ins$mid,ins[,cond2],col=inscols[2])
	abline(h=0)
	abline(v=bor$mid,lty=2)
	axis(2,at=pretty(range(tot)),labels=as.character(pretty(range(tot))))
	par(mar=c(1,3,1,1))
	tri(tab2,id=id,vmax=vmax,main=paste(cond1,id),rev=T)
	abline(v=bor$mid,lty=2)
}
