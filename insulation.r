#libraries
library(psych)

tab2matrices=function(tab) {
	mats=list()
	subs=unique(tab$ID)
	for (i in 1:length(subs)) {
		sub=tab[tab$ID==subs[i],]
		from=min(sub$coord_bin1)
		to=max(sub$coord_bin1)
		mat=matrix(0,ncol=to-from+1,nrow=to-from+1)
		colnames(mat)=from:to
		rownames(mat)=from:to
		cor=from-1
		for (j in 1:dim(sub)[1]) {
			mat[sub[j,"coord_bin1"]-cor,sub[j,"coord_bin2"]-cor]=sub[j,"val"]
		}
		mats[[subs[i]]]=mat
	}
	return(mats)
}


insulation=function(mat,w) {
	ins=rep(NA,ncol(mat))
	names(ins)=colnames(mat)
	for (d in (w+1):(length(ins)-w-1)) {
		sub=mat[d:(d+w-1),(d-w):(d-1)]
		ins[d]=sum(sub)
	}
	return(ins)
}
#for a 5kb matrix and a 25kb window, w=5
wholeinsulation=function(mats,w) {
	ins=list()
	for (i in 1:length(mats)) {
		ins[[i]]=insulation(mats[[i]],w)
		names(ins)[i]=names(mats)[i]
	}
	tot=NULL
	for (i in 1:length(ins)) {
		tot=c(tot,ins[[i]])
	}
	gm=geometric.mean(tot)
	for (i in 1:length(ins)) {
		v=log2(ins[[i]]/gm)
		ins[[i]]=v
	}
	return(ins)
}

instobed=function(tab,ins,ofn) {
	d=unique(tab$coord_bin1)
	bed=NULL
	for (i in 1:length(d)) {
		sel=tab[tab$coord_bin1==d[i],]
		id=sel[1,"ID"]
		insa=ins[[id]]
		insa=insa[names(insa)==d[i]]
		tabl=data.frame("chr"=sel[1,"chr"],"ID"=id,"from"=sel[1,"from1"],"to"=sel[1,"to1"],"ins"=insa)
		bed=rbind(bed,tabl)
	}
	write.table(bed,ofn,col.names=T,row.names=F,quote=F,sep="\t")
}

inspipe=function(ifn,ofn,w) {
	tab=read.table(ifn,header=T,stringsAsFactors=F)
	names(tab)[9]="val"
	mats=tab2matrices(tab)
	ins=wholeinsulation(mats,w)
	instobed(tab,ins,ofn)
}


collateinsulation=function(folder) {
	files=list.files(folder)
	files=files[grep("insulation",files)]
	table=data.frame()
	for (i in 1:length(files)) {
		fn=paste0(folder,files[i])
		tab=read.table(fn,header=T,stringsAsFactors=F)
		win=unlist(strsplit(fn,"\\.|_"))
		win=win[length(win)-1]
		names(tab)[5]=as.character(win)
		if (i==1) {
			table = tab
		} else {
			table = merge(table,tab)
		}
	}
	ofn=paste0(folder,"allins.txt")
	write.table(table,ofn,col.names=T,row.names=F,quote=F,sep="\t")
	return(table)	
}


