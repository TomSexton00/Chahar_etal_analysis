#!/bin/env Rscript

options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

#Arguments
args = commandArgs(T)
if (length(args)==0) {
	cat(sprintf("usage: %s <capture table> <o_contact file> <bin file> <field name> <output file prefix>\n", script.name))
	q(status=1)
}

cap.fn = args[1]
ifn = args[2]
bin.fn = args[3]
field.name = paste0(args[4],"_normalised")
ofn.prefix = args[5]

#Libraries
library(Matrix)


#Functions

loss.function = function(mat) {
	tmp = rowSums(mat)
	tmp = tmp[tmp > 0]
	rs = abs(tmp-1)
	tmp = colSums(mat)
	tmp = tmp[tmp > 0]
	cs = abs(tmp-1)
	return( 0.5*sum(c(rs,cs)))
}

create.complete.matrix = function(x) {
	t.sum = as.data.frame(summary(x), stringsAsFactors=FALSE)
	test = sparseMatrix(i = t.sum$i, j = t.sum$j, x = t.sum$x, dims = dim(x), symmetric = TRUE)
	rownames(test)=rownames(x)
	colnames(test)=colnames(x)
	return (test)
}

IPF.alg = function(lfm, numberOfIterations) {
	d = dim(lfm)[2]
	x = matrix(1, nrow=d, ncol=numberOfIterations)
	y = matrix(1, nrow=d, ncol=numberOfIterations)
	lf = numeric(numberOfIterations)
	counter = 1
	while (counter <= numberOfIterations) {
		s.temp = rowSums(lfm)
		x[,counter] = ifelse(s.temp>0, s.temp, 1)
		lfm = lfm/x[,counter]
		s.temp = colSums(lfm)
		y[,counter] = ifelse(s.temp>0, s.temp, 1)
		lfm = t(t(lfm)/y[,counter])
		lf = c(lf, loss.function(lfm))
		print(paste("loss function: ", lf[length(lf)], "counter: ",counter))
		counter = counter + 1
	}
	return(list(lfm=lfm, x=x, y=y, lf=lf))
}

#Main
cap = read.table(cap.fn,header=TRUE)
obs = read.table(ifn,header=TRUE)

for (i in 1:dim(cap)[1]) {
	cat(sprintf("Normalising sub-matrix at chr %s : %s \n", cap[i,"chr"], cap[i, "gene"]))
	tab = obs[obs$chr1==cap[i,"chr"] & obs$from1 > cap[i,"start"] & obs$to1 < cap[i,"end"] & obs$chr2==cap[i,"chr"] & obs$from2 > cap[i,"start"] & obs$to2 < cap[i,"end"],]
	#resMatrix = Matrix(0, nrow=max(tab$coord_bin1), ncol = max(tab$coord_bin1), sparse=TRUE)
	resMatrix = Matrix(0,nrow=pmax(max(tab$coord_bin1),max(tab$coord_bin2)),ncol=pmax(max(tab$coord_bin1),max(tab$coord_bin2)),sparse=TRUE)
	resMatrix[ cbind(tab$coord_bin1,tab$coord_bin2) ] = tab$observed_count
	cm = create.complete.matrix(resMatrix)
	tp = IPF.alg(cm, numberOfIterations=20)
	x = tp$x
	y = tp$y
	tp.x = data.frame(x)
	test.x = exp(Reduce('+',log(tp.x)))
	tp.y = data.frame(y)
	test.y = exp(Reduce('+',log(tp.y)))
	fac = (test.x/test.y)
	van = rowSums(cm)==0
	test.x[van]=fac
	test.x=test.x/sqrt(fac)
	tp = cm/test.x
	tp = t(t(tp)/test.x)
	idx = which(tp!=0,arr.ind=TRUE)
	c=cbind(idx,tp[idx])
	c=as.data.frame(c)
	names(c) = c("coord_bin1","coord_bin2",field.name)
	c = c[order(c$coord_bin1),]
	tmp.fn = paste0(ofn.prefix,"tmp.txt")
	write.table(c,tmp.fn,col.names=T,row.names=F,quote=F,sep="\t")
	ofn = paste0(ofn.prefix,"_",cap[i,"gene"],"_normalised.txt")
	
	command = paste0("lscripts/cbin2coords.pl ",tmp.fn, " ",bin.fn, " ",cap[i,"chr"], " ", cap[i,"gene"], " ", field.name, " ", ofn)
	cat(sprintf("%s\n", command))
	if (system(command) != 0) {
		stop("error in cbins2coord script\n")
	}
	system(paste("rm",tmp.fn)) 
}
cat("Combining sub-matrices...\n")
comb = NULL
for (i in 1:dim(cap)[1]) {
	ifn = paste0(ofn.prefix,"_",cap[i,"gene"],"_normalised.txt")
	c = read.table(ifn,header=T)
	if (i == 1) {
		comb = c
	}
	else {
		comb = rbind(comb,c)
	}
}
ofn = paste0(ofn.prefix,"_all_normalised.txt")
write.table(comb,ofn,col.names=T,row.names=F,quote=F,sep="\t")
