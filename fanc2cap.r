#!/shared/software/miniconda3/envs/r_3.3.2/bin/Rscript

options(warn=-1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <chr> <ID> <mat> <reg> <ofn>\n", script.name))
  q(status=1) 
}

chr = args[1]
id = args[2]
mat.fn = args[3]
reg.fn = args[4]
ofn = args[5]

mat=read.table(mat.fn,header=F)
reg=read.table(reg.fn,header=F)

if ((dim(mat)[1] != dim(mat)[2]) | (dim(mat)[1] != dim(reg)[1])) {
	print("Must be square matrix with all rows annotated\n")
	q(status=1)
}

tab=NULL
for (i in 1:dim(reg)[1]) {
	tmp=data.frame("chr"=chr,"ID"=id,"coord_bin1"=i,"coord_bin2"=1:dim(reg)[1],"from1"=reg[i,2],"to1"=reg[i,3],"from2"=NA,"to2"=NA,"score"=NA)
	for (j in 1:dim(reg)[1]) {
		tmp[j,"from2"]=reg[j,2]
		tmp[j,"to2"]=reg[j,3]
		tmp[j,"score"]=mat[i,j]
	}
	tab=rbind(tab,tmp)
}

write.table(tab,ofn,col.names=T,row.names=F,quote=F,sep="\t")


