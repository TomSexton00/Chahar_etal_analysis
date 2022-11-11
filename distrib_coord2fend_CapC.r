#!/shared/software/miniconda3/envs/r_3.3.2/bin/Rscript

options(warn=-1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <dataset> <input dir> <probe.fends table> <read length> <facing threshold> <max segment length> <split dir> <output dir> <qsub directory> <working directory> <required memory>\n", script.name))
  q(status=1) 
}

dataset = args[1]
input.dir = args[2]
fends.fn = args[3]
read.length = as.numeric(args[4])
facing.threshold = as.numeric(args[5])
max.seg.len = as.numeric(args[6])
split.dir = args[7]
output.dir = args[8]
qsub.dir = args[9]
wd = args[10]
req.mem=as.numeric(args[11])

qsub.dir=paste0(qsub.dir,"_coord2fend")
cat(sprintf("cleaning qsub dir: %s\n", qsub.dir))
system(paste("rm -rf", qsub.dir))
system(paste("mkdir -p", qsub.dir))

cat(sprintf("cleaning split dir: %s\n", split.dir))
system(paste("rm -rf", split.dir))
system(paste("mkdir -p", split.dir))

cat(sprintf("cleaning output dir: %s\n", output.dir))
system(paste("rm -rf", output.dir))
system(paste("mkdir -p", output.dir))

ifiles = list.files(input.dir)
if (length(ifiles) == 0) {
  stop(sprintf("no files in input directory: %s\n", input.dir))
} else {
  cat(sprintf("found %d input files in input directory: %s\n", length(ifiles), input.dir))
}

commands = NULL
for (i in seq_along(ifiles)) {
  oprefix.job = paste(split.dir, "/", i, sep="")
  file = paste(input.dir, ifiles[i], sep="/")
  command = sprintf("%s/lscripts/coords2fends.pl %s %d %d %d %s %s", wd, fends.fn, read.length, facing.threshold, max.seg.len, oprefix.job, file)
  commands = c(commands, command)
}

time = proc.time()
source(paste(wd, "/R/distrib.r", sep=""))
distrib(commands, working.dir=wd, qsub.dir=qsub.dir, jobname=paste("coord2fends", dataset, sep="_"), batch.max.jobs=400, total.max.jobs.fn="/home/eitany/maxjobs", sleep.time=10, req.mem=req.mem)
cat(sprintf("distributed computation took %.1f seconds\n", (proc.time() - time)[3]))

# give time for nfs to update files
max.sleeps = 0
while (1) {
  command = sprintf("grep done %s/*.stdout | wc -l", qsub.dir)
  done = as.numeric(system(command, intern=T))
  if (done == length(ifiles))
    break
  cat("All jobs not done yet, sleeping for 30 secs...\n")
  Sys.sleep(30)
  max.sleeps = max.sleeps + 1
  if (max.sleeps > 10)
    stop("not all distributed tasks finished successfully, stopping")
}

# unite general stats
all.fns = NULL
for (i in seq_along(ifiles)) {
  all.fns = c(all.fns, paste(split.dir, "/", i, "_all.mat.stats", sep=""))
}

unite.stats = function(ifns, ofn) {
  if (length(ifns) == 0)
    stop("no files to unite")
  for (i in seq_along(ifns)) {
    data = read.delim(ifns[i])
    if (i > 1)
      result = result + data
    else
      result = data
  }
  write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

unite.stats(all.fns, paste(output.dir, "/all.stats", sep=""))
            
# unite type specific stats
types = c("s0", "s1", "s2")
for (type in types) {
  read.stats.fns = NULL
  for (i in seq_along(ifiles)) {
    oprefix.job = paste(split.dir, "/", i, "_", type, sep="")
    
    read.stats.fns = c(read.stats.fns, paste(oprefix.job, ".read.stats", sep=""))
  }
  unite.stats(read.stats.fns, paste(output.dir, "/", type, ".read.stats", sep=""))

  if (type == "s0") {
	oprefix = paste(output.dir,"/",type,sep="")
	command = sprintf("%s/lscripts/merge_mat_files2.pl %s %s %s %s", wd, fends.fn, oprefix, split.dir, type)
	cat(sprintf("command: %s\n", command))
	if (system(command)!=0)
	   stop(sprintf("error in command: %s",command))

	command = sprintf("%s/lscripts/merge_mat_probe_files.pl %s %s %s", wd, fends.fn, oprefix, type)
	cat(sprintf("command: %s\n",command))
	system(command)
	}
}
