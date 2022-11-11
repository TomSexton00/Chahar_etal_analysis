#!/shared/software/miniconda3/envs/r_3.3.2/bin/Rscript

options(warn=-1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <file with batch commands> <odir> <qsub tmp dir> <jobname> <max batch jobs> <total max jobs fn> <required memory (Mb)> <working dir>\n", script.name))
  q(status=1) 
}

commands.fn = args[1]
odir = args[2]
qsub.dir = args[3]
jobname = args[4]
batch.max.jobs = as.numeric(args[5])
total.max.jobs.fn = args[6]
req.mem = args[7]
wd = args[8]

commands = read.delim(commands.fn, header=F)
commands = as.vector(commands[,1])

source(paste(wd, "/R/distrib.r", sep=""))
cat(sprintf("distrib(working.dir=%s, qsub.dir=%s, jobname=%s, batch.max.jobs=%d, total.max.jobs.fn=%s, req.mem=%s)\n",
    odir, qsub.dir, jobname, batch.max.jobs,total.max.jobs.fn, req.mem))
distrib(commands=commands, working.dir=odir, qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs, total.max.jobs.fn=total.max.jobs.fn, req.mem=req.mem)
