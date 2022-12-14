options(warn=-1)
get.job.return.codes=function(files)
{
  codes = rep(0, length(files))
  for (i in 1:length(files)) {
    file = paste(files[i], ".code", sep="")
    rep = 300
    while (!file.exists(file)) {
      cat(sprintf("could not find file: %s, sleeping a minute to let the NFS recover\n", file))
      Sys.sleep(60)
      rep = rep - 1
      if (rep == 0)
        stop("Return code file %s not found")
    }
    codes[i] = read.delim(file, header=F)[1,1]
  }
  codes
}

batch.job.count=function(jobids)
{
  if (is.null(jobids))
    return (0)
  str = as.numeric(system("squeue -h -o %A",intern=T))
  str = intersect(jobids, str)
    njobs = length(str)
  njobs
}

# small helper functions
total.job.count=function()
{
    as.numeric(system("squeue -h | wc -l", intern=T))
}

total.max.jobs=function(total.max.jobs.fn, default=300)
{
  if (file.exists(total.max.jobs.fn))
    return (read.delim(total.max.jobs.fn, head=F)[,1])
  else
    return (default)
}

get.hosts=function(min.total = 12)
{  
  command = "qstat -f"
  result = system(command, intern=T)
  result = result[!grepl("---", result) & result != "queuename" & grepl("@", result)]

  # get hosts
  hosts = sapply(strsplit(result, " "), function(x) x[1])

  # get number of tasks on host
  counts = sapply(strsplit(result, '\\s+', perl=T), function(x) x[3])
  used = sapply(strsplit(counts, '/'), function(x) as.numeric(x[2]))
  total = sapply(strsplit(counts, '/'), function(x) as.numeric(x[3]))
  ind = total >= min.total
  host.map = data.frame(host=hosts[ind], key=rep(NA, sum(ind)), count=used[ind])
  host.map
}

assign.keys.to.hosts=function(host.map, keys, tasks.per.host=12)
{
  N = length(keys)

  hosts = rep(NA, N)
  unique.keys = unique(keys)

  if (length(unique.keys) > dim(host.map)[1])
    stop(sprintf("more unique keys (%d) than hosts (%s)", length(unique.keys) > dim(host.map)[1]))

  # assign a key to all hosts
  hosts.per.key = floor(dim(host.map)[1] / length(unique.keys))
  for (i in seq_along(unique.keys))
    host.map$key[1:hosts.per.key + (i-1)*hosts.per.key] = unique.keys[i]
  
  for (i in 1:N) {
    # find hosts which handle this key
    indices = (host.map$key == keys[i]) & !is.na(host.map$key)
    
    if (!any(indices))
      stop("assertion: no host is handling key")

    # can we add to a host
    index = which((host.map$count < tasks.per.host) & (host.map$count > 0) & indices)[1]

    # if non found we add to host with minimal load
    if (is.na(index))
      index = which(indices & host.map$count == min(host.map$count[indices]))[1]

    if (is.na(index))
      stop("assertion: no host found")
    
    # update host count
    host.map$count[index] = host.map$count[index] + 1
    
    # assign job to host
    hosts[i] = host.map$host[index]
  }

  for (key in unique.keys) {
    msg = gsub("all.q@", "", gsub(".mcl4.weizmann.ac.il", "", unique(hosts[keys == key])))
    cat(sprintf("key %s assigned to: %s\n", key, paste(msg, collapse=", ")))
  }
  
  hosts
}

# files: list of bash files to run
# host.keys: if supplied, we send all files with same key to same host
distrib.files = function(
  files,
  host.keys=NULL, tasks.per.host=12, min.host.total=12,
  batch.max.jobs=400, total.max.jobs.fn="/home/eitany/maxjobs", sleep.time=10, retries=3, req.mem=NA)
{
  options(stringsAsFactors=F)

  if (!is.na(retries) && retries < 0)
    stop("not all jobs were completed and number of retries was exhausted")
  
  N = length(files)
  
  if (!is.null(host.keys)) {
    if (length(host.keys) != length(files))
      stop("files and host.keys must be the same length")
    host.map = get.hosts(min.total=min.host.total)
    hosts = assign.keys.to.hosts(host.map, host.keys, tasks.per.host)
  }

  cat(sprintf("currently %d jobs are running on cluster (all users)\n", total.job.count()))
  cat(sprintf("starting to submit %d job wrapper scripts, sleeping %s seconds between each qsub command (%d retries left)\n", N, sleep.time, retries))
  
  script.names = NULL
  jobids = NULL
  jobid.to.index = list()
  for (i in 1:N)
  {
    while (total.job.count() >= total.max.jobs(total.max.jobs.fn) || batch.job.count(jobids) >= batch.max.jobs) {
      cat(sprintf("[%s] running in full capacity (batch_pending=%d, batch_running=%d, total_user_running=%d, total_running=%d)\n",
                  date(), N-i+1, batch.job.count(jobids), total.job.count(),total.job.count()))
      # wait a bit
      Sys.sleep(sleep.time)
    }
    
    # submit job
    # cat(sprintf("submitting job %d\n", i))
    script.file = files[i]

    if(is.na(req.mem)) {
    #command = "sbatch --mem=3000"
    command = "sbatch --mem=3000 --exclude=phantom-node29"
    }

    else {
    #command = paste0("sbatch --mem=",req.mem)
    command = paste0("sbatch --mem=",req.mem," --exclude=phantom-node29")
    }

    # if host key is supplied we send all jobs with same key to same host
    if (!is.null(host.keys))
      command = paste(command, "-q", hosts[i])
    
    command = paste(command, script.file)
    # cat(sprintf("command: %s\n", command))
    jobid = as.numeric(system(command, intern=T))
    
    if (length(jobid) != 1)
      stop("error in qsub command")
    
    # cat(sprintf("submitted jobid: %s\n", jobid))
    jobids = c(jobids, jobid)
    script.names = c(script.names, script.file)
    jobid.to.index[[as.character(jobid)]] = i

    # sleep a bit, in attempt to solve the auto mount bug in the cluster
    Sys.sleep(1)
  }
  
  cat(sprintf("all jobs submitted successfully\n"))
  njobs = batch.job.count(jobids)
  
  prev.njobs = 0
  # wait for jobs to finish running
  while (njobs > 0) {
    # give user a message if new jobs finished
    if (prev.njobs != njobs)
      cat(sprintf("[%s] running in full capacity (batch_running=%d, total_running=%d)\n",
                  date(), batch.job.count(jobids), total.job.count()))
    
    # wait a bit
    Sys.sleep(sleep.time)
    
    # recount running batch jobs 
    prev.njobs = njobs
    njobs = batch.job.count(jobids)
  }
  
  # finally check exit status of all jobs
  Sys.sleep(60)
  
  codes = get.job.return.codes(files)
  failed = which(codes != 0)
  if (length(failed)>0)
  {
    cat(sprintf("Jobs failed:\n"))
    for (i in seq_along(failed)) {
      cat(sprintf("  script: %s, return_code: %s\n", script.names[failed[i]], codes[failed[i]]))
    }

    if (length(failed) == N)
      stop("All jobs failed")
    
    if (!is.na(retries)) {
      cat(sprintf("sending again %d failed jobs (%d retries left)...\n", length(failed), retries))
      distrib.files(files[failed], host.keys=host.keys[failed], tasks.per.host=tasks.per.host, min.host.total=min.host.total,
                    batch.max.jobs=batch.max.jobs, total.max.jobs.fn=total.max.jobs.fn, sleep.time=sleep.time, retries=(retries-1), req.mem=req.mem)
      
    }
    else {
      stop("All jobs failed and out of retries")
    }
  }
  cat("all jobs complete\n")
}

generate.distrib.script.file=function(commands, prefix, jobname, taskname, working.dir, path=NULL, rprofile=NULL)
{  
  script.fn = paste(prefix, ".sh", sep="")
  stdout.log.file = paste(prefix, ".stdout", sep="")
  cl.log.file = paste(prefix, ".command_line", sep="")
  rc.file = paste(script.fn, ".code", sep="")

  system(paste("rm -rf", script.fn, stdout.log.file, cl.log.file))
    
  # add header to script
  write("#!/bin/sh", file=script.fn, append=T)
  write(paste("#$ -N ", taskname, sep=""), file=script.fn, append=T)
  write(paste("#SBATCH -o ", stdout.log.file, sep=""), file=script.fn, append=T)
  write("#$ -S /bin/sh", file=script.fn, append=T)

  write("echo hostname: `hostname`", file=script.fn, append=T)
  
  write(paste("cd", working.dir), file=script.fn, append=T)
  
  assert = "if [ $? -ne 0 ]; then exit 1; fi"
  for (i in 1:length(commands)) {
    command = commands[i]
    
    # add command linee to script
    write(command, file=script.fn, append=T)
    if (i < length(commands))
      write(assert, file=script.fn, append=T)
  
    # save command line
    write(command, file=cl.log.file, append=T)
  }
  write(sprintf("echo $? > %s", rc.file), file=script.fn, append=T)
  
  return (script.fn)
}

# command: a list of commands
# working directory: cd there before running commands
# qsub dir: tmp directory with sh script, stdout, stderr and the command line
distrib = function(
  commands, working.dir=getwd(), qsub.dir, jobname="test", 
  batch.max.jobs=400, total.max.jobs.fn="/home/eitany/maxjobs", sleep.time=10, path=NULL, rprofile=NULL, host.keys=NULL, retries=3, req.mem=NA)
{
  
  # generate scripts
  cat(sprintf("preparing bash wrapper scripts in directory: %s\n", qsub.dir))
  files = NULL
  for (i in 1:length(commands))
  {
    prefix = paste(qsub.dir, "/", jobname, "_", i, sep="")
    command = commands[i]
    taskname = paste(jobname, "_", i, sep="")
    files[i] = generate.distrib.script.file(commands=command, prefix=prefix ,jobname=jobname, taskname=taskname, working.dir=working.dir)
  }

  distrib.files(files=files, batch.max.jobs=batch.max.jobs, total.max.jobs.fn=total.max.jobs.fn, sleep.time=sleep.time, host.keys=host.keys, retries=retries, req.mem=req.mem)
}
