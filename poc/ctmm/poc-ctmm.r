#!/usr/bin/env Rscript --vanilla

# ==== Input parameters ====

#TODO: segment output files should have a seg_name column
'
Usage:
poc_ctmm.r <dat> <out> [--db=<db>] [--sesid=<sesid>] [--seed=<seed>] [--parMethod=<parMethod>] [--cores=<cores>] [--mpilogs=<mpilogs>] [-t] 
poc_ctmm.r (-h | --help)

Control files:
ctfs/individual.csv

Parameters:
dat: path to the csv file containing segment info. Required fields:
  individual_id,seg_start,seg_end
out: output directory for all ctmm objects

Options:
-h --help     Show this screen.
-v --version     Show version.
-r --sesid=<sesid>  Id that uniquely identifies a script run
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-t --test         Indicates script is a test run, will not save output parameters or commit to git
-p --parMethod=<parMethod>  Either <mpi | mc>. If not passed in, script will run sequentially.
-c --cores=<cores>  The number of cores
-m --mpilogs=<mpilogs> Directory for the mpi log files
' -> doc

if(interactive()) {
  #library(here)

  #.sesid <- 'test4'
  
  .wd <- '~/projects/ms3/analysis/full_workflow_poc/test4'
  .seed <- NULL
  .test <- TRUE

  rd <- here::here
  
  .datPF <- file.path(.wd,'data/seg_dates.csv')
  .dbPF <- '~/projects/ms3/analysis/full_workflow_poc/data/move.db'
  .outP <- file.path(.wd,'ctmm')
  
  .parMethod <- NULL
  # .parMethod <- 'mc'
  # .cores <- 3
} else {
  library(docopt)
  library(rprojroot)

  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .script <-  thisfile()
  .seed <- ag$seed
  .test <- as.logical(ag$test)
  rd <- is_rstudio_project$make_fix_file(.script)
  .sesid <- ag$sesid
  .parMethod <- ag$parMethod
  .cores <- ag$cores

  source(rd('src/funs/input_parse.r'))
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datPF <- makePath(ag$dat)
  .outP <- makePath(ag$out)
  .dbPF <- makePath(ifelse(length(ag$db)!=0, ag$db,'data/move.db'))
  .mpiLogP <- makePath(ifelse(is.null(ag$mpilogs),'mpilogs',ag$mpilogs))

}

# ==== Setup ====

#---- Initialize Environment ----#

.seed <- ifelse(is.null(.seed),5326,as.numeric(.seed)) 

set.seed(.seed)
t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(ctmm)
    library(iterators)
    library(foreach)
    library(DBI)
    library(RSQLite)
  }))

source(rd('src/funs/auto/breezy_funs.r'))

invisible(assert_that(file.exists(.dbPF)))

#---- Local parameters ----#
.timingPF <- file.path(.outP,'timing.csv')

#---- Functions ----#
ctmm.select_q <- quietly(ctmm::ctmm.select)
akde_q <- quietly(ctmm::akde)

#---- Files and directories ----#

dir.create(file.path(.outP,'tel'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(.outP,'vg'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(.outP,'mod'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(.outP,'akde'),recursive=TRUE,showWarnings=FALSE)

#---- Load control files ----#

inds <- read_csv(file.path(.wd,'ctfs/individual.csv'),col_types=cols()) %>%
  filter(as.logical(run)) %>% select(-run) 

#---- Load data ----#

segs <- read_csv(.datPF,col_types=cols()) %>%
  inner_join(inds %>% select(individual_id,local_identifier),by='individual_id')

# ==== Start cluster and register backend ====
if(is.null(.parMethod)) {
  message('No parallel method defined, running sequentially.')
  #foreach package as %do% so it is loaded even if the parallel packages are not
  `%mypar%` <- `%do%`
} else if(.parMethod=='mpi') {
  message('Registering backend doMPI')
  library(doMPI)
  
  dir.create(.mpiLogP,showWarnings=FALSE,recursive=TRUE)
  #start the cluster. number of tasks, etc. are defined by slurm in the init script.
  message('Starting mpi cluster.')
  cl <- startMPIcluster(verbose=TRUE,logdir=.mpiLogP)
  registerDoMPI(cl)
  setRngDoMPI(cl) #set each worker to receive a different stream of random numbers
  
  `%mypar%` <- `%dopar%`
  
} else if(.parMethod=='mc') {
  #.cores <- strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1)) #for testing on hpc
  message(glue('Registering backend doMC with {.cores} cores'))
  library(doMC)
  RNGkind("L'Ecuyer-CMRG")
  
  registerDoMC(.cores)
  
  `%mypar%` <- `%dopar%`
  
} else {
  stop('Invalid parallel method')
}

# ==== Perform analysis ====
foreach(i=icount(nrow(segs)),.combine='rbind') %mypar% { #

  #i <- 1
  tic()
  seg <- segs[i,]
  segName <- glue('{seg$individual_id}_{seg$year}') #Construct seg name. In the future, generalize to have seg_name column in segment file
  message(glue('Starting individual {i}, {seg$local_identifier} ({seg$individual_id}), {seg$year}'))
  
  #---- initialize database ----#
  db <- dbConnect(RSQLite::SQLite(), .dbPF)
  invisible(assert_that(length(dbListTables(db))>0))
  indtb <- tbl(db,'individual')
  
  #---- Load data ----#
  
  sql <- 'select event_id,individual_id,lon,lat,timestamp
  from event
  where individual_id = {seg$individual_id}
    and timestamp >= {seg$seg_start}
    and timestamp < {seg$seg_end + 1}
    and event_id not in (select event_id from outlier where individual_id = {seg$individual_id})
  order by individual_id,timestamp' %>% glue_sql(.con=db)
  
  tic()
  dat <- dbGetQuery(db, sql) %>% as_tibble %>%
    mutate(timestamp=as.POSIXct(timestamp,tz='UTC'))
  toc()
  
  dbDisconnect(db)
  
  message(glue('Segment has {nrow(dat)} records'))
  
  #Create telemetry object. ctmm wants movebank column names.
  #Note using individual_id instead of local.identifier as names
  tel <- dat %>% 
    #sample_n(100) %>% #!!!!!!!!!!!!!!!! REMOVE !!!!!!!!!!!!!!!!
    select(location.long=lon,location.lat=lat,timestamp) %>%  
    mutate(individual.local.identifier=segName) %>%
    arrange(timestamp) %>%
    as.data.frame %>% 
    as.telemetry
  
  message('Estimating variogram')
  tic()
  vg <- variogram(tel)
  toc()
  
  
  message('Estimating model')
  tic()
  GUESS <- variogram.fit(vg, interactive=FALSE)
  
  mod <- ctmm.select_q(tel,
    CTMM = GUESS,
    method = "pHREML",
    control=list(method="pNewton"))
  toc()

  message('Estimating akde')
  tic()
  akde <- akde_q(tel, mod$result, res = 50)
  toc()
  
  segF <- glue('{segName}.rds')
  
  message('Saving ctmm objects')
  saveRDS(tel,file.path(.outP,'tel',segF))
  saveRDS(vg,file.path(.outP,'vg',segF))
  saveRDS(mod$result,file.path(.outP,'mod',segF))
  saveRDS(akde$result,file.path(.outP,'akde',segF))
  
  tot <- toc()
  totm <- round((tot$toc - tot$tic)/60,2)
  
  message(glue('{seg$local_identifier} ({seg$individual_id}), {seg$year} complete in {totm} minutes'))
  
  #In the future, seg_id will be available in the segment file, and will be used as FK here, so no seg_name, ind id req.
  tibble(individual_id=seg$individual_id,
         seg_id=NA,
         seg_name=segName,
         npts=nrow(tel), #!!!!! TEST THIS !!!!
         minutes=totm) %>%
    write_csv(.timingPF,append=file.exists(.timingPF),na="")
  
  return(TRUE)
} -> status

# status %>% 
#   as_tibble(.name_repair = 'minimal') %>%
#   rename(status=1) %>% 
#   write_csv(.statusPF)

# ==== Finalize script ====
if(!.test) {
  suppressWarnings(
    suppressPackageStartupMessages({
      library(git2r)
      library(uuid)
    }))
  
  .runid <- UUIDgenerate()
  .parPF <- file.path(.wd,"run_params.csv")
  
  #Update repo and pull out commit sha
  repo <- repository(rd('src'))
  
  rstat <- status(repo)
  if(length(rstat$staged) + 
     length(rstat$unstaged) + 
     length(rstat$untracked) > 0) {
    add(repo,'.')
    commit(repo, glue('script auto update. runid: {.runid}'))
  }
  
  
  .git_sha <- sha(repository_head(repo))
  
  #Save all parameters to csv for reproducibility
  #TODO: write this to a workflow database instead
  saveParams(.parPF)
}

message(glue('Script complete in {diffmin(t0)} minutes'))

if(!is.null(.parMethod) && .parMethod=='mpi') { #seems nothing after mpi.quit() is executed, so make sure this is the last code
  closeCluster(cl)
  mpi.quit()
}

