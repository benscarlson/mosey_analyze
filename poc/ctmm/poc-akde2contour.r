#!/usr/bin/env Rscript --vanilla
# chmod 744 script_template.r #Use to make executable

# This script implements the breezy philosophy: github.com/benscarlson/breezy

#TODO: deal with issue writing Polygon & Multi Polygon to the gpkg
#TODO: make this into a script

# ==== Breezy setup ====

'
Load AKDE objects and saves contours to a geopackage.

Usage:
akde2contour <dat> <segs> <out> [--contour=<contour>] [--seed=<seed>] [-t]
akde2contour (-h | --help)

Control files:
  ctfs/individual.csv

Parameters:
  dat: path to folder containing akde objects
  out: output geopackage.
  segs: csv file containing segment info. Fields: individual_id, year

Options:
-h --help     Show this screen.
-v --version     Show version.
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-t --test         Indicates script is a test run, will not save output parameters or commit to git
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)

  .sesid <- 'test2'
  
  .wd <- file.path('~/projects/ms3/analysis/full_workflow_poc',.sesid)
  .seed <- NULL
  .test <- TRUE
  rd <- here::here
  
  .datP <- file.path(.wd,'ctmm/akde')
  .outPF <- file.path(.wd,'ctmm/contour.gpkg')
  .segPF <- file.path(.wd,'data/seg_dates.csv')
  
  .contour <- 0.95
  
} else {
  library(docopt)
  library(rprojroot)

  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  thisfile()
  .seed <- ag$seed
  .test <- as.logical(ag$test)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('src/funs/input_parse.r'))
  
  .datP <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
}

#---- Initialize Environment ----#
.seed <- ifelse(is.null(.seed),5326,as.numeric(.seed))

set.seed(.seed)
t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(ctmm)
    library(sf)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Local parameters ----#

#---- Load control files ----#
inds <- read_csv(file.path(.wd,'ctfs/individual.csv'),col_types=cols()) %>%
  filter(as.logical(run)) %>% select(-run)

#---- Load data ----#
message('Loading data...')
segs <- read_csv(.segPF,col_types=cols()) %>%
  inner_join(inds %>% select(individual_id,local_identifier),by='individual_id')

#====

#---- Perform analysis ----#

#Read in akde object, convert to sf, calculate area, save to geo package db
#Note, just keep the estimate for the contour. Not using ci at this time.
message(glue('Saving to {.outPF}'))

dir.create(dirname(.outPF),recursive=TRUE,showWarnings=FALSE)

#TODO: error handling if akde rds is not available

result <- segs %>%
  mutate(status=map2_lgl(individual_id,year,~{
    readRDS(file.path(.datP,glue('{.x}_{.y}.rds'))) %>%
      SpatialPolygonsDataFrame.UD(level.UD=.contour,level=0.95) %>%
      st_as_sf %>%
      mutate(individual_id=.x, year=.y, name=as.character(name)) %>%
      separate(name,c('seg_name','per','est'),sep=' ') %>%
      filter(est=='est') %>%
      mutate(area_km2=as.numeric(st_area(.))/1e3^2) %>%
      st_transform(4326) %>%
      select(individual_id,year,seg_name,per,area_km2) %>% 
      st_write(.outPF,layer='contour',append=file.exists(.outPF))
      
    return(TRUE)
  }))

#TODO: deal with this warning. Maybe convert everything to MULTIPOLYGON?
# In CPL_write_ogr(obj, dsn, layer, driver, as.character(dataset_options),  :
#                    GDAL Message 1: A geometry of type POLYGON is inserted into layer contour of geometry type MULTIPOLYGON, which is not normally allowed by the GeoPackage specification, but the driver will however do it. To create a conformant GeoPackage, if using ogr2ogr, the -nlt option can be used to override the layer geometry type. This warning will no longer be emitted for this combination of layer and feature geometry type.

# st_read(.outPF)
# #https://gis.stackexchange.com/questions/341718/how-do-i-read-a-layer-from-a-gpkg-file-whilst-selecting-on-an-attribute
# #ogrinfo contour.gpkg #This will tell you the name of the table/layer
# st_read(.outPF,query='SELECT * FROM contour WHERE individual_id=8863778')

#---- Finalize script ----#

if(!.test) {
  library(git2r)
  library(uuid)
  
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