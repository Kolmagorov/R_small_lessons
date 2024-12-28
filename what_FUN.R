

import_cdf <- function(file_path
                       , rngRT = c(10,90)
                       , rng_relative = TRUE
                       , extension = "cdf"){
  
  if(any(rngRT < 0)) {stop("\nRT limits cannot be negative")}
  
  if(rng_relative){
    
    if(any(rngRT > 100)) {stop("\nRelative RT limits must be within (0-100) %")}
  }
  
  
  
  ext <- paste("\\.", extension, sep = "")
  file_path <- file_path[grepl(x = file_path, pattern = ext)]
  
  cat("\nNumber of Files provided: ", length(file_path),"\n")
  
  require(RNetCDF)
  
  
  META_VARS <- c("sample_name",
                 "sample_type",
                 "sample_injection_volume",
                 "detector_name",
                 "injection_date_time_stamp")
  
  meta <- NULL
  dataSet <- NULL
  
  
  
  
  for (i in file_path){
   
    nc <- open.nc(con = i, write = FALSE)
    
    signal_value <- var.get.nc(nc, variable = "ordinate_values")
    
    points_number <- length(signal_value)
    samling_rate <- 1/var.get.nc(nc, variable = "actual_sampling_interval")
    run_time <- var.get.nc(nc, variable = "actual_run_time_length")
    rt_scale <- (0:(points_number-1))/samling_rate
    
    if (!is.null(downsamle)){
      
      if(downsamle < 100){
        down_rate <- 0.01*last_point*downsamle
      } else {
        down_rate <- downsamle
      }
      
      new_wav <- seq(0, max(run_time*samling_rate),length.out = down_rate)
      
      as<-prospectr::resample(X = signal_value
                              , wav = rt_scale
                              , new.wav = new_wav)
    }
    
    
    
    
    
    
    meta_row<-sapply(META_VARS, 
                function(x) att.get.nc(ncfile = nc
                                       , variable = "NC_GLOBAL"
                                       , attribute = x))
    
    
    close.nc(nc)
    
    
    if(rng_relative){
      
      idx_lims <- floor(0.01*rngRT*points_number)
      
    } else {
      
      idx_lims <- floor(rngRT*60*samling_rate)
    }
    
    
    rn <- idx_lims[1]:idx_lims[2]
    local_max <- max(signal_value[rn])
    apex_time <- (0:points_number)[signal_value == local_max]/samling_rate/60
    score <- ((sum(signal_value[rn]**2)**0.5)/points_number)*1e6


    
    meta_row<-c(meta_row
                , points_number
                , round(samling_rate, digits = 0)
                , run_time/60
                , round(apex_time, digits = 2)
                , round(local_max, digits = 4)
                , round(score, digits = 2)
                , i)
    
    names(meta_row)<-c(META_VARS
                  , "points_number"
                  , "samling_rate"
                  , "run_time"
                  , "apex_time"
                  , "local_max"
                  , "score"
                  , "file")
    
    
    meta<-rbind(meta, meta_row)
    
    tmp <- data.frame(AU = signal_value
                      , norm_AU = signal_value/local_max
                      , retention_time = rt_scale
                      , injection_date_time_stamp = meta_row[["injection_date_time_stamp"]])
    
    dataSet <- rbind(dataSet, tmp)
  }
  
 
  meta <- as.data.frame(meta)
  row.names(meta) <- NULL
  
  dataSet <- merge(x = dataSet, y = meta[c("sample_name", "injection_date_time_stamp")])
  
  meta$injection_date_time_stamp <- strptime(meta$injection_date_time_stamp, format = "%Y%m%d%H%M%S")
  
  list(meta = meta, dataSet = dataSet)
}







 



