# ================== Signal Simulation =======================
# n - length of image sequence;
# points_number - a numeric scalar, number of points to generate per signal instance;
# sig_mag - numeric scalar, signal magnitude;
# m - baseline noise magnitude, i.e. low frequency noise;
# dur_time - numeric, animation duration;
# file_name - output file name w/o extension;
# path - a  path to a temp folder. Default is a current working directory;
# spread - width of the simulated peak signal as a percentage of points_number

signalSim <- function(n = 50, sig_mag = 3, m = 0, ext = 1.5
                      , dur_time = 10
                      , points_number = 1000
                      , spread = 10
                      , file_name
                      , seed = 1988
                      , path = getwd()
                      ){
  require(magick)
  require(magrittr)
  
  # CONSTANTS
  MAX_NPC <- 50
  MAX_TIME <- 10
  SPREAD_LIMS <- c(1, 20)
  
  file_name <- paste(file_name, "gif", sep = ".")
  
  # Create a temporary directory to store an image sequence
  if (file.exists(file.path(path, "out"))){
    dir_out <- file.path(path, "out")
  } else {
    dir.create(file.path(path, "out"))
    dir_out <- file.path(path, "out")
  }
  
  # Checking arguments
  if (spread < SPREAD_LIMS[1] | spread > SPREAD_LIMS[2]) {
    warning("\n spread value must be within ", SPREAD_LIMS[1], " and ", SPREAD_LIMS[2]
            , "\n automatically set to ", SPREAD_LIMS[2], "%", immediate. = TRUE)
    spread <- SPREAD_LIMS[2]
  }
  
  if (dur_time > MAX_TIME) {
    dur_time <- MAX_TIME
    warning("\ndur_time exceeds maximum duration limit and was set to ", MAX_TIME, " (sec)\n", immediate. = TRUE)
    }
  
  # Image ID to generate
  image_id <- MAX_NPC %>% 
    seq(from = 1, to = n, length.out = .) %>% 
    round(., digits = 0) %>% 
    unique()
  
  # Compute fps value for the GIF render
  fps <- ceiling(length(image_id)/dur_time)
  fps <- c(fps, fps)
  
  
  # Work Around FPS Limitation, image_animate requires fps as a factor of 100 !!!
  # Searching in both direction if fps is not a factor of 100
  while (!any(100 %% fps == 0)){
    fps <- fps + c(1, -1)
  }
  
  # Getting right fps value
  fps <- fps[100 %% fps == 0]
  fps <- max(fps)
  
  signal <- 0
  noise <- 0
  pic <- NULL
  s2n <- NULL
  
  # Converts spread from % to number of points
  spread <- round(0.01*points_number*spread, digits = 0)
  if (spread < 1) stop("\n Number of points to compute a signal width is too small, 
                       either spread or points_number need to adjust")
  
  # Max value for the Barplot
  ymax <- sig_mag*n*ext
  
  # Point Sequence per Plot
  point_seq <- 1:points_number
  
  # PURE Signal
  pure <- sig_mag*exp(-(0.5*length(point_seq) - point_seq)^2/spread^2)
  
  # Baseline low frequency noise
  lfnoise <- sin((point_seq)/150)
  
  # Image GENERATION
  for (i in 1:length(image_id)){
    
    # Stacking image names
    pic <- c(pic, paste0(image_id[i], ".png"))
    
    # Generate noise components 
    set.seed(seed = seed)
    noise <- sapply(1:image_id[i]
                    , function(x) rnorm(n = points_number, mean = 0, sd = 1))
    
    # Summation noise over sequence
    noise <- apply(noise, 1, sum)
    
    # Noise Amplitude
    h <- (noise + image_id[i]*m*lfnoise) %>% 
      range() %>% 
      diff()
    
    # Signal Summation
    signal <- noise + image_id[i]*pure
    
    # Compute Signal-to-Noise Ratio 
    s2n <- c(s2n, round(2*max(signal)/h, digits = 1))
    
    # Setup image Render
    png(filename = file.path(dir_out, pic[i])
        , width = 1200
        , height = 680
        , res = 160)
    
    # Layout setup
    layout(matrix(c(1,3,2,3), ncol=2, byrow=TRUE), widths = (c(2,1)))
    
    # Setting margins and style
    par(mar=c(1, 3, 3, 3), 
        bg = "black",
        col.axis ="green",
        col.lab = "green", 
        fg = "green")
    
    # MAIN UPPER PLOT 
    plot(signal, type="l", xlab = NA, ylab = 'Intensity', bty="n", col = "cyan")
    title(main = list(paste0("N = ", image_id[i]), col="red", font=3))
    
    # MAIN BOTTOM PLOT
    par(mar=c(3, 3, 1, 3))
    
    plot(x = image_id[1:i], s2n, type = 'l', col ="red", xlab = NA, ylab = "Apex",
         bty = "n",
         xlim = c(0, n),
         ylim = c(0, max(s2n) + 1),
         xaxs = "i"
         )
    lines(x = image_id[1:i]
          , y = s2n[1]*image_id[1:i]**0.5 # Theoretic SNR
          , col = "lightblue"
          , lty = 2)
    points(x = image_id[i], y = s2n[i], col = "red", pch = 19)
    abline(h = s2n[i], col = "magenta", lty = 2)
    abline(v = image_id[i], col = "magenta", lty = 2)
    
    
    # SIDE BARPLOT
    par(mar = c(3, 1, 3, 3))
    
    barplot(c(max(signal), h), col = c("cyan","yellow"),
            border = "violet",
            names.arg = c("Signal","Noise"),
            ylim = c(0, ymax)
            )
    title(main = paste("S/N = ", s2n[i], sep = ""), col.main = "red")
    
    dev.off()
  }
  
  # Setup a GIF Render
  file.path(dir_out, pic) %>% 
    lapply(., image_read) %>% 
    image_join(.) %>% 
    image_animate(., fps = fps) %>% 
    image_write(., path = file.path(path, file_name))
  
  # NOTIFICATION
  cat("\nDuration: ", dur_time,"(sec)\n")
  cat("\nFrames per second: ", fps,"\n")
  cat("\nInitial SNR: ", s2n[1],"\n")
  cat("\nFinal SNR: ", s2n[length(s2n)],"\n")
  cat("\nMaximum SNR: ", max(s2n)," is achived at ",image_id[s2n == max(s2n)], "th out of ", n," iters\n" , sep = "")
  cat("\nTemporary directory created: ", dir_out,"\n")
  cat("\nGenerated GIF is saved to ", file.path(path, file_name))
  
  # Deleting temporary files
  gc()
  rmvd <- file.remove(list.files(dir_out, full.names = TRUE))
}


signalSim(n = 50, sig_mag = 8, file_name = "TEST", dur_time = 10, seed = 1988, spread = 3)


