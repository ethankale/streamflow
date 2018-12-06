#' Estimating Offsets for Log-Log Stage Discharge Relationships
#'
#' Given paired stage and discharge values, calculate the best stage offset
#' to create a linear relationship in log-log space.
#'
#' USGS Water Supply Paper 2175 gives a process for creating stage-discharge
#' relationships in streams.  Basically, you plot out the measured discharges
#' and their associated stages in log-log space.  Then you do it again,
#' but this time you subtract an offset from all the stage values.  You continue
#' this successive approximation until you get a straight line in log-log space.
#' 
#' This function does that, but faster and better.
#' 
#' At each interval the function creates a log-log relationship, then plots it.
#' Estimates happen at 0.01 intervals from 0 to just below the minimum measured
#' stage value.  
#' 
#' @param stages A numeric vector of stage values, associated with discharges.
#' @param discharges A numeric vector of discharges, associated with stages.
#' @param era.name What time interval?
#' @param stream.name What stream?
#' @param export.dir Directory to put the graphs into.
#' @param plots One of "all" (default), "offsets", "none".
#' @keywords stream stage discharge stage-discharge flow cfs
#' @export
#' @examples 
#' stages <- c(1.2, 1.5, 1.3, 1.0)
#' discharges <- c(22, 55, 31, 21)
#' era.name <- '2008-2012'
#' stream.name <- 'Demo Stream'
#' export.dir <- './results/offsets/'
#' find_best_offset(stages, discharges, era.name, stream.name, export.dir)

find_best_offset <- function(stages, discharges, era.name, stream.name, export.dir,
                             plots = "all") {
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(dplyr))
  
  #################
  # Validation
  #################
  
  if (!dir.exists(export.dir)) { 
    message(paste0(export.dir, " does not exist!"))
    return(NULL); 
  }

  era.name <- gsub("[^[:alnum:][:space:]]", "", era.name)
  
  #############
  # Plots
  #############
  
  message(paste("Calculating offsets for", stream.name, era.name ))
  
  d.new <- tibble(stage = stages, discharge = discharges) %>%
    mutate(d.cfs.log = log(discharge))
  
  # Earlier era; find the best offset/gzf.
  gage.minimum <- min(d.new$stage, na.rm = TRUE)
  gage.offsets <- seq(-0, gage.minimum-0.01, 0.01)
  
  offset.stats <- tibble(offsets = gage.offsets,
                         r2 = rep(0, length(offsets)))
  
  pb <- txtProgressBar(min = 0, max = length(gage.offsets))
  
  i <- 0
  while(i < length(gage.offsets)) {
    i <- i+1
    e <- gage.offsets[i]
    
    setTxtProgressBar(pb, i)
    
    d.new.offset <- d.new %>% 
      mutate(stage.shifted = stage - e,
             stage.shifted.log = log(stage.shifted))
  
    # Build a relationship; make the line & calculate stats
    d.lm <- with(d.new.offset, lm(d.cfs.log ~ stage.shifted.log))
    d.new.offset$d.cfs.pred <- exp(predict(d.lm, newdata = d.new.offset))
    d.r2 <- summary(d.lm)$adj.r.squared
  
    offset.stats$r2[i] <- d.r2
    
    if (plots == "all") {
      ggplot(d.new.offset, aes(x = discharge,
                               y = stage.shifted)) +
        geom_line(aes(x = d.cfs.pred,
                      y = stage.shifted)) +
        geom_point(size = 3) +
        labs(title = paste0(stream.name, " Stage-Discharge"),
             subtitle = paste0(era.name, " Measurements; ",
               "E = ", e, 
               "; R2 = ", round(d.r2, 3)),
             x = "Discharge (cfs)",
             y = "Stage - Offset (feet)") +
        theme_minimal()
      
      ggsave(file.path(export.dir,
                       paste0("stage-discharge-",
                       era.name, "-offset-", 
                       i, ".png")),
             width = 7,
             height = 4)
    }
  }
  
  close(pb)
  
  e.best <- offset.stats$offsets[which.max(offset.stats$r2)]
  
  if(!plots == "none") {
    ggplot(offset.stats %>% filter(r2 >= median(offset.stats$r2)), aes(x = offsets,
                             y = r2)) +
      geom_line() +
      geom_point() +
      geom_point(data = offset.stats %>% filter(offsets==e.best),
                 color="red",
                 size = 3)+
      theme_minimal() +
      labs(title=paste0(stream.name, " Stage Offset Estimate"),
           subtitle=paste0(era.name, "; Best e = ", e.best),
           x = "Offsets (feet)",
           y = "R2 Values (for log-log regression)") +
      scale_x_continuous(breaks = seq(0, max(d.new$stage), 0.1))
    
    ggsave(paste0(export.dir, "offsets_", era.name, ".png"),
           width = 7,
           height = 4)
  }
  
  return(offset.stats[offset.stats$offsets == e.best,])
}
