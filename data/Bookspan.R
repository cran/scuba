# data from Jolie Bookspan's book

# halftimes from Bookspan p 16
Bookspan <-
  list(Halftimes = list(
         Haldane   = c(5, 10, 20, 40, 75),
         USN	   = c(5, 10, 20, 40, 80, 120),
         DSAT	   = c(5, 10, 20, 30, 40, 60, 80, 120),
         OrcaEdge  = c(5, 11, 17, 24, 37, 61, 87, 125, 197, 271, 392, 480),
         MicroBrain= c(4, 11, 31, 86, 238, 396),
         Aladin	   = c(4, 12, 26, 54, 108, 304),
         ZHL12     = c(4, 7.94, 12.2, 18.5, 26.5, 37, 53, 79, 114,
                       146, 185, 238, 304, 503, 635)))

# surfacing M-values from Bookspan p 23

Bookspan$Mvalues.ata <- list(Haldane = rep(2 * 0.79, 5))

Bookspan$Mvalues.fsw <-
  list(USN  = c(104, 88, 72, 56, 54, 52),
       DSAT = c(99.08, 82.68, 66.89, 59.74, 55.73, 51.44, 49.21, 46.93))

Bookspan$Mvalues.ata$USN  <- Bookspan$Mvalues.fsw$USN/32.646
Bookspan$Mvalues.ata$DSAT <- Bookspan$Mvalues.fsw$DSAT/32.646

Bookspan$Mvalues.fsw$Haldane <- Bookspan$Mvalues.ata$Haldane* 32.646

