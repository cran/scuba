# Mvalues.R
# data from Jolie Bookspan's book

.ATA <-
  list(
       Haldane = rep(2 * 0.79, 5)
       )

.FSW <-
  list(
       USN  = c(104, 88, 72, 56, 54, 52),
       DSAT = c(99.08, 82.68, 66.89, 59.74, 55.73, 51.44, 49.21, 46.93)
       )

Mvalues <- 
  list(ata=append(.ATA,
         lapply(.FSW, function(x) { x/32.646 })),
       Mvalues.fsw=append(.FSW,
              lapply(.ATA, function(x) { x * 32.646 }))
       )

rm(.ATA, .FSW)

