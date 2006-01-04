Mvalues.fsw <- list(
	USN = c(104, 88, 72, 56, 54, 52),
	DSAT = c(99.08, 82.68, 66.89, 59.74, 55.73, 51.44, 49.21, 46.93)
)
Mvalues.ata <- lapply(Mvalues.fsw, function(x) { x/32.646})

Mvalues <- list(ata=Mvalues.ata, fsw=Mvalues.fsw)
