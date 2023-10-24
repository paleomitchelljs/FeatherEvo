###Function: add spindle aspects to a plot
spindlePlot <- function(y, at, widthFactor=2, absWidth=100, col='black', type="both", breaks=NULL, breakPts=NULL, counts=NULL, Border='black')	{
	if (!is.null(breaks) && !is.null(counts))	{
		Values <- list()
		Values$breaks <- breaks
		Values$counts <- counts
	}
	else		{
		if (is.null(breakPts))	{
			Values <- hist(y, plot=FALSE)
		}
		else		{
			Values <- hist(y, breaks=breakPts, plot=FALSE)
		}
	}
	if (is.null(absWidth))	{
		absWidth <- max(Values$counts)
	}
		
	#Compute Y values
	Yvec <- as.numeric(sapply(Values$breaks, function(x) rep(x, 2)))
	Yvec <- Yvec[-1]
	Yvec <- Yvec[-length(Yvec)]
	
	Total <- (sum(Values$counts) / absWidth) * (widthFactor / 2)

	if (type == "both")	{
		Widths <- (Values$counts / absWidth) * (widthFactor / 2)

		#Compute X values
		Right <- Widths + at
		Left <- at - Widths
	
		xR <- as.numeric(sapply(Right, function(x) rep(x,2)))
		xL <- as.numeric(sapply(Left, function(x) rep(x,2)))
		polygon(c(at, xR, at, rep(at, length(xR))), c(Yvec[1], Yvec, Yvec[length(Yvec)], rev(Yvec)), col=col, border=Border)
		polygon(c(at, xL, at, rep(at, length(xL))), c(Yvec[1], Yvec, Yvec[length(Yvec)], rev(Yvec)), col=col, border=Border)
	}
	if (type == "left")	{
		Widths <- (Values$counts / absWidth) * (widthFactor)

		#Compute X values
		Left <- at - Widths
	
		xL <- as.numeric(sapply(Left, function(x) rep(x,2)))
		polygon(c(at, xL, at, rep(at, length(xL))), c(Yvec[1], Yvec, Yvec[length(Yvec)], rev(Yvec)), col=col, border=Border)
	}
	if (type == "right")		{
		Widths <- (Values$counts / absWidth) * (widthFactor)

		#Compute X values
		Right <- Widths + at
			
		xR <- as.numeric(sapply(Right, function(x) rep(x,2)))
		polygon(c(at, xR, at, rep(at, length(xR))), c(Yvec[1], Yvec, Yvec[length(Yvec)], rev(Yvec)), col=col, border=Border)		
	}
}