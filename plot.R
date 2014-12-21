#!/usr/bin/Rscript

# parse the command line arguments
args = commandArgs(T)
datafilename = args[1]
plotfilename = args[2]

# read in the csv file
d = read.csv(datafilename)

# use the smallest and largest values in the first column as the x-range
xmin = min(d[,1])
xmax = max(d[,1])

# use the remaining columns for the y-range, but include 0
ymin = min(d[,2:(dim(d)[2])], 0)
ymax = max(d[,2:(dim(d)[2])], 0)

# use points iff there are 50 or fewer values, otherwise use lines
if (length(d[,1]) <= 50) {
    plottype = "p"
} else {
    plottype = "l"
}

# create a new plot
pdf(plotfilename)
plot(d[,1], d[,2], xlim=c(xmin,xmax), ylim=c(ymin,ymax), type=plottype)

for (i in 3:(dim(d)[2])) {
    lines(d[,1], d[,i], type=plottype)
}
dev.off()
