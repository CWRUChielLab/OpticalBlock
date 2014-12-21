#!/usr/bin/Rscript

library(yaml)

strip_yaml_extension = function (filename) {
    return(strsplit(filename, "\\.yaml")[[1]][1]);
}

# parse the command line arguments
args = commandArgs(T)
if (length(args) < 1 || length(args) > 3) {
    print("usage: ", quote=F)
    print("    plot.R <yaml file> [<data file> [<plot file>]]", quote=F)
    quit(status=1)
}
yamlfilename = args[1]

if (length(args) >= 2) {
    datafilename = args[2]
} else {
    datafilename = paste0(strip_yaml_extension(yamlfilename), ".csv")
}

if (length(args) >= 3) {
    plotfilename = args[3]
} else {
    plotfilename = paste0(strip_yaml_extension(yamlfilename), ".pdf")
}

# read in the yaml file
y = yaml.load_file(yamlfilename)

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
plot(d[,1], d[,2], xlim=c(xmin,xmax), ylim=c(ymin,ymax), type=plottype,
     main=y$plot_title, xlab=y$plot_x_axis_label, ylab=y$plot_y_axis_label)

for (i in 3:(dim(d)[2])) {
    lines(d[,1], d[,i], type=plottype)
}
dev.off()
