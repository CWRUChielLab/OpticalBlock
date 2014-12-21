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

# choose appropriate x and y columns
x_cols = y$plot_x_variable
if (is.null(x_cols)) {
    x_cols = colnames(d)[1]
}

y_cols = y$plot_y_variable
if (is.null(y_cols)) {
    y_cols = colnames(d)[2:(dim(d)[2])]
}

# use the smallest and largest values in the first column as the x-range
xmin = min(d[,x_cols])
xmax = max(d[,x_cols])

# use the remaining columns for the y-rang
ymin = min(d[,y_cols])
ymax = max(d[,y_cols])

# use points iff there are 50 or fewer values, otherwise use lines
print(" one")
if (length(d[,x_cols]) <= 50) {
    plottype = "p"
} else {
    plottype = "l"
}

# create a new plot
pdf(plotfilename)
plot(d[,x_cols], d[y_cols][,1], xlim=c(xmin,xmax), ylim=c(ymin,ymax), type=plottype,
     main=y$plot_title, xlab=y$plot_x_axis_label, ylab=y$plot_y_axis_label)

if (length(y_cols) > 1) {
    for (i in 2:length(y_cols)) {
        lines(d[,x_cols], d[y_cols][,i], type=plottype)
    }
}

dev.off()
