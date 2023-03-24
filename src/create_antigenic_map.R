# to install Racmacs:
# install.packages("devtools")
# devtools::install_github("acorg/Racmacs")

# Load the Racmacs package
library(Racmacs)

# Read in the titer table
path_to_titer_file = "/home/jack/Documents/School/Spring 2023/CS523/cs523_project2/results/edge_titer_table.csv"
# path_to_titer_file <- system.file("extdata/h3map2004_hitable.csv", package = "Racmacs")
titer_table        <- read.titerTable(path_to_titer_file)
map <- acmap(
  titer_table = titer_table
)

set.seed(42) # for reproducibility
map <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 2,
  number_of_optimizations = 500,
  minimum_column_basis    = "10"
)

# don't show labels for sites which clustered close together
old_names = agNames(map)
# old_names = c("514" "517" "390" "518" "519" "522" "525" "528" "529" "530" "415" "427"
#               "441" "447" "449" "459" "460" "339" "469" "475" "349" "478" "479" "485"
#               "487" "366" "373" "376")
new_names = c("514", "", "390", "", "", "", "", "", "", "", "", "", 
							"441", "447", "449", "", "460", "339", "", "475", "349", "", "", "485",
							"487", "366", "", "376")
map = edit_agNames(map, old_names, new_names)
map = edit_srNames(map, old_names, new_names)

pdf("antigenic_map.pdf")
plot(map, plot_labels=TRUE, label.offset=0.3, plot_sr=FALSE)
