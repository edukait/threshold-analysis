tsvreader <- function(filename){
  ## This function reads tsv file outputs from Angicart C++.  It was written to allows for easier extraction of the first seven columns of the tsv files, a task that apparently is not possible with the standard R file reading functions.  This functions reads line by line the tsv file and extracts the first seven columns of data corresponding to: (1) Vessel name "nodeid", (2) vessel volume "vol", (3) vessel length "length", (4) vessel radius as backcalculated from volume assuming cylindrical shape "radius_vol", (5) vessel radius as calculated by averaging along lenght of vessel "radius", (6) parent vessel to current vessel "parent", (7) number of child branches, or bifurcation value "n".
  
  # Find number of rows.
  data <- readLines(filename)
  rowcount <- length(data)
  
  # Initialize column vectors
  nodeid <- mat.or.vec(rowcount, 1)
  vol <- mat.or.vec(rowcount, 1)
  length <- mat.or.vec(rowcount, 1)
  radius_vol <- mat.or.vec(rowcount, 1)
  radius <- mat.or.vec(rowcount, 1)
  parent <- mat.or.vec(rowcount, 1)
  n <- mat.or.vec(rowcount, 1)
  
  # Loop through readLines output and extract row information, store in column vectors.
  for(i in 2:rowcount){
    # index starts at 2 because the first row is the header, and we don't need it
    datarow <- strsplit(x = data[i], split = "\t")[[1]]
    nodeid[i] <- datarow[1]
    vol[i] <- datarow[2]
    length[i] <- datarow[3]
    radius_vol[i] <- datarow[4]
    radius[i] <- datarow[5]
    parent[i] <- datarow[6]
    n[i] <- datarow[7]
  }
  
  # Convert characters to numerics for all but nodeid and parent.  For nodeid and parent, use auto convert to set to either numeric or character (whichever is defaulted to), but also fixes automatically the N/A values.
  nodeid <- type.convert(x = nodeid, na.strings = "N/A", as.is = TRUE)
  vol <- as.numeric(vol)
  length <- as.numeric(length)
  radius_vol <- as.numeric(radius_vol)
  radius <- as.numeric(radius)
  parent <- type.convert(x = parent, na.strings = "N/A", as.is = TRUE)
  n <- as.numeric(n)
  
  treeDF <- data.frame("nodeid" = nodeid, "vol" = vol, "length" = length, "radius_vol"= radius_vol, "radius" = radius, "parent" = parent, "n" = n)
  
  # Remove first row of zeroes
  treeDF <- treeDF[-1, ]
  
  return(treeDF)
}

