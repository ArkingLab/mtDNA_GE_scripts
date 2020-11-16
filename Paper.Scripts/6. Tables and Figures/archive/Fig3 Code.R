

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0044451","nucleoplasm part", 0.895, 4.569, 3.924, 4.945,-47.5255,0.447,0.000),
c("GO:0098805","whole membrane", 0.888,-2.408,-4.504, 4.941,-8.6677,0.934,0.000),
c("GO:1902494","catalytic complex", 3.734,-4.886, 3.212, 5.565,-58.1393,0.762,0.000),
c("GO:0097458","neuron part", 0.320,-4.208,-3.014, 4.498,-7.0540,0.852,0.059),
c("GO:0031975","envelope", 2.324,-0.126,-5.703, 5.359,-15.2879,0.836,0.073),
c("GO:0005929","cilium", 0.221, 6.699, 1.014, 4.338,-8.0937,0.621,0.230),
c("GO:0005694","chromosome", 1.505, 3.395, 5.554, 5.170,-16.8077,0.538,0.290),
c("GO:0005739","mitochondrion", 2.156, 3.491, 2.633, 5.327,-22.1099,0.532,0.400),
c("GO:0005794","Golgi apparatus", 0.969, 3.801, 1.102, 4.979,-14.7382,0.473,0.426),
c("GO:0035770","ribonucleoprotein granule", 0.131, 0.449, 4.702, 4.109,-8.8682,0.514,0.430),
c("GO:1990234","transferase complex", 1.223,-4.697, 4.920, 5.080,-37.0791,0.735,0.439),
c("GO:0098687","chromosomal region", 0.306, 4.664, 4.870, 4.479,-7.4754,0.513,0.467),
c("GO:0044429","mitochondrial part", 1.203, 4.417, 3.065, 5.073,-14.5545,0.464,0.481),
c("GO:0098796","membrane protein complex", 2.473,-4.957, 3.573, 5.386,-7.0498,0.762,0.484),
c("GO:0015630","microtubule cytoskeleton", 0.900, 3.840, 5.946, 4.947,-13.3919,0.507,0.523),
c("GO:1990904","ribonucleoprotein complex", 5.291,-4.877, 4.062, 5.717,-17.3005,0.757,0.543),
c("GO:0032838","cell projection cytoplasm", 0.014, 5.995,-2.409, 3.136,-7.9594,0.678,0.601),
c("GO:0000151","ubiquitin ligase complex", 0.232,-3.386, 5.040, 4.358,-17.1811,0.616,0.604),
c("GO:0005681","spliceosomal complex", 0.250, 1.287, 3.948, 4.392,-13.3503,0.415,0.628),
c("GO:0005768","endosome", 0.319, 4.073, 0.568, 4.497,-9.2019,0.502,0.669),
c("GO:0031248","protein acetyltransferase complex", 0.152,-3.634, 5.450, 4.175,-7.6483,0.625,0.693),
c("GO:0031461","cullin-RING ubiquitin ligase complex", 0.159,-3.711, 4.837, 4.195,-9.0096,0.624,0.695),
c("GO:0044441","ciliary part", 0.139, 6.193, 2.393, 4.137,-7.8004,0.532,0.711),
c("GO:0044431","Golgi apparatus part", 0.608, 4.218, 1.796, 4.777,-12.2954,0.421,0.713),
c("GO:0016607","nuclear speck", 0.091, 5.789, 4.217, 3.953,-18.0175,0.515,0.722),
c("GO:0042175","nuclear outer membrane-endoplasmic reticulum membrane network", 0.771, 3.000,-4.031, 4.880,-9.4714,0.680,0.731),
c("GO:0005813","centrosome", 0.185, 4.677, 5.418, 4.261,-10.0874,0.483,0.747),
c("GO:0005730","nucleolus", 0.664, 4.033, 4.209, 4.815,-7.3250,0.436,0.759),
c("GO:0016604","nuclear body", 0.189, 5.408, 4.071, 4.269,-29.9712,0.494,0.770),
c("GO:0044440","endosomal part", 0.124, 4.634, 1.152, 4.087,-8.8982,0.473,0.785),
c("GO:0005815","microtubule organizing center", 0.350, 4.227, 5.138, 4.537,-13.0272,0.463,0.792),
c("GO:0005740","mitochondrial envelope", 0.901, 4.737, 3.010, 4.948,-10.1093,0.474,0.822),
c("GO:0000139","Golgi membrane", 0.403, 4.414, 1.656, 4.598,-12.3780,0.432,0.865));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
# changed plot size to be smaller!
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = (plot_size/2)), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
# p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 

# 4.20.2020: selected terms to include based on interest (first 3) and based on dispensability (last 5)
terms.include <- c('mitochondrion', 'ubiquitin ligase complex', 'spliceosomal complex', 'nucleoplasm part', 'catalytic complex', 'whole membrane', 'neuron part', 'envelope')
ex <- subset(one.data, description %in% terms.include)

p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

# added to ex for readability:
terms.include2 <- c('nuclear outer membrane-endoplasmic reticulum membrane network', 'cell projection cytoplasm', 'endosome', 'cilium')
ex2 <- subset(one.data, description %in% terms.include2)
p1 <- p1 + geom_text( data = ex2, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );

# remove the legend size and change color title....
p1 + guides(size = F) + labs(col = '-log 10 P-value') + theme_classic()

