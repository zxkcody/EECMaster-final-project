library(Momocs)
library(tidyr)

# load data
data <- read.csv('522.csv')
data.s <- data %>%
  select(filename, avg_blue, avg_green, avg_red, masks_length, i, contour_coordinates)

# subset background ones
subset <- data[data$avg_blue > 20 | data$avg_green > 20 | data$avg_red > 20, ] %>%
  select(filename, avg_blue, avg_green, avg_red, masks_length, i, contour_coordinates)
# exclude from data
data.s <- data[!(rownames(data) %in% rownames(subset)), ]


# create a blank list for coordinates
matrix_list <- list()
# run through the contour_coordinates column
for (i in 1:length(data.s$contour_coordinates)){
  # each row is coordinates of a mask
  coordinates <- data.s$contour_coordinates[i]
  # remove square brackets
  coordinates <- gsub("\\[|\\]", "", coordinates)
  # split the string by commas and parentheses
  co_lst <- strsplit(coordinates, ", |\\(|\\)")
  # remove empty elements from the list
  co_lst <- co_lst[[1]][co_lst[[1]] != ""]
  # create a matrix with two columns
  matrix <- matrix(co_lst, ncol = 2, byrow = TRUE)
  # turn the character matrix into numeric one
  class(matrix) <- "numeric"
  # combine all object coordinates into list
  matrix_list <- c(matrix_list, list(matrix))
}

# plot all shapes to overview
coo_listpanel(matrix_list) # shapes are divergent

########################
# split matrix_list into subset for shorter processing time
matrix_list1 <- list()
matrix_list2 <- list()
matrix_list3 <- list()
matrix_list4 <- list()
for (x in seq_along(matrix_list)){
  matx <- matrix_list[[x]]
  nrow <- nrow(matx)
  if (!is.na(nrow) && nrow > 0 && nrow %in% 5:50){
    matrix_list1 <- c(matrix_list1, list(matx))
  }
  if (!is.na(nrow) && nrow > 50 && nrow <= 156) {
    matrix_list2 <- c(matrix_list2, list(matx))
  }
  if (!is.na(nrow) && nrow > 156 && nrow <= 339) {
    matrix_list3 <- c(matrix_list3, list(matx))
  }
  if (!is.na(nrow) && nrow > 339) {
    matrix_list4 <- c(matrix_list4, list(matx))
  }
}

# plot to overview
coo_listpanel(matrix_list1)

# convert to Out format
matrix_list1 <- Out(matrix_list1)

# obtain harmonic number
#calibrate_harmonicpower_efourier(matrix_list1)

matrix_list1.f <- matrix_list1[1:7]
#stack(matrix_list1.f, title = 'without fgProcrustes')
# interpolate & fgProcrustes
matrix_list1.f <- Out(matrix_list1[1:7])
matrix_list1.f <- matrix_list1.f %>%
  coo_interpolate(50) %>%
  #coo_sample(50) %>%
  Ldk() %>%
  fgProcrustes() %>% 
  coo_untiltx() %>%
  coo_slidedirection("right")
stack(matrix_list1.f, title = 'with fgProcrustes')

# EF analysis
ef1 <- efourier(matrix_list1, norm = FALSE)
pca <- PCA(ef1, center = TRUE)
plot_PCA(pca, labelpoints = TRUE, zoom = 1)
plot_PCA(pca, center_origin=T, zoom=1,
         labelpoints = F, chull=F, palette=pal_seq_magma, legend = F)%>%
  layer_points(transp=0.1, cex=0.01)%>%
  layer_legend(probs=seq(0,1,0.5), cex=0.5)


matrix_list1.f <- coo_close(matrix_list1.f)
ef1.f <- efourier(matrix_list1.f, norm = FALSE)
pca <- PCA(ef1.f, center = TRUE)
plot_PCA(pca, center_origin=T, zoom=1,
         labelpoints = F, chull=F, palette=pal_seq_magma, legend = F) %>%
  layer_points(transp=0.1, cex=0.01) %>%
  layer_legend(probs=seq(0,1,0.5), cex=0.5)
########################  

# EFA on the whole dataset
# convert into coo format for closed outlines
coo_list <- Out(matrix_list) %>%
  coo_interpolate(300) %>%
  Ldk() %>%
  # full generalized procrustes alignment
  fgProcrustes()
  
stack(coo_list) # shapes are divergent

# calculate EF coefficients
coo_list2 <- Out(coo_list)
ef1 <- efourier(coo_list2, norm = FALSE)
pca <- PCA(ef1, center = TRUE)
plot_PCA(pca, center_origin=T, zoom=1,
         labelpoints = F, chull=F, palette=pal_seq_magma, legend = F) %>%
  layer_points(transp=0.1, cex=0.01) %>%
  layer_legend(probs=seq(0,1,0.5), cex=0.5)

# hierarchical clustering
dist_matrix <- dist(coo_list)  # Compute the distance matrix
hclust_result <- hclust(dist_matrix)  # Perform hierarchical clustering




