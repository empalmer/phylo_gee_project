library(tidyverse)
vec <- 1:10
id <- rep(1:5, each = 2)

test <- data.frame(id,vec) %>% 
  group_by(id) %>% 
  group_map(~.x$vec %*% t(.x$vec))
data.frame(id,vec) %>% 
  group_by(id) %>% 
  group_map(~tcrossprod(.x$vec)) %>% 
  reduce(cbind)

mtcars %>%
  group_by(cyl) %>%
  group_map(~ head(.x, 2L))

test[1,3][[1]]

matrix(vec, ncol = 5) %>% 
  as.data.frame() %>% 
  group_by %>% 
  
  summarize(mean(vec))
  
