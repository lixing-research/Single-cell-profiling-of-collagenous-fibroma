left_axes = function(scRNA){

pc12 <- Embeddings(object = scRNA,reduction = 'umap') %>%
  data.frame()

head(pc12,3)

lower <- floor(min(min(pc12$umap_1),min(pc12$umap_2))) - 2

linelen <- abs(0.3*lower) + lower

mid <- abs(0.3*lower)/2 + lower

axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('umap_2','umap_1'),each = 2))

label <- data.frame(lab = c('umap_2','umap_1'),angle = c(90,0),
                    x = c(lower - 3,mid),y = c(mid,lower - 2.5))

result <- list(axes = axes, label = label)

return(result)
}
