
library(igraph)
vertices$names <- gsub("_", " ", vertices$names)
net <- graph_from_data_frame(d=edgetab, vertices=vertices$names, directed=T) 


plot(net, edge.arrow.size=.2,vertex.label=vertices$names, 
     vertex.color = vertices$colors,
     edge.width = edgetab$width,
     edge.color=edgetab$edge.color, 
     vertex.size=5,
     #edge.curved=1,
     layout=as.matrix(vertices[, c("xpos1", "ypos1")])) # layout_randomly
# Otro layout por si acaso
L<- layout.star(net) #layout.davidson.harel fruchterman.reingold gem star
vertices$x2 <- L[, 1]
vertices$y2 <- L[, 2]

library(plotly)
edgelist <- list()
for(i in 1:nrow(edgetab)){
  edge_shape <- list(
    type="line",
    showarrow=FALSE,
    arrowhead=1,
    arrowsize=2,
    line = list(color = edgetab$edge.color[i], width=edgetab$width[i]*edge_width_factor),
    x0 = vertices$xpos1[vertices$names==edgetab$from[i]],
    y0 =vertices$ypos1[vertices$names==edgetab$from[i]],
    x1 = vertices$xpos1[vertices$names==edgetab$to[i]],
    y1 =vertices$ypos1[vertices$names==edgetab$to[i]]
  )
  edgelist[[i]]<- edge_shape
}

axis <- list(title="", showgrid=FALSE, showticklabels=FALSE, zeroline=FALSE)
p <- plot_ly(type="scatter", 
             x=vertices$xpos1, 
             y=vertices$ypos1, 
             mode="markers+text", 
             textfont = list(color=C_NS, size=12),
             marker= list(color = vertices$colors, size=20),
             text=vertices$names, 
             textposition = "center left",
             hoverinfo="text") #%>% 
#add_annotations(edgelist, axref="x", ayref="y")
fig <- layout(p,
              shapes = edgelist, 
              xaxis=axis, yaxis=axis,
              title = "Mediation Analysis - IMC")
fig