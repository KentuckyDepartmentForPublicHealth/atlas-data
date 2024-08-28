setwd("/Users/amistry/Desktop")
um<-read.csv("tsne_74_5000.csv",header = T)

Sys.setlocale(locale="C")

library(plotly)
library(ggplot2)
library(ggpubr)

ggplotly(ggplot(um[um$classifier!="#N/A",], aes(V1,V2,color=classifier, label=subgroup, label2=filename, label3=loc, label4=dataset))+
           geom_point(size=1)+theme_classic())

ggplot(um[um$classifier!="#N/A",], aes(V1,V2,color=factor(classifier)))+geom_point(size=0.75)+theme_pubr()+
  theme(legend.position = "none")+rremove("axis")+rremove("axis.text")+rremove("ticks")+rremove("xylab")

summary(reorder(factor(um$classifier),factor(um$classifier), FUN=length))

library(shiny)

# UI  ----
ui <- fluidPage(plotlyOutput("plot"),
                tableOutput("click"))

# server  ----
server <- function(input, output) {
  output$plot <- renderPlotly({
    p <- um[um$hist2=="AT",] %>% 
      ggplot(aes(V1,V2,color=subgroup, label1=dataset, label2=loc, customdata = filename))+geom_point(size = 0.25)+theme_classic()
    fig <- ggplotly(p)
    event_register(fig, "plotly_selected")
    fig
  })
  
  output$click <- renderTable({
    plotly_event_data <- event_data(event = "plotly_selected", priority = "event")
    req(plotly_event_data)
    filter(um, filename %in% plotly_event_data$customdata)
  })
}

# app ----
shinyApp(ui = ui, server = server)

library(dbscan)
res <- optics(um[,c(2:3)], eps = 10, minPts = 3) #10
res <- extractXi(res, xi = 0.03) #0.025
hullplot(um[,c(2:3)],res)
res$cluster

#not accurate
library(factoextra)
r<-dbscan(um[,c(2:3)], 0.4, 4)
r$cluster<-res$cluster
fviz_cluster(r, um[,c(2:3)], geom = "point")+theme_bw()+geom_point(aes(colour=cluster))+
  theme(legend.position = "none")

library(dbscan)
kNNdistplot(um[,c(2:3)], k = 2)
abline(h=1.75, col = "red", lty=2)

res <- dbscan(um[,c(2:3)], eps = 1.75, minPts = 3)
hullplot(um[,c(2:3)], res)
