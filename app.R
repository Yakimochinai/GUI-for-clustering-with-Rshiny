  library(shiny)
  library(fclust)
  library(plotrix)
  
  reactiveConsole(TRUE)
  make.crisp <- function(u){
    return(apply(u,1,which.max))
  }
  
  
  compute.points.circular <- function(u,central.cluster=1,m=2){
    
    cluster.angles <- seq(0,(ncol(u)-2),length.out=ncol(u)-1)
    
    cluster2angle <- function(clu.ind){
      if (clu.ind<=central.cluster){
        return((clu.ind-1)*2*pi/(ncol(u)-1))
      }else{
        return((clu.ind-2)*2*pi/(ncol(u)-1))
      }
    }
    
    
    
    cluster.centres <- matrix(0,nrow=ncol(u),ncol=2)
    for (i in 1:nrow(cluster.centres)){
      if (i!=central.cluster){
        clu.angle <- cluster2angle(i)
        cluster.centres[i,1] <- cos(clu.angle)
        cluster.centres[i,2] <- sin(clu.angle)
      }
    }
    
    data.points <- matrix(1,nrow=nrow(u),ncol=2)
    for (i in 1:nrow(u)){
      ordi <- order(u[i,],decreasing=T)
      second.clu <- ordi[1]
      if (second.clu==central.cluster){
        second.clu <- ordi[2]
      }
      
      consti <- 0.5
      if (u[i,second.clu]==1){
        consti <- 1
      }else{
        if (u[i,central.cluster]!=0){
          consti <- (u[i,second.clu]/u[i,central.cluster])^((m-1)/2)
        }
      }
      r <- consti/(1+consti)
      data.points[i,1] <- r*cos(cluster2angle(second.clu))
      data.points[i,2] <- r*sin(cluster2angle(second.clu))
    }
    
    
    return(list(cluster.centres=cluster.centres,data.points=data.points))
  }
  
  ui <- fluidPage(sidebarLayout(
    sidebarPanel( numericInput("nClusters", 'total No. of Clusters', 5),
                  numericInput("nowClusters", 'Cluster number now', 5),
                  h1(textOutput("text")),
                  tags$head(tags$style("#text {color: red;
                                    font-size: 20px;
                                    font-style: italic;
                                   }"
                  )
                  )
                  
    ),
    mainPanel(navbarPage(
      title = "Plots",
      tabPanel(title = "the clusted data",
               plotOutput("plot1",click = "plot_click",width = "500px",height = "500px"),
               verbatimTextOutput("info1")
      ),
      tabPanel(title = " clusterPage ",
               plotOutput("plot2", click = "plot_click",width = "500px",height = "500px"),
               verbatimTextOutput("info2"),
               actionButton("mark","mark"),
               actionButton("Cancel","Cancel")
      ),
      tabPanel(title = "angle",
               plotOutput("plot3",click = "plot_click",width = "500px",height = "500px"),
               verbatimTextOutput("info3"),
               actionButton("mark1","mark"),
               actionButton("Cancel1","Cancel")
      )
      
    )
    
    )
  ))
  
  
  
  server <- function(input,output){
    
    
    n.per.cluster <- 100
    mu1 <- c(0,0)
    sigma1 <- 1
    mu2 <- c(0,7)
    sigma2 <- 1
    mu3 <- c(6.5,3.5)
    sigma3 <- 1
    mu4 <- c(4.5,6)
    sigma4 <- 1
    x <- c(rnorm(n.per.cluster,mean=mu1[1],sd=sigma1),rnorm(n.per.cluster,mean=mu2[1],sd=sigma2),rnorm(n.per.cluster,mean=mu3[1],sd=sigma3),rnorm(n.per.cluster,mean=mu4[1],sd=sigma4))
    y <- c(rnorm(n.per.cluster,mean=mu1[2],sd=sigma1),rnorm(n.per.cluster,mean=mu2[2],sd=sigma2),rnorm(n.per.cluster,mean=mu3[2],sd=sigma3),rnorm(n.per.cluster,mean=mu4[2],sd=sigma4))
    print("Init OK")
    
    getClures <- function(nCl) {
      clures <- FKM(data.frame(x,y),nCl,2)
      return(clures)
    }
    
    getClures1 <- function(nCl) {
      clures1<- FKM(iris[,-5],nCl,2)
      return(clures1)
    }
    
    getClus <- function(clures1) {
      clus <- make.crisp(clures1$U)
      return(clus)
    }
    
    getCluvi <- function(clures, nowCl) {
      cluvi = compute.points.circular(clures$U,central.cluster=nowCl)
      return(cluvi)
    }
    
    getXycoord <- function(clures1, nowCl) {
      xycoord = angle.mapping(t(clures1$H[nowCl,]),iris[,-5])
      return(xycoord)
    }
    
    outputRenderPlot1 <- function(clures) {
      #print("OutputRP1 OK1")
      output$plot1 <- renderPlot({
        plot(x,y,col=make.crisp(clures$U),xlim=c(-3,11),ylim=c(-3,11),pch=16)
      })
      #print("OutputRP1 OK2")
      output$info1 <- renderText({
        paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
      })
      #print("OutputRP1 OK3")
    }
    
    outputRenderPlot2 <- function(clures,cluvi,axis1) {
      #print("OutputRP2 OK1")
      output$plot2 <- renderPlot({
        #print("OutputRP2 OK2")
        plot(cluvi$data.points,xlim=c(-1,1),ylim=c(-1,1),xlab="",ylab="",col=make.crisp(clures$U))
        points(cluvi$cluster.centres,pch=15,col=1:ncol(clures$U))
        draw.circle(0,0,1)
        paint1(axis1,cluvi)
      })
      #print("OutputRP2 OK3")
      output$info2 <- renderText({
        paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
      })
      #print("OutputRP2 OK4")
    }
    
    
    
    
    outputRenderPlot3 <- function(cluvi,axis2,xycoord) {
     output$plot3 <- renderPlot({
       plot(xycoord,col=clus,xlab="",ylab="")
       paint2(axis2,xycoord)
       
     })
     output$info3 <- renderText({
       paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
     })
     
    }
    
    nClusters = 5
    nowClusters = 5
    
    gClures <<- NULL
    gClus <<- NULL
    gCluvi <<- NULL
    gClures1 <<- NULL
    gReady <<- FALSE
    
    axis1 <- reactiveValues(kValue= list(),selectedRoundTwo= list(),distance= list())
    axis2 <- reactiveValues(kValue= list(),selectedRoundTwo= list(),distance= list())
    gAxis1 <<- axis1
    gAxis2 <<- axis2
    
    observeEvent(input$nClusters,{
      nClusters = input$nClusters
      nowClusters = input$nowClusters
      if(nClusters < nowClusters || nClusters <= 1 || nowClusters < 1) {
        output$text <- renderText({
          print("Invalid")
        })
      } else {
        output$text <- renderText({
        })
        
        clures = getClures(nClusters)
        gClures <<- clures
        
        clures1 = getClures1(nClusters)
        gClures1 <<- clures1
        
        clus = getClus(gClures1)
        gClus <<- clus
        
        cluvi = getCluvi(gClures, nowClusters)
        gCluvi <<- cluvi
        
        xycoord = getXycoord(gClures1, nowClusters)
        gXycoord <<- xycoord
  
        outputRenderPlot1(gClures)
        outputRenderPlot2(gClures,gCluvi,gAxis1)
        outputRenderPlot3(gCluvi,gAxis1,gXycoord)
        gReady <<- TRUE
      }
    })
    
    observeEvent(input$mark1,{mark1(gAxis2,gXycoord)})
    observeEvent(input$mark,{mark(gAxis1,gCluvi)})
    observeEvent(input$Cancel,{delete(gAxis1)})
    observeEvent(input$Cancel1,{delete(gAxis2)})
    
    observeEvent(input$nowClusters,{
      nClusters = input$nClusters
      nowClusters = input$nowClusters
      if(nClusters < nowClusters || nClusters <= 1 || nowClusters < 1) {
        output$text <- renderText({
          print("Invalid")
        })
      } else {
        output$text <- renderText({
          print("")
        })
        
        if(gReady && !is.null(gClures)) {
          gCluvi <<- getCluvi(gClures, nowClusters)
          
          clures = getClures(nClusters)
          gClures <<- clures
          
          clures1 = getClures1(nClusters)
          gClures1 <<- clures1
          
          clus = getClus(gClures1)
          gClus <<- clus
          
          cluvi = getCluvi(gClures, nowClusters)
          gCluvi <<- cluvi
          
          xycoord = getXycoord(gClures1, nowClusters)
          gXycoord <<- xycoord
          
          outputRenderPlot1(gClures)
          outputRenderPlot2(gClures,gCluvi,gAxis1)
          outputRenderPlot3(gCluvi,gAxis2,gXycoord)
        }
      }
    })
    
    obs <- observe({
      #print("OBS OK")
      print(paste("observer: ", length(gAxis1$kValue)))
      #print("TEST!!!!!!!!!!")
    })    
    
    paint1 <- function(axis1,cluvi){
      #print("Paint1 OK")
      if(length(axis1$kValue) > 0){
        for(i in 1:length(axis1$kValue)){
          index = axis1$kValue[i]
          index = as.integer(index)
          points(cluvi$data.points[index,1],cluvi$data.points[index,2],pch=15,col="blue4")
        }
      }
    }
    
    paint2 <- function(axis1,xycoord){
      #print("Paint1 OK")
      if(length(axis1$kValue) > 0){
        for(i in 1:length(axis1$kValue)){
          index = axis1$kValue[i]
          index = as.integer(index)
          points(xycoord[index,1],xycoord[index,2],pch=15,col="blue4")
        }
      }
    }
    
    round2 <- function(m,n){
      #print("Round2 OK1")
      if(is.numeric(m) && is.numeric(n)) {
        posneg = sign(m)
        z = abs(m)*10^n
        z = z + 0.5
        z = trunc(z)
        z = z/10^n
        z*posneg
      }
      #print("Round2 OK2")
    }
    
    mark <- function(axis1,cluvi){
      index <- -1
      for (k in 1:length(x)) {
        if(round2(input$plot_click$x,1) == round2(cluvi$data.points[k,1],1) 
           && round2(input$plot_click$y,1) == round2(cluvi$data.points[k,2],1)){
          axis1$selectedRoundTwo[length(axis1$selectedRoundTwo)+1] <- k
        }
      }
      if(length(axis1$selectedRoundTwo) > 0) {
        for (i in 1:length(axis1$selectedRoundTwo)){
          ki = as.integer(axis1$selectedRoundTwo[i])
          k1Point = (cluvi$data.points[ki,1])
          k2Point = (cluvi$data.points[ki,2])
          dXPowTwo <-(k1Point - input$plot_click$x)^2
          dyPowTwo <-(k2Point - input$plot_click$y)^2   
          axis1$distance[length(axis1$distance)+1] <- (dXPowTwo+dyPowTwo)
        }
        index = which.min(axis1$distance)
      }
      if(index > 0) {
        axis1$kValue[length(axis1$kValue)+1] <- axis1$selectedRoundTwo[index]
      }
      axis1$selectedRoundTwo <- list()
      axis1$distance <- list()
      print(axis1)
      return(axis1)
    }
    
    mark1 <- function(axis2,xycoord){
      #print("Mark OK1")
      index <- -1
      for (k in 1:(length(xycoord)/2)) {
        if(round2(input$plot_click$x,1) == round2(xycoord[k,1],1) 
           && round2(input$plot_click$y,1) == round2(xycoord[k,2],1)){
          axis2$selectedRoundTwo[length(axis2$selectedRoundTwo)+1] <- k
        }
      }
    
      #print("Mark OK2")
      if(length(axis2$selectedRoundTwo) > 0) {
        for (i in 1:length(axis2$selectedRoundTwo)){
          ki = as.integer(axis2$selectedRoundTwo[i])
          k1Point = (xycoord[ki,1])
          k2Point = (xycoord[ki,2])
          dXPowTwo <-(k1Point - input$plot_click$x)^2
          dyPowTwo <-(k2Point - input$plot_click$y)^2   
          axis2$distance[length(axis2$distance)+1] <- (dXPowTwo+dyPowTwo)
        }
        index = which.min(axis2$distance)
      }
      #print("Mark OK3")
      if(index > 0) {
        axis2$kValue[length(axis2$kValue)+1] <- axis2$selectedRoundTwo[index]
      }
      axis2$selectedRoundTwo <- list()
      axis2$distance <- list()
      #print(axis2)
      #print("Mark OK4")
      return(axis2)
    }
    
    
    delete <- function(axis1){
      #print("Delete OK")
      axis1$kValue[length(axis1$kValue)] <- NULL
      return(axis1)
    }
    
    
    
    
    make.crisp <- function(u){
      return(apply(u,1,which.max))
    }
    
    compute.radius <- function(x){
      return(sqrt(sum(x^2)))
    }
    
    
    fuci <- 2*pi
    
    cut.angle <- function(alpha){  
      alpha.mod <- alpha
      while(alpha.mod<0 | alpha.mod>fuci){
        if (alpha.mod<0){
          alpha.mod <- alpha.mod + fuci
        }else{
          alpha.mod <- alpha.mod - fuci
        }
      }
      
      return(alpha.mod)
    }
    
    
    angle.diff <- function(alpha1,alpha2){
      andi <- abs(cut.angle(alpha1)-cut.angle(alpha2))
      if (andi>pi){
        andi <- fuci - andi
      }
      return(andi)
    }
    
    angle.mapping <- function(cluster.centre,data.mat){
      xy <- matrix(0,nrow=nrow(data.mat),ncol=2)
      
      data.matrix <- data.mat - matrix(rep(cluster.centre,nrow(data.mat)),nrow=nrow(data.mat),ncol=ncol(data.mat),byrow=T)  
      if (nrow(data.matrix)==1){
        xy[1,1] <- compute.radius(data.matrix[1,])
        xy[1,2] <- 0
        return(xy)
      }
      
      r.values <- apply(data.matrix,1,compute.radius)
      rd.mat <- matrix(Inf,nrow=nrow(data.matrix),ncol=nrow(data.matrix))
      ang.mat <- matrix(0,nrow=nrow(data.matrix),ncol=nrow(data.matrix))
      for (i in 1:(nrow(data.matrix)-1)){
        for (j in (i+1):nrow(data.matrix)){
          if (r.values[i]>0 & r.values[j]>0){
            rd.mat[i,j] <- 2*abs(r.values[i] - r.values[j])/(r.values[i] + r.values[j])
            rd.mat[j,i] <- rd.mat[i,j]
            
            ang.mat[i,j] <- acos(min(c(max(c(t(data.matrix[i,])*t(data.matrix[j,])/(r.values[i]*r.values[j]),-1)),1)))
            ang.mat[j,i] <- ang.mat[i,j]
          } 
        }
      }
      
      
      
      angles <- rep(0,nrow(data.matrix))
      for (i in 1:nrow(data.matrix)){
        if (r.values[i]>0){
          no.pos.i <- sum(r.values[1:i]>0)
          if (no.pos.i==1){
            angles[i] <- 0
          }else{
            if (no.pos.i==2){
              ind <- subset(1:(i-1),r.values[1:(i-1)]>0)
              angles[i] <- ang.mat[i,ind]
            }else{
              ord <- order(rd.mat[1:(i-1),i])
              angplus <- angles[ord[1]] + ang.mat[i,ord[1]]
              angdifp <- angle.diff(angplus,angles[ord[2]])
              
              angminus <- angles[ord[1]] - ang.mat[i,ord[1]]
              angdifm <- angle.diff(angminus,angles[ord[2]])
              
              angles[i] <- min(c(angplus,angminus)) 		  
            }
          }
          xy[i,1] <- r.values[i]*cos(angles[i])	  
          xy[i,2] <- r.values[i]*sin(angles[i])	  
        }
      }
      
      return(xy)
    }
    
  
    
    
  }
  
  shinyApp(server = server,ui = ui)
  
