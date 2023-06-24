library(shiny)
library(shinydashboard)
library(tidyverse)
library(HDInterval)
library(ggpubr)

samples <- readRDS("samples_no_loglik.rds")
#to standardize based on included data
meanheartw <-442.8382
meanheight <- 174.1788
input <- list()

meanline <- rep(NA, 140)
upperline <- rep(NA,140)
lowerline <- rep(NA, 140)

input$height <- 178
input$heart <- 500

server <- function(input, output) {
  
  # Simulate heart weights
  mydata <- eventReactive(input$goButton, {
    iter = 8000
    

    
    mu_hyper <- log(samples$beta[,2]+samples$beta_height[,2]*(input$height/meanheight))
    
    mu_norm <- log(samples$beta[,1]+samples$beta_height[,1]*(input$height/meanheight))
    # simulate heart weights
    hypertrophy <- rlnorm(iter,
                         mu_hyper,
                        samples$sigma
    )*meanheartw
    
    nonhypertrophy <- rlnorm(iter,
                              mu_norm,
                              samples$sigma_norm
    )*meanheartw
    
    data <- list(
      hypertrophy  = hypertrophy,
      nonhypertrophy = nonhypertrophy,
      iter = iter,
      heart = input$heart
    )
    data
  })

  
  mydata3 <- eventReactive(input$goButton, {

    heartweight <- input$heart/meanheartw
    height <- input$height/meanheight

    uncond_lik <- matrix(nrow=8000,ncol=2)
    uncond_lik[,1] <-  dlnorm(heartweight, log(samples$beta[,1]+samples$beta_height[,1]*height), samples$sigma_norm)
    uncond_lik[,2] <-  dlnorm(heartweight, log(samples$beta[,2]+samples$beta_height[,2]*height), samples$sigma)
    
    bayesfactor <- sum(uncond_lik[,2])/sum(uncond_lik[,1])
    
    
    for(i in 1:140){
      
      uncond_lik2 <- matrix(nrow=8000,ncol=2)
      
      uncond_lik2[,1] <- dlnorm(((i+60)/100),
                               log(samples$beta[,1]+samples$beta_height[,1]*height),
                                           samples$sigma_norm)
      
      uncond_lik2[,2] <- dlnorm(((i+60)/100),
                              log(samples$beta[,2]+samples$beta_height[,2]*height),
                                     samples$sigma)
      
      
      prob2 <-  sum(uncond_lik2[,2])/sum(uncond_lik2[,1])
      
      meanline[i] <- prob2 #dian(prob2)
      # intermed <- hdi(prob2)
      
      # lowerline[i] <- intermed[1]
      # upperline[i] <- intermed[2]
    }
    
    
    #Wrap it in a dataframe for export
    data <- list(
      prob = bayesfactor,
      mean = (meanline),
      # meanlow = (lowerline),
      # meanupp = (upperline),
      hw =(((((1:140)+60)/100))) * meanheartw,
      heart = input$heart
    )
    data
  })
  
  #plot heart weight histograms
  output$plot <- renderPlot({
    
    data <- mydata()
    heart <- data$heart
    
    ggdat <- data.frame(heartweight = c( data$nonhypertrophy, data$hypertrophy),
                        Group = c(rep("Normal", 8000), rep("Hypertrophy", 8000)))
    ggdat$Group <- factor(ggdat$Group, levels = c("Normal", "Hypertrophy"))
    hwplot <- ggplot(data = ggdat, aes(x = heartweight, group = Group, fill = Group))+
      geom_histogram(position = "identity", color = "black", alpha = 0.5, bins = 60) +
      geom_vline(xintercept = heart, linetype = "dashed", size = 1, color = "black") +
      xlab("Heart weight (g)\n Dashed line at measured heart weight") +
      ylab("")+
      scale_fill_manual(values = c("Normal" = "#0073B7", "Hypertrophy"  = "#F39C12"))+
      scale_color_manual(values = c("Normal" = "#1E3247", "Hypertrophy"  = "#C68F43"))+
      scale_x_continuous(limits = c(0, 1000), breaks = c(seq(0,1000,100)))+
      scale_y_continuous(expand = c(0,0)) +
      theme_bw()+
      theme(legend.position="bottom", legend.title = element_blank())
    hwplot
    
  
    })
  
  #plot bayesfactor of hypertrophic probability
  output$plot3 <- renderPlot({
    data <- mydata3()
    heart <- data$heart
    
    data <- data.frame(mean = data$mean,
                       # meanlow = data$meanlow,
                       # meanupp = data$meanupp,
                       hw = data$hw)
    
    scientific_10 <- function(x) {   parse(text=gsub("1e\\+*", "10^", scales::scientific_format()(x))) }
    
    hwplot <- ggplot(data=data) +
      geom_line(data = data, aes(y=mean, x=hw), size = 1.2)+
      # geom_line(data = data, aes(y=meanlow, x = hw)) +
      # geom_line(data = data, aes(y=meanupp, x = hw)) +
      geom_vline(xintercept = heart, linetype = "dashed", size = 1, color = "black") +
      # geom_ribbon(aes(ymin=meanlow, ymax=meanupp, x = hw), fill = "#001F3F", alpha = 0.3)+
      xlab("Heart weight (g) \n Dashed line at measured heart weight")+#, shaded area covers 95% HPDI") +
      ylab("") +
      scale_x_continuous(limits = c(200,700), breaks = seq(200,700,100))+
      # scale_y_continuous(limit = c(0,100), breaks =seq(0,100,10), labels = function(x) paste0(x, "%"))+
      scale_y_continuous(limits = c(NA, 1500), breaks = c(0, 1, 3.2, 10, 10^2, 10^3, 10^4),
                         # labels = scientific_10)+                         
                         labels = c(0, 1, 3.2,10, 10^2, 10^3, 10^4))+
      geom_hline(yintercept = 1, linetype = "dashed", size = 0.5, color = "black") +
      annotate(geom = "text", x = 200, y = 1, hjust = 0, vjust = 0, size = 2.5,
               label=expression(paste(1, " Support barely worth mentioning")))+
      geom_hline(yintercept = 3.2, linetype = "dashed", size = 0.5, color = "black") +
      annotate(geom = "text", x = 200, y = 3.2, hjust = 0, vjust = 0, size = 2.5,
               label=expression(paste(3.2," Substantial support"))) +
      
      geom_hline(yintercept = 10, linetype = "dashed", size = 0.5, color = "black") +
      annotate(geom = "text", x = 200, y = 10, hjust = 0, vjust = 0, size = 2.5,
               label=expression(paste(10," Strong support"))) +
      
      geom_hline(yintercept = 100, linetype = "dashed", size = 0.5, color = "black") +
      annotate(geom = "text", x = 200, y = 100, hjust = 0, vjust = 0, size = 2.5,
               label=expression(paste(100," Decisive support"))) +
      coord_trans(y = "log10")+
      theme_bw() 
    hwplot
    
  })

  # 
  output$text <- renderText({
    data <- mydata()
    
    HDIodds_norm <- signif(hdi(data$nonhypertrophy, credMass = 0.95), digits = 3)
    
    HTML(paste("mean =<b>  ", signif(mean(data$nonhypertrophy), digits = 3),
               " g </b><br>",HDIodds_norm[1]," - ",HDIodds_norm[2]," g, 95% HPDI"))
  })
  
  output$text2 <- renderText({
    data <- mydata()
    
    HDIodds_card <- signif(hdi(data$hypertrophy, credMass = 0.95), digits = 3)
    
    HTML(paste("mean =<b>", signif(mean(data$hypertrophy), digits = 3),
               " g </b><br>",HDIodds_card[1]," - ",HDIodds_card[2]," g, 95% HPDI"))
  })
  
  output$text3 <- renderText({
    data <- mydata3()
    prob <- data$prob
    # HDIodds <- round(hdi(prob, credMass = 0.95), digits = 0)
    HTML(paste("<b>", round(mean(data$prob), digits = 1), "</b>"))#,
               # "</b><br>",HDIodds[1]," - ",HDIodds[2],", 95% HPDI"))
  })
  
  
}

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "calculator",
            fluidRow(
              column(width = 8,
                     box(
                       title = "Histogram of 8,000 simulations",
                       width = NULL,
                       plotOutput("plot")
                     ),
                     box(
                       title = uiOutput("plot2_title"),
                       width = NULL,
                       plotOutput("plot3")
                     )
              ),
              column(width = 4,
                     box(
                       title = "Modeled heart weight",
                       width = NULL,
                       background = "blue",
                       htmlOutput("text")
                     ),
                     box(
                       title = "Modeled hypertrophic heart weight",
                       width = NULL,
                       background = "yellow"
                       ,
                       htmlOutput("text2")
                     ),
                     box(
                       title = "Bayes factor in favour of hypertrophy hypothesis",
                       width = NULL,
                       background = "navy",
                       htmlOutput("text3")
                     ),
                     helpText(HTML(
                       'The source code of this app is <a href="https://github.com/formedum/HeartWeightCalc">on Github</a>.'
                     )
                     ),
                     helpText(HTML(
                       'The results are presented using 95% highest posterior density intervals (HPDI).'
                     )
                     )
                     
                     
              )
            )
    )
  )
)

ui <-dashboardPage(
  dashboardHeader(title = "Heart weight calculator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      numericInput("height", "Height (cm)",180,min=150, max = 200),
      numericInput("heart", "Heart weight (g)", 388 ,min=0, max = 4000),
      actionButton("goButton", "Simulate heart weights†")
    ),
    h6(em("†To improve computation time the app allows for minor simulation error"))
  ),
  body
  
)




shinyApp(ui = ui, server = server)


