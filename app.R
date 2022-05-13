library(shiny)
library(shinydashboard)
library(tidyverse)
library(HDInterval)
library(emg)
library(ggpubr)


samples <- readRDS("samples.rds")
#to standardize based on included data
meanheartw <-445.059
sdheartw <- 122.5996
meanage <- 54.45597
sdage <- 17.12676
meanweight <- 82.19719
sdweight <- 22.00117
meanheight <- 174.0564
sdheight <-9.648763

meanbmi <- 26.97812
sdbmi <- 6.327531

meanline <- rep(NA, 100)
upperline <- rep(NA,100)
lowerline <- rep(NA, 100)

# Create inverse logit function for convenience
inv_logit <- function( x ) {
  p <- 1 / ( 1 + exp( -x ) )
  p <- ifelse( x==Inf , 1 , p )
  p
}

server <- function(input, output) {

  # Simulate heart weights
  mydata <- eventReactive(input$goButton, {
    iter = 50000

    ##Create mu
    mu_prep <- rep(0, 24000)
    mu_prep <- (samples$gamma_sex*as.numeric(input$sex)+
                  samples$gamma_age*((input$age-meanage)/sdage)+
                  samples$gamma_weight*((input$weight-meanweight)/sdweight)+
                  samples$gamma_height*((input$height-meanheight)/sdheight))
    mu_prep <- sample(mu_prep, iter, replace = TRUE )

    beta <- sample(samples$beta, iter, replace = TRUE)
    sigma <- sample(samples$sigma, iter, replace = TRUE)
    lambda <- sample(samples$lambda, iter, replace = TRUE)
    
    # simulate heart weights
    cardiomegaly <- remg(iter,
                              mu = beta + mu_prep,
                              sigma = sigma,
                              lambda = lambda
                              )*sdheartw + meanheartw
    noncardiomegaly <- rnorm(
      iter,
      mean = beta+mu_prep,
      sd = sigma)*sdheartw + meanheartw

    data <- list(
      cardiomegaly  = cardiomegaly,
      noncardiomegaly = noncardiomegaly,
      iter = iter,
      heart = input$heart

    )
    data
  })

  # Calculate the probability of cardiomegaly
  mydata2 <- eventReactive(input$goButton, {
    iter = 1000
    bmi <- input$weight/((input$height/100)^2)
    
    mu_prep2 <- rep(NA, 24000)
    mu_prep2 <- (samples$gamma_sex*as.numeric(input$sex)+
                  samples$gamma_age*((input$age-meanage)/sdage)+
                  samples$gamma_weight*((input$weight-meanweight)/sdweight)+
                  samples$gamma_height*((input$height-meanheight)/sdheight))

  lik_norm <- dnorm((input$heart - meanheartw)/sdheartw,
                  mean = samples$beta + mu_prep2,
                  sd = samples$sigma)

  lik_card <- demg((input$heart - meanheartw)/sdheartw,
                   mu = (samples$beta + mu_prep2),
                   sigma = samples$sigma,
                   lambda = samples$lambda)

  theta_rep <- rep(0, 24000)
  theta_rep <-inv_logit(samples$theta_alpha +
                          samples$theta_beta_age*((input$age-meanage)/sdage)+
                          samples$theta_beta_sex*as.numeric(input$sex) +
                         samples$theta_beta_bmi*((bmi-meanbmi)/sdbmi))

  bf <- lik_card/lik_norm
  podds <- theta_rep/(1-theta_rep)
  
  postprob <- bf*podds
  
  postprob <- postprob/(1+  postprob)
  


  beta <- sample(samples$beta, iter, replace = TRUE)
  sigma <- sample(samples$sigma, iter, replace = TRUE)
  lambda <- sample(samples$lambda, iter, replace = TRUE)
  mu_prep2 <- sample(mu_prep2, iter, replace = TRUE)
  theta_rep <- sample(theta_rep, iter, replace = TRUE)

  for(i in 1:100){

  noncardioest <- dnorm((i/10)-5,
                        mean= beta + mu_prep2,
                        sd = sigma)

    cardioest <- demg((i/10)-5,
                      mu = beta + mu_prep2,
                      sigma = sigma,
                      lambda = lambda)
 
      
      bf <- cardioest/noncardioest
      podds <- theta_rep/(1-theta_rep)
      
      prob2 <- bf*podds
      
      prob2 <- prob2/(1+prob2)

    meanline[i] <- mean(prob2)
    intermed <- hdi(prob2)

    lowerline[i] <- intermed[1]
    upperline[i] <- intermed[2]
  }


  #Wrap it in a dataframe for export
  data <- list(
    prob = postprob*100,
    mean = (meanline)*100,
    meanlow = (lowerline)*100,
    meanupp = (upperline)*100,
    hw =((((1:100)/10)-5) * sdheartw) + meanheartw,
    iter = iter,
    heart = input$heart
  )
  data
})

  #plot heart weight histograms
output$plot <- renderPlot({

  data <- mydata()
  iter <- data$iter
  heart <- data$heart

  ggdat <- data.frame(heartweight = c( data$noncardiomegaly, data$cardiomegaly),
                      Group = c(rep("Normal", iter), rep("Cardiomegaly", iter)))
  ggdat$Group <- factor(ggdat$Group, levels = c("Normal", "Cardiomegaly"))
  hwplot <- ggplot(data = ggdat, aes(x = heartweight, group = Group, fill = Group))+
    geom_histogram(position = "identity", color = "black", alpha = 0.5) +
    geom_vline(xintercept = heart, linetype = "dashed", size = 1, color = "black") +
    xlab("Heart weight (g)\n Dashed line at measured heart weight") +
    ylab("")+
    scale_fill_manual(values = c("Normal" = "#0073B7", "Cardiomegaly"  = "#F39C12"))+
    scale_color_manual(values = c("Normal" = "#1E3247", "Cardiomegaly"  = "#C68F43"))+
    scale_x_continuous(limits = c(0, 1000), breaks = c(seq(0,1000,100)))+
    scale_y_continuous(expand = c(0,0)) +
    theme_bw()+
    theme(legend.position="bottom", legend.title = element_blank())
  hwplot

})
#plot the probability of cardiomegaly
  output$plot2 <- renderPlot({
    data <- mydata2()
    iter <- data$iter
    heart <- data$heart

    data <- as.data.frame(data)
    hwplot <- ggplot(data=data) +
      geom_line(data = data, aes(y=mean, x=hw), size = 0.5)+
      geom_line(data = data, aes(y=meanlow, x = hw)) +
      geom_line(data = data, aes(y=meanupp, x = hw)) +
      geom_vline(xintercept = heart, linetype = "dashed", size = 1, color = "black") +
      geom_hline(yintercept = 50, linetype = "dotted", size = 0.5, color = "black") +
      geom_ribbon(aes(ymin=meanlow, ymax=meanupp, x = hw), fill = "#001F3F", alpha = 0.3)+
      xlab("Heart weight (g) \n Dashed line at measured heart weight, shaded area covers 95% HPDI") +
      ylab("") +
      scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%"))+
      scale_x_continuous(limits = c(200,700), breaks = c(seq(200,700,100)))+
      annotate(geom = "text", x = 200, y = 55, hjust = 0, vjust = 0, size = 3,
               label=paste("Evidence for cardiomegaly")) +
      annotate(geom = "text", x = 200, y = 45, hjust = 0, vjust = 1, size = 3,
               label=paste("Evidence against cardiomegaly")) +
      theme_bw()


    hwplot

  })
  
  # Plot exponentially modified gaussian examples
  output$plot3 <- renderPlot({

    
    hwplot1 <- ggplot() +
      geom_histogram(aes(x = rnorm(100000,0,1)), fill = "#1E3247", alpha =.5, bins = 50) +
      xlim(-10, 10) +
      labs( x = expression(paste(mu==0~","~sigma==1)),
            title = "Normal distribution")+
      theme(plot.title = element_text(size=10), 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank())
    
    hwplot2 <- ggplot() +
      geom_histogram(aes(x = rexp(100000, .4)), fill = "#C68F43", alpha =.5, bins = 50) +
      xlim(-10, 10) +
      labs( x = expression(paste(lambda==0.4)),
            title = "Exponential distribution")+
      theme(plot.title = element_text(size=10), 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank())
    
    hwplot3 <- ggplot() +
      geom_histogram(aes(x = remg(100000,0,1,.4)), fill ="#001F3F", alpha =.5, bins = 50) +
      xlim(-10, 10) +
      labs( x = expression(paste(mu==0~","~sigma==1~","~lambda==0.4)),
            title = "Exponentially modified normal distribution")+
      theme(plot.title = element_text(size=10), 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank())
    
    hwplot4 <- ggarrange(hwplot1, hwplot2, hwplot3, ncol =1) 
      
    hwplot4
    
  })
  output$text <- renderText({
    data <- mydata()

    HDIodds_norm <- signif(hdi(data$noncardiomegaly, credMass = 0.95), digits = 3)

    HTML(paste("mean =<b>  ", signif(mean(data$noncardiomegaly), digits = 3),
               " g </b><br>",HDIodds_norm[1]," - ",HDIodds_norm[2]," g, 95% HPDI"))
  })

  output$text2 <- renderText({
    data <- mydata()

    HDIodds_card <- signif(hdi(data$cardiomegaly, credMass = 0.95), digits = 3)

    HTML(paste("mean =<b>", signif(mean(data$cardiomegaly), digits = 3),
               " g </b><br>",HDIodds_card[1]," - ",HDIodds_card[2]," g, 95% HPDI"))
  })

  output$text3 <- renderText({
    data <- mydata2()

    HDIodds <- signif(hdi(data$prob, credMass = 0.95), digits = 2)
    HTML(paste("mean = <b>", signif(mean(data$prob), digits = 2),
               "% </b><br>",HDIodds[1]," - ",HDIodds[2],"%, 95% HPDI"))
})

}

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "calculator",
            fluidRow(
              column(width = 8, 
                     box(
                       title = "Histogram of 50,000 simulations",
                       width = NULL
                       ,
                       plotOutput("plot")
                     ),
                     box(
                       title = "Probability of cardiomegaly",
                       width = NULL
                       ,
                       plotOutput("plot2")
                     )
              ),
              column(width = 4,
                     box(
                       title = "Modeled heart weight",
                       width = NULL,
                       background = "blue"
                       ,
                       htmlOutput("text")
                     ),
                     box(
                       title = "Modeled cardiomegalic heart weight",
                       width = NULL,
                       background = "yellow"
                       ,
                       htmlOutput("text2")
                     ),
                     box(
                       title = "Probability of cardiomegaly",
                       width = NULL,
                       background = "navy"
                       ,
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
            ), style = "width:650px;align:center;"
    )
    ,
    tabItem(tabName = "information",
            withMathJax(),
            helpText("This heart weight calculator models heart weight as a mixture of normal and cardiomegalic hearts.
                     Normal hearts are modeled as
                     $$Heart\\ weight \\sim normal(\\mu, \\sigma)$$"),
            helpText("while cardiomegalic hearts are modeled as an exponentially modified normal distribution
                     $$Cardiomegalic\\ heart\\ weight \\sim ExModNormal(\\mu, \\sigma, \\lambda)$$"),
            helpText("Since both distributions share the same mean and standard deviation this equivalent
                     to cardiomegalic hearts being modeled as
                     $$Cardiomegalic\\ heart weight = normal\\ heart\\ weight + exponential(\\lambda)$$
                     as seen in the graphs below."),
            plotOutput("plot3")

    )
  )
)

ui <-dashboardPage(
  dashboardHeader(title = "Heart weight calculator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      menuItem("Information", tabName = "information", icon = icon("info")),
      numericInput("age", "Age",60,min=18, max = 150),
      radioButtons("sex", "Sex",
                   c("Female" = 0,
                     "Male" = 1), selected = 1),
      numericInput("height", "Height (cm)",180,min=0, max = 300),
      numericInput("weight", "Weight (kg)",80,min=0, max = 300),
      numericInput("heart", "Heart weight (g)", 540 ,min=0, max = 4000),
      actionButton("goButton", "Simulate heart weights*")

      ),
    h6(em("*To improve computation time the app allows for minor simulation error"))
    ),
  body

  )



shinyApp(ui = ui, server = server)


