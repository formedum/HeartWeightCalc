library(shiny)
library(shinydashboard)
library(tidyverse)
library(HDInterval)
library(emg)
library(ggpubr)

samples <- readRDS("samples.rds")
#to standardize based on included data
meanheartw <-445.059
meanheight <- 174.0564
input <- list()

meanline <- rep(NA, 100)
upperline <- rep(NA,100)
lowerline <- rep(NA, 100)

# Create inverse logit function for convenience
inv_logit <- function( x ) {
  p <- 1 / ( 1 + exp( -x ) )
  p <- ifelse( x==Inf , 1 , p )
  p
}

input$weight <- 60
input$height <- 178
input$heart <- 500
input$sex <- 1
input$age <- 30

server <- function(input, output) {
  
  # Simulate heart weights
  mydata <- eventReactive(input$goButton, {
    iter = 240000
    
    
    sigma <- sample(samples$sigma, iter, replace = TRUE)
    alpha <- sample(samples$alpha, iter, replace = TRUE)
    lambda <- sample(samples$lambda, iter, replace = TRUE)
    hw_incr <- sample(samples$HW_incr, iter, replace = TRUE)
    
    mu <- ((input$height/meanheight)^3)*alpha
    # simulate heart weights
    hypertrophy <- remg(iter,
                        mu*hw_incr,
                        sigma,
                        (lambda/((input$height/meanheight)^3))
    )*meanheartw
    
    nonhypertrophy <- rnorm(iter,
                            mu,
                            sigma
    )*meanheartw
    
    data <- list(
      hypertrophy  = hypertrophy,
      nonhypertrophy = nonhypertrophy,
      iter = iter,
      heart = input$heart
      
    )
    data
  })
  
  # Calculate the probability of hypertrophy
  mydata2 <- eventReactive(input$goButton, {
    inv_logit <- function( x ) {
      p <- 1 / ( 1 + exp( -x ) )
      p <- ifelse( x==Inf , 1 , p )
      p
    }
    
    iter = 1000
    BMI <- ((input$weight/((input$height/100)^2)))
    heartweight <- input$heart/meanheartw
    height <- input$height/meanheight
    age <- input$age
    
    iter2 <- sample(1:24000, iter)
    sigma <- samples$sigma[iter2]
    alpha <- samples$alpha[iter2]
    hw_incr <- samples$HW_incr[iter2]
    lambda <- samples$lambda[iter2]
    
    theta <- samples$theta[iter2]
    theta_beta_BMI <- samples$theta_beta_BMI[iter2,]
    theta_beta_age <- samples$theta_beta_age[iter2,]
    
    if(age < 40){
      agecat <- 1
    }
    if(age>= 40 & age< 50){
      agecat <- 2
    }
    if(age>= 50 & age< 60){
      agecat <- 3
    }
    if(age>= 60 & age< 70){
      agecat <- 4
    }
    if(age>= 70){
      agecat <- 5
    }
    
    if(BMI < 18.5){
      BMIcat <- 1
    }
    if(BMI>= 18.5 & BMI< 25){
      BMIcat <- 2
    }
    if(BMI>= 25 & BMI< 30){
      BMIcat <- 3
    }
    if(BMI>= 30 & BMI< 35){
      BMIcat <- 4
    }
    if(BMI>= 35){
      BMIcat <- 5
    }
    
    baserate <- inv_logit(theta+theta_beta_BMI[,BMIcat]+theta_beta_age[,agecat])
    mu <- (height^3)*(alpha)
    cond_lik <- matrix(nrow=iter,ncol=2)
    cond_lik[,1] <-  (1-baserate) * dnorm(heartweight, mu, sigma)
    cond_lik[,2] <-   baserate* demg(heartweight, mu*hw_incr, sigma, lambda/(height^3))
    
    marginal.prob <- rep(NA, iter)
    for(i in 1:iter)
      marginal.prob[i] <- sum(cond_lik[i,1:2])
    
    postprob <- cond_lik/marginal.prob
    
    
    for(i in 1:300){
      
      # i = 100
      cond_lik2 <- matrix(nrow=iter,ncol=2)
      
      cond_lik2[,1] <- (1-baserate)*dnorm((i/100),
                                          mean= mu,
                                          sd = sigma)
      
      cond_lik2[,2] <- baserate*demg((i/100),
                                     mu,
                                     sigma,
                                     (lambda/((input$height/meanheight)^3)))
      
      marginal.prob2 <- rep(NA, iter)
      for(n in 1:iter){
        marginal.prob2[n] <- sum(cond_lik2[n,])
      }
      
      prob2 <- cond_lik2/marginal.prob2
      
      meanline[i] <- mean(prob2[,2])
      intermed <- hdi(prob2[,2])
      
      lowerline[i] <- intermed[1]
      upperline[i] <- intermed[2]
    }
    
    
    #Wrap it in a dataframe for export
    data <- list(
      prob = postprob[,2]*100,
      mean = (meanline)*100,
      meanlow = (lowerline)*100,
      meanupp = (upperline)*100,
      hw =((((1:300)/100))) * meanheartw,
      iter = iter,
      heart = input$heart
    )
    data
  })
  
  mydata3 <- eventReactive(input$goButton, {
    inv_logit <- function( x ) {
      p <- 1 / ( 1 + exp( -x ) )
      p <- ifelse( x==Inf , 1 , p )
      p
    }
    
    iter = 1000
    
    heartweight <- input$heart/meanheartw
    height <- input$height/meanheight
    
    iter2 <- sample(1:24000, iter)
    sigma <- samples$sigma[iter2]
    alpha <- samples$alpha[iter2]
    hw_incr <- samples$HW_incr[iter2]
    lambda <- samples$lambda[iter2]
    
    
    mu <- (height^3)*(alpha)
    uncond_lik <- matrix(nrow=iter,ncol=2)
    uncond_lik[,1] <-  dnorm(heartweight, mu, sigma)
    uncond_lik[,2] <- demg(heartweight, mu*hw_incr, sigma, lambda/(height^3))
    
    bayesfactor <- uncond_lik[,2]/uncond_lik[,1]
    
    
    for(i in 1:300){
      
      uncond_lik2 <- matrix(nrow=iter,ncol=2)
      
      uncond_lik2[,1] <- dnorm((i/100),
                               mean= mu,
                               sd = sigma)
      
      uncond_lik2[,2] <- demg((i/100),
                              mu,
                              sigma,
                              (lambda/((input$height/meanheight)^3)))
      
      
      prob2 <- uncond_lik2[,2]/uncond_lik2[,1]
      
      meanline[i] <- mean(prob2)
      intermed <- hdi(prob2)
      
      lowerline[i] <- intermed[1]
      upperline[i] <- intermed[2]
    }
    
    
    #Wrap it in a dataframe for export
    data <- list(
      prob = bayesfactor,
      mean = (meanline),
      meanlow = (lowerline),
      meanupp = (upperline),
      hw =((((1:300)/100))) * meanheartw,
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
    
    ggdat <- data.frame(heartweight = c( data$nonhypertrophy, data$hypertrophy),
                        Group = c(rep("Normal", iter), rep("Hypertrophy", iter)))
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
  #plot the probability of Hypertrophy
  output$plot2 <- renderPlot({
    data <- mydata2()
    iter <- data$iter
    heart <- data$heart
    
    data <- data.frame(mean = data$mean,
                       meanlow = data$meanlow,
                       meanupp = data$meanupp,
                       hw = data$hw)
    
    hwplot <- ggplot(data=data) +
      geom_line(data = data, aes(y=mean, x=hw), size = 0.5)+
      geom_line(data = data, aes(y=meanlow, x = hw)) +
      geom_line(data = data, aes(y=meanupp, x = hw)) +
      geom_vline(xintercept = heart, linetype = "dashed", size = 1, color = "black") +
      geom_ribbon(aes(ymin=meanlow, ymax=meanupp, x = hw), fill = "#001F3F", alpha = 0.3)+
      xlab("Heart weight (g) \n Dashed line at measured heart weight, shaded area covers 95% HPDI") +
      ylab("") +
      scale_x_continuous(limits = c(200,700), breaks = seq(200,700,100))+
      scale_y_continuous(limit = c(0,100), breaks =seq(0,100,10), labels = function(x) paste0(x, "%"))+
      theme_bw() 
    
    hwplot
    
  })
  #plot bayesfactor of hypertrophic probability
  output$plot3 <- renderPlot({
    data <- mydata3()
    iter <- data$iter
    heart <- data$heart
    
    data <- data.frame(mean = data$mean,
                       meanlow = data$meanlow,
                       meanupp = data$meanupp,
                       hw = data$hw)
    
    scientific_10 <- function(x) {   parse(text=gsub("1e\\+*", "10^", scales::scientific_format()(x))) }
    
    hwplot <- ggplot(data=data) +
      geom_line(data = data, aes(y=mean, x=hw), size = 0.5)+
      geom_line(data = data, aes(y=meanlow, x = hw)) +
      geom_line(data = data, aes(y=meanupp, x = hw)) +
      geom_vline(xintercept = heart, linetype = "dashed", size = 1, color = "black") +
      geom_ribbon(aes(ymin=meanlow, ymax=meanupp, x = hw), fill = "#001F3F", alpha = 0.3)+
      xlab("Heart weight (g) \n Dashed line at measured heart weight, shaded area covers 95% HPDI") +
      ylab("") +
      scale_y_continuous(limits = c(NA, 10^8), breaks = c(0, 1, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8),
                         # labels = scientific_10)+                         
                         labels = c(0, 1, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8))+
      scale_x_continuous(limits = c(200,700), breaks = seq(200,700,100))+
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
      # ymax(1000)+
      theme_bw()
    #theme( axis.text.y = element_blank())
    
    
    hwplot
    
  })
  
  output$plot_choose = renderUI({
    if(input$type==0){
      # paste("Probability of Hypertrophy")
      plotOutput("plot2")
    }else{
      plotOutput("plot3")
    }
  })
  
  output$plot2_title = renderUI({ #title of BF or prob plot
    
    if(input$type==0){
      paste("Probability of hypertrophy")
      
    }else{
      paste("Bayes factor of hypertrophy")
    }
    
  })
  
  output$box3_title = renderUI({ #title of BF or probbox
    
    if(input$type==0){
      paste("Probability of hypertrophy")
      
    }else{
      paste("Bayes factor of hypertrophy")
    }
    
  })
  
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
    data <- mydata2()
    prob <- data$prob
    HDIodds <- signif(hdi(prob, credMass = 0.95), digits = 2)
    HTML(paste("mean = <b>", signif(mean(data$prob), digits = 2),
               "% </b><br>",HDIodds[1]," - ",HDIodds[2],"%, 95% HPDI"))
  })
  
  output$text4 <- renderText({
    data <- mydata3()
    prob <- data$prob
    HDIodds <- round(hdi(prob, credMass = 0.95), digits = 0)
    HTML(paste("mean = <b>", round(mean(data$prob), digits = 0),
               "</b><br>",HDIodds[1]," - ",HDIodds[2],", 95% HPDI"))
  })
  
  output$text_choose = renderUI({
    if(input$type==0){
      htmlOutput("text3")
    }else{
      htmlOutput("text4")
    }
  })
  
}

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "calculator",
            fluidRow(
              column(width = 8,
                     box(
                       title = "Histogram of 240,000 simulations",
                       width = NULL,
                       plotOutput("plot")
                     ),
                     box(
                       title = uiOutput("plot2_title"),
                       width = NULL,
                       uiOutput("plot_choose")
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
                     # uiOutput("box2"),
                     box(
                       # title = "Bayes factor or probability of hypertrophy hypothesis",
                       title = uiOutput("box3_title"),
                       width = NULL,
                       background = "navy",
                       uiOutput("text_choose")
                       # htmlOutput("text3")
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
    ,
    tabItem(tabName = "information",
            withMathJax(),
            helpText("The probability of hypertrophy is calculated as 
                     $$Pr(hypertrophy|heart\\ weight)=$$"),
            helpText("$$\\frac{Pr(hypertrophy)\\times Pr(heart\\ weight|hypertrophy)}{Pr(hypertrophy)\\times Pr(heart\\ weight|hypertrophy)+Pr(normal)\\times Pr(heart\\ weight|normal)}
$$"),
            helpText(div(HTML("The probability is a function both of the likelihoods (i.e., <em>Pr(heart weight|hypertrophy)</em> and <em>Pr(heart weight|normal)</em> ) and the
                              baserate (<em>(Pr(hypertrophy)</em>). As the baserate is a function of age and BMI 
                                       there are subpopulations where the baserate is very high, 
                                       meaning that the a priori there is a very slim probability of that heart not being hypertrophic"))),
            helpText(div(HTML("If the decedent in question belongs to such a population (e.g., if BMI > 35 and
                    age > 70) the baserate of hypertrophy is so high that even if a normal or low heart weight is measured (i.e., the likelihood of hypertrophy [<em>Pr(heart weight|hypertrophy)</em>] 
                    is very low) the posterior probability of hypertrophy would still be high. As such, when baserates are very large or small 
                                      it is instead recommended to instead use the Bayes factor (BF)."))),
            helpText("$$BF=\\frac{Pr(heart\\ weight|hypertrophy)}{Pr(heart\\ weight|normal)}$$"),
            helpText("The BF is the ratio between the two likelihoods and is a continous value between 0 and infinity. 
                     It is harder to interpret than a pure probability but it has the advantage of being invariant to the baserate.
                     The BF reflects the degree to which only the measured  heart weight itself supports the diagnosis of hypertrophy."),
            helpText("The interpretation guidelines presented in the graph are adapted from"),
            helpText(div(HTML("Kass RE, Raftery AE.  Bayes Factors. <em>Journal of the American Statistical Association</em> 1995;90(430):773–95. <a href =https://doi.org/10.1080/01621459.1995.10476572>https://doi.org/10.1080/01621459.1995.10476572</a>.
"))),
    )
  )
)

ui <-dashboardPage(
  dashboardHeader(title = "Heart weight calculator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      menuItem("Information*", tabName = "information", icon = icon("info")),
      numericInput("age", "Age",60,min=18, max = 150),
      radioButtons("type", "Presentation type*",
                   c("Posterior probability" = 0,
                     "Bayes factor" = 1), selected = 1),
      numericInput("height", "Height (cm)",180,min=150, max = 220),
      numericInput("weight", "Weight (kg)",80,min=0, max = 120),
      numericInput("heart", "Heart weight (g)", 388 ,min=0, max = 4000),
      actionButton("goButton", "Simulate heart weights†")
      
    ),
    h6(em("†To improve computation time the app allows for minor simulation error"))
  ),
  body
  
)

if(input$heart==1){
  
}



shinyApp(ui = ui, server = server)





