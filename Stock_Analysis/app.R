library(shiny)
library(quantmod)
library(TTR)
library(dplyr)
library(ggplot2)
library(shinythemes)
library(prophet)
library(tseries)
library(PerformanceAnalytics)
library(quadprog)

ui <- fluidPage(
  theme = shinytheme("darkly"),
  titlePanel(
    div(
      img(src = "logo.png", height = "100px"),
      "Securities Analysis"
    )
  ),
  sidebarLayout(
    sidebarPanel(
      textInput("stock1", "Enter first stock symbol:", value = "AAPL"),
      textInput("stock2", "Enter second stock symbol:", value = "GOOG"),
      dateRangeInput("dates", "Date range:", start = Sys.Date() - 365, end = Sys.Date())
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Comparison",
                 plotOutput("stockPlot1"),
                 plotOutput("stockPlot2"),
                 tableOutput("stockMetrics")
        ),
        tabPanel("Portfolio Optimization",
                 plotOutput("portfolioPlot"),
                 plotOutput("efficientFrontierPlot")
        ),
        tabPanel("30-Day Forecast",
                 plotOutput("forecastPlot1"),
                 plotOutput("forecastPlot2")
        )
      )
    )
  )
)

server <- function(input, output) {
  
  stockData <- reactive({
    symbols <- c(input$stock1, input$stock2)
    data_list <- list()
    for (sym in symbols) {
      tryCatch({
        data_list[[sym]] <- getSymbols(sym, src = 'yahoo', from = input$dates[1], to = input$dates[2], auto.assign = FALSE)
      }, error = function(e) {
        data_list[[sym]] <- NULL
      })
    }
    data_list
  })
  
  stockMetrics <- reactive({
    stocks <- c(input$stock1, input$stock2)
    
    metric_template <- data.frame(
      Stock = character(),
      Avg_Close = numeric(),
      Avg_Volume = numeric(),
      Std_Dev = numeric(),
      stringsAsFactors = FALSE
    )
    
    metrics <- lapply(stocks, function(stock) {
      data <- stockData()[[stock]]
      if (is.null(data)) return(metric_template[0, ])
      
      avg_close <- mean(Cl(data), na.rm = TRUE)
      avg_volume <- mean(Vo(data), na.rm = TRUE)
      std_dev <- sd(Cl(data), na.rm = TRUE)
      
      result <- data.frame(
        Stock = stock,
        Avg_Close = avg_close,
        Avg_Volume = avg_volume,
        Std_Dev = std_dev,
        stringsAsFactors = FALSE
      )
      
      return(result)
    })
    
    metrics <- lapply(metrics, function(df) {
      colnames(df) <- colnames(metric_template)
      df
    })
    
    metrics <- do.call(rbind, metrics)
    metrics
  })
  
  output$stockPlot1 <- renderPlot({
    stock1 <- stockData()[[input$stock1]]
    if (!is.null(stock1)) {
      chartSeries(stock1, theme = chartTheme("white"), TA = "addSMA(n = 20, col = 'blue')")
      title(main = input$stock1, line = -1)
    }
  })
  
  output$stockPlot2 <- renderPlot({
    stock2 <- stockData()[[input$stock2]]
    if (!is.null(stock2)) {
      chartSeries(stock2, theme = chartTheme("white"), TA = "addSMA(n = 20, col = 'blue')")
      title(main = input$stock2, line = -1)
    }
  })
  
  output$stockMetrics <- renderTable({
    stockMetrics()
  })
  
  output$portfolioPlot <- renderPlot({
    stock1 <- stockData()[[input$stock1]]
    stock2 <- stockData()[[input$stock2]]
    
    if (!is.null(stock1) && !is.null(stock2)) {
      returns1 <- na.omit(ROC(Cl(stock1)))
      returns2 <- na.omit(ROC(Cl(stock2)))
      
      returns <- cbind(returns1, returns2)
      colnames(returns) <- c(input$stock1, input$stock2)
      
      cov_matrix <- cov(returns)
      exp_returns <- colMeans(returns)
      
      Dmat <- 2 * cov_matrix
      dvec <- rep(0, 2)
      Amat <- cbind(rep(1, 2), diag(2))
      bvec <- c(1, rep(0, 2))
      
      sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=1)
      weights <- sol$solution
      
      barplot(weights, names.arg = c(input$stock1, input$stock2), col = c("blue", "red"), main = "Optimized Portfolio Weights")
    }
  })
  
  output$efficientFrontierPlot <- renderPlot({
    stock1 <- stockData()[[input$stock1]]
    stock2 <- stockData()[[input$stock2]]
    
    if (!is.null(stock1) && !is.null(stock2)) {
      returns1 <- na.omit(ROC(Cl(stock1)))
      returns2 <- na.omit(ROC(Cl(stock2)))
      
      returns <- cbind(returns1, returns2)
      colnames(returns) <- c(input$stock1, input$stock2)
      
      portfolio <- tseries::portfolio.optim(returns)
      
      frontier <- function(returns, n_portfolios = 100) {
        results <- matrix(0, nrow = n_portfolios, ncol = 3)
        colnames(results) <- c("return", "risk", "sharpe")
        
        for (i in seq_len(n_portfolios)) {
          weights <- runif(ncol(returns))
          weights <- weights / sum(weights)
          
          port_return <- sum(weights * colMeans(returns))
          port_risk <- sqrt(t(weights) %*% cov(returns) %*% weights)
          sharpe_ratio <- port_return / port_risk
          
          results[i, ] <- c(port_return, port_risk, sharpe_ratio)
        }
        return(as.data.frame(results))
      }
      
      frontier_data <- frontier(returns)
      
      ggplot() +
        geom_point(data = frontier_data, aes(x = risk, y = return), color = 'blue') +
        geom_point(aes(x = portfolio$ps, y = portfolio$pm), color = 'red', size = 3) +
        labs(title = "Markowitz Efficient Frontier", x = "Risk (Standard Deviation)", y = "Return") +
        theme_minimal(base_family = "Arial") +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, color = "white"),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white")
        )
    }
  })
  
  output$forecastPlot1 <- renderPlot({
    stock1 <- stockData()[[input$stock1]]
    if (!is.null(stock1)) {
      stock1_df <- data.frame(ds = index(stock1), y = as.numeric(Cl(stock1)))
      
      m <- prophet(stock1_df)
      future <- make_future_dataframe(m, periods = 30)
      forecast <- predict(m, future)
      
      ggplot(forecast) +
        geom_line(aes(x = ds, y = yhat), color = 'blue') +
        geom_ribbon(aes(x = ds, ymin = yhat_lower, ymax = yhat_upper), alpha = 0.2) +
        ggtitle(paste("30-Day Forecast for", input$stock1))
    }
  })
  
  output$forecastPlot2 <- renderPlot({
    stock2 <- stockData()[[input$stock2]]
    if (!is.null(stock2)) {
      stock2_df <- data.frame(ds = index(stock2), y = as.numeric(Cl(stock2)))
      
      m <- prophet(stock2_df)
      future <- make_future_dataframe(m, periods = 30)
      forecast <- predict(m, future)
      
      ggplot(forecast) +
        geom_line(aes(x = ds, y = yhat), color = 'red') +
        geom_ribbon(aes(x = ds, ymin = yhat_lower, ymax = yhat_upper), alpha = 0.2) +
        ggtitle(paste("30-Day Forecast for", input$stock2))
    }
  })
}

shinyApp(ui, server)
