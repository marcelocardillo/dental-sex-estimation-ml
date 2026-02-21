library(shiny)
library(tidymodels)

# Cargar modelos
xgb_fit_CS <- readRDS("models/xgb_fit_CS.rds")
svm_rad_fit_CI <- readRDS("models/svm_rad_fit_CI.rds")

ui <- fluidPage(
  
  titlePanel("Sex Estimation from Dental Metrics"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("tooth",
                  "Select Tooth:",
                  choices = c("Maxillary Canine (CS)" = "CS",
                              "Mandibular Canine (CI)" = "CI")),
      
      numericInput("MD", "MD", value = NA),
      numericInput("BL", "BL", value = NA),
      numericInput("MDCu", "MDCu", value = NA),
      numericInput("BLCu", "BLCu", value = NA),
      
      actionButton("predict_btn", "Predict")
    ),
    
    mainPanel(
      h3("Prediction"),
      verbatimTextOutput("prediction"),
      verbatimTextOutput("probability"),
      verbatimTextOutput("confidence")
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$predict_btn, {
    
    new_data <- data.frame(
      MD = input$MD,
      BL = input$BL,
      MDCu = input$MDCu,
      BLCu = input$BLCu
    )
    
    if (input$tooth == "CS") {
      model <- xgb_fit_CS
    } else {
      model <- svm_rad_fit_CI
    }
    
    pred_class <- predict(model, new_data)
    pred_prob  <- predict(model, new_data, type = "prob")
    
    output$prediction <- renderText({
      paste("Predicted sex:", pred_class$.pred_class)
    })
    
    output$probability <- renderText({
      paste("Probability (F):", round(pred_prob$.pred_F, 3))
    })
    
    output$confidence <- renderText({
      if (pred_prob$.pred_F >= 0.75 | pred_prob$.pred_F <= 0.25) {
        "Confidence level: High"
      } else if (pred_prob$.pred_F >= 0.65 | pred_prob$.pred_F <= 0.35) {
        "Confidence level: Moderate"
      } else {
        "Confidence level: Low"
      }
    })
    
  })
}

shinyApp(ui = ui, server = server)