library(shiny)
library(tidymodels)
library(ggplot2)
library(scales)
library(xgboost)

# =========================================
# Load trained models
# =========================================

xgb_fit_CS <- readRDS("models/xgb_fit_CS.rds")
svm_rad_fit_CI <- readRDS("models/svm_rad_fit_CI.rds")

# =========================================
# USER INTERFACE
# =========================================

ui <- fluidPage(
  
  # ---------------------------------------
  # Title
  # ---------------------------------------
  
  titlePanel(
    
    div(
      
      style = "text-align:center;",
      
      h2(
        "Sex Estimation from Dental Metrics",
        style = "font-size:28px; margin-bottom:5px;"
      ),
      
      h4(
        "Estimación de sexo a partir de métricas dentales",
        style = "font-size:18px; color:gray; margin-top:0px;"
      )
    )
  ),
  
  sidebarLayout(
    
    # =====================================
    # Sidebar inputs
    # =====================================
    
    sidebarPanel(
      
      selectInput(
        "tooth",
        "Select Tooth:",
        choices = c(
          "Maxillary Canine (CS)" = "CS",
          "Mandibular Canine (CI)" = "CI"
        )
      ),
      
      numericInput("MDCo", "MDCo (mm)", value = NA, step = 0.01),
      numericInput("BLCo", "BLCo (mm)", value = NA, step = 0.01),
      numericInput("MDCu", "MDCu (mm)", value = NA, step = 0.01),
      numericInput("BLCu", "BLCu (mm)", value = NA, step = 0.01),
      
      br(),
      
      actionButton("predict_btn", "Predict")
    ),
    
    # =====================================
    # Main panel
    # =====================================
    
    mainPanel(
      
      h3("Prediction"),
      
      verbatimTextOutput("prediction"),
      verbatimTextOutput("probability"),
      
      plotOutput("prob_plot", height = "250px"),
      
      verbatimTextOutput("confidence"),
      
      # -----------------------------------
      # Separator
      # -----------------------------------
      
      tags$hr(),
      
      # ===================================
      # About section
      # ===================================
      
      h4("About / Acerca de"),
      
      p(
        "Machine-learning based application for sex estimation ",
        "using dental metrics from maxillary and mandibular canines."
      ),
      
      p(
        "Aplicación basada en aprendizaje automático para ",
        "estimación de sexo a partir de métricas dentales ",
        "de caninos superiores e inferiores."
      ),
      
      # -----------------------------------
      # Models used
      # -----------------------------------
      
      strong("Implemented models / Modelos implementados:"),
      
      tags$ul(
        tags$li("XGBoost — Maxillary Canine (CS)"),
        tags$li("SVM Radial — Mandibular Canine (CI)")
      ),
      
      p(
        strong("Predictor variables / Variables predictoras:")
      ),
      
      tags$ul(
        tags$li("MDCo: Mesiodistal crown diameter / Diámetro mesiodistal de la corona"),
        tags$li("BLCo: Buccolingual crown diameter / Diámetro bucolingual de la corona"),
        tags$li("MDCu: Mesiodistal cervical diameter / Diámetro mesiodistal del cuello dental"),
        tags$li("BLCu: Buccolingual cervical diameter / Diámetro bucolingual del cuello dental")
      ),
      # -----------------------------------
      # Methodological warning
      # -----------------------------------
      
      p(
        strong("Important / Importante: "),
        "Predictions are probabilistic and should be interpreted ",
        "together with osteological and archaeological evidence."
      ),
      
      # -----------------------------------
      # Repository and DOI
      # -----------------------------------
      
      br(),
      
      strong("Source code and reproducible workflow:"),
      
      br(),
      
      tags$a(
        href = "https://github.com/TU_USUARIO/dental-sex-estimation-ml",
        target = "_blank",
        "GitHub Repository"
      ),
      
      br(),
      
      tags$a(
        href = "https://doi.org/10.5281/zenodo.20380057",
        target = "_blank",
        "Zenodo DOI (archived repository)"
      )
    )
  )
)

# =========================================
# SERVER
# =========================================

server <- function(input, output) {
  
  observeEvent(input$predict_btn, {
    
    # -------------------------------------
    # Create new dataframe
    # Allows decimal commas
    # -------------------------------------
    
    new_data <- data.frame(
      MDCo = as.numeric(gsub(",", ".", input$MDCo)),
      BLCo = as.numeric(gsub(",", ".", input$BLCo)),
      MDCu = as.numeric(gsub(",", ".", input$MDCu)),
      BLCu = as.numeric(gsub(",", ".", input$BLCu))
    )
    
    # -------------------------------------
    # Validation
    # At least one measurement required
    # -------------------------------------
    
    if (all(is.na(new_data))) {
      
      output$prediction <- renderText({
        "Please enter at least one measurement."
      })
      
      output$probability <- renderText("")
      
      output$prob_plot <- renderPlot(NULL)
      
      output$confidence <- renderText("")
      
      return(NULL)
    }
    
    # -------------------------------------
    # Select model according to tooth
    # -------------------------------------
    
    if (input$tooth == "CS") {
      model <- xgb_fit_CS
    } else {
      model <- svm_rad_fit_CI
    }
    
    # -------------------------------------
    # Predictions
    # -------------------------------------
    
    pred_class <- predict(model, new_data)
    
    pred_prob <- predict(model, new_data, type = "prob")
    
    # =====================================
    # Predicted class
    # =====================================
    
    output$prediction <- renderText({
      
      paste("Predicted sex:", pred_class$.pred_class)
    })
    
    # =====================================
    # Probabilities (text)
    # =====================================
    
    output$probability <- renderText({
      
      prob_F <- pred_prob$.pred_F
      
      if (is.na(prob_F)) {
        return("Probability could not be computed.")
      }
      
      prob_M <- 1 - prob_F
      
      paste0(
        "Probability Female: ",
        sprintf("%.1f%%", prob_F * 100),
        "\n",
        "Probability Male: ",
        sprintf("%.1f%%", prob_M * 100)
      )
    })
    
    # =====================================
    # Probability plot
    # =====================================
    
    output$prob_plot <- renderPlot({
      
      prob_F <- pred_prob$.pred_F
      
      if (is.na(prob_F)) return(NULL)
      
      prob_M <- 1 - prob_F
      
      df_plot <- data.frame(
        Sex = c("Female", "Male"),
        Probability = c(prob_F, prob_M)
      )
      
      ggplot(df_plot,
             aes(x = Sex,
                 y = Probability,
                 fill = Sex)) +
        
        geom_col(width = 0.6) +
        
        scale_y_continuous(
          labels = percent_format(),
          limits = c(0, 1)
        ) +
        
        labs(
          title = "Predicted Probabilities",
          y = "Probability",
          x = NULL
        ) +
        
        theme_minimal(base_size = 14) +
        
        theme(
          legend.position = "none"
        )
    })
    
    # =====================================
    # Confidence + Missing data assessment
    # =====================================
    
    output$confidence <- renderText({
      
      prob_F <- pred_prob$.pred_F
      
      if (is.na(prob_F)) {
        return("Confidence could not be computed.")
      }
      
      # -----------------------------------
      # Confidence according to probability
      # -----------------------------------
      
      prob_level <- if (
        prob_F >= 0.75 | prob_F <= 0.25
      ) {
        "High"
      } else if (
        prob_F >= 0.65 | prob_F <= 0.35
      ) {
        "Moderate"
      } else {
        "Low"
      }
      
      # -----------------------------------
      # Missing variables assessment
      # -----------------------------------
      
      n_na <- sum(is.na(new_data))
      
      missing_note <- if (n_na == 0) {
        
        "No missing variables."
        
      } else if (n_na == 1) {
        
        "One variable missing – minimal expected degradation."
        
      } else if (n_na == 2) {
        
        "Two variables missing – moderate expected degradation."
        
      } else {
        
        "Three or more variables missing – substantial expected degradation."
      }
      
      paste0(
        "Confidence level: ", prob_level, "\n",
        "Missing variables: ", n_na, "\n",
        missing_note
      )
    })
  })
}

# =========================================
# Run app
# =========================================

shinyApp(ui = ui, server = server)