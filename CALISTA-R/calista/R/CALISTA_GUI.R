#' A calista Function
#'
#'  Usage:
#'  calista(INPUTS)
#' @param INPUTS
#' @keywords calista
#' @export
#' @examples
#' CALISTA_GUI()




calista <- function(INPUTS) {
  #source("./R/initialization.R")
  ## Declare global variables
  #Results <<- Results
  #DATA <<-DATA
  INPUTS <<- INPUTS
  Results <<- list()
  aa<<-1
  #DATA<<-list()

  ## Execute the first part of CALISTA
  #source("./R/initialization.R")
  # % Upload and pre-process data
  DATA <<- import_data(INPUTS)

  # %% *** 2-SINGLE-CELL CLUSTERING ***
  # %
  # % Please check comments in 'CALISTA_clustering_main' for more information.
  # %


  CALISTA_clustering_main_results=CALISTA_clustering_main(DATA,INPUTS)
  Results <<- CALISTA_clustering_main_results$Results
  DATA <<- CALISTA_clustering_main_results$DATA
  INPUTS <<- CALISTA_clustering_main_results$INPUTS

  ###cluster cutting
  writeLines('Press 1 if you want to remove one cell cluster, 0 otherwise:')
  cluster_cut=readLines(n=1)
  cluster_cut=as.integer(cluster_cut)

  if(cluster_cut==1){
    writeLines('key cluster number(e.g 1 or 5 3)')
    which_cut=strsplit(readLines(n=1), " ")
    which_cut=as.integer(unlist(which_cut))
    cells_2_cut2=which(Results$final_groups %in% which_cut)
    cells_2_cut2=matrix(cells_2_cut2,nrow = 1)
    write.table(cells_2_cut2,"Cells_2_remove.csv",row.names = FALSE,col.names = FALSE)
    writeLines("cell's indices to remove are saved in 'cells 2 remove.csv',\n
               please rerun MAIN script again \n")
    stopQuietly_calista()
  }


  writeLines(' Press 1 if you want to perform additional analysis (e.g. lineage inference, cell ordering) , 0 otherwise: ')
  Proceed=readLines(n=1)
  if(as.integer(Proceed)!=1){
    stopQuietly_calista()
  }

  # %% *** 3-RECONSTRUCTION OF LINEAGE PROGRESSION ***
  # %
  # % Please check comments in 'CALISTA_transition_main' for more information.
  #Results=CALISTA_transition_main(DATA,Results)
  Results <<- cluster_distance_GUI(DATA,INPUTS,Results)
  # saveRDS(INPUTS,"INPUTS.rds")
  # saveRDS(DATA,"DATA.rds")
  # saveRDS(Results,"Results.rds")
  #
  #require(shiny)
  app <-shinyApp(
    ui <- fluidPage(
      useShinyjs(),
      # App title ----
      titlePanel("CALISTA"),
      tags$button(
        id = 'close',
        type = "button",
        class = "btn action-button",
        onclick = "setTimeout(function(){window.close();},500);",  # close browser
        "Close window"
      ),
      br(),
      br(),
      splitLayout(cellWidths = c("80%", "20%"),
                  plotlyOutput("p_time",height = "400px")
      ),

      br(),

      splitLayout(cellWidths = c("80%", "20%"),
                  plotlyOutput("p_cluster",height = "400px")
      ),

      br(),

      conditionalPanel(
        condition = "output$plot_p_transition == 'Yes'",
        splitLayout(cellWidths = c("80%", "20%"),
                    plotlyOutput("p_transition",height = "400px"),
                    verticalLayout(

                      wellPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 260px",
                                #h5('Cluster distance:'),
                                checkboxGroupInput("Clusterdistance", "Cluster distance:",
                                                   choiceNames =Results$TRANSITION$check_boxes,
                                                   choiceValues =1:length(Results$TRANSITION$check_boxes),
                                                   selected = 1:nrow(Results$TRANSITION$h_edge))

                      ),
                      wellPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 120px",
                                # checkboxInput("select", "Select all/Default", FALSE),
                                #verbatimTextOutput("value"),
                                actionButton("selectall", "Select all",width = "80px",style='height:30px;'),
                                actionButton("default", "Default",width = "65px",style='height:30px;'),
                                br(),
                                br(),

                                actionButton("ok", "OK",width = "150px",style='height:30px;color: #fff; background-color: #337ab7; border-color: #2e6da4'),
                                tags$style(type='text/css', "#run_report { width:50%; margin-left: 5px;}")
                      )
                    )
        )
      ),
      br(),
      splitLayout(cellWidths = c("80%", "20%"),
                  plotlyOutput("p_transition_genes",height = "400px"),
                  shinyjs::hidden(wellPanel(id = "transition_genes_panel", style = "overflow-y:scroll; max-height: 400px",
                                            radioButtons("transition_genes_radio_buttons", "Transition edges:",
                                                         choiceValues =   1,
                                                         choiceNames =  "Initializing"
                                            )
                  ))),

      br(),

      splitLayout(cellWidths = c("80%", "20%"),
                  plotlyOutput("p_ordering",height = "400px"),
                  shinyjs::hidden(wellPanel(id = "ordering_panel", style = "overflow-y:scroll; max-height: 400px",
                                            radioButtons("ordering_radio_buttons", "Color info:",
                                                         c("Clusters" = "button_cluster",
                                                           "Time/Cell stage" = "button_time",
                                                           "Pseudotime" = "button_pseudotime"
                                                         ))
                  ))),

      br()
    ),
    server <- function(input, output, session) {
      #source("./R/initialization.R")
      shinyjs::hide(id="ordering_panel")
      observe({
        if (input$close > 0) stopApp()                             # stop shiny
      })

      output$p_time <- renderPlotly({
        Results$p_time
      })


      output$p_cluster <- renderPlotly({
        Results$p_cluster
      })
      if (is.null(Results$TRANSITION[["p"]])) {
        output$plot_p_transition <- renderPrint({
          "No"
        })
      }else {

        output$plot_p_transition <- renderPrint({
          "Yes"
        })

        # Check first the select all checkbox
        observeEvent(input$selectall, {
          updateCheckboxGroupInput(session, "Clusterdistance",
                                   choiceNames =Results$TRANSITION$check_boxes,
                                   choiceValues =1:length(Results$TRANSITION$check_boxes),
                                   selected = 1:nrow(Results$TRANSITION$cluster_distance_selected)
          )
        })
        observeEvent(input$default, {
          updateCheckboxGroupInput(session, "Clusterdistance",
                                   choiceNames =Results$TRANSITION$check_boxes,
                                   choiceValues =1:length(Results$TRANSITION$check_boxes),
                                   selected = 1:nrow(Results$TRANSITION$h_edge)
          )
        })
        ################################################################################################

        observeEvent(input$ok, {
          if(length(unique(DATA$timeline))==1 && unique(DATA$timeline)==0){
            Results <<- reorder_clusters_GUI(DATA,INPUTS,Results)
            output$p_cluster <- renderPlotly({
              Results$p_cluster
            })
            output$p_transition <- renderPlotly({
              Results$TRANSITION$p
            })
          }



          # %% *** 5-PSEUDOTEMPORAL ORDERING OF CELLS ***
          # %
          # % Please check comments in 'CALISTA_ordering_main' for more information.
          #
          Results<<-CALISTA_ordering_main(DATA,Results)
          shinyjs::show(id="ordering_panel")

          # # %% *** 4-DETERMINATION OF TRANSITION GENES ***
          # %
          # % Please check comments in 'CALISTA_transition_genes_main' for more information.
          Results<<-CALISTA_transition_genes_main(DATA,INPUTS,Results)

          # Generate radio buttons labels
          Results$GENES$button_labels=paste(Results$TRANSITION$h_edge[,1],'-',Results$TRANSITION$h_edge[,2])

          updateRadioButtons(session, "transition_genes_radio_buttons",
                             choiceValues =   1: length(Results$GENES$button_labels),
                             choiceNames =  Results$GENES$button_labels,
                             selected = 1

          )
          output$p_ordering <- renderPlotly({
               p <- switch(input$ordering_radio_buttons,
                         button_cluster = Results$ORDERING$p_cluster,
                         button_time = Results$ORDERING$p_time,
                         button_pseudotime = Results$ORDERING$p_pseudotime
            )

            p
          })

          shinyjs::show(id="transition_genes_panel")
          output$p_transition_genes <- renderPlotly({
            edge_chosen <- as.numeric(input$transition_genes_radio_buttons)
            p <- plot_ly(
              x =  Results$GENES$final_transition_genes[[edge_chosen]],
              y = Results$GENES$sorted_prob[[edge_chosen]],
              name = "Transition Genes",
              type = "bar"
            ) %>%
              layout(
                title = "Transition genes",
                xaxis = list(title = "",
                             categoryorder = "array",
                             categoryarray = ~Results$GENES$final_transition_genes[[edge_chosen]]),
                yaxis = list(title = "logP")
              )

            p
          })


        })



        ################################################################################################

        output$p_transition <- renderPlotly({

          edges_chosen <- as.numeric(input$Clusterdistance)
          print(edges_chosen)
          nodes=Results$TRANSITION$cluster_distance_selected[edges_chosen,]
          # check if the graph is conected otherwise set to the default graph
          observeEvent(input$Clusterdistance, {
            if (is.vector(nodes) || length(unique(c(nodes[,1],nodes[,2])))<Results$expected_clusters ){
              print('hello')
              # Can also set the label and select items
              updateCheckboxGroupInput(session, "Clusterdistance",
                                       choiceNames =Results$TRANSITION$check_boxes,
                                       choiceValues =1:length(Results$TRANSITION$check_boxes),
                                       selected = 1:nrow(Results$TRANSITION$h_edge)
              )
              nodes=Results$TRANSITION$h_edge

            }
          })
          Results <<- update_transition_GUI(DATA,INPUTS,Results,nodes,edges_chosen)
          Results$TRANSITION$p

        })
        output$value <- renderText({ input$Clusterdistance })
      }

    }
  )

  runApp(app)

}
