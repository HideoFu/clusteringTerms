# Clustering Terms


library(shiny)
library(readr)
library(dplyr)
library(stringr)

library(tm)
library(dendextend)
library(slam)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Clustering Terms"),
   
   # Sidebar
   sidebarLayout(
      sidebarPanel(
        fileInput("file", "Select file:"),
        fileInput("sw", "Stop Words File:"),  # Stop words
        
        sliderInput("word_n", "Word number:", min = 0, max = 160, value = 40),
        sliderInput("freq_word_n", "Rows:", min = 0, max = 20, value = 5),
        textInput("kw", "Key word:"),
        
        actionButton("go", "Draw")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("dendro"),
         tableOutput("freq"),
         plotOutput("asoc")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Read Data
  raw_data <- reactive({
    data <- read_csv(input$file$datapath)
    if("clust" %in% names(data)){
      data <- data[,-which(names(data) == "clust")]
    }
    
    return(data)
  })
  
  # Make Clean Corpus
  clean_corp <- reactive({
    ## clean special chars
    data <- raw_data()
    
    data$Abstract <- gsub("\\r*\\n*", "", data$Abstract)
    data$PubDate <- data$PubDate %>% 
      str_replace_all("\\r\\n", "/") %>%
      str_replace_all("^/|/$", "")
    data$AuthorList <- gsub("\\r\\n", ", ", data$AuthorList)
    data$KeywordList <- data$KeywordList %>%
      str_replace_all("\\r\\n", ", ") %>%
      str_replace("^, ", "")
    
    
    ## rename / sort columns
    data <- data %>%
      mutate(doc_id = seq(1, nrow(data))) %>%
      select(doc_id, Abstract)
    
    colnames(data)[2] <- "text"
    
    
    ## convert to corpus
    source <- DataframeSource(data)
    corp <- VCorpus(source)
    
    
    
    ### Stemming function
    stemming <- function(text){
      char_vec <- unlist(strsplit(text, split = " "))
      stem_doc <- stemDocument(char_vec)
      new_text <- paste(stem_doc, collapse = " ")
      return(new_text)
    }
    
    ### Function for remove "[number]%", "[num]/[num]"
    rm_res_num <- function(text){
      # new_text <- gsub("([[:digit:]]|\\.)+%", "", text)
      # new_text <- gsub("([[:digit:]]|/)", "", new_text)
      new_text <- text %>% str_remove_all("([[:digit:]]|\\.)+%") %>% 
        str_remove_all("[[:digit:]]/[[:digit:]]")
      
      return(new_text)
    }
    
    
    clean_corp <- tm_map(corp, content_transformer(tolower)) %>%
      tm_map(content_transformer(rm_res_num)) %>%
      tm_map(removePunctuation, preserve_intra_word_dashes = T) %>%
      tm_map(stripWhitespace) %>%
      tm_map(content_transformer(stemming))  # must place after remove punctuations
    
    ## clean corpus
    
    if(!is.null(input$sw)){
      sw <- scan(input$sw$datapath, what = character())
      sw <- unique(stemDocument(sw))
      clean_corp <- tm_map(clean_corp, removeWords, sw)
      
      # sw <- paste0(" ", sw, " ")
      # sw <- paste(sw, collapse = "|")
      # rm_sw <- function(text){
      #   #new_text <- gsub(sw, "", text)
      #   new_text <- str_replace_all(text, sw, " ")
      #   return(new_text)
      # }
      # clean_corp <- tm_map(clean_corp, content_transformer(rm_sw))
    }
    
    return(clean_corp)
  })
  
  # Make TDM
  tdm <- reactive({
    TermDocumentMatrix(clean_corp(), control = list(weighting = weightBin, wordLengths = c(2, Inf)))
  })
  
  # Word Clustering
  ## key word
  stem_kw <- reactive({
    input$go
    
    isolate({stemDocument(input$kw)})
  })
  
  ## Reduce word number
  sub_tdm <- reactive({  
        sub_tdm <- tdm()
        sparse = 1
        
        while(sub_tdm$nrow > input$word_n){
          sparse = sparse - 0.005
          sub_tdm <- removeSparseTerms(tdm(), sparse = sparse)
        }
        
        
        ## Check if key word remains in terms
        if(input$kw != ""){
          df <- as.data.frame(as.matrix(sub_tdm))
          if(!(stem_kw() %in% rownames(df))){
            tdm <- as.data.frame(as.matrix(tdm()))
            df <- rbind(df, tdm[stem_kw(),])
            sub_tdm <- as.TermDocumentMatrix(df, weighting = weightBin, wordLengths = c(2, Inf))
          }
        }
        
        return(sub_tdm)
    })
      
  
  ## Plot dendrogram
  output$dendro <- renderPlot({
    if(is.null(input$file)){
      return(NULL)
    }
    
    input$go
    
    isolate({
      
    sub_dist_mat <- dist(as.matrix(sub_tdm()))
    sub_hc <- hclust(sub_dist_mat)
    sub_hcd <- as.dendrogram(sub_hc, hang = 0.1)
    
    if(input$kw != ""){
      sub_hcd <- branches_attr_by_labels(sub_hcd, stem_kw(), "red")
    }
    
    sub_hcd %>% set("labels_cex", 10 / sqrt(input$word_n)) %>%
    plot()
    })
  })
  
  ## Plot frequentry terms
  tdm_mat <- reactive({
    input$go
    
    isolate({
      tdm_vec <- row_sums(sub_tdm())
      names(tdm_vec) <- sub_tdm()$dimnames$Terms
      tdm_vec <- tdm_vec[order(tdm_vec, decreasing =T)] / dim(sub_tdm())[2]
      
      tdm_mat <- data.frame(Term = names(tdm_vec), Ratio = tdm_vec)
    })
  })
  
  output$freq <- renderTable({
    if(is.null(input$file)){
      return(NULL)
    }
    
    input$go
    
    isolate({
      head_tdm <- tdm_mat()[1:input$freq_word_n,]
      if(input$kw != "" & !(stem_kw() %in% head_tdm$Term)){
        head_tdm <- rbind(tdm_mat()[stem_kw(),], head_tdm)
      }
      
      head_tdm
    })
  })
  
  # Plot co occurence
  output$asoc <- renderPlot({
    if(is.null(input$file)){
      return(NULL)
    } 
    
    asoc <- findAssocs(tdm(), stem_kw(), 0.2)
    
    df <- as.data.frame(unlist(asoc))
    
    if(dim(df)[1] >0){
      colnames(df) <- "corr"
      rownames(df) <- str_remove(rownames(df), paste0(stem_kw(), "."))
      
      df$names <- rownames(df)
      df$id <- as.factor(seq(1:dim(df)[1]))
      
      
      p <- ggplot(df[1:15,], aes(x=id, y=corr)) + 
        geom_bar(stat="identity") + 
        scale_x_discrete(breaks=df$id,labels=df$names) +
        theme(plot.margin = unit(c(1,1,2,1), "cm")) +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(2.0)),
              axis.title.x = element_blank(),
              plot.title = element_text(size = rel(3.0))) +
        ggtitle(paste0(input$kw, " (",round(tdm_mat()[stem_kw(), "Ratio"],2), ")"))
      
      print(p)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

