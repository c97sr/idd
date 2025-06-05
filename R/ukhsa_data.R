if (FALSE) {
  
  # Install the 'httr' package if not already installed
  if (!requireNamespace("httr", quietly = TRUE)) {
    install.packages("httr")
  }
  
  library(httr)
  library(jsonlite)
  library(jsonlite)
  library(ggplot2)
  
  # Define the API endpoint
  # Next thing is to refactor below
  url <- "https://api.ukhsa-dashboard.data.gov.uk/themes/infectious_disease/sub_themes/vaccine_preventable/topics/Measles/geography_types/Nation/geographies/England/metrics/measles_cases_casesByOnsetWeek?page_size=365&age=all"
  
  # Make the GET request to the API with the API key
  response <- GET(url, add_headers(Authorization = paste("Bearer")))
  
  # Check the status code of the response
  if (status_code(response) == 200) {
    # If the request was successful, parse the JSON content
    content <- content(response, "text")
    json_data <- fromJSON(content)
    
    # Print the parsed JSON data
    print(json_data)
  } else {
    # If the request was not successful, print the status code
    print(paste("Error:", status_code(response)))
  }
  
  
  df <- json_data$results
  
  # Ensure the date column is in the correct date format
  df$date <- as.Date(df$date, format = "%Y-%m-%d")
  
  # Plot the data
  ggplot(data = df, aes(x = date, y = metric_value)) +
    geom_line() +
    labs(title = "Metric Value Over Time", x = "Date", y = "metric_value") +
    theme_minimal()
  
  
}