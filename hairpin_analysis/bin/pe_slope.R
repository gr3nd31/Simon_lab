pe_slope <- function(df = datum){
  df$pe_slope <- 0
  df$pe_intecept <- 0
  for (i in c(1:nrow(df))){
    j <- df[i,]$PE
    j <- str_remove(string = j, pattern = "\\(")
    j <- str_remove(string = j, pattern = "\\)")
    j <- str_split(j, " ")[[1]]
    j <- j[2:length(j)]
    j <- as.numeric(j)
    interim <- tibble("Position" = c(1:length(j)),
                      "PE" = j)
    midpoint <- as.integer(length(j)/2)
    interim$relative_Position <- interim$Position
    interim[interim$Position > midpoint,]$relative_Position <- 1+abs(interim[interim$Position > midpoint,]$Position-max(interim$Position))
    pe_slope <- lm(PE~relative_Position, data = interim)
    df[i,]$pe_slope <- pe_slope$coefficients[2]
    df[i,]$pe_intecept <- pe_slope$coefficients[1]
  }
  return(df)
}