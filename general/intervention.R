library(ggpubr)
library(tidyverse)

incidence_color <- function(unprotected_incidence = 5.6,
                            protected_incidence = 3.9,
                            intervention_incidence = 1,
                            protection_rate = 0,
                            population_size = 1000000){
  for (i in c(0:100)){
    interim <- data.frame("intervention_percent" = i,
                          "infection_percent" = c(0:100))
    if (!exists("datum")){
      datum <- interim
    } else {
      datum <- rbind(datum, interim)
    }
  }
  datum$vaccinated <- (datum$intervention_percent/100)*population_size
  datum$unvaccinated <- population_size-datum$vaccinated
  datum$theoretically_infected <- (datum$infection_percent/100)*population_size
  # Number of vaccinated people who get infected
  datum$infected_vaccinated <- ((datum$infection_percent/100)*datum$vaccinated)*((100-protection_rate)/100)
  datum$infected_unvaccinated <- (datum$infection_percent/100)*datum$unvaccinated
  
  # Number of people getting hurt by vaccination
  datum$intervention_incidence <- datum$vaccinated*(intervention_incidence/100)
  # Number of vaccinated people getting hurt by infection
  datum$protected_incidence <- datum$infected_vaccinated*(protected_incidence/100)
  datum$unprotected_incidence <- datum$infected_unvaccinated*(unprotected_incidence/100)
  datum$total_incidence <- datum$intervention_incidence+datum$protected_incidence+datum$unprotected_incidence
  
  #print(head(datum))
  draft <- ggplot(data = datum, aes(x = intervention_percent, y = infection_percent, fill = total_incidence))+
    geom_tile()+
    scale_fill_viridis_b(option="inferno")+
    theme_bw()+
    labs(x = "Percent of Population Vaccinated",
         y = "Percent of Population Exposed to Infectious Dose",
         fill = "Total Incidence",
         title = paste0("Population: ", population_size, " (PR: ", protection_rate, "%, VacInd: ",intervention_incidence,"%, UnInc: ", unprotected_incidence,"%, VacInc: ",protected_incidence,"%)"))+
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 20))+
    geom_hline(yintercept = 10, color = "white", alpha = 0.5, linewidth = 2, linetype = 2)+
    geom_vline(xintercept = 20, color = "white", alpha = 0.5, linewidth = 2, linetype = 2)
  print(draft)
  write_csv(datum, "datum.csv")
}
incidence_color(protection_rate = 29)

