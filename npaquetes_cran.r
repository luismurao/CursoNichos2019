library(Ecdat)
library(dygraphs)
library(rvest)
library(magrittr)
library(dplyr)
library(zoo)
data("CRANpackages")
str(CRANpackages)



CRANpackages$Version <- as.character(CRANpackages$Version)
CRANpackages <- rbind(CRANpackages, 
                      data.frame(Version = "3.2", 
                                 Date = as.Date("2016-04-27"), 
                                 Packages = 8329, 
                                 Source = "Andrie de Vries"))







url <- "https://cran.r-project.org/web/packages/available_packages_by_date.html"

page <- read_html(url)
pkgs <- page %>%
  html_node("table") %>%
  html_table() %>%
  mutate(count = rev(1:nrow(.))) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(Month = format(Date, format="%Y-%m")) %>%
  group_by(Month) %>%
  summarise(published = min(count)) %>%
  mutate(Date = as.Date(as.yearmon(Month)))




pkgs %>%
  filter(Date > as.Date("2012-12-01")) %>%
  mutate(publishedGrowth = c(tail(.$published, -1), NA) / published) %>%
  mutate(counter = 1:nrow(.)) -> new_pkgs




tail(new_pkgs)


write.csv(pkgs,"Arichivos_01Intro/paquetesTS.csv",row.names = F)

pkgs_old <- CRANpackages[1:8,c(3,2)]
names(pkgs_old) <- names(pkgs[,c(2,3)])
n_pkgs <- rbind(pkgs_old,pkgs[,c(3,2)])
#cumsum(n_pkgs$published)

#diff(n_pkgs$published)

#data_pkgs <- CRANpackages[,2:3]
#pkgs <- data.frame(pkgs[,c(1,2)])
row.names(new_pkgs) <- as.Date(new_pkgs$Date)
dygraph(new_pkgs) %>% dyRangeSelector() 


