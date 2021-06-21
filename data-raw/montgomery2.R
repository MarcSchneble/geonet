library(spatstat)
library(lubridate)
library(dplyr)
library(readxl)
library(igraph)
library(Matrix)

path <- "~/NetworkSplines/"
source(file = paste0(path, "Functions.R"))

data <- read.csv(file = paste0(path, "Data/CarCrashesMontgomery_Incidents.csv")) %>% as_tibble() %>%
  mutate(lat = as.numeric(substring(sub("\\,.*", "", Location), first = 2)),
         lon = as.numeric(gsub("[)]" ,"" , sub("^\\S+", "", Location))),
         date = as.POSIXct(substring(Crash.Date.Time, 1, 19), format = "%m/%d/%Y %H:%M:%OS", tz = "America/New_York"),
         date = date + hours(12)*(substring(Crash.Date.Time, 21, 22) == "PM")) %>%
  filter(year(date) != 2020,
         Route.Type %in% c("Maryland (State)", "Interstate (State)", "US (State)"),
         hour(date) >= 6, hour(date) <= 21)

V <- read_excel(paste0(path, "Data/VerticesMontgomery.xlsx"))%>% filter(!is.na(lat), is.na(remove)) %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>% arrange(Id)
E <- read_excel(paste0(path, "Data/EdgesMontgomery.xlsx")) %>% filter(!is.na(from))

E$from.lat <- V$lat[match(E$from, V$Id)]
E$from.lon <- V$lon[match(E$from, V$Id)]
E$to.lat <- V$lat[match(E$to, V$Id)]
E$to.lon <- V$lon[match(E$to, V$Id)]

min.lon <- min(V$lon)
min.lat <- min(V$lat)
V$lon <- V$lon - min.lon
V$lat <- V$lat - min.lat

P <- ppp(x = V$lon, y = V$lat,
         window = owin(xrange = c(min(V$lon), max(V$lon)), yrange = c(min(V$lat), max(V$lat))))
L <- linnet(vertices = P, edges = as.matrix(E[, 2:3]))
s <- L$dpath[191, 193]/2.2
L <- spatstat.geom::rescale(L, s = s, unitname = c("kilometer", "kilometers"))

data <- data %>%
  mutate(lon.net = (lon - min.lon)/s, lat.net = (lat - min.lat)/s)

L.ppp <- ppp(x = data$lon.net, y = data$lat.net, L$window, marks = data$Report.Number)
L.psp <- as.psp(L)
projection <- project2segment(L.ppp, L.psp)

retain <- which(projection$d < 0.1)
seg <- projection$mapXY[retain]
tp <- projection$tp[retain]
marks <- projection$Xproj$marks[retain]
L.lpp <- as.lpp(seg = seg, tp = tp, L = L, marks = marks)

W <- owin(xrange = c(15, 33), yrange = c(0, 12), unitname = c("kilometer", "kilometers"))
L.old <- L
L <- L[W, snip = FALSE]
L$ind.edges <-  which((E$from.lat - min.lat)/s <= 12 & (E$from.lat - min.lat)/s >= 0 &
                        (E$to.lat - min.lat)/s <= 12 & (E$to.lat - min.lat)/s >= 0 &
                        (E$from.lon - min.lon)/s <= 33 & (E$from.lon - min.lon)/s >= 15 &
                        (E$to.lon - min.lon)/s <= 33 & (E$to.lon - min.lon)/s >= 15)
X <- as.ppp(L.lpp[W, snip = FALSE])
data <- data %>% mutate(on = factor(1*(Report.Number %in% X$marks)))
plot(L, box = TRUE)

delta <- 0.05
h <- 0.025
r <- 2
L <- augment.linnet(L, delta, h, r)
L.lpp <- lpp(X, L)
plot(unmark(L.lpp))

L.lpp <- unmark(L.lpp)
L.lpp$data$hour <- hour(data %>% filter(on == 1) %>% pull(date)) + floor(minute(data %>% filter(on == 1) %>% pull(date))/15)/4

internal.covariates <- tibble(id = 1:length(L$ind.edges),
                              type = factor(E$type[L$ind.edges], levels = c("state", "interstate", "US")),
                              direction = factor(E$direction[L$ind.edges], levels = c("SN", "EW", "SENW", "SWNE")))

montgomery <- geonet::as_gnpp(L.lpp)
montgomery$network$lins <- left_join(montgomery$network$lins, internal.covariates, by = "id")

plot(montgomery)

usethis::use_data(montgomery, overwrite = TRUE)

