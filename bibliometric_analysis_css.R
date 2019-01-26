##---------------------------------------------------------------------------------------------##
##               GENDER & AUTHORSHIP IN COMPUTATIONAL SOCIAL SCIENCE PUBLICATIONS              ##
##---------------------------------------------------------------------------------------------##


## R version 3.4.3 (2017-11-30)

## Author: Lisa Hehnke (dataplanes.org | @DataPlanes)


#-------#
# Setup #
#-------#

# Set working directory
setwd("~/Documents/GitHub/bibliometric-analysis")

# Install and load pacman if not already installed
if (!require("pacman")) install.packages("pacman")
library(pacman)

# Load packages
p_load(bibliometrix, genderizeR, ggthemes, igraph, reshape2, splitstackshape, tidyverse)


#-------------#
# Import data #
#-------------#

# Import bibtex citations
bib <- readFiles("CSS_WoS.bib")

# Convert to data frame
bib_df <- convert2df(bib, dbsource = "isi", format = "bibtex")


#-----------------------------#
# Extract authors' full names #
#-----------------------------#

# Modify isibib2df() function to obtain authors' full names by replacing AU with AF
isibib2df_mod <- function (D) 
{
  D <- gsub("\\@", "Manuscript =", D)
  Papers <- which(regexpr("Manuscript =", D) == 1)
  Papers <- c(Papers, length(D))
  nP <- length(Papers) - 1
  uniqueTag <- c("AF", "TI", "SO", "JI", "VO", "NU", "PP", 
                 "BO", "DT", "DT2", "DE", "ID", "AB", "C1", "CR", "TC", 
                 "PY", "SC", "UT", "DI", "RP")
  DATA <- data.frame(matrix(NA, nP, length(uniqueTag)))
  names(DATA) <- uniqueTag
  Tag <- c("Author =", "Title =", "Journal =", "Journal-ISO =", 
           "Volume =", "Number =", "Pages =", "Booktitle =", "Manuscript =", 
           "Type =", "Keywords =", "Keywords-Plus =", "Abstract =", 
           "Affiliation =", "Cited-References =", "Times-Cited =", 
           "Year =", "Web-of-Science-Categories  =", "Unique-ID =", 
           "DOI =")
  for (i in 1:nP) {
    iP <- Papers[i]
    iPs <- Papers[i + 1] - 1
    if (i%%100 == 0 | i == nP) 
      cat("Articles extracted  ", i, "\n")
    iPiPs <- seq(iP, iPs)
    for (j in 1:length(Tag)) {
      POS <- which(regexpr(Tag[j], D[iPiPs]) == 1) + iP - 
        1
      if (length(POS) == 1) {
        Seq <- seq(POS, iPs)
        END <- which(regexpr(".*\\}", D[Seq]) == 1)[1]
        POSEND <- seq(POS, (POS + END - 1))
        if (uniqueTag[j] != "C1") {
          DATA[[uniqueTag[j]]][i] <- paste0(D[POSEND], 
                                            collapse = "")
        }
        else {
          DATA[[uniqueTag[j]]][i] <- paste0(gsub(";", 
                                                 ",", D[POSEND]), collapse = ";")
        }
        if (uniqueTag[j] == "DI") {
          DOI <- gsub("DOI = ", "", D[POS])
          DATA[[uniqueTag[j]]][i] <- gsub(",", "", DOI)
        }
      }
    }
  }
  DATA <- as.data.frame(apply(DATA, 2, function(d) gsub("\\{", 
                                                        "", d)), stringsAsFactors = FALSE)
  DATA <- as.data.frame(apply(DATA, 2, function(d) gsub("\\},", 
                                                        "", d)), stringsAsFactors = FALSE)
  DATA <- as.data.frame(apply(DATA, 2, function(d) gsub("\\}", 
                                                        "", d)), stringsAsFactors = FALSE)
  for (i in 1:length(Tag)) {
    DATA[[uniqueTag[i]]] <- gsub(paste0(Tag[i], " ", collapse = ""), 
                                 "", DATA[[uniqueTag[i]]], fixed = TRUE)
  }
  DATA$DT <- gsub("Manuscript =", "", unlist(lapply(strsplit(DATA$DT, 
                                                             ","), function(l) l[1])))
  DATA$DT <- gsub("ISI:.*", "", DATA$DT)
  listAF <- strsplit(DATA$AF, " and ")
  AF <- lapply(listAF, function(l) {
    lastname <- trim(gsub(",.*", ",", l)) # add comma as separator
    firstname <- trim(gsub(".*,", "", l))
    AF <- paste(lastname, unlist(firstname), sep = " ", collapse = ";")
    return(AF)
  })
  DATA$AF <- unlist(AF)
  DATA$TC <- as.numeric(sub("\\D*(\\d+).*", "\\1", DATA$TC))
  DATA$PY <- as.numeric(sub("\\D*(\\d+).*", "\\1", DATA$PY))
  DATA$UT <- gsub(":", "", DATA$UT, fixed = TRUE)
  DATA$RP <- unlist(lapply(strsplit(DATA$C1, "\\."), function(l) l[1]))
  DATA <- data.frame(lapply(DATA, toupper), stringsAsFactors = FALSE)
  DATA$ID <- gsub("   ", ";", DATA$ID)
  DATA$DE <- gsub("   ", ";", DATA$DE)
  DATA$DB = "ISI"
  ind <- which(is.na(DATA$SO))
  DATA$SO[ind] <- DATA$BO[ind]
  DATA <- DATA[, -(which(names(DATA) == "BO"))]
  DATA$NU <- as.numeric(gsub("[^0-9]", "", DATA$NU))
  return(DATA)
}

# Get full names
authors_full <- isibib2df_mod(bib)[, 1]

# Replace AU with AF
bib_df$AU <- authors_full


#-----------------------#
# Bibliometric analysis #
#-----------------------#

# Descriptive analysis of bibliographic df
results <- biblioAnalysis(bib_df, sep = ";")

# Summarize and plot results
results_sum <- summary(results, 10, pause = FALSE) # top 10 results for each metric
plot(results, 10, pause = FALSE)


#--------------------------#
# Genderize authors' names #
#--------------------------#

# Extract authors' names and split strings
authors_total <- data.frame(results$Authors)
authors_total <- select(authors_total, "AU")
authors_split <- colsplit(authors_total$AU, ",", c("lastname", "firstname"))

# Predict sex and genderize the original character vector (via genderizeR)
givenNames <- findGivenNames(authors_split$firstname, progress = FALSE)
authors_gender <- genderize(authors_split$firstname, genderDB = givenNames, progress = FALSE)

# Merge genederized data with authors' full names
authors_master <- cbind(authors_split, authors_gender[, -c("text")])

# Concatenate authors' names
authors_master$AU <- paste(authors_master$lastname, authors_master$firstname, sep = ",")
authors_gen <- authors_master[, c("AU", "gender")]

# Export genederized data
saveRDS(authors_gen, "authors_genderized.rds")

# Import data
authors_gen <- readRDS("authors_genderized.rds")


#------------------------#
# Recode genderized data #
#------------------------#

# Recode NA
authors_gen$gender[is.na(authors_gen$gender)] <- "unknown"

# Manually recode gender
authors_gen$gender_rec <- authors_gen$gender

authors_gen$gender_rec[authors_gen$AU == "HELBING, D."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SERESINHE, CHANUKI ILLUSHKA"] <- "female"
authors_gen$gender_rec[authors_gen$AU == "SHAH, DHAVAN V."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SQUAZZONI, FLAMINIO"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "TONG, HANGHANG"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "AGGARWAL, CHARU"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "ARENAS, A."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BAINBRIDGE, WS"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BALAN, G"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BALIETTI, S."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BANKES, S"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BAR-YAM, YANEER"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BARONCHELLI, ANDREA"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BARTUMEUS, F."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "BONELLI, G."] <- "female"
authors_gen$gender_rec[authors_gen$AU == "BORGE-HOLTHOEFER, J."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CAO, YANPENG"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CARAVIELLO, MICHELE"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CENTELLEGHER, SIMONE"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CHEN, KWANG-CHENG"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CHEN, SICONG"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CHENG, YUNSHENG"] <- "male"
# CHENG, ZILONG ?
authors_gen$gender_rec[authors_gen$AU == "CHUAH, CHEN-NEE"] <- "female"
authors_gen$gender_rec[authors_gen$AU == "CIOFFI-REVILLA, C"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "COELHO, H"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "CONTE, R."] <- "female"
# DAVID, N ?
authors_gen$gender_rec[authors_gen$AU == "DE HAAN, FJALAR"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "DEFFUANT, G."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "DIAZ-GUILERA, A."] <- "male"
# DONG, LEI ?
authors_gen$gender_rec[authors_gen$AU == "DONG, YUXIAO"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "FLACHE, A."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "FUJITA, HAMIDO"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "GENG, TIEMING"] <- "male" # ?
authors_gen$gender_rec[authors_gen$AU == "GILBERT, N."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "GOMEZ, VICENC"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "GUTIERREZ-ROIG, M."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "HAILEGIORGIS, ATESMACHEW"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "HALL, MARGERET"] <- "female"
authors_gen$gender_rec[authors_gen$AU == "HSIAO, JEN-HAO"] <- "unknown"
authors_gen$gender_rec[authors_gen$AU == "HU, YUENING"] <- "female"
authors_gen$gender_rec[authors_gen$AU == "KERTESZ, J."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "KIM, HWALBIN"] <- "male"
# KIM, KWANSOO ?
authors_gen$gender_rec[authors_gen$AU == "KUZNAR, LA"] <- "unknown" # Lawrence Kuznar?
authors_gen$gender_rec[authors_gen$AU == "LEMPERT, R"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "LIN, CHING-YUNG"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "LIU, HUAN"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "LORETO, V."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "LUKE, S"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "MARTINEZ-TORRES, M. R."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "MASON, WINTER."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "MOAT, S."] <- "female"
authors_gen$gender_rec[authors_gen$AU == "MORENO, Y."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "NADAL, J. -P."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "NOWAK, A."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "OLMEDILLA, M."] <- "female"
# OLTRA, A. ?
authors_gen$gender_rec[authors_gen$AU == "PALMER, J. R. B."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "PANAIT, L"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "PARMAR, NIKI JITENDRA"] <- "unknown"
authors_gen$gender_rec[authors_gen$AU == "PERELLO, J."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "PIEDRAHITA, P."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "POPPER, S"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SAGARRA, O."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SAGARRA, OLEGUER"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SALLACH, DL"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SAN MIGUEL, M."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SANCHEZ, A."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SAUNDERS-NEWTON, D"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SHI, LEI"] <- "unknown" # male?
authors_gen$gender_rec[authors_gen$AU == "SONG, HYUNJIN"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "SUDHAHAR, SAATVIGA"] <- "female"
# SULLIVAN, K?
authors_gen$gender_rec[authors_gen$AU == "SUNDSOY, PAL ROE"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "TANG, JIE"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "TANG, JIJUN"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "TORAL, S. L."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "TUMASJAN, ANDRANIK"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "VARSHNEY, LAV R."] <- "male"
authors_gen$gender_rec[authors_gen$AU == "VESCOVI, MICHELE"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "WU YOUYOU WU YOUYOU,"] <- "female" # WU YOUYOU
authors_gen$gender_rec[authors_gen$AU == "WU, ZHENGWEI"] <- "female"
authors_gen$gender_rec[authors_gen$AU == "XIA, RUOFAN"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "XU, ZESHUI"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "YOO, YOUNGJIN"] <- "unknown" # not sure if male or female (I'd say male...)
authors_gen$gender_rec[authors_gen$AU == "ZHANG, QINGPENG"] <- "male"
authors_gen$gender_rec[authors_gen$AU == "ZHANG, YULIN"] <- "male" # ?


#----------------------------------#
# Genderize first and last authors #
#----------------------------------#

# Extract first authors' names and merge with genderized df
authors_first <- data.frame(results$FirstAuthors)
colnames(authors_first) <- "AU"
authors_first_gen <- merge(authors_first, authors_gen, key = "AU")

# Extract last authors' names
authors_raw <- data.frame(bib_df$AU)
authors_last <- data.frame(sub(".*;", "", authors_raw$bib_df.AU)) # remove all characters up to last ;
colnames(authors_last) <- "AU"
authors_last_gen <- merge(authors_last, authors_gen, key = "AU")


#--------------------------#
# Theme for visualizations #
#--------------------------#

# Set theme for visualizations
viz_theme <- theme(
  strip.background = element_rect(colour = "transparent", fill = "grey90"),
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  legend.key = element_rect(fill = "white"),
  strip.text = element_text(size = rel(1), face = "bold"),
  plot.caption = element_text(colour = "grey50"),
  text = element_text(family = "Avenir"))


#-------------------#
# Descriptive stats #
#-------------------#

# Plot number of authors by gender
authors_gen %>%  
  count(gender_rec, sort = TRUE) %>%
  ungroup() %>%
  mutate(gender_rec = reorder(gender_rec, n)) %>%
  ggplot(aes(gender_rec, n, label = n)) +
  #geom_col(fill = c("#00BFC4", "#F8766D", "grey"), alpha = 0.9) + # instead of geom_segment & geom_point
  geom_segment(aes(x = gender_rec, 
                   xend = gender_rec, 
                   y = 0, 
                   yend = n), col = "grey50", alpha = 0.9) +
  geom_point(col = c("#00BFC4", "#F8766D", "grey"), size = 4, alpha = 0.9) + 
  geom_text(hjust = -0.5, size = 4) +
  theme(text = element_text(size = 20, color = "#1f232e")) + 
  xlab("") + ylab("") + ggtitle("Number of authors by gender", subtitle = " ") +
  ylim(0, 500) + coord_flip() + viz_theme

ggsave("authors_gender.png", width = 12, height = 8, units = "in", dpi = 100)

# Plot number of first authors by gender
authors_first_gen %>%  
  count(gender_rec, sort = TRUE) %>%
  ungroup() %>%
  mutate(gender_rec = reorder(gender_rec, n)) %>%
  ggplot(aes(gender_rec, n, label = n)) +
  #geom_col(fill = c("#00BFC4", "#F8766D", "grey"), alpha = 0.9) + # instead of geom_segment & geom_point
  geom_segment(aes(x = gender_rec, 
                   xend = gender_rec, 
                   y = 0, 
                   yend = n), col = "grey50", alpha = 0.9) +
  geom_point(col = c("#00BFC4", "#F8766D", "grey"), size = 4, alpha = 0.9) + 
  geom_text(hjust = -0.5, size = 4) +
  theme(text = element_text(size = 20, color = "#1f232e")) + 
  xlab("") + ylab("") + ggtitle("Number of first authors by gender", subtitle = " ") +
  ylim(0, 200) + coord_flip() + viz_theme

ggsave("firstauthor_gender.png", width = 12, height = 8, units = "in", dpi = 100)

# Plot number of last authors by gender
authors_last_gen %>%  
  count(gender_rec, sort = TRUE) %>%
  ungroup() %>%
  mutate(gender_rec = reorder(gender_rec, n)) %>%
  ggplot(aes(gender_rec, n, label = n)) +
  #geom_col(fill = c("#00BFC4", "#F8766D", "grey"), alpha = 0.9) + # instead of geom_segment & geom_point
  geom_segment(aes(x = gender_rec, 
                   xend = gender_rec, 
                   y = 0, 
                   yend = n), col = "grey50", alpha = 0.9) +
  geom_point(col = c("#00BFC4", "#F8766D", "grey"), size = 4, alpha = 0.9) +  
  geom_text(hjust = -0.5, size = 4) +
  theme(text = element_text(size = 20, color = "#1f232e")) + 
  xlab("") + ylab("") + ggtitle("Number of last authors by gender", subtitle = " ") +
  ylim(0, 200) + coord_flip() + viz_theme

ggsave("lastauthor_gender.png", width = 12, height = 8, units = "in", dpi = 100)


#------------------------#
# Publications over time #
#------------------------#

# Count publications by year
pub_ts <- bib_df %>% 
  dplyr::count(PY)

# Plot timeline
ggplot(pub_ts, aes(PY, n)) + geom_point(col = "grey") + 
  geom_line(col = "grey", size = 1) +
  #facet_wrap(~channelTitle, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Count", title = "Number of publications over time", subtitle = " ") +
  scale_x_continuous(breaks = seq(1998, 2018)) + 
  theme(text = element_text(size = 20)) + 
  viz_theme + ylim(0, 50) + theme(axis.text.x = element_text(angle = 65, vjust = 0.5))

ggsave("citations_timeline.png", width = 12, height = 8, units = "in", dpi = 100)

# Count publications by year and gender of first author
bib_df$FA <- gsub(";.*", "", bib_df$AU)
authors_first_merge <- unique(authors_first_gen)
colnames(authors_first_merge) <- c("FA", "FA_gen", "FA_gen_rec")
bib_df <- merge(bib_df, authors_first_merge, by = "FA")

pub_gen_ts <- bib_df %>% 
  group_by(FA_gen_rec) %>%
  dplyr::count(PY)

# Plot timeline
ggplot(pub_gen_ts, aes(PY, n, color = FA_gen_rec)) + geom_point() + 
  geom_line(size = 1) +
  scale_color_manual(values = c(female = "#F8766D", male = "#00BFC4", unknown = "grey"), name = "Gender") +
  #facet_wrap(~FA_gen, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Count", title = "Number of publications over time by first author's gender", subtitle = " ") +
  scale_x_continuous(breaks = seq(1998, 2018)) + 
  theme(text = element_text(size = 20)) + 
  viz_theme + ylim(0, 50) + theme(axis.text.x = element_text(angle = 65, vjust = 0.5))

ggsave("publications_timeline_firstauthor.png", width = 12, height = 8, units = "in", dpi = 100)

# Count publications by year and gender of last author
bib_df$LA <- gsub(".*;", "", bib_df$AU)
authors_last_merge <- unique(authors_last_gen)
colnames(authors_last_merge) <- c("LA", "LA_gen", "LA_gen_rec")
bib_df <- merge(bib_df, authors_last_merge, by = "LA")

pub_gen_ts <- bib_df %>% 
  group_by(LA_gen_rec) %>%
  dplyr::count(PY)

# Plot timeline
ggplot(pub_gen_ts, aes(PY, n, color = LA_gen_rec)) + geom_point() + 
  geom_line(size = 1) +
  scale_color_manual(values = c(female = "#F8766D", male = "#00BFC4", unknown = "grey"), name = "Gender") + # , guide = "none"
  #facet_wrap(~LA_gen, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Count", title = "Number of publications over time by last author's gender", subtitle = " ") +
  scale_x_continuous(breaks = seq(1998, 2018)) + 
  theme(text = element_text(size = 20)) + 
  viz_theme + ylim(0, 50) + theme(axis.text.x = element_text(angle = 65, vjust = 0.5))

ggsave("publications_timeline_lastauthor.png", width = 12, height = 8, units = "in", dpi = 100)


#---------------------#
# Citations over time #
#---------------------#

# Count citations by year
cit_ts <- bib_df %>% 
  group_by(PY) %>%
  summarise(Freq = sum(TC))

# Plot timeline
ggplot(cit_ts, aes(PY, Freq)) + geom_point(col = "grey") + 
  geom_line(col = "grey", size = 1) +
  #facet_wrap(~channelTitle, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Count", title = "Number of citations over time", subtitle = " ") +
  scale_x_continuous(breaks = seq(1998, 2018)) + 
  theme(text = element_text(size = 20)) + 
  viz_theme + ylim(0, 1000) + theme(axis.text.x = element_text(angle = 65, vjust = 0.5))

ggsave("citations_timeline.png", width = 12, height = 8, units = "in", dpi = 100)

# Count citations by year and gender of first author
cit_gen_ts <- bib_df %>% 
  group_by(FA_gen_rec, PY) %>%
  summarise(Freq = sum(TC))

# Plot timeline
ggplot(cit_gen_ts, aes(PY, Freq, color = FA_gen_rec)) + geom_point() + 
  geom_line(size = 1) +
  scale_color_manual(values = c(female = "#F8766D", male = "#00BFC4", unknown = "grey"), name = "Gender") +
  #facet_wrap(~FA_gen, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Count", title = "Number of citations over time by first author's gender", subtitle = " ") +
  scale_x_continuous(breaks = seq(1998, 2018)) + 
  theme(text = element_text(size = 20)) + 
  viz_theme + ylim(0, 1000) + theme(axis.text.x = element_text(angle = 65, vjust = 0.5))

ggsave("citations_timeline_firstauthor.png", width = 12, height = 8, units = "in", dpi = 100)

# Count citations by year and gender of last author
cit_gen_ts <- bib_df %>% 
  group_by(LA_gen_rec, PY) %>%
  summarise(Freq = sum(TC))

# Plot timeline
ggplot(cit_gen_ts, aes(PY, Freq, color = LA_gen_rec)) + geom_point() + 
  geom_line(size = 1) +
  scale_color_manual(values = c(female = "#F8766D", male = "#00BFC4", unknown = "grey"), name = "Gender") + # , guide = "none"
  #facet_wrap(~LA_gen, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Count", title = "Number of citations over time by last author's gender", subtitle = " ") +
  scale_x_continuous(breaks = seq(1998, 2018)) + 
  theme(text = element_text(size = 20)) + 
  viz_theme + ylim(0, 1000) + theme(axis.text.x = element_text(angle = 65, vjust = 0.5))

ggsave("citations_timeline_lastauthor.png", width = 12, height = 8, units = "in", dpi = 100)


#-----------------------#
# Bibliometric networks #
#-----------------------#

## (1) Co-citation network

# Create co-citation network
bib_net_mat_cit <- biblioNetwork(bib_df, analysis = "co-citation", network = "references", sep = ".  ")

# Network statistics
bib_net_stats_cit <- networkStat(bib_net_mat_cit)

# Plot co-citation network
bib_net_cit <- networkPlot(bib_net_mat_cit, n = bib_net_stats_cit$network$networkSize, Title = "Co-citation network", type = "fruchterman", size.cex = T, remove.multiple = F, labelsize = 0.7, edgesize = 5)


## (2) Collaboration network

# Create collaboration network
bib_net_mat_collab <- biblioNetwork(bib_df, analysis = "collaboration", network = "authors", sep = ";")

# Network statistics
bib_net_stats_collab <- networkStat(bib_net_mat_collab)

# Plot collaboration network
bib_net_collab <- networkPlot(bib_net_mat_collab, n = bib_net_stats_collab$network$networkSize, Title = "Collaboration network", type = "fruchterman", size.cex = T, remove.multiple = F, labelsize = 0.7, edgesize = 5)


#------------------------------------#
# Visualize collab network by gender #
#------------------------------------#

# Transform matrix to undirected graph object
g <- graph_from_adjacency_matrix(bib_net_mat_collab, mode = "undirected")

# Remove loops
g <- igraph::simplify(g, remove.loops = TRUE, remove.multiple = FALSE)

# Descriptive stats
V(g)
E(g)

# Specify layouts
layout <- layout_nicely(g)
layout_kk <- layout_with_kk(g)
layout_mds <- layout_with_mds(g)

# Assign color of vertices by gender
authors_net <- data.frame(V(g)$name)
colnames(authors_net) <- "AU"
library(plyr)
authors_gen_net <- join(authors_net, authors_gen, by = "AU") # additional/missing observation?
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")

V(g)$gender <- authors_gen_net$gender_rec
V(g)$color <- ifelse(V(g)$gender %in% c("male", "female"), "#00BFC4", "grey")
V(g)$color <- ifelse(V(g)$gender %in% c("female"), "#F8766D", V(g)$color)

# Compute degree centrality
V(g)$degree <- degree(g, mode = "total")/2

# Plot network (degree centralities as vertex sizes for aesthetical reasons)
## Layout 1 (nicely)
plot(g, edge.arrow.size = 0.2, vertex.color = V(g)$color, vertex.label = NA, vertex.frame.color = "slategray", layout = layout, 
     vertex.size = (degree(g, mode = "total")/2))

## Layout 2 (kk)
plot(g, edge.arrow.size = 0.2, vertex.color = V(g)$color, vertex.frame.color = "slategray", layout = layout_kk, 
     vertex.size = (degree(g, mode = "total")/2), vertex.label.cex = 2.2 * V(g)$degree/max(V(g)$degree) + 0.2, 
     vertex.label.family = "sans", vertex.label.color = "black", edge.color = "lightgray")

plot(g, edge.arrow.size = 0.2, vertex.color = V(g)$color, vertex.label = NA, vertex.frame.color = "slategray", layout = layout_kk, 
     vertex.size = (degree(g, mode = "total")/2), edge.color = "lightgray")

# Layout 3 (mds)
plot(g, edge.arrow.size = 0.2, vertex.color = V(g)$color, vertex.label = NA, vertex.frame.color = "slategray", layout = layout_mds, 
     vertex.size = (degree(g, mode = "total")/2), edge.color = "lightgray")

# Save plot
png("collaboration_network.png", width = 8, height = 8, units = 'in', res = 500)
plot(g, edge.arrow.size = 0.2, vertex.color = V(g)$color, vertex.label = NA, vertex.frame.color = "slategray", layout = layout, 
     vertex.size = (degree(g, mode = "total")/2))
dev.off()


#-----------------#
# Inspect network #
#-----------------#

# Extract edges and split string
#coop_gen <- data.frame(as_ids(E(g)[inc(V(g)[gender == "female"])]))
coop_gen <- data.frame(as_ids(E(g)))
colnames(coop_gen) <- "AUs"
coop_gen <- colsplit(coop_gen$AUs, "\\|", c("AU_1", "AU_2"))

# Add data on gender
library(plyr)
coop_gen <- merge(coop_gen, authors_gen, by.x = c("AU_1"), by.y = c("AU"))
colnames(coop_gen) <- c("AU_1", "AU_2", "gender", "gender_1")
coop_gen <- merge(coop_gen, authors_gen, by.x = c("AU_2"), by.y = c("AU"))
coop_gen <- coop_gen[, c("AU_1", "AU_2", "gender_1", "gender_rec")]
colnames(coop_gen) <- c("AU_1", "AU_2", "gender_1", "gender_2")
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")

# Code type of cooperation (male only, female only, mixed, unknown)
coop_gen$type <- ifelse(coop_gen$gender_1 == "male" & coop_gen$gender_2 == "male", "male only", 
                        ifelse(coop_gen$gender_1 == "female" & coop_gen$gender_2 == "female", "female only", 
                               ifelse(coop_gen$gender_1 == "unknown" | coop_gen$gender_2 == "unknown", "unknown", 
                                      "mixed")))

# Plot number of collaborations by type
coop_gen %>%    
  count(type, sort = TRUE) %>%
  ungroup() %>%
  mutate(type = reorder(type, n)) %>%
  ggplot(aes(type, n, label = n)) +
  #geom_col(fill = c("#00BFC4", "black", "grey", "#F8766D"), alpha = 0.9) + # instead of geom_segment & geom_point
  geom_segment(aes(x = type, 
                   xend = type, 
                   y = 0, 
                   yend = n), col = "grey50", alpha = 0.9) +
  geom_point(col = c("#00BFC4", "black", "grey", "#F8766D"), size = 4, alpha = 0.9) + 
  geom_text(hjust = -0.5, size = 4) +
  theme(text = element_text(size = 20, color = "#1f232e")) + 
  xlab("") + ylab("") + ggtitle("Number of collaborations by type", subtitle = " ") +
  ylim(0, 800) + coord_flip() + viz_theme

ggsave("collaborations_type.png", width = 12, height = 8, units = "in", dpi = 100)


