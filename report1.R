library(ggplot2)
library(gt)
library(reshape)
library(dplyr)

sample1 <- '1-407-3-20-21.pnBase.txt'
UKVOC <- "20I/501Y.V1"
UKpos <- factor(c(8782,14408,14676,15279,23063,23403,28881,28882))

mydata1 <- read.table('1-407-3-20-21.pnBase.txt', sep="\t", header=TRUE)

#remove insertion type column
mydata1 <- mydata1[,-10]
Coverage <- rowSums(mydata1[,4:9])
mydata1 <- cbind(mydata1,Coverage)

mydata2 <- mydata1[mydata1$Pos %in% UKpos,-1]
mydata2 %>% gt()

cleanCoverage <- mydata2[,-10]
mydataMelt <- melt(cleanCoverage, id=c("Pos","Nucl"))

TopVariants <- mydataMelt %>% 
  group_by(Pos) %>% 
  mutate(max_score = max(value)) %>% 
  ungroup() %>% 
  filter(value == max_score)

TopVariants <- TopVariants[order(TopVariants$Pos),]

#make a table of top variants for report
TopVariants %>% gt()

#make a plotting theme
mytheme <- theme(
  text=element_text(family="Avenir Next"),
  axis.text.x = element_text(angle = 90,size = 5),
  axis.text.y = element_text(size = 5),
  panel.background = element_rect(fill = NA),
  panel.border = element_rect(fill = NA, color = "grey75"),
  axis.ticks = element_line(color = "grey85"),
  panel.grid.major = element_line(color = "grey95", size = 0.2),
  panel.grid.minor = element_line(color = "grey95", size = 0.2)
)

#plot coverage
p1 <-ggplot(data=mydata1, aes(x=Pos, y=Coverage)) +
  geom_area(color="lightgray",fill="lightgray",alpha=0.5) +
  geom_point(aes(color = Coverage)) +
  scale_color_gradientn(colours=c('red','orange',"cyan","blue",'blueviolet')) +
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 50)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ggtitle(sample1) +
  mytheme


#plot coverage per variant
mydata3 <- mydata2[,-10]
mydata4 <- cbind(mydata3,Var = TopVariants$variable)
mydata4 <- melt(mydata4, id=c("Pos","Nucl", "Var"))

p2 <-ggplot(mydata4, aes(x=as.factor(Pos)), y=value) +
  geom_col(aes(y=value, fill=variable, group=Nucl), width = 0.7) +
  ggtitle(sample1) +
  geom_text(aes(y = -100, label = Nucl), color = "black") +
  geom_text(aes(y = max(value)+200, label = Var), color = "gray") +
  mytheme

## draw the plots
p1
p2


