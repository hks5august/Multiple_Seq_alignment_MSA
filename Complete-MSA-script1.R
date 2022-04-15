library(gplots) #library that contains functions for drawing plots like heatmap.2

quartzFonts(avenir = c("Avenir Book", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"))

#load and prepare data
msa1 <- read.csv("https://raw.githubusercontent.com/pine-bio-support/COVID-19-origin/main/SARS_CoV_2_Origin_project_FullTable_1.csv", header=TRUE)

#example 1: visualize 50 rows
msa1 <- msa1[,-1] #remove the "POS" column
msa1m <- data.matrix(msa1) #prepare a numeric matrix for msa1 data frame

#select only 50 rows
msa1sel <- msa1[21000:21100,] #select only 100 rows from the letter data frame
msa1msel <- data.matrix(msa1sel) #select only 100 rows from the numeric matrix

#transpose
msa1selT <- t(msa1sel)
msa1mselT <- t(msa1msel)

#set colors and margins for plot
TCGAcolors <- c("white", "lightgreen", "pink", "lightblue", "yellow") 
names(TCGAcolors) = c(1,2,3,4,5) #labels = c("-","A","T","C","G")

par(cex.main=0.9, family="avenir") #set plot margins

#draw heatmap for first 50 positions of alignment
heatmap.2(msa1mselT, #data source matrix
          #main settings
          cexRow = 0.7, #row name font size
          col = TCGAcolors, #set colors
          dendrogram = "none", #remove dendrogram
          Rowv = FALSE, #no reordering for rows
          Colv = FALSE, #no reordering for columns
          density.info="none", #remove density info
          trace="none", #remove row and column lines
          
          offsetRow=0.1, #change position of the row names
          offsetCol=0.1, #change position of the column names
          
          #add gray borders between cells
          sepwidth=c(0.05,0.05), #sets separation width and height
          sepcolor="gray", #color for border
          colsep=1:ncol(msa1mselT), #add separation for number of columns in source data
          rowsep=1:nrow(msa1mselT), #add separation for number of rows in source data
          
          #plot title
          main = "SARS_CoV_2_Origin_project_FullTable_1.csv", #heat map title
          
          #plot margins
          margins = c(5,10), #set margins
          lwid=c(0.2,4), 
          lhei=c(0.9,3),
          
          #adding letters inside the heatmap
          notecex=1.0, #size of font inside each cell
          cellnote = msa1selT, #data to use in cells
          notecol="black", #font color for cells
          
          #legend
          key = FALSE
)

#select only data where any sample has a variant
new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9] & msa1sel[2] == msa1sel[10] & msa1sel[2] == msa1sel[11] & msa1sel[2] == msa1sel[12] & msa1sel[2] == msa1sel[13] & msa1sel[2] == msa1sel[14]),] 
new_dfT <- t(new_df)

new_matrix <- msa1msel[!(msa1msel[,2] == msa1msel[,3] & msa1msel[,2] == msa1msel[,4] & msa1msel[,2] == msa1msel[,5] & msa1msel[,2] == msa1msel[,6] & msa1msel[,2] == msa1msel[,7] & msa1msel[,2] == msa1msel[,8] & msa1msel[,2] == msa1msel[,9] & msa1msel[,2] == msa1msel[,10] & msa1msel[,2] == msa1msel[,11] & msa1msel[,2] == msa1msel[,12] & msa1msel[,2] == msa1msel[,13] & msa1msel[,2] == msa1msel[,14]),]
new_matrixT <- t(new_matrix)

#draw heatmap for only variants
heatmap.2(new_matrixT, #data source
          
          #main settings
          cexRow = 0.7, #row name font size
          col = TCGAcolors, #set colors
          dendrogram = "none", #remove dendrogram
          Rowv = FALSE, #no reordering for rows
          Colv = FALSE, #no reordering for columns
          density.info="none", #remove density info
          trace="none", #remove row and column lines
          
          offsetRow=0.1, #change position of the row names
          offsetCol=0.1, #change position of the column names
          
          #add gray borders between cells
          sepwidth=c(0.05,0.05), #sets separation width and height
          sepcolor="gray", #color for border
          colsep=1:ncol(new_matrixT), #add separation for number of columns in source data
          rowsep=1:nrow(new_matrixT), #add separation for number of rows in source data
          
          #plot title
          main = cbind("Number of Variants Found: ",ncol(new_matrixT)), #heat map title
          
          #plot margins
          margins = c(5,10), #set margins
          lwid=c(0.2,4), 
          lhei=c(0.9,3),
          
          #adding letters inside the heatmap
          notecex=1.0, #size of font inside each cell
          cellnote = new_dfT, #data to use in cells
          notecol="black", #font color for cells
          
          #legend
          key = FALSE
)

#dev.off()  #used to reset plotting