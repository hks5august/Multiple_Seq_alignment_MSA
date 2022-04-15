#Load data
#msa <- read.csv("https://raw.githubusercontent.com/pine-bio-support/COVID-19-origin/main/SARS_CoV_2_1.csv", header=TRUE)
msa <- read.table("input_data - Copy.txt", sep="\t",  header=TRUE)


#Define selected positions (Here, we are defining position for spike protein)

startnt <- 21563
endnt <- 25384


#select specified rows from the letter data frame
msa_sel1 <- msa[startnt:endnt,] 

#Extract position column
pos <- msa_sel1[1]

msa_sel <- msa_sel1[,-1]  #remove the "POS" column
#Report on the number of mutations in each sample:
x <- ncol(msa_sel)

#create an empty data frame to store the new data
daf2 <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))
colnames(daf2) <- colnames(msa_sel)

#define i value
i=2

for (i in 2:x) {
  daf2[,i] <- ifelse(msa_sel[,1]== msa_sel[,i],0,1)
}

#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
daf3 <- select(daf2, -ncol(daf2))


#count the number of mutations in each sample
summ1 <- colSums(daf3)

summ1
#make a barplot
par(mar=c(11,4,4,4))
barplot(summ1, las=2, main= "ORF1ab mutations", cex.axis = 0.6, cex.lab=0.5, cex.names=0.5,  ylab="Number of mutations", col="blue", space=0.2)
b<- barplot(summ1, las=2, main= "ORF1ab mutations", cex.axis = 0.6, cex.lab=0.5, cex.names=0.5,  ylab="Number of mutations", col="blue", space=0.2)


#Check dimensions
dim(daf3)



########### Variant report######



#select specified rows from the letter data frame
msa_sel1 <- msa[startnt:endnt,] 

#Extract position column
pos <- msa_sel1[1]

msa_sel <- msa_sel1[,-1]  #remove the "POS" column
#Report on the number of mutations in each sample:
x <- ncol(msa_sel)

#create an empty data frame to store the new data
df_new1 <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))
colnames(df_new1) <- colnames(msa_sel)

#define i value
i=2

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_new1[,i] <- ifelse(msa_sel[,1]== msa_sel[,i],"-", paste0(msa_sel[,1],">",msa_sel[,i]))
}


#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
df_new2 <- select(df_new1,-1, -ncol(df_new1))

#Add position column
final_df <- cbind(pos, df_new2)

#Check whether output is correct
head(final_df)

#Write into a file
write.table(final_df, file="ORF1ab_variant_info.txt", row.names = F, sep ="\t")




################# MSA for variants ############
#load libraries
library(gplots)

#load and prepare data
#msa1 <- read.csv("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/SARS_CoV_2_1.csv", header=TRUE)
msa <- read.table("input_data - Copy.txt", sep="\t",  header=TRUE)

#prepare a numeric matrix for msa1 data frame
msa1m <- data.matrix(msa) 

#remove the "POS" column
msa1 <- msa[,-1] 



#select from the letter data frame
msa1sel <- msa1[startnt:endnt,]

#select from the numeric matrix
msa1msel <- msa1m[startnt:endnt,]

#for numeric matrix, make sure to add the position to row names
row.names(msa1msel) <- msa1msel[,1] 

#now, remove the POS column from the input data
msa1msel <- msa1msel[,-1]

#select only data where any sample has a variant to create a report
new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9]),] 
new_dfT <- t(new_df)

#select only data where any sample has a variant in the matrix
new_matrix <- msa1msel[!(msa1msel[,2] == msa1msel[,3] & msa1msel[,2] == msa1msel[,4] & msa1msel[,2] == msa1msel[,5] & msa1msel[,2] == msa1msel[,6] & msa1msel[,2] == msa1msel[,7] & msa1msel[,2] == msa1msel[,8] & msa1msel[,2] == msa1msel[,9]),]

#transpose
new_matrixT <- t(new_matrix)

#set colors and margins for plot
if("-" %in% new_matrixT){ #first, check if we have 4 or 5 characters (I am assuming that TCGA will be present, so only looking for "-")
  TCGAcolors <- c("white", "lightgreen", "pink", "lightblue", "yellow") 
  names(TCGAcolors) = c(1,2,3,4,5) #labels = c("-","A","T","C","G")
} else {TCGAcolors <- c("lightgreen", "pink", "lightblue", "yellow") }




#draw heat map for only variants
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
          main = paste(" No. of Variants Found in ORF1ab:",ncol(new_matrixT)), #heat map title
          cex.main=0.6,
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

dim(new_matrixT)



