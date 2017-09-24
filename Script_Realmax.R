#1: View and extract dataset from the packages
#We need to download several package and file which contain the statistic and map function
#There are many packages need dependency packages, and need update some data. When we rerun the code,
#the error always happens, so we dicide to move all install and load packages in the top to remove error.
install.packages("adegenet", dependencies = T)
install.packages("cape")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")# We will need to use the Barplot and Maps function in this R script for visualization
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
library(LEA)#install and load LEA for snmf function use
library(stringi)#get the function stri_list2matrix
library(cape)#install and load packages which contains read geno function
library(adegenet)#install and load the package we need to use H3N2 file

data(package="adegenet")#View what dataset in this package
#We need to go through the information about H3N2 dataset, H3N2 is a genind object contain several data frame,factor, list,etc. 
#So we need to abstract the data from it; we view all the data and do some important tranformatioin first 
data(H3N2)#load H3N2 data
H3N2#View H3N2
#Frist we can know that the ploidy of all individual are 1, thus we can define ploidy =1 in snmf function
#type:tell us values are numbers of alleles, summing up to the individuals’ ploidies
#The loc.fac is the factor of matrix tab which we already know the level of each number.
#call:the matched call; loc.n.all: integer vector giving the number of alleles per locus
#based on the analysis of each data, we find we need tab, other
View(H3N2$other$x);View(H3N2$other$xy)

#we need to view the all.names data to know the component of the locus
H3N2$all.names
H3N2$loc.fac# it is a list in different length, we need to convert it to data frame for better viewer
allelesname <- as.data.frame(t(stri_list2matrix(H3N2$all.names)))
colnames(allelesname) <- unique(unlist(sapply(H3N2$all.names,names)))
View(allelesname)

#2: data preprocessing
#First we need to test if the SNP of alleles in tab data is pure noise.
#The test is basic code, so we decide to write in into function, thus we can make our code clear
l=ncol(H3N2$tab)
#get the probability of all the alleles in tab data. and combined all the larger probability and the number into one data frame
df=data.frame()
for (p in 2:l)
{
  table<-table(H3N2$tab[,p])
  pop=prop.table(table);pop
  pop=as.data.frame(pop)
  u=max(pop$Freq)
  pop=subset(pop,Freq>=u)
  cn<-p
  proba<-cbind(cn,pop)
  df=rbind(df,proba)
}#get the maximum probability of each coloumn 
SNP1=subset(df,Var1==1)
View(SNP1)
p=min(SNP1$Freq);p
#choose the minimum one and test the randomness, if we can make sure the minimum one is not random distributed, we have more confident that the larger one are not
pro=function()
{
  v=c("0","1");p=c(0.99,0.01)
  s=sample(v,1903,replace = T,prob=p)
  x=table(s)
  prob=prop.table(x)
  return(prob[2])
}
t=replicate(10000,pro())
plot(density(t),main="Sample Distribution")
polygon(density(t),col="green")
gap=abs(p-mean(t));gap
t1=mean(t)-gap;t2=mean(t)+gap
t1;t2
abline(v=t1);abline(v=t2)
y=t[t<t1|t>t2]
pv=length(y)/length(t)
pv
label=paste("P Value=",pv)
text(0.01,75,label,srt=0,cex=2.8,col="brown")
arrows(0.01,90,x1=0.015,y1=120,length=0.2,code=2,col="brown",lwd=3)
text(0.017,120,"Reject!",cex=3,col="brown")#For bette observation in larger version, we set cex to adjust the size. 
#Based on the p value, we can almost sure that the data we get are not pure noise.
#For using snmf function from LEA we need to get geno format data first
View(H3N2$tab) # Check the tab dataframe, which includes the 
write.geno(H3N2$tab[,29], "genotype.geno")#convert matrix to geno type
G = read.geno("genotype.geno")#read the geno data we just get

#3: Get the admixture coefficient as the measurement, Then visualize the findings
obj.snmf = snmf("genotypes.geno", K = 1:10, ploidy = 1, entropy = T,alpha = 100, project = "new")  #evaluate population structure using the snmf function for K = 1-10.
# Snmf function will get the cross-entropy of different K(cluster numbers)
# alpha indicates the regularization parameter of snmf
# entropy = T indicates we want to calculate the cross-entropy criterion.
# project = "new" meansthe current project is removed and a new one is created to store the result.

# Plot the cross-entropy curve on diffrent value of K, then select the best K
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)  #Looking at the results, the cross-entropy criterion does not exhibit a plateau when K=4
qmatrix = Q(obj.snmf, K = 4) # use Q to get the admixture coefficient for 4 clusters, qmatrix object contains the matrix of ancestry coefficients for each individual and for K = 4 clusters.
View(qmatrix)  # Take a look of the admixture coefficient of qmatrix
barplot(t(qmatrix), col = c("orange","violet","lightgreen","lightblue"), border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients") # Use Barplot to display the ancestry coefficient in qmatrix

#Get the coordinate for every sample in H3N2 data
coord = H3N2$other$x[,c(15,16)] # Extract the lontitude and latitude of every sample from data frame x 
plot(x=coord$lon,y=coord$lat, xlab = "Longitude", ylab = "Latitude", col="red",type = "p", pch=19) # plot all the sample locations using their lon and lat number

opar = par() 
par(bg="white",mfrow=c(2,2),las=2)

asc.data1="C:/Users/Donglin Jia/Downloads/Project/Europe.asc"  # Get the map data in asc format and assign it to asc.raster object
grid=createGridFromAsciiRaster(asc.data1)  # usecreateGridFromAsciiRaster() function to get all the lon and lat in the asc.raster object and put them into data frame grid
maps(matrix = qmatrix, coord, method = "max", grid, main = "Ancestry coefficients Europe", xlab = "Longitude", ylab = "Latitude", cex = .3) #Use the maps() function to plot all the admixture coefficient, using the coordinates data and limited to the "grid" data range from the map we provided
map(add = T, fill = FALSE, interior = T) # Add the national borders.

asc.data2="C:/Users/Donglin Jia/Downloads/Project/Asia.asc"
grid=createGridFromAsciiRaster(asc.data2)
maps(matrix = qmatrix, coord, method = "max", grid, main = "Ancestry coefficients Asia", xlab = "Longitude", ylab = "Latitude", cex = .3)
map(add = T, fill = FALSE, interior = T)

asc.data3="C:/Users/Donglin Jia/Downloads/Project/Central_Asia.asc"
grid=createGridFromAsciiRaster(asc.data3)
maps(matrix = qmatrix, coord, method = "max", grid, main = "Ancestry coefficients Central Asia", xlab = "Longitude", ylab = "Latitude", cex = .3)
map(add = T, fill = FALSE, interior = T)

asc.data4="C:/Users/Donglin Jia/Downloads/Project/North_America.asc"
grid=createGridFromAsciiRaster(asc.data4)
maps(matrix = qmatrix, coord, method = "max", grid, main = "Ancestry coefficients North America", xlab = "Longitude", ylab = "Latitude", cex = .3)
map(add = T, fill = FALSE, interior = T)

par(opar)
