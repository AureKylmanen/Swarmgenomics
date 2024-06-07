files<-read.table("files_list.txt",h=F)

ins_length<-200 ##### Minimum NUMT length to be retained
no_numts<-c()
nntm<-c()
for (inp in 1:length(files[,1])){ # loop through the diferent species
print(inp)
data<-read.table(paste(files[inp,1]),h=F)

dat200<-data[data[,3]>=ins_length,] # Keep only insertions at least 'ins_length' (defined above) long

##################################################################################################################
#### Sanity check 3: Introduce an if statement so that only not-empty data frames (i.e.: species with at least) ##
#### one numt longer than 200 (or whatever else specified in 'ins_length') bp) will be processed. Preventing    ##
#### the R script from stopping. The name os the files that have no numts longer than 200 bp will be output in  ##
#### a text file named "No_long_numts.txt".                                                                     ##
##################################################################################################################

if (dim(dat200)[1]!=0){

for (i in 1:length(dat200[,1])){ # this will switch values of column 4 and 5 (start and end positions of the numt in the nuclear genome) to make sure the start pos is always smaller than the end pos (i.e.: "correct reverse blast hits")
if (dat200[i,4]>dat200[i,5]){
bla<-dat200[i,4]
dat200[i,4]<-dat200[i,5]
dat200[i,5]<-bla
}}


dat200_sorted <- dat200[order(dat200[,1], dat200[,4]),] #Sort remaining blast hits and output a bed file to extract all numt fragments separately (output also the entire file as backup)
dat200_bed<-dat200_sorted[,c(1,4,5)]
# write.table(dat200_sorted,paste(files[inp,1],".filtered.sorted",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
write.table(dat200_bed,paste(files[inp,1],"_AllFrags.bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")

bla<-c() # Now, merge together blast hits for numts that are less then 10kb apart from each other. Store merged numts info in "bla"

##################################################################################################################
#### Sanity check 4: Introduce if statement so that numts are merged only in species that has at least 2 NUMTs, ##
#### preventing the R script from exiting. All species' names for which merged NUMTs will not be created are    ##
#### stored in the text file "No_numts_to_merge                                                                 ##
################################################################################################################## 

#if (dim(dat200)[1]>1){
#for (j in 2:length(dat200_bed[,1])){
#if (dat200_bed[j,1]==dat200_bed[j-1,1] && (dat200_bed[j,2]-dat200_bed[j-1,3])<10000){
#dat200_bed[j,2]<-dat200_bed[j-1,2]
#bla[j]<-j-1
#}}} else {
#nntm[inp]<-paste(files[inp,1])}

#if (length(bla)!=0){
#bla<-bla[!is.na(bla)]
#dat200_bed_10k<-dat200_bed[-bla,]
# and including the nuclear sequence between them)
#write.table(dat200_bed_10k,paste(files[inp,1],"_AllFrags_merged10k.bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")} # Output bed file to extract numts (fragments close to each other extracted as part of the same numt
} else {
no_numts[inp]<-paste(files[inp,1])
}
#}
#if (length(no_numts!=0)){
#no_numts<-no_numts[!is.na(no_numts)]
#write.table(no_numts,"No_long_numts.txt",quote=F,col.names=F,row.names=F,sep="\t")}
#if (length(nntm!=0)){
#nntm<-nntm[!is.na(nntm)]
#write.table(nntm,"No_numts_to_merge.txt",quote=F,col.names=F,row.names=F,sep="\t")
}