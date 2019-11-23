# 生物信息学科研训练笔记

时间：11月23日   

作者：rabbit



## 1-RSQLite

- Read online documents and practice
  - [How to Use SQLite with R](https://www.bioconductor.org/help/course-materials/2006/rforbioinformatics/labs/thurs/SQLite-R-howto.pdf)
  - [RSQLite Tutorial](https://github.com/ysquared2/RSQLiteTutorial)
  - [SQLite in R](https://www.datacamp.com/community/tutorials/sqlite-in-r)
  - [A Gentle Introduction to SQL Using SQLite](https://a-gentle-introduction-to-sql.readthedocs.io/en/latest/)

- Implement class schedule query using SQLite  

```R
#first: load the data
df<-read.csv("classtable.csv",header = T)
df<-df[,-1]

#second: load the RSQLite Library
library(RSQLite)

#third: create a connection to our new database, classDB.db, you can check that the .db file has been created on your working directory
conn <- dbConnect(RSQLite::SQLite(), "classDB.db")

#fourth: write the df dataset into a table names classtable_data
dbWriteTable(conn, "classtable_data", df)

#fifth: list all the tables available in the database
dbListTables(conn)

#sixth: gather the first 10 rows in the classtable_data table
dbGetQuery(conn, "SELECT * FROM classtable_data LIMIT 10")

#seventh: we also can do
dbGetQuery(conn, "SELECT day,week,secondclass FROM classtable_data LIMIT 10")
dbGetQuery(conn,"SELECT * FROM classtable_data
                 WHERE day = 903")
dbGetQuery(conn,"SELECT day,week,secondclass FROM classtable_data
                 WHERE day = 903")
dbGetQuery(conn,"SELECT day,week,secondclass FROM classtable_data
                 WHERE day LIKE '9%' AND week IN (5,6,7)")

#finally: Close the database connection to classDB
dbDisconnect(conn)
```

代码和课表的下载链接在这里：  

* 代码：[classtable.R](https://github.com/PPeachi/RSQLite/blob/master/classtable.R)
* 表格：[classtable.csv](https://github.com/PPeachi/RSQLite/blob/master/classtable.csv)



## 2-fasta

- Write an R function to parse FASTA file
  - [seq1](https://github.com/PPeachi/fasta/blob/master/flu_seq.fas)
  - [seq2](https://github.com/PPeachi/fasta/blob/master/flu_seq_v2.fas)

- Calculate base frequency table

```R
#读入fasta格式的序列数据
read_fas<-function(file){
  line<-readLines(file)
  i<-grep(">",line)
  id<-sub(">","",line[i])
  start<-(i+1)
  end<-c(i[-1]-1,length(line))
  seq<-sapply(seq_along(start), function(i) paste0(line[start[i]:end[i]],collapse = ""))
  df<-data.frame(id=id,seq=seq)
  df$id<-as.character(df$id)
  df$seq<-as.character(df$seq)
  return (df)
}
x<-read_fas("flu_seq.fas")

#计算碱基比例        
base_freq<-function(fasta){
  seq_num <- length(fasta$id)
  cat(seq_num,"sequences in total","\n","\n")
  cat("#Labels:","\n")
  id<-fasta$id
  writeLines(id)
  ids<-c()
  len<-c()
  for (i in 1:seq_num){
    ids[i]<-unlist(strsplit(id[i],split = "[|]"))[1]
    len[i]<-nchar(fasta$seq[i])
  }
  df<-data.frame(id=ids,length=len)
  for (i in 1:seq_num){
    count_A<-0
    count_T<-0
    count_G<-0
    count_C<-0
    s<-toupper(unlist(strsplit(fasta$seq[i],"")))
    for (j in 1:length(s)){
      if(s[j]=='A'){count_A=count_A+1}
      if(s[j]=='T'){count_T=count_T+1}
      if(s[j]=='G'){count_G=count_G+1}
      if(s[j]=='C'){count_C=count_C+1}
    }
    m<-length(s)
    df$A[i]<-count_A/m
    df$T[i]<-count_T/m
    df$G[i]<-count_G/m
    df$C[i]<-count_C/m
  }
  rownames(df)<-df[,1]
  tb<-df[,-1]
  res<-df[,-(1:2)]
  writeLines(paste0("\n","#Base composition:"))
  print(tb)
  return (res)
}
y<-base_freq(x)
y

#画热图
pheatmap::pheatmap(y)
```

做出的热图如下：  

<img src="C:\Users\18307\Documents\大四\课件\科研训练\R_y叔\2-fasta\plot.png" alt="heatmap" style="zoom:75%;" />

代码下载链接在这里：  

- 代码：[fasta](https://github.com/PPeachi/fasta/blob/master/fasta_base_freq.R)



## 3-genbank2fasta

- Write an R function to extract accession number and sequence data from genbank file
- Write an R function export the output to a FASTA file

```R
#读入genbank格式的序列数据
read_genbank<-function(file){
  seq<-readLines(file)
  a_n<-unlist(strsplit(seq[1],split = "\\s+"))[2]
  st<-grep("ORIGIN",seq)
  ed<-grep("^//",seq)
  gb_seq<-gsub(" ","",gsub("\\d","",paste(seq[(st+1):(ed-1)],collapse = "")))
  count<-nchar(gb_seq)
  writeLines(paste0("accession number:","\n",a_n,"\n","\n","sequence:","\n",gb_seq,"\n","\n","sequence length:","\n",count,"\n","\n","Base compositon:","\n"))
  s<-tolower(unlist(strsplit(gb_seq,"")))
  count_A<-0
  count_T<-0
  count_G<-0
  count_C<-0
  for (j in 1:length(s)){
      if(s[j]=='a'){count_A=count_A+1}
      if(s[j]=='t'){count_T=count_T+1}
      if(s[j]=='g'){count_G=count_G+1}
      if(s[j]=='c'){count_C=count_C+1}
  }
  m<-length(s)
  freq_a<-count_A/m
  freq_t<-count_T/m
  freq_g<-count_G/m
  freq_c<-count_C/m
  freq<-data.frame(a_n=a_n,a=freq_a,t=freq_t,g=freq_g,c=freq_c)
  print(freq)
}
x<-read_genbank("AB115403.gb")
x

#将genbank格式转为fasta格式
genbank2fasta<-function(file){
  seq<-readLines(file)
  a_n<-paste0(">",unlist(strsplit(seq[1],split = "\\s+"))[2])
  st<-grep("ORIGIN",seq)
  ed<-grep("^//",seq)
  gb_seq<-gsub(" ","",gsub("\\d","",paste(seq[(st+1):(ed-1)],collapse = "")))
  res<-paste0(a_n,"\n",gb_seq)
  return(res)
}
y<-genbank2fasta("AB115403.gb")
y
writeLines(y)
cat(y)
write.table(y,
            file = "AB115403.fas",
            row.names = F,
            col.names = F,
            quote = F)
```

代码下载链接在这里：  

- 代码：[genbank](https://github.com/PPeachi/genbank/blob/master/genbank2fasta.R)



## 4-批量下载genbank序列

- read the API of NCBI Entrez system
- write an R function to download genbank sequences

```R
dir.create("gb")
dir.create("fasta")
getwd()
wd<-"C:/Users/18307/Documents/大四/课件/科研训练/R_y叔/4-download_genbankORfasta"
setwd(wd)
setwd("gb")
setwd("fasta")
accn<-paste("AJ5345", 26:49, sep="")
accn
list.files()
#acc<-accn
#database<-"nucleotide"
#returntype<-"gb"
#returnmode<-"text"
#base<- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
download_GBorFASTA<-function(acc,database,returntype,returnmode){
  base<- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  for(i in 1:length(acc)){
    url=paste(base,"?db=",database,"&rettype=",returntype,"&retmode=",returnmode,"&id=",acc[i],sep = "")
    if(returntype=="gb"){files<-paste(acc[i],".gb",sep = "")}
    if(returntype=="fasta"){files<-paste(acc[i],".fas",sep = "")}
    download.file(url,destfile = files)
  }
}
download_GBorFASTA(accn,"nucleotide","gb","text")
download_GBorFASTA(accn,"nucleotide","fasta","text")
list.files()
```

代码下载链接在这里：  

- 代码：[download](https://github.com/PPeachi/download_gbORfasta/blob/master/download_gbORfasta.R)



## 5-使用动态规划实现双序列比对

话不多说，上代码：  

```R
x<-'TTCATA'
y<-'TGCTCGTA'
global_align<-function(x,y,match,mismatch,gap){
  m<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  n<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  l<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  for(j in 2:ncol(m)){
    m[1,j]<-(m[1,j-1]+gap)
    n[1,j]<-'←'
    l[1,j]<-'left'
  }
  for(i in 2:nrow(m)){
    m[i,1]<-(m[i-1,1]+gap)
    n[i,1]<-'↑'
    l[i,1]<-'up'
  }
  for(i in 2:nrow(m)){
    for(j in 2:ncol(m)){
      xl<-unlist(strsplit(x,split = ""))
      yl<-unlist(strsplit(y,split = ""))
      if(xl[i-1]==yl[j-1]){s1<-m[i-1,j-1]+match} else {s1<-m[i-1,j-1]+mismatch}
      s2<-(m[i-1,j]+gap)
      s3<-(m[i,j-1]+gap)
      m[i,j]<-max(s1,s2,s3)
      if(s1==m[i,j]){n[i,j]<-'↖';l[i,j]<-'diag'} 
      if(s2==m[i,j]){n[i,j]<-'↑';l[i,j]<-'up'}
      if(s3==m[i,j]){n[i,j]<-'←';l[i,j]<-'left'}
      if(m[i,j]==s1&&m[i,j]==s2){n[i,j]<-paste0('↖','↑');l[i,j]<-'diag_up'}
      if(m[i,j]==s1&&m[i,j]==s3){n[i,j]<-paste0('↖','←');l[i,j]<-'diag_left'}
      if(m[i,j]==s2&&m[i,j]==s3){n[i,j]<-paste0('↑','←');l[i,j]<-'up_left'}
      if(m[i,j]==s1&&m[i,j]==s2&&m[i,j]==s3){n[i,j]<-paste0('↖','↑','←');l[i,j]<-'diag_up_left'}
    }
  }
  cat("Sequence x:",x,"\n")
  cat("Sequence y:",y,"\n")
  cat("Scoring system:",match,"for match,",mismatch,"for mismatch,",gap,"for gap","\n","\n")
  cat("Dynamic programming matrix 1:","\n")
  print(m)
  cat("Dynamic programming matrix 2:","\n")
  print(n)
  cat("Dynamic programming matrix 3:","\n")
  print(l)
  writeLines(paste0("\n","Alignment:"))
  lx<-unlist(strsplit(x,split = ""))
  ly<-unlist(strsplit(y,split = ""))
  L2=length(lx)
  L1=length(ly)
  D<-c('left','up','diag')
  if(l[L2+1,L1+1]==D[1]){ 
    RX2<-c(ly[L1]);RX1<-c("-");i=L2+1;j=L1 
  } else {
    if(l[L2+1,L1+1]==D[2]){
      RX2<-c("-");RX1<-c(lx[L2]);i=L2;j=L1+1
    } else { 
      RX2<-c(ly[L1]);RX1<-c(lx[L2]);i=L2;j=L1 
    } 
  } 
  while((i>1)&&(j>1)){ 
    # browser() 
    if(l[i,j]==D[1]){ 
      RX2<-c(ly[j-1],RX2);RX1<-c("-",RX1);j=j-1 
    } 
    else if(l[i,j]==D[2]){ 
      RX2<-c("-",RX2);RX1<-c(lx[i-1],RX1);i=i-1 
    } 
    else {RX2<-c(ly[j-1],RX2);RX1<-c(lx[i-1],RX1);j=j-1;i=i-1} 
  }
  RX3<-c()
  for(r in 1:length(RX1)){
    if(RX1[r]==RX2[r]){RX3<-c(RX3,"|")}
    else {RX3<-c(RX3," ")}
  }
  #hamming_distance<-(max(nchar(x),nchar(y))-length(grep("[|]",RX3)))
  if(length(lx)>length(ly)){lmer<-ly;s<-lx} else {lmer<-lx;s<-ly}
  hamming_distance<-numeric(length(s)-length(lmer)+1)
  for(i in 1:(length(s)-length(lmer)+1)){
    ss<-s[i:(i+length(lmer)-1)]
    for(j in 1:length(lmer)){
      if(lmer[j]!=ss[j]){hamming_distance[i]<-(hamming_distance[i]+1)}
    }
  }
  hamming_distance<-min(hamming_distance)
  cat(" x: ",RX1,"\n","   ",RX3,"\n","y: ",RX2,"\n")
  writeLines(paste0("\n","#1 score: ",m[L2+1,L1+1],"\n","#2 hamming-distance: ",hamming_distance))
}
#writeLines(paste0(paste0(unlist(strsplit(ss,split = "\n"))[1],unlist(strsplit(ss2,split = "\n"))[1]),"\n",paste0(unlist(strsplit(ss,split = "\n"))[2],unlist(strsplit(ss2,split = "\n"))[2]),"\n",paste0(unlist(strsplit(ss,split = "\n"))[3],unlist(strsplit(ss2,split = "\n"))[3])))
z<-global_align(x,y,5,-2,-6)
```

得到的比对结果为：  

```R
#Sequence x: TTCATA 
#Sequence y: TGCTCGTA 
#Scoring system: 5 for match, -2 for mismatch, -6 for gap 
#
#Dynamic programming matrix 1: 
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#[1,]    0   -6  -12  -18  -24  -30  -36  -42  -48
#[2,]   -6    5   -1   -7  -13  -19  -25  -31  -37
#[3,]  -12   -1    3   -3   -2   -8  -14  -20  -26
#[4,]  -18   -7   -3    8    2    3   -3   -9  -15
#[5,]  -24  -13   -9    2    6    0    1   -5   -4
#[6,]  -30  -19  -15   -4    7    4   -2    6    0
#[7,]  -36  -25  -21  -10    1    5    2    0   11
#Dynamic programming matrix 2: 
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#[1,] "0"  "←"  "←"  "←"  "←"  "←"  "←"  "←"  "←" 
#[2,] "↑"  "↖"  "←"  "←"  "↖←" "←"  "←"  "↖←" "←" 
#[3,] "↑"  "↖↑" "↖"  "↖←" "↖"  "←"  "←"  "↖←" "←" 
#[4,] "↑"  "↑"  "↖↑" "↖"  "←"  "↖"  "←"  "←"  "←" 
#[5,] "↑"  "↑"  "↖↑" "↑"  "↖"  "↖←" "↖"  "↖←" "↖" 
#[6,] "↑"  "↖↑" "↖↑" "↑"  "↖"  "↖"  "↖←" "↖"  "←" 
#[7,] "↑"  "↑"  "↖↑" "↑"  "↑"  "↖"  "↖"  "↑"  "↖" 
#Dynamic programming matrix 3: 
#     [,1] [,2]      [,3]      [,4]        [,5]        [,6]       
#[1,] "0"  "left"    "left"    "left"      "left"      "left"     
#[2,] "up" "diag"    "left"    "left"      "diag_left" "left"     
#[3,] "up" "diag_up" "diag"    "diag_left" "diag"      "left"     
#[4,] "up" "up"      "diag_up" "diag"      "left"      "diag"     
#[5,] "up" "up"      "diag_up" "up"        "diag"      "diag_left"
#     [,7]        [,8]        [,9]  
#[1,] "left"      "left"      "left"
#[2,] "left"      "diag_left" "left"
#[3,] "left"      "diag_left" "left"
#[4,] "left"      "left"      "left"
#[5,] "diag"      "diag_left" "diag"
#[6,] "diag_left" "diag"      "left"
#[7,] "diag"      "up"        "diag"
#
#Alignment:
# x:  T - - T C A T A 
#     |     | |   | | 
# y:  T G C T C G T A 
#
#1 score: 11
#2 hamming-distance: 2
```

代码下载链接在这里：  

- 代码：[program](https://github.com/PPeachi/dynamic_programming/blob/master/program.R)



## 6-make构建流程

- Write makefile
  - download genbank
  - convert it to fasta
  - do alignment
  - calculate hamming distance

写了一个makefile：

```bash
all: download_gb.Rout gb2fasta.Rout align.Rout

download_gb.Rout: accn.txt download_gb.R
	R CMD BATCH download_gb.R

gb2fasta.Rout: accn.txt download_gb.Rout gb2fasta.R
	R CMD BATCH gb2fasta.R

align.Rout: accn.txt download_gb.Rout gb2fasta.Rout align.R
	R CMD BATCH align.R
```

把需要的文件搞到服务器上，输入命令`make`就可以了  

文件和代码的下载链接在这里：

- 整个跑出结果的文件夹：[make_out](https://github.com/PPeachi/make/tree/master/make2_out)



## 7-git和github

收藏了一个很好用的教程：  

- 教程：[git bash](https://blog.csdn.net/xiangwanpeng/article/details/54178653)

主要就是下面这几条命令：

```bash
#初始化
git init

#本地仓库和远程仓库建立连接
git remote add origin [远程仓库链接]

#把文件推送到远程仓库
git add [文件]
git commit -m'first commit'

#从远程仓库扒拉一个readme文件
git fetch origin master

#合并
git merge origin/master --allow-unrelated-histories

#最后上传
git push -u origin master
```



## 8-rmarkdown

我最后选择用typora，确实是好用！  

教程很多，比如菜鸟教程的markdown（[菜鸟教程](https://www.runoob.com/markdown/md-tutorial.html)）



# 谢谢你们看我表演