# include libraries
library("igraph")
library("stringr")
library("plotrix")

# parameter setting
year_start = 2005
year_end = 2010
stepsize = 1
threshold_size = 0.05
threshold_strength = 0.30
scale_cir = 800

time = seq(year_start, year_end, by=stepsize)

col_age = read.csv("C:/Rwork/col_age.csv",header=TRUE)
col_normalizedsize = read.csv("C:/Rwork/col_normalizedsize.csv",header=TRUE)

y_lowerlim = max(floor(min(col_age)*10)/10 - 0.1, 0)
y_upperlim = ceiling(max(col_age)*10)/10 + 0.1

# depict timeline plot of first year

year_ini = as.character(year_start)
string0_1 = paste("C:/Rwork/cluster",year_ini,".csv", sep="")
string0_2 = paste("C:/Rwork/cluster_sta",year_ini,".csv", sep="")
string0_3 = paste("C:/Rwork/collector_new",year_ini,".csv", sep="")
x_start_cluster = read.csv(string0_1, header=TRUE)
x_start = read.csv(string0_2, header=TRUE)
eval(parse(text=paste("data_sta",year_ini," = x_start",sep="")))
collector2 = read.csv(string0_3, header=FALSE)

# Initialize the timeline plot
plot(1, 1, xlab="Year", ylab="Average age", xlim=c(year_start-0.5, year_end+0.9), ylim=c(y_lowerlim,y_upperlim), type="n")
for (j in 1:nrow(x_start)){
    if (eval(parse(text=paste("x_start[",j,",3]",sep=""))) >= threshold_size){
        Immediate = collector2[j,j]/x_start[j,1]
        floating.pie(year_start, x_start[j,2], c(1-Immediate+0.001,Immediate+0.001), radius=x_start[j,1]/scale_cir, col=c("gray","black"))
    }
}

# depict timeline plot of subsequent years
col_link = 0
col_link_info = 0
for (i in time[2:length(time)]){
   
    year_new = as.character(i)
    year_past = as.character(i - stepsize)
    string1 <- paste("C:/Rwork/cluster_sta",year_past,".csv", sep="")
    string2 <- paste("C:/Rwork/cluster_sta",year_new,".csv", sep="")
    string3 <- paste("C:/Rwork/data_",year_past,".csv", sep="")
    string4 <- paste("C:/Rwork/data_",year_new,".csv", sep="")
    x_past = read.csv(string1,header=TRUE)
    x_new = read.csv(string2,header=TRUE)
    eval(parse(text=paste("data_sta",year_new," = x_new",sep="")))
    x_past_cluster = read.csv(string3,header=TRUE)
    x_new_cluster = read.csv(string4,header=TRUE)
    # load transistion matrix
    string5 <- paste("C:/Rwork/collector_old",year_new,".csv", sep="")
    string6 <- paste("C:/Rwork/collector_new",year_new,".csv", sep="")

    collector = read.csv(string5,header=FALSE)
    collector2 = read.csv(string6,header=FALSE)
    
    # collect string strength of each inter-year
    for (j in 1:nrow(x_new)){
        if (eval(parse(text=paste("x_new[",j,",3]",sep=""))) >= threshold_size){
            for (k in 1:nrow(x_past)){
                if (eval(parse(text=paste("x_past[",k,",3]",sep=""))) >= threshold_size) {
                    col_link = c(col_link, eval(parse(text=paste("collector[",k,",",j,"]/x_past[",k,",1]",sep=""))))
                }
            }
        }
    }
    col_link = col_link[2:length(col_link)]
    write.table(col_link, file = "C:/Rwork/col_link.csv", sep = ",", row.names = FALSE, col.names = TRUE, qmethod = "double")

    # Begin to depict latter fronts and links
    for (j in 1:nrow(x_new)){
        if (eval(parse(text=paste("x_new[",j,",3]",sep=""))) >= threshold_size){
            Inheritable = 0
            for (k in 1:nrow(x_past)){
                if (eval(parse(text=paste("x_past[",k,",3]",sep=""))) >= threshold_size && eval(parse(text=paste("collector[",k,",",j,"]/x_past[",k,",1]",sep=""))) > threshold_strength){
                    segments(i-stepsize+x_past[k,1]/scale_cir, x_past[k,2], i-x_new[j,1]/scale_cir, x_new[j,2], col= 'red', lwd = 0.5)
                    Inheritable = Inheritable + eval(parse(text=paste("collector[",k,",",j,"]/x_new[",j,",1]",sep="")))
                    col_link_info = c(col_link_info, paste(i-stepsize,"_", k, sep=""), paste(i,"_", j, sep=""))
                }
            }
            Immediate = eval(parse(text=paste("collector2[",j,",",j,"]/x_new[",j,",1]",sep="")))
            floating.pie(i, x_new[j,2], c(Inheritable+0.001,1-Inheritable-Immediate+0.001,Immediate+0.001), radius=x_new[j,1]/scale_cir, col=c("white","gray","black"))
        }
    }
}

if (length(col_link_info) > 1){
  col_link_info = col_link_info[2:length(col_link_info)]
  el2 = t(matrix(col_link_info,2,byrow=F))
  g2 = graph.edgelist(el2,directed=FALSE)
  c2 <- clusters(g2)

  col_year = "T1"
  col_membership = "-----"
  col_size = "-----"
  col_average_age = "-----"
  col_normalized_size = "-----"
  col_keyword = "-----"
  col_newkeyword = "-----"
  col_uspc = "-----" #
  col_smybol_year = 0
  col_smybol_average_age = 0
  colbar = rainbow(c2$no)
  for (j_1 in 0:(c2$no-1)){
      trajectory_temp1 = sort(V(g2)$name[which(c2$membership == j_1)])
      for (j_2 in 1:(length(trajectory_temp1))){
          trajectory_temp2 = str_split(trajectory_temp1[j_2],"_")[[1]]
          eval(parse(text=paste("trajectory_temp3 = data_sta",trajectory_temp2[1],"[",trajectory_temp2[2],",]",sep="")))
          col_membership = c(col_membership, trajectory_temp2[2])
          col_year = c(col_year, trajectory_temp2[1])
          col_size = c(col_size, trajectory_temp3[[1]])
          col_average_age = c(col_average_age, trajectory_temp3[[2]])
          col_normalized_size = c(col_normalized_size, trajectory_temp3[[3]])
          col_keyword = c(col_keyword, paste(trajectory_temp3[[4]]))
          col_newkeyword = c(col_newkeyword, paste(trajectory_temp3[[5]]))
          col_uspc = c(col_uspc, paste(trajectory_temp3[[6]])) #
          if (j_2 == (length(trajectory_temp1))){
              col_smybol_year = c(col_smybol_year, as.numeric(trajectory_temp2[1]))
              col_smybol_average_age = c(col_smybol_average_age, trajectory_temp3[[2]])
              countsimilar = length(which(col_smybol_average_age[which(col_smybol_year == as.numeric(trajectory_temp2[1]))]==trajectory_temp3[[2]]))
              text(as.numeric(trajectory_temp2[1])+(countsimilar)*0.25, trajectory_temp3[[2]], paste("T",j_1+1,sep=""), col= colbar[j_1+1])
          }
      }
      if (j_1 < c2$no-1){
          col_year = c(col_year, paste("T",j_1+2,sep=""))
          col_membership = c(col_membership, "-----")
          col_size = c(col_size, "-----")
          col_average_age = c(col_average_age, "-----")
          col_normalized_size = c(col_normalized_size, "-----")
          col_keyword = c(col_keyword, "-----")
          col_newkeyword = c(col_newkeyword, "-----")
          col_uspc = c(col_uspc, "-----") #
      }
  }

  for (j_1 in 0:(c2$no-1)){
      trajectory_temp1 = sort(V(g2)$name[which(c2$membership == j_1)])
      for (j_2 in 1:(length(trajectory_temp1)-1)){
          trajectory_temp2 = str_split(trajectory_temp1[j_2],"_")[[1]]
          eval(parse(text=paste("trajectory_temp3 = data_sta",trajectory_temp2[1],"[",trajectory_temp2[2],",]",sep="")))
          for (j_3 in 1:(length(trajectory_temp1))){
              trajectory_temp4 = str_split(trajectory_temp1[j_3],"_")[[1]]
              trajectory_temp4
              eval(parse(text=paste("trajectory_temp5 = data_sta",trajectory_temp4[1],"[",trajectory_temp4[2],",]",sep="")))
              check = 0
              len = length(which(trajectory_temp1[j_3]==el2[,2]))
              if (len > 0){
                  for (j_4 in 1:len){
                      check = check + length(which((which(trajectory_temp1[j_2]==el2[,1]) == which(trajectory_temp1[j_3]==el2[,2])[j_4])==TRUE)==1)
                  }
                  check
                  if (check==1){
                      segments(as.numeric(trajectory_temp4[1])-as.numeric(trajectory_temp5[[1]])/scale_cir, as.numeric(trajectory_temp5[[2]]), as.numeric(trajectory_temp2[1])+as.numeric(trajectory_temp3[[1]])/scale_cir, as.numeric(trajectory_temp3[[2]]), col= colbar[j_1+1], lwd = 0.7)
                  }
              }
          }
      }
  }

  x4 = data.frame(col_year,col_membership, col_size,col_average_age,col_normalized_size,col_keyword, col_newkeyword, col_uspc)#

  write.table(x4, file = "C:/Rwork/col_trajectory.csv", sep = ",", row.names = FALSE, col.names = TRUE, qmethod = "double")
}
 