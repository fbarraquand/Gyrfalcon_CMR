### Reconstructing average weather for NE Iceland with data up to end 2019/beginning 2020
### FBarraquand 21/04/2020 -- from older code ExploWeatherIceland 

rm(list=ls())
graphics.off()

DB=read.table("Stod_495_Grimsstadir.ManMedal.txt",header=T,encoding = "latin1")
head(DB)
ndb=names(DB)
names(DB)[1:3]=c("site","year","month")
names(DB)
### Other stations
DB2=read.table("Stod_422_Akureyri.ManMedal.txt",header=T,encoding = "latin1")
DB3=read.table("Stod_479_Manarbakki.ManMedal.txt",header=T,encoding = "latin1")
names(DB2)[1:3]=c("site","year","month")
names(DB3)[1:3]=c("site","year","month")
DB$date=paste(DB$year,DB$month,15,sep="-")
DB2$date=paste(DB2$year,DB2$month,15,sep="-")
DB3$date=paste(DB3$year,DB3$month,15,sep="-")

### Adding NAs to Manarbakki and Grimsstadir -> Akureyri =  most up to date file
dummy = DB2[DB2$year%in%c(2017,2018,2019,2020),]
dummy$site = 479
dummy[,4:18] = NA
dummy
DB3=rbind(DB3,dummy)
DB3[DB3$year%in%c(2016,2017,2018,2019,2020),]
### Do that for Grimsstadir
dummy = DB2[(DB2$year%in%c(2020))&(DB2$month%in%c(2,3)),]
dummy$site = 495
dummy[,4:18] = NA
dummy
DB=rbind(DB,dummy)
DB[DB$year%in%c(2019,2020),]


######## More data on rainfall
DB4=read.table("Stod_468_Reykjahlid.ManMedal.txt",header=T,encoding = "latin1")
DB5=read.table("Stod_462_Myri.ManMedal.txt",header=T,encoding = "latin1")
DB6=read.table("Stod_473_Stadarholl.ManMedal.txt",header=T,encoding = "latin1")
DB7=read.table("Stod_502_Raufarhofn.ManMedal.txt",header=T,encoding = "latin1") # very coastal, a bit outside. Not used
DB8=read.table("Stod_502_Raufarhofn.ManMedal.txt",header=T,encoding = "latin1")
DB9=read.table("Stod_448_Lerkihlid.ManMedal.txt",header=T,encoding = "latin1")
### 
names(DB4)[1:3]=c("site","year","month")
names(DB5)[1:3]=c("site","year","month")
names(DB6)[1:3]=c("site","year","month")
names(DB7)[1:3]=c("site","year","month")
names(DB8)[1:3]=c("site","year","month")
names(DB9)[1:3]=c("site","year","month")
DB4$date=paste(DB4$year,DB4$month,15,sep="-")
DB5$date=paste(DB5$year,DB5$month,15,sep="-")
DB6$date=paste(DB6$year,DB6$month,15,sep="-")
DB7$date=paste(DB7$year,DB7$month,15,sep="-")
DB8$date=paste(DB8$year,DB8$month,15,sep="-")
DB9$date=paste(DB9$year,DB9$month,15,sep="-")

#### DB4 til 2019/12-- ### We're missing 2018 and 2019 1 to 7
dummy = DB2[(DB2$year%in%c(2020)),]
dummy$site = 468
dummy[,4:18] = NA
dummy
DB4=rbind(DB4,dummy)
DB4[DB4$year%in%c(2019,2020),]


# Loop through years and then months
for (year in 1975:2019){
  for (month in 1:12){
  ### Take averages of temp and log Rainfall
    
  # Temp
  t1 = DB$t[(DB$year==year)& (DB$month==month)]
  if(length(t1)==0){t1=NA}
  t2 = DB2$t[(DB2$year==year)& (DB2$month==month)]
  if(length(t2)==0){t2=NA}
  t3 = DB3$t[(DB3$year==year)& (DB3$month==month)]
  if(length(t3)==0){t3=NA}
  t=c(t1,t2,t3)
  av.temp = mean(t,na.rm = T)
  
  # Rainfall 
  lr1 = log(DB$r[(DB$year==year)& (DB$month==month)])
  if(length(lr1)==0){lr1=NA}
  lr2 = log(DB3$r[(DB3$year==year)& (DB3$month==month)])
  if(length(lr2)==0){lr2=NA}
  lr3 = log(DB4$r[(DB4$year==year)& (DB4$month==month)])
  if(length(lr3)==0){lr3=NA}
  lr4 = log(DB5$r[(DB5$year==year)& (DB5$month==month)])
  if(length(lr4)==0){lr4=NA}
  lr5 = log(DB6$r[(DB6$year==year)& (DB6$month==month)])
  if(length(lr5)==0){lr5=NA}
  lr6 = log(DB9$r[(DB9$year==year)& (DB9$month==month)])
  if(length(lr6)==0){lr6=NA}
  # Removes numeric(0) and transform in NAs
  lr=c(lr1,lr2,lr3,lr4,lr5,lr6)
  av.logRainfall = mean(lr,na.rm = T)
  date=paste(year,month,15,sep="-")
  ### Construct data structure here
  if ((year== 1975)&(month==1)){
  av.weather=data.frame(date,year,month,av.temp,av.logRainfall)
  }
  else {# add to data structure
  new_row = data.frame(date,year,month,av.temp,av.logRainfall)
  av.weather=rbind(av.weather,new_row)  
    }
  }
}

#av.weather$av.temp[is.nan(av.weather$av.temp)]=NA
#av.weather$av.logRainfall[is.nan(av.weather$av.logRainfall)]=NA

# Verif
year=1987
month=4
lr1 = log(DB$r[(DB$year==year)& (DB$month==month)])
if(length(lr1)==0){lr1=NA}
lr2 = log(DB3$r[(DB3$year==year)& (DB3$month==month)])
if(length(lr2)==0){lr2=NA}
lr3 = log(DB4$r[(DB4$year==year)& (DB4$month==month)])
if(length(lr3)==0){lr3=NA}
lr4 = log(DB5$r[(DB5$year==year)& (DB5$month==month)])
if(length(lr4)==0){lr4=NA}
lr5 = log(DB6$r[(DB6$year==year)& (DB6$month==month)])
if(length(lr5)==0){lr5=NA}
lr6 = log(DB9$r[(DB9$year==year)& (DB9$month==month)])
if(length(lr6)==0){lr6=NA}
# Removes numeric(0) and transform in NAs
lr=c(lr1,lr2,lr3,lr4,lr5,lr6)
av.logRainfall = mean(lr,na.rm = T)
# seems OK

names(av.weather)[4:5]=c("temp","logRainfall")
head(av.weather)
write.csv(av.weather,"average_weatherNEIceland.csv")

### Plots
pdf("logRainfall_NEIceland_May_withAverage.pdf",width = 14,height=8)
plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on May)",col="cyan",)
legend("bottomright",c("Grimsstadir","Manarbakki","Reykjahlid","Myri","Stadarholl","Lerkihlid","Average"),pch=21,pt.bg=c("cyan","blue","red","pink","orange","green","black"))
#lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)
lines(as.Date(DB4$date[DB4$year>1980]),log(DB4$r[DB4$year>1980]),type="o",col="red",lty=2)
lines(as.Date(DB5$date[DB5$year>1980]),log(DB5$r[DB5$year>1980]),type="o",col="pink",lty=2)
lines(as.Date(DB6$date[DB6$year>1980]),log(DB6$r[DB6$year>1980]),type="o",col="orange",lty=2)
lines(as.Date(DB9$date[DB9$year>1980]),log(DB9$r[DB9$year>1980]),type="o",col="green",lty=2)
lines(as.Date(av.weather$date),av.weather$logRainfall,type="o",col="black",lty=2)

### Add points for May weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),log(DB$r[DB$year>1980& (DB$month==5)]),col="cyan",pch=16,lwd=3)
#lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),log(DB2$r[DB2$year>1980 & (DB2$month==5)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),log(DB3$r[DB3$year>1980 & (DB3$month==5)]),col="blue",pch=16,lwd=3)
lines(as.Date(DB4$date[DB4$year>1980 & (DB4$month==5)]),log(DB4$r[DB4$year>1980 & (DB4$month==5)]),col="red",pch=16,lwd=3)
lines(as.Date(DB5$date[DB5$year>1980 & (DB5$month==5)]),log(DB5$r[DB5$year>1980 & (DB5$month==5)]),col="pink",pch=16,lwd=3)
lines(as.Date(DB6$date[DB6$year>1980 & (DB6$month==5)]),log(DB6$r[DB6$year>1980 & (DB6$month==5)]),col="orange",pch=16,lwd=3)
lines(as.Date(DB9$date[DB9$year>1980 & (DB9$month==5)]),log(DB9$r[DB9$year>1980 & (DB9$month==5)]),col="green",pch=16,lwd=3)
lines(as.Date(av.weather$date[av.weather$month==5]),av.weather$logRainfall[av.weather$month==5],type="o",col="black",lwd=3)
dev.off()


pdf("Temp_NEIceland_May_withAverage.pdf",width = 14,height=8)
### Same with temp
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on May)",col="cyan",)
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki","Average"),pch=21,pt.bg=c("cyan","green","blue","black"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)
lines(as.Date(av.weather$date),av.weather$logRainfall,type="o",col="black",lty=2)

### Add points for May weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),DB$t[DB$year>1980& (DB$month==5)],col="cyan",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),DB2$t[DB2$year>1980 & (DB2$month==5)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),DB3$t[DB3$year>1980 & (DB3$month==5)],col="blue",pch=16,lwd=3)
lines(as.Date(av.weather$date[av.weather$month==5]),av.weather$temp[av.weather$month==5],type="o",col="black",lwd=3)
dev.off()
