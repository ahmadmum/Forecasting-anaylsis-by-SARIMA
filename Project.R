
##########    Project for " forcasting of Air Temperature by Sarima "


##########    library

library(forecast)
library(zoo)
library(chron)
library(hydroTSM)
library("tseries")


##########     Setting working directory

setwd("C:/Users/mumta/Desktop/Time Series Analysis/Assignment 2/Input data")


##########     Setting the class


setClass("Water_Ts",representation(name="character",unit="character", rawdata="data.frame",timeseries="zoo"))

tl = new("Water_Ts"); tl@name="air temperature"      ; tl@unit = "[°C]"

##########     reading csv files

tl@rawdata <- read.csv(file = "tl.txt",sep="\t", header=F,dec=",",col.names=c("date","value","no data"), stringsAsFactors=F) 

##########     eliminate duplicated values

tl@rawdata = tl@rawdata[!duplicated(tl@rawdata[,1]),]



##########     convert to zoo object incl. converting first column string to time values 

tl@timeseries <- read.zoo(tl@rawdata[,1:2],format="%d.%m.%Y %H:%M",tz="UTC")



##########     plotting whole data

plot(tl@timeseries, main="    Air Temp for 6 years (1996-2001)")





##########     extract year

extract_year <- function(timeseries,year)
{
  # extract year
  ts=window(timeseries,start=as.POSIXct(sprintf("%4d-01-01 00:00:00",year),"UTC"),
            end  =as.POSIXct(sprintf("%4d-12-31 23:59:59",year),"UTC")) 
  return(ts)
}


##########     print MetaData of a time series 


printMetaData <- function(water_ts)
{
  print(paste(water_ts@name,": values",length(water_ts@timeseries),"NaN:",sum(is.na(water_ts@timeseries)),
              sprintf("%.1f%%",(sum(is.na(water_ts@timeseries)/length(water_ts@timeseries)*100)))))
  print(sprintf("                      min  %6.2f %s",min(water_ts@timeseries,na.rm=TRUE),water_ts@unit))
  print(sprintf("                      max  %6.2f %s",max(water_ts@timeseries,na.rm=TRUE),water_ts@unit))
  print(sprintf("                      mean %6.2f %s",mean(water_ts@timeseries,na.rm=TRUE),water_ts@unit))
  print(sprintf("                      sd   %6.2f",sd(water_ts@timeseries,na.rm=TRUE))) 
  
  #
  # annual values
  #
  title_text="                      year "
  na_text   ="                      NA   "
  na2_text  ="                      NA % "
  min_text  ="                      min  "
  max_text  ="                      max  "
  mean_text ="                      mean "
  sd_text   ="                      sd   "
  #
  # loop on all years
  #
  for(i in 1:6) 
  {
    year_time_series = extract_year(water_ts@timeseries,1995+i)
    title_text = paste(title_text,sprintf("  %4d",(1995+i)))
    na_text   = paste(na_text,sprintf("%6d",sum(is.na(year_time_series))))
    na2_text  = paste(na2_text,sprintf("%6.2f",(100.-sum(!is.na(year_time_series))*100/(365*24*6))))
    min_text  = paste(min_text,sprintf("%6.2f",min(year_time_series,na.rm=TRUE)))
    max_text  = paste(max_text,sprintf("%6.2f",max(year_time_series,na.rm=TRUE)))
    mean_text = paste(mean_text,sprintf("%6.2f",mean(year_time_series,na.rm=TRUE)))
    sd_text   = paste(sd_text,sprintf("%6.2f",sd(year_time_series,na.rm=TRUE)))
  }
  # print the results line by line 
  print(title_text)
  print(na_text)
  print(na2_text)
  print(min_text)
  print(max_text)
  print(mean_text)
  print(sd_text)
}

printMetaData(tl)


#( mising value is 1.8 % only ,so it will be ignore )



##########     filter values in a time window for min max


filterTimeWindow <- function(time_series,start,end,minimum,maximum)
{
  is = which(time(time_series)==start)
  ie = which(time(time_series)==end)
  print(paste("start",is," ",start," ",time(time_series)[is]))
  print(paste("end",ie," ",end," ",time(time_series)[ie]))
  for(it in is:ie)
  {
    value = coredata(time_series)[it]
    print(paste(it,value,time(time_series)[ie]))
    if(is.na(value))next;
    if(value >= maximum) {coredata(time_series)[it] <- NA}
    if(value <= minimum) {coredata(time_series)[it] <- NA}  
  }
  return(time_series)  
}

##########     function for preprocessing


tl_pre_processing <- function(year)
{
  # extract year 
  tl=window(tl@timeseries,start=as.POSIXct(sprintf("%4d-01-01 00:00:00",year),"UTC"),
            end  =as.POSIXct(sprintf("%4d-12-31 23:59:59",year),"UTC")) 
  # transform to regular time series zoo object
  tl = izoo2rzoo(tl,date.fmt="%Y-%m-%d %H:%M:%S",tstep="10 min") # hydroTSM !!
 
  return(tl)
}


##########     tl air temperature for diff year


tl_1996=tl_pre_processing(1996)
tl_1997=tl_pre_processing(1997)
tl_1998=tl_pre_processing(1998)
tl_1999=tl_pre_processing(1999)
tl_2000=tl_pre_processing(2000)
tl_2001=tl_pre_processing(2001)


tl6 = rbind(tl_1996,tl_1997,tl_1998,tl_1999,tl_2000,tl_2001)

#plot(tl6)

##########     daily data

#tl6d = aggregate(tl6,as.Date(time(tl6)),FUN=mean,na.rm=TRUE)
#plot(tl6d)

##########     monthly data


tl6m = aggregate(tl6,as.yearmon(time(tl6)),FUN=mean,na.rm=TRUE)

#plot(tl6m)

##########     Forcasting of air temp

##auto-regressive (AR) and moving average (MA) series using ACF and PACF plots.

#AR I MA
#p  d  q

##########     test for stationarity-first visual check for stationary variance

#augmented dickey--Fuller test

adf.test(tl6m,alternative  = "stationary")

##kpss test
#kpss.test(tl6m)

##########     ACF and PACF plot

#ACF plots show correlation b/w a series and its lag
#PACF plots display correlation b/w and its lags that explained by previous lags

par(mfcol=c(1,2))
acf(tl6m) 
pacf(tl6m)

##########     decomposition of the data - take seasonality , trend and cycle into account


decom<-decompose(ts(tl6m, frequency = 12))
autoplot(decom)

##other method
#decom<- stl(tl6m,s.window="periodic")
#autoplot(decom<- stl(tl6m,s.window="periodic"))


##########     removeing seasonality

deseasonal_tl6m <-seasadj(decom)
autoplot(deseasonal_tl6m)



##########     Difference of 1 is sufficent(d=1)


diffdesea_tl6m = diff(deseasonal_tl6m, differences = 1)
plot(diffdesea_tl6m)

#adf testdifferenced series
adf.test(diffdesea_tl6m,alternative  = "stationary")


#ACF and PACF for differenced series
par(mfcol=c(1,2))
acf(diffdesea_tl6m, main= 'ACF for Differenced Series')
pacf(diffdesea_tl6m, main= 'PACF for Differenced Series')


tsdisplay(diffdesea_tl6m,lag.max = 40,main= 'ACF and PACF for Differenced Series' )



##########     fitting an arima model

#get auto fit p,d,q values

auto.arima(deseasonal_tl6m, seasonal = FALSE)


##########     further testing analysis of other models


arima(deseasonal_tl6m,order= c(1,1,1))
arima(deseasonal_tl6m,order= c(1,1,1))
arima(deseasonal_tl6m,order= c(2,1,1))
arima(deseasonal_tl6m,order= c(0,1,4))


arima(deseasonal_tl6m,order= c(0,1,1))
fittlseas<-arima(deseasonal_tl6m,order= c(0,1,1))
seas_fcast<- forecast(fittlseas, h= 24)
plot(seas_fcast)




##########     checking for best model

##smallest AICc value. That is so our best model.

fit_arima<-auto.arima(tl6m,d=1,D=1,stepwise = FALSE,approximation = FALSE,trace = TRUE)

print(summary(fit_arima))

checkresiduals(fit_arima)

forecast(fit_arima)
autoplot(forecast(fit_arima),h=24)







##########     evaluate and iterate 

fit<-auto.arima(deseasonal_tl6m, seasonal = FALSE)
tsdisplay(residuals(fit),lag.max = 40, main = '(1,1,1) Model Residuals')


#graph show lag at 4, so modify  model for p or q=4

arima(deseasonal_tl6m, order = c(0,1,4))
fit2 = arima(deseasonal_tl6m, order = c(0,1,4))
tsdisplay(residuals(fit2),lag.max = 20, main = 'Seasonal Model Residual')

fcast <- forecast(fit2, h=24)
autoplot(fcast)

##########     test model performance with a holdout set

##hold<- window(ts(deseasonal_tl6m),start= 150)
#fit_no_holdout = arima(ts(deseasonal_tl6m[-c(150:172)]),order=c(0,1,4))
#fcast_no_holdout<- forecast(fit_no_holdout, h=24)
#plot(fcast_no_holdout, main="")





##########     models needs seasonality


auto.arima(deseasonal_tl6m, seasonal = TRUE)
fit_tl_seasonality = auto.arima(deseasonal_tl6m, seasonal = TRUE)
seas_fcast<- forecast(fit_tl_seasonality, h= 24)
autoplot(seas_fcast, main = "                                                      Forecasts from ARIMA(0,1,1)(0,1,1)[12] " )


##########     further testing analysis of our models

arima(deseasonal_tl6m,order= c(0,1,1))
fittlseas<-arima(deseasonal_tl6m,order= c(0,1,1))
seas_fcast<- forecast(fittlseas, h= 24)
plot(seas_fcast)


arima(deseasonal_tl6m,order= c(1,1,1))
arima(deseasonal_tl6m,order= c(2,1,1))
arima(deseasonal_tl6m,order= c(0,1,4))

fittlseas<-arima(deseasonal_tl6m,order= c(0,1,1))
seas_fcast<- forecast(fittlseas, h= 24)
plot(seas_fcast)


##2nd method autoarima






##auto forecasting 

##The forecast package includes the auto.arima() function. The function uses a variation of the
##Hyndman and Khandakar algorithm which combines unit root tests, minimization of the AICc and 
##maximum likelihood estimation (MLE) to obtain an optimized ARIMA model. The auto.arima() generates
#a set of optimal (p, d, q) that optimizes model fit criteria by searching through combinations of order parameters.

auto.arima1<-auto.arima(tl6m)
forecast1<-forecast(auto.arima1,h=24)
autoplot(forecast1)
summary(auto.arima1)

## plot the residuals over time to see congrience or variance

plot(forecast1$residuals)

acf(forecast1$residuals)
pacf(forecast1$residuals)





