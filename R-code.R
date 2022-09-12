setwd("")
data<- read.csv("OTE.CZ.price.csv")

id.select<- 2 
price<- data[,id.select] 
price.names<- names(data)[id.select]
time.utc<- as.POSIXct( strptime(data[,1], format="%Y-%m-%d %H:%M:%S", tz="UTC"))
time.numeric<- as.numeric(time.utc) 
local.time.zone<- "CET"
time.lt<- as.POSIXct( time.numeric, origin = "1970-01-01", tz=local.time.zone) 

#date/time definition
S<- 24*60*60/ as.numeric(names(which.max(table(diff(time.numeric)))))
start.end.time.S<- strptime(format(c(time.lt[1], time.lt[length(time.lt)])),format="%Y-%m-%d %H:%M:%S", tz="UTC")
time.S.numeric<- seq(as.numeric(start.end.time.S[1]), as.numeric(start.end.time.S[length(start.end.time.S)]), by=24*60*60/S)
time.S<- as.POSIXct( time.S.numeric, origin = "1970-01-01", tz="UTC") 
dates.S<- unique(as.Date(time.S, tz=local.time.zone))

#this loop converts time zone; particularly cirtical for electricity price forecasting 
source("DST.trafo.R")

price.S.out<- DST.trafo(X=price, Xtime=time.utc, Xtz=local.time.zone) 
price.S<- price.S.out[,,1]



####SECTION 1: NAIVE AND EXPERT MODELS###


#
forecast.naive<- function(Y, days){
	weekday <- format(days[length(days)] +1, "%a") ##weekday of next day (day to forecast)
	if(weekday %in% c("Mon", "Sat", "Sun") ) {
		forecast<- Y[dim(Y)[1]-7+1,]## values a week ago
	} else {
		forecast<- Y[dim(Y)[1]-1+1,] ## last day values
	}
	## compare with strptime
	forecast
}

## The extended naive model is constructed according to the Task given on the page 12 of the script
forecast.naive.extended =function(Y,days){
	weekday <- format(days[length(days)] +1, "%a") ##weekday of next day (day to forecast)
	if(weekday %in% c("Sat", "Sun") ) {
		forecast<- Y[dim(Y)[1]-7+1,]## values a week ago
	} else if (weekday %in% c("Fri")){
		forecast = 0.5*(Y[dim(Y)[1]-1+1,] + Y[dim(Y)[1]-7+1,])
	} else if (weekday %in% c("Mon")){
		forecast = 0.5*(Y[dim(Y)[1]-3+1,] + Y[dim(Y)[1]-7+1,])
	} else {
		forecast<- Y[dim(Y)[1]-1+1,] 
	}
	forecast
}


model_estimation=function(){   # This function looks for optimal lags (i.e. autoregressive days) 
	D=365*2+1
	index=1:D

	combinations=expand.grid(as.data.frame(matrix(0:1,2,7))) # define all possible combinations of a vector (1:7) 
	combinations=combinations[-1,] #Neglect a first row of the matrix, which is a bunch of zeros 
	daysvec=c(1,2,3,4,5,6,7)
	for (i in 1:nrow(combinations)){ #replace 1's obtained by the expand.grid command with a respective day number 
		for (j in 1:ncol(combinations)){
			if (combinations[i,j]!=0){
				combinations[i,j]= daysvec[j]
			} else {
			combinations[i,j]=NA}
			}
		}

	## model_order_estimation
	days=dates.S[index]
	BIC_current=Inf #Initial setting of the BIC criterion 
	Y=price.S[index,]
	S<- dim(Y)[2]
	days.ext<- c(days, days[length(days)]+1) 
	weekdays.num<- as.numeric(format(days.ext, "%w")) 
	WD<- t( sapply(weekdays.num , "==", c(1:6,0)) ) + 0 
	dimnames(WD)<- list(NULL, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) #
	get.lagged<- function(lag, Z) c( rep(NA, lag),  Z[(1+lag):(length(Z)+1) -lag ]  ) 
	## expert specification
	for (i in 1:nrow(combinations)){ #each regression model conducted in the function will thus have a unique expert.wd and expert.lags combination. 
		for (j in 1:nrow(combinations)){ # 
			expert.wd=combinations[i,]
			expert.wd=expert.wd[!is.na(expert.wd)]
			#
			expert.lags=combinations[j,]
			expert.lags=expert.lags[!is.na(expert.lags)]
			#		
			for(s in 1:S){
				XLAG<- sapply(expert.lags, get.lagged, Z=price.S[index,][,s] ) 
				dimnames(XLAG)<- list(NULL, paste("lag", expert.lags))			
				XREG<- cbind(1, WD[,expert.wd], XLAG)
				YREG<- Y[,s]
				act.index<- which( !apply(is.na(XREG),1, any)[-dim(XREG)[1]] ) # without NA's and last
				model<- lm.fit(XREG[act.index,], YREG[act.index])
				BIC_model= length(model$res)*(log(2*pi)+log((t(model$res)%*%(model$res)) /length(model$res))+1) + log(D)*(ncol(XREG)+1)
																		# Since 'lm.fit' coommand is recalled, the BIC criterion should be 
																		# calculated manually. The formula aboves is equivalent to the one  
																		# used by the build-in command 'BIC(model)'. 
				if (BIC_model[1,1] < BIC_current){
					current_model<<-model #Note, that this variable and two below are global
					current_expert_wd<<-expert.wd 
					current_expert_lags<<-expert.lags 
					BIC_current=BIC_model 
				} 
			}
		}
	}
}


model_estimation() #The function is turned off for the time being. Its execution requires my PC around 15 minutes for processing.  

current_expert_wd=c(1,5,6,7)
current_expert_lags=c(1,3,5,7) 


###The function below combines two expert forecasts: the first is made under the conventional combination of the lags and the week day dummies, the second
###is made under the combination of current_expert_wd and current_expert_lags vectors.

## forecast function
forecast.expert<- function(Y, days){
	forecast<- numeric()
	forecast_list = list()
	S<- dim(Y)[2]
	days.ext<- c(days, days[length(days)]+1) 
	weekdays.num<- as.numeric(format(days.ext, "%w")) 
	WD<- t( sapply(weekdays.num , "==", c(1:6,0)) ) + 0 
	dimnames(WD)<- list(NULL, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) #
	get.lagged<- function(lag, Z) c( rep(NA, lag),  Z[(1+lag):(length(Z)+1) -lag ]  ) 
	## expert specification
		wd.spec = list(c(1,6,7), current_expert_wd) #a list of two possible options 
		lags.spec = list(c(1,2,7), current_expert_lags) 
		for (i in 1:2){ 		#Clearly, this "selection" loop must be exectuted before the loop responsible for running regressions 
			expert.wd=wd.spec[[i]][1:length(wd.spec[[i]])]
			expert.lags=lags.spec[[i]][1:length(lags.spec[[i]])]		
				for(s in 1:S){
					XLAG<- sapply(expert.lags, get.lagged, Z=Y[,s] ) 
					dimnames(XLAG)<- list(NULL, paste("lag", expert.lags))
					XREG<- cbind(1, WD[,expert.wd], XLAG)
					YREG<- Y[,s]
					act.index<- which( !apply(is.na(XREG),1, any)[-dim(XREG)[1]] ) 
					model<- lm.fit(XREG[act.index, ], YREG[act.index] )
					forecast[s]<- model$coef %*% XREG[dim(XREG)[1],]
					}
			forecast_list[[i]]=forecast # As the result, two forecasts are stored under the list type variable
			}
		forecast_list
	}

D<- 365*2+1
N<- 365*4 #
oos.dates<- dates.S[D+1:N] #
model.names<- c("naive", "naive_extended", "expert", "expert_par_estimated")
M<- length(model.names) ## the number of models under consideration is set to adjust automatically

ERROR<- array(, dim=c(N, S, M))
dimnames(ERROR)<- list(oos.dates, 1:S, model.names)
for(i.N in 1:N){
	index<- 1:D + i.N -1  ## defining index set for input data
	## model 1
	forecast_naive<- forecast.naive(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 1] <- forecast_naive - price.S[index[length(index)]+1,]
	##model 2
	forecast_naive_extended<- forecast.naive.extended(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 2] <- forecast_naive_extended - price.S[index[length(index)]+1,]	
	## model 3
	forecast_expert<- forecast.expert(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 3] <- forecast_expert[[1]] - price.S[index[length(index)]+1,] #Remember, that the function forecast.expert produces a list of two forecasts
	ERROR[i.N,, 4] <- forecast_expert[[2]] - price.S[index[length(index)]+1,]
}#i.N

## forecasting:
MAE<- apply( abs(ERROR) ,3, mean)
RMSE<- sqrt(apply( ERROR^2 ,3, mean))
MAE
RMSE

## DM-test function
DM.test<- function(E1, E2, hmax=1, power=1){
## as dm.test with alternative == "less"
	L1<- apply( abs(as.matrix(E1))^power, 1, sum)^(1/power)
	L2<- apply( abs(as.matrix(E2))^power, 1, sum)^(1/power)
	delta <- L1 - L2
	delta.var <- var(delta)/length(delta)## estimation of the variance 
	STATISTIC <- mean(delta, na.rm = TRUE)/sqrt(delta.var)
	delta.len <- length(delta)
	k <- ((delta.len + 1 - 2 * hmax + (hmax/delta.len) * (hmax - 1))/delta.len)^(1/2)
	STATISTIC <- STATISTIC * k
	PVAL <- pt(STATISTIC, df = delta.len - 1)
	list(stat=STATISTIC, p.val=PVAL)
}


#The function below simplifies the application of the DM-test. In doing so, the function uses seemingly non-trivial method and works well for 
#	the pairwise comparison of the obtained models, among which it selects the fittest one. 

compare_models=function(){
index_vec=(1:dim(ERROR)[3]) 
i=0
while (i < (length(index_vec)*10)){	# The number of iterations is set to be excessive to ensure that no models remained aside of the comparison  
	k=sample(index_vec,1) # Randomly select one model  
	j=sample(index_vec,1) # Randomly select another model
		if (j!=k & k>j){ # Conditions are necessary for clear interpretation of the outcome 
			if (DM.test(ERROR[,,k], ERROR[,,j])$p.val<0.05){ # Since k is surely greater than j, the model k is significantly better than model j if p-value is small
				index_vec=index_vec[! index_vec %in% j] # The less fitting model j is thus excluded from the comparison pool, whereas the model k is retained
				i= i+ 1
			}else if (DM.test(ERROR[,,k], ERROR[,,j])$p.val>0.05 & DM.test(ERROR[,,j], ERROR[,,k])$p.val>0.05) {
				print(paste("The model", model.names[j], "is equivalent to the model", model.names[k]))	
				index_vec=index_vec[! index_vec %in% k] # if one model is not significantly better than another, the more complex model is eliminated
				i= i+ 1
			}else if (DM.test(ERROR[,,k], ERROR[,,j])$p.val>0.05){
				print(paste("The model", model.names[k], "is not significantly better than the model", model.names[j])) 
				index_vec=index_vec[! index_vec %in% k]
				i= i+ 1
				}
			}
		}
print(paste("The best model is", model.names[index_vec]))} #The final conclusion as to the model selection is printed in the console

compare_models()


####SECTION 2: EVALUATING THE CORRELATION STRUCTURE OF THE DATA####

get.cpacf<- function(Y, k=1){ 
	CPACF<- matrix(, S,S) 
	for(s in 1:S){
		for(l in 1:S){
			Y.s<- Y[(k+1):length(Y[,s])  , s] 
			Y.l.lagged<- Y[(k+1):length(Y[,s]) -k , l] 
			CPACF[s, l]<- cor( Y.s, Y.l.lagged )
		}
	}
	CPACF
}

k<- 1
cpacf<- get.cpacf( tail( price.S, 5*D ), k=1)

###FIGURE 1 ###

library(plotrix) 

#pdf(paste("cpacff_k=",k,".pdf", sep=""), width=8, height=6)
par(mar=c(4.4, 4.4, 1, 4), family="Times") ## plotting area and font
cred<- 		c(.5,0,0,0,1,1,1) ## red
cgreen<- 	c(.5,0,1,1,1,0,0) ## green
cblue<- 	c(.5,1,1,0,0,0,1) ## blue
crange<- c(-1,1) ## range
cpacftmp= cpacf
cpacftmp2= cpacf
cpacftmp[1,1]<- 0

COLplot<- color.scale(cpacf,cred, cgreen, cblue , xrange=crange ) 
COLleg<- color.scale( seq(crange[1],crange[2], length=100), cred, cgreen, cblue , alpha=.5)
color2D.matplot(cpacf, cellcolors=COLplot, show.values=2, vcol=rgb(0,0,0), vcex=.6, axes=FALSE, xlab="l", ylab="s")
axis(1,1:S-.5,1:S, cex.axis=.8) ## draw x-axis
axis(2, S:1-.5,1:S, cex.axis=.8, las=2) ## draw y-axis
lableg <- formatC(crange[1] + (crange[2]-crange[1])*seq(0,1,length=11), format="f", digits=2) 
pardat <- par()
color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), COLleg, align="rb", gradient="y")
#dev.off() ## closes/saves pdf 


pcor<- function(X,Y, Z) cor(lm.fit(Z,Y)$res,lm.fit(Z, X)$res)#


MAT<- matrix(, S, S)
for(s in 1:S){
	for(l in 1:S){
		X<- price.S[index,s]
		Y<- price.S[index-1,l] 
		Z<- cbind(price.S[index-1,s]) 
		MAT[s, l] <- pcor(X,Y, Z)  
	}
}
	

###FIGURE 2###

#pdf( paste("MAT_k=",k,".pdf", sep=""), width=8, height=6)
par(mar=c(4.4, 4.4, 1, 4), family="Times") ## plotting area and font
cred<- 		c(.5,0,0,0,1,1,1) ## red
cgreen<- 	c(.5,0,1,1,1,0,0) ## green
cblue<- 	c(.5,1,1,0,0,0,1) ## blue
crange<- c(-1,1.000001) ## range
COLplot<- color.scale(MAT,cred, cgreen, cblue , xrange=crange ) 
COLleg<- color.scale( seq(crange[1],crange[2], length=100), cred, cgreen, cblue , alpha=.5)
color2D.matplot(MAT, cellcolors=COLplot, show.values=2, vcol=rgb(0,0,0), vcex=.6, axes=FALSE, xlab="l", ylab="s")
axis(1,1:S-.5,1:S, cex.axis=.8) ## draw x-axis
axis(2, S:1-.5,1:S, cex.axis=.8, las=2) ## draw y-axis
lableg <- formatC(crange[1] + (crange[2]-crange[1])*seq(0,1,length=11), format="f", digits=2) 
pardat <- par()
color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), COLleg, align="rb", gradient="y")

#dev.off() ## closes/saves pdf 


 
MAT2<- matrix(, S, S)
for(s in 1:S){
	for(l in 1:S){
		X<- price.S[index,s]
		Y<- price.S[index-1,l]
		Z<- cbind(price.S[index-1,S]) 
		MAT2[s, l] <- pcor(X,Y, Z)
	}
}

###FIGURE 3###


#pdf( paste("MAT2_k=",k,".pdf", sep=""), width=8, height=6)

par(mar=c(4.4, 4.4, 1, 4), family="Times") ## plotting area and font
cred<- 		c(.5,0,0,0,1,1,1) ## red
cgreen<- 	c(.5,0,1,1,1,0,0) ## green
cblue<- 	c(.5,1,1,0,0,0,1) ## blue
crange<- c(-1,1.000001) ## range
COLplot<- color.scale(MAT2,cred, cgreen, cblue , xrange=crange ) 
COLleg<- color.scale( seq(crange[1],crange[2], length=100), cred, cgreen, cblue , alpha=.5)
color2D.matplot(MAT2, cellcolors=COLplot, show.values=2, vcol=rgb(0,0,0), vcex=.6, axes=FALSE, xlab="l", ylab="s")
axis(1,1:S-.5,1:S, cex.axis=.8) ## draw x-axis
axis(2, S:1-.5,1:S, cex.axis=.8, las=2) ## draw y-axis
lableg <- formatC(crange[1] + (crange[2]-crange[1])*seq(0,1,length=11), format="f", digits=2) 
pardat <- par()
color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), COLleg, align="rb", gradient="y")

#dev.off() ## closes/saves pdf 


MAT3<- matrix(, S, S)
for(s in 1:S){
	for(l in 1:S){
		X<- price.S[index,s]
		Y<- price.S[index-1,l]
		Z<- cbind(price.S[index-1,S], price.S[index-1,s]) 
		MAT3[s, l] <- pcor(X,Y, Z)
	}
}


###FIGURE 4###


#pdf( paste("MAT3_k=",k,".pdf", sep=""), width=8, height=6)

par(mar=c(4.4, 4.4, 1, 4), family="Times") ## plotting area and font
cred<- 		c(.5,0,0,0,1,1,1) ## red
cgreen<- 	c(.5,0,1,1,1,0,0) ## green
cblue<- 	c(.5,1,1,0,0,0,1) ## blue
crange<- c(-1,1.000001) ## range
COLplot<- color.scale(MAT3,cred, cgreen, cblue , xrange=crange ) 
COLleg<- color.scale( seq(crange[1],crange[2], length=100), cred, cgreen, cblue , alpha=.5)
color2D.matplot(MAT3, cellcolors=COLplot, show.values=2, vcol=rgb(0,0,0), vcex=.6, axes=FALSE, xlab="l", ylab="s")
axis(1,1:S-.5,1:S, cex.axis=.8) ## draw x-axis
axis(2, S:1-.5,1:S, cex.axis=.8, las=2) ## draw y-axis
lableg <- formatC(crange[1] + (crange[2]-crange[1])*seq(0,1,length=11), format="f", digits=2) 
pardat <- par()
color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), COLleg, align="rb", gradient="y")

#dev.off() ## closes/saves pdf 


MAT4<- matrix(, S, S)
for(s in 1:S){
	for(l in 1:S){
		X<- price.S[index,s]
		Y<- price.S[index-1,l]
		Z<- cbind(price.S[index-1,S], price.S[index-1,7], price.S[index-1,21], price.S[index-1,s]) 
		MAT4[s, l] <- pcor(X,Y, Z)
	}
}


###FIGURE 5###


#pdf( paste("MAT4_k=",k,".pdf", sep=""), width=8, height=6)

par(mar=c(4.4, 4.4, 1, 4), family="Times") ## plotting area and font
cred<- 		c(.5,0,0,0,1,1,1) ## red
cgreen<- 	c(.5,0,1,1,1,0,0) ## green
cblue<- 	c(.5,1,1,0,0,0,1) ## blue
crange<- c(-1,1.000001) ## range
COLplot<- color.scale(MAT4,cred, cgreen, cblue , xrange=crange ) 
COLleg<- color.scale( seq(crange[1],crange[2], length=100), cred, cgreen, cblue , alpha=.5)
color2D.matplot(MAT4, cellcolors=COLplot, show.values=2, vcol=rgb(0,0,0), vcex=.6, axes=FALSE, xlab="l", ylab="s")
axis(1,1:S-.5,1:S, cex.axis=.8) ## draw x-axis
axis(2, S:1-.5,1:S, cex.axis=.8, las=2) ## draw y-axis
lableg <- formatC(crange[1] + (crange[2]-crange[1])*seq(0,1,length=11), format="f", digits=2) 
pardat <- par()
color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), COLleg, align="rb", gradient="y")

#dev.off() ## closes/saves pdf 


#Having observed the graphs, the two regression models are conceived: 

###expert_par_estimated model extended with the previous day's hour 24
forecast.expert_24<- function(Y, days){
	forecast<- numeric()
	forecast_list=list()
	S<- dim(Y)[2]
	days.ext<- c(days, days[length(days)]+1) 
	weekdays.num<- as.numeric(format(days.ext, "%w")) 
	WD<- t( sapply(weekdays.num , "==", c(1:6,0)) ) + 0 
	dimnames(WD)<- list(NULL, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) 
	get.lagged<- function(lag, Z) c( rep(NA, lag),  Z[(1+lag):(length(Z)+1) -lag ]  ) 
	expert.wd<-  current_expert_wd 
	expert.lags<- current_expert_lags 	
	XLAGS <- get.lagged(lag=1, Z=Y[, S]) 
	for(s in 1:S){
		XLAG<- sapply(expert.lags, get.lagged, Z=Y[,s] )
		dimnames(XLAG)<- list(NULL, paste("lag", expert.lags))
		if (s==S){ 		# One out of the two similar regressors is dropped and the matrix X'X is thus almost certain to be invertible  
			XREG<- cbind(1, XLAG, WD[,expert.wd])
		} else {
			XREG<- cbind(1, XLAG, XLAGS, WD[,expert.wd])}
		YREG<- Y[,s]
		act.index<- !apply(is.na(XREG),1, any)[-dim(XREG)[1]] 
		model<- lm.fit(XREG[act.index, ], YREG[act.index] )    
		forecast[s]<- model$coef %*% XREG[dim(XREG)[1],] 
		}
	forecast
}


###expert_par_estimated model extended with the previous day's hours 7, 21, 24
forecast.expert_7_21_24<- function(Y, days){
	forecast<- numeric()
	forecast_list=list()
	S<- dim(Y)[2]
	days.ext<- c(days, days[length(days)]+1) 
	weekdays.num<- as.numeric(format(days.ext, "%w")) #
	WD<- t( sapply(weekdays.num , "==", c(1:6,0)) ) + 0 
	dimnames(WD)<- list(NULL, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) 
	get.lagged<- function(lag, Z) c( rep(NA, lag),  Z[(1+lag):(length(Z)+1) -lag ]  ) 
	expert.wd<-  current_expert_wd #
	expert.lags<- current_expert_wd 
	LAGS=c(7,21,24)
		XLAGS=matrix(,nrow(Y)+1, length(LAGS))	
	for (i in 1:length(LAGS)){
		XLAGS[,i] <- get.lagged(lag=1, Z=Y[, LAGS[i]])} ## 
	for(s in 1:S){
		XLAG<- sapply(expert.lags, get.lagged, Z=Y[,s] )
		dimnames(XLAG)<- list(NULL, paste("lag", expert.lags))
		if (s %in% LAGS){	
			XREG<- cbind(1, XLAG, WD[,expert.wd])
		} else {
			XREG<- cbind(1, XLAG, XLAGS, WD[,expert.wd])}
		YREG<- Y[,s]
		act.index<- !apply(is.na(XREG),1, any)[-dim(XREG)[1]] # without NA's and last
		model<- lm.fit(XREG[act.index, ], YREG[act.index] )
		forecast[s]<- model$coef %*% XREG[dim(XREG)[1],]
		}
	forecast
}


#Forecasting study: 

D<- 365*2+1 
N<- 365*4 #
oos.dates<- dates.S[D+1:N] #

new_models = c("expert_24", "expert_7_21_24)") # For the time being, the algorighm of the forecasting study was modified: 
							     #          the errors are computed only for the models elaborated in the section. 
							     # The modification thus requires the code to be executed in a consecutive order.
model.names=c(model.names, new_models)
M<- length(model.names) #


ERROR<- array(ERROR, dim=c(N, S, M)) # The redefined ERROR array is thus not constructed from scratch, but is spanned from its predecessor 
dimnames(ERROR)<- list(oos.dates, 1:S, model.names)
for(i.N in 1:N){
	index<- 1:D + i.N -1  ## defining index set for input data
	#Model 4
	forecast_expert_LAG24<- forecast.expert_24(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 5] <- forecast_expert_LAG24 - price.S[index[length(index)]+1,]
	#Model 5
	forecast_expert_LAGS<- forecast.expert_7_21_24(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 6] <- forecast_expert_LAGS - price.S[index[length(index)]+1,]
}#i.N


MAE<- apply( abs(ERROR) ,3, mean)
RMSE<- sqrt(apply( ERROR^2 ,3, mean))
MAE
RMSE
 
compare_models() # As was noted, the function allows a rapid and convenient comparison of the elaborated models.




### SECTION 3: SEASONAL EFFECTS ##		

A<- 365.24 ## annual periodicity

Aindex<- 0:floor(2*A+1)
M<- 3
Fourier.basis<- matrix(, length(Aindex), 2*M+1)
Fourier.basis[,1]<- 1
for(i in 1:M) {
	Fourier.basis[, 2*i ]<- sin( Aindex*2*pi / A*i)
	Fourier.basis[, 2*i+1 ]<- cos( Aindex*2*pi / A*i)	
}

library(splines)
##periodic-B-spline:
get.pbas<- function( Bindex , period=24, dK = period/4, ord=4){
	## ord=4 --> cubic splines 
	## dK = equidistance distance
	## S= season
	## support will be 1:n
	n<- length(Bindex)
	stp<- dK
	x<- 1:period ## must be sorted!
	lb<- x[1]
	ub<- x[length(x)]
	knots<- seq(lb-(0)*stp,ub+(1)*stp, by=stp)
	derivs<- numeric(length(x))
	## some stuff adjusted from pbc-package to be faster
	degree<- ord-1
	nKnots = length(knots)
	Aknots = c( knots[1] - knots[nKnots] + knots[nKnots - degree:1] , knots,  knots[nKnots] + knots[1:degree + 1] - knots[1] )
	basisInterior <- splineDesign(Aknots, x, ord, derivs) 
	basisInteriorLeft <- basisInterior[, 1:(ord-1), drop = FALSE]
	basisInteriorRight <- basisInterior[, (ncol(basisInterior) - ord+2):ncol(basisInterior), drop = FALSE]
	basis <- cbind(basisInterior[, -c(1:(ord-1), (ncol(basisInterior) - ord+2):ncol(basisInterior)), drop = FALSE], basisInteriorLeft + basisInteriorRight)
	t(array(t(basis), dim= c(dim(basis)[2] ,n)))
}
K<- 6 ## 
Bspline.basis<- get.pbas(Aindex, period=A, dK= A/K, ord=4)

#The function below is the forecast.expert_7_21_24 extended with the Fourier.basis and Bspline basis regressors. 
#Note, that the regressors are not simultaneously included into the XREG matrix. 

forecast.expert_fourier_bspline<- function(Y, days){
	forecast<- numeric()
	forecast_list=list()
	S<- dim(Y)[2]
	days.ext<- c(days, days[length(days)]+1) 
	weekdays.num<- as.numeric(format(days.ext, "%w")) 
	WD<- t( sapply(weekdays.num , "==", c(1:6,0)) ) + 0 
	dimnames(WD)<- list(NULL, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) 
	get.lagged<- function(lag, Z) c( rep(NA, lag),  Z[(1+lag):(length(Z)+1) -lag ]  ) 
	## expert specification
	expert.wd<-  current_expert_wd
	expert.lags<- current_expert_lags
	LAGS=c(7,21,24)
	XLAGS=matrix(,nrow(Y)+1, length(LAGS))	
	for (i in 1:length(LAGS)){
		XLAGS[,i] <- get.lagged(lag=1, Z=Y[, LAGS[i]])}
	forecast_list=list() ## 
	for (j in 1:2){ # Analogously to the "forecast_expert" function, this one produces a list that contains two different forecasts
		for(s in 1:S){
			XLAG<- sapply(expert.lags, get.lagged, Z=Y[,s] )
			dimnames(XLAG)<- list(NULL, paste("lag", expert.lags))
			if (s %in% LAGS){
				if (j==1){
					XREG<- cbind(1, XLAG, Fourier.basis[,2:ncol(Fourier.basis)], WD[,expert.wd]) #The unitary column of the Fourier.basis is ignored
				} else if (j==2){
					XREG<- cbind(1, XLAG, Bspline.basis[,2:ncol(Bspline.basis)], WD[,expert.wd])} #The first column of the Bsline.basis is neglected as well
			} else {
				if (j==1){
					XREG<- cbind(1, XLAG, XLAGS, Fourier.basis[,2:ncol(Fourier.basis)], WD[,expert.wd])
				} else if (j==2){
					XREG<- cbind(1, XLAG, XLAGS, Bspline.basis[,2:ncol(Bspline.basis)], WD[,expert.wd])}}
			YREG<- Y[,s]
			act.index<- !apply(is.na(XREG),1, any)[-dim(XREG)[1]] # without NA's and last
			model<- lm.fit(XREG[act.index, ], YREG[act.index] )
			forecast[s]<- model$coef %*% XREG[dim(XREG)[1],]
			}
		forecast_list[[j]]=forecast
	}
	forecast_list
}


D<- 365*2+1
N<- 365*4 
oos.dates<- dates.S[D+1:N] 

new_models = c("expert_fourier", "expert_bspline")
model.names=c(model.names, new_models)
M<- length(model.names) 

ERROR<- array(ERROR, dim=c(N, S, M))
dimnames(ERROR)<- list(oos.dates, 1:S, model.names)
for(i.N in 1:N){
	index<- 1:D + i.N -1  ## defining index set for input data
	##Model 6
	forecast_expert_fourier<- forecast.expert_fourier_bspline(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 7] <- forecast_expert_fourier[[1]] - price.S[index[length(index)]+1,]
	forecast_expert_fourier<- forecast.expert_fourier_bspline(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 8] <- forecast_expert_fourier[[2]] - price.S[index[length(index)]+1,]
}#i.N

## forecasting:
MAE<- apply( abs(ERROR) ,3, mean)
RMSE<- sqrt(apply( ERROR^2 ,3, mean))
MAE
RMSE

compare_models()



###SECTION 4: NON-LINEAR EFFECTS###



forecast.expert.minmax<- function(Y, days){ 
	forecast<- numeric()
	forecast_list=list()
	S<- dim(Y)[2]
	days.ext<- c(days, days[length(days)]+1) 
	weekdays.num<- as.numeric(format(days.ext, "%w")) 
	WD<- t( sapply(weekdays.num , "==", c(1:6,0)) ) + 0 
	dimnames(WD)<- list(NULL, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) #
	get.lagged<- function(lag, Z) c( rep(NA, lag),  Z[(1+lag):(length(Z)+1) -lag ]  )
	expert.wd<-  current_expert_wd
	expert.lags<- current_expert_lags
	LAGS=c(7,21,24)
		XLAGS=matrix(,nrow(Y)+1, length(LAGS))	
	for (i in 1:length(LAGS)){
		XLAGS[,i] <- get.lagged(lag=1, Z=Y[, LAGS[i]])} ## 
	Ymin<- apply(Y, 1, min)
	Ymax<- apply(Y, 1, max)
	YminLAG<- sapply(1, get.lagged, Z=Ymax) 
	YmaxLAG<- sapply(1, get.lagged, Z=Ymin) 
	Xnonlin<- cbind(YminLAG, YmaxLAG)
	dimnames(Xnonlin)<- list(NULL, c("min", "max"))
	for(s in 1:S){
		XLAG<- sapply(expert.lags, get.lagged, Z=Y[,s] )
		dimnames(XLAG)<- list(NULL, paste("lag", expert.lags))
		if (s %in% LAGS){
			XREG<- cbind(1, XLAG, Xnonlin, WD[,expert.wd])
		} else {
			XREG<- cbind(1, XLAG, XLAGS,Xnonlin, WD[,expert.wd])}
		XREG
		YREG<- Y[,s]
		act.index<- !apply(is.na(XREG),1, any)[-dim(XREG)[1]] # without NA's and last
		model<- lm.fit(XREG[act.index, ], YREG[act.index] )
		forecast[s]<- model$coef %*% XREG[dim(XREG)[1],]
		}
	forecast
}

D<- 365*2+1 
N<- 365*4 
oos.dates<- dates.S[D+1:N] #

new_models = c("expert_nonlin")
model.names=c(model.names, new_models)
M<- length(model.names) #

ERROR<- array(ERROR, dim=c(N, S, M))
dimnames(ERROR)<- list(oos.dates, 1:S, model.names)
for(i.N in 1:N){
	index<- 1:D + i.N -1  ## defining index set for input data
	##Model 7
	forecast_expert_minmax<- forecast.expert.minmax(price.S[index,], days=dates.S[index])
	ERROR[i.N,, 9] <- forecast_expert_minmax - price.S[index[length(index)]+1,]
}#i.N

## forecasting:
MAE<- apply( abs(ERROR) ,3, mean)
RMSE<- sqrt(apply( ERROR^2 ,3, mean))
MAE
RMSE

compare_models()


####
