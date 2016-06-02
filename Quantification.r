#uncomment the following lines and install the packages first
#install.packages('quadprog')
#install.packages('LiblineaR')
#install.packages('NlcOptim')
#install.packages('e1071')
#install.packages('linprog')



#Optimization Functions
quad.constrain<-function(Y,X){ 
	require(quadprog)
	K<-ncol(X)
	Amat<-t(rbind(rep(1,K),diag(1,K),diag(-1,K)))
	solve.QP(t(X)%*%X,matrix(Y,nrow=1)%*%X, Amat, c(1, rep(0,K), rep(-1,K)), meq=1)$solution
}
	
hellinger<-function(Y,X, start=NULL){
	require(NlcOptim)
	K<-ncol(X)
	if (is.null(start)) start = rep(1/K,K)
	Aeq=matrix(1,ncol=K,nrow=1)
	Beq=matrix(1,1,1)
	conf<-function(x) return(list(ceq=NULL,c=NULL))	
	objfun<-function(beta){Yh<-(X%*%beta);Yh[Yh<0]<-0;1-sum(sqrt(Y*Yh))}
	as.numeric(NlcOptim(X=start, objfun=objfun,Aeq=Aeq,Beq=Beq, lb=matrix(0,K,1), ub=matrix(1,K,1),confun=conf)$p)
}

L1<-function(Y,X,start=NULL){
	require(NlcOptim)
	K<-ncol(X)
	if (is.null(start)) start = rep(1/K,K)
	Aeq=matrix(1,ncol=K,nrow=1)
	Beq=matrix(1,1,1)
	conf<-function(x) return(list(ceq=NULL,c=NULL))
	objfun<-function(beta){q<-X%*%beta; q[q<0]<-0; p<-Y; sum(abs(p-q))}
	as.numeric(NlcOptim(X=matrix(start,K,1), objfun=objfun,Aeq=Aeq,Beq=Beq, lb=matrix(0,K,1), ub=matrix(1,K,1),confun=conf)$p)
}

LAD<-function(Y,X){
	require(linprog)
	M<-nrow(X)
	K<-ncol(X)
	I<-diag(1,M,M)
	Amat<-rbind(cbind(I,X),cbind(I,-X),c(rep(0,M), rep(1,K)))
	bvec<-c(Y,-Y,1)
	cvec<-c(rep(1,M), rep(0,K))
	l<-solveLP(cvec,bvec,Amat,const.dir=c(rep(">=",2*M),"=="),lpSolve=T)
	as.numeric(tail(l$solution[],K))
}

#Base classifiers
svmModelFunction<-function(xTrain){
	require(e1071)
	svm(xTrain[,-1],xTrain[,1],probability=TRUE, kernel="linear", type="C-classification")
}

svmPredFunction<-function(PLT, testData){
	pred<-predict(PLT,testData,probability=TRUE); 
	probs<-attr(pred,"probabilities")
	probs[,order(as.numeric(colnames(probs)))]
	
}

liblinearModelFunction<-function(xTrain){
	require(LiblineaR)
	LiblineaR(data=xTrain[,-1],target=xTrain[,1], type=7,cost=1, bias=TRUE)
}

liblinearPredFunction<-function(PLT, testData){
	probs<-predict(PLT,testData, proba=T)$probabilities
	probs[,order(as.numeric(colnames(probs)))]
}

# methods

va1<-function(trainMatrix, testMatrix, trials=300, n=15){
	bintodec<-function(x){
		 sum(x * 2^ ((length(x)-1):0))
	}

	getProb <- function(testSet, trainingSet) {
        temp <- rbind(testSet, trainingSet)
        prob.wt <- apply(temp[, -1], 2, function(x) prod(prop.table(table(x))))
        prob.wt[colSums(testSet[, -1]) == 0] <- 0
        prob.wt[colSums(trainingSet[, -1]) == 0] <- 0
        prob.wt[colSums(testSet[, -1]) == nrow(testSet)] <- 0
        prob.wt[colSums(trainingSet[, -1]) == nrow(trainingSet)] <- 0
        prop.table(prob.wt)
	}
	K<-length(table(trainMatrix[,1]))
    resList<-NULL
    prob=getProb(testMatrix,trainMatrix)
    
    XL<-NULL
    YL<-NULL
    
    for (i in 1:trials){
			cols<-sample(2:ncol(trainMatrix),n,prob=prob)
			XR<-NULL
			profilesTest<-apply(testMatrix[,cols],1,bintodec)
			profilesTrain<-apply(trainMatrix[,cols],1,bintodec)		
			profiles<-unique(union(profilesTest,profilesTrain))
			Y<-prop.table(table(factor(profilesTest, levels=profiles)))
			for (j in 1:K)
				XR<-cbind(XR,prop.table(table(factor(profilesTrain[trainMatrix[,1]==j], levels=profiles))))
			XL<-rbind(XL,XR)
			YL<-c(YL,Y)
	}
	tryCatch(quad.constrain(YL,XL),error = function(e) cat(paste(e, collapse = ","), "\n"))
}

friedman<-function(testPred,trainingFractions,cvPred, cvPredY){
		K<-ncol(testPred)
		
		fF<-matrix(0,nrow(testPred),K)
		for (i in 1:K) fF[,i]<-testPred[,i]>=trainingFractions[i]
	
		flF<-rep(0,K)
		for (i in 1:K) flF[i]<-mean(fF[,i],na.rm=T)			

		fA<-matrix(0,nrow(cvPred),K)
		for (i in 1:K) fA[,i]<-cvPred[,i]>=trainingFractions[i]

		flA<-matrix(0,K,K)
		for (l in 1:K) 
			for (k in 1:K) 
				flA[l,k]<-mean(fA[cvPredY==k,l],na.rm=T)
		flA<-flA+matrix(rnorm(nrow(flA)*ncol(flA),sd=0.0001),nrow(flA),ncol(flA))
		 		
		tryCatch(quad.constrain(flF,flA),error = function(e) cat("friedman  warning\n"))
			 
}

probAveraging<-function(testPred,cvPred, cvPredY){
		K<-ncol(testPred)
		
		Y<-rep(0,K)
		for (i in 1:K) Y[i]<-mean(testPred[,i])		

		X<-matrix(0,K,K)
		for (l in 1:K) 
			for (k in 1:K) 
				X[l,k]<-mean(cvPred[which(cvPredY==k),l],na.rm=T)
		 X<-X+matrix(rnorm(K*K,sd=0.0001),K,K)
		 		
		 tryCatch(quad.constrain(Y,X),error = function(e)  cat("probAveraging warning\n"))
}



mixtureL2<-function(testPred,cvPred, cvPredY,numberOfBins=15){
		K<-ncol(testPred)
		
		binProbs=seq(0, 1, by=1/numberOfBins)
		X<-matrix(0,(numberOfBins)*K,K)
		Y<-matrix(0,(numberOfBins)*K,1)
		
		for (k in 1:K){
				binNumsTrain<-findInterval(cvPred[,k],binProbs,rightmost.closed=T)
				binNumsTest<-findInterval(testPred[,k], binProbs,rightmost.closed=T)
				for (ki in 1:numberOfBins){
					Y[(k-1)*numberOfBins+ki]<-mean(binNumsTest<=ki,na.rm=T)
					for (j in 1:K) X[(k-1)*numberOfBins+ki,j]<-mean(binNumsTrain[cvPredY==j]<=ki,na.rm=T)
				}
		}
		X<-X+matrix(rnorm(numberOfBins*K*K,sd=0.0001),numberOfBins*K,K)
		tryCatch(quad.constrain(Y,X),error = function(e)  cat("mixtureLS warning\n"))
}

rdirichlet <- function(alpha) prop.table(rgamma(length(alpha), alpha))

hellingerPMF<-function(testPred,cvPred, cvPredY,start=NULL,numberOfBins=15){
		K<-ncol(testPred)
		
		binProbs=seq(0, 1, by=1/numberOfBins)
		X<-matrix(0,(numberOfBins)*K,K)
		Y<-matrix(0,(numberOfBins)*K,1)
		
		for (k in 1:K){
				binNumsTrain<-findInterval(cvPred[,k],binProbs,rightmost.closed=T)
				binNumsTest<-findInterval(testPred[,k], binProbs,rightmost.closed=T)
				for (ki in 1:numberOfBins){
					Y[(k-1)*numberOfBins+ki]<-mean(binNumsTest==ki,na.rm=T)
					for (j in 1:K)
						X[(k-1)*numberOfBins+ki,j]<-mean(binNumsTrain[cvPredY==j]==ki,na.rm=T)
				}
		}
		X<-X+matrix(rnorm(numberOfBins*K*K,sd=0.0001),numberOfBins*K,K)
		
		res<-tryCatch(hellinger(Y,X, start),error = function(e) cat("hellingerPMF warning\n"))
		if (!is.null(res)) return(res)
		res<-tryCatch(hellinger(Y,X),error = function(e) cat("hellingerPMF warning\n"))
		if (!is.null(res)) return(res)
		for (i in 1:30){
			res<-tryCatch(hellinger(Y,X, start=rdirichlet(rep(1,K))),error = function(e) {})
			if (!is.null(res)) {
				cat("Success!\n")
				return(res)
			}
		}
}

hdx<-function(trainMatrix, testMatrix,start=NULL){
		K<-length(table(trainMatrix[,1]))
		
		X<-NULL
		Y<-NULL

		for (i in 1:(ncol(trainMatrix)-1)){
		    Y<-c(Y,prop.table(table(factor(testMatrix[,i+1],levels=c(0,1)))))
		    XRow<-NULL
		    for (j in 1:K)
		    	XRow<-cbind(XRow,prop.table(table(factor(trainMatrix[trainMatrix[,1]==j,i+1],levels=c(0,1)))))
		    X<-rbind(X,XRow)
		 }
		 
		 res<-tryCatch(hellinger(Y,X,start),error = function(e) cat("hdx warning\n"))
		 if (!is.null(res)) return(res)
		 res<-tryCatch(hellinger(Y,X),error = function(e) cat("hdx warning\n"))
		 if (!is.null(res)) return(res)
		 for (i in 1:30){
			res<-tryCatch(hellinger(Y,X, start=rdirichlet(rep(1,K))),error = function(e) {})
			if (!is.null(res)) {
				cat("Success!\n")
				return(res)
			}
		}
}

adjCount<-function(testPred,cvPred, cvPredY){
		K<-ncol(testPred)
		
		fF<-matrix(0,nrow(testPred),K)
		for (i in 1:K) fF[,i]<-testPred[,i]>=apply(testPred,1,max)
	
		flF<-rep(0,K)
		for (i in 1:K) flF[i]<-mean(fF[,i],na.rm=T)			

		fA<-matrix(0,nrow(cvPred),K)
		for (i in 1:K) fA[,i]<-cvPred[,i]>=apply(cvPred,1,max)

		flA<-matrix(0,K,K)
		for (l in 1:K) 
			for (k in 1:K) 
				flA[l,k]<-mean(fA[cvPredY==k,l],na.rm=T)
		flA<-flA+matrix(rnorm(nrow(flA)*ncol(flA),sd=0.0001),nrow(flA),ncol(flA))
		 		
		tryCatch(quad.constrain(flF,flA),error = function(e) cat("adjCount  warning\n"))	
}

medianSweep<-function(testPred,cvPred, cvPredY){
		K<-ncol(testPred)
			
		binSize<-100
		betas<-array(0,c(K,binSize))
		thresholds<-seq(0,1-1/binSize,1/binSize)
		counter<-1
		
		for (threshold in thresholds){
			fF<-matrix(0,nrow(testPred),K)
			for (i in 1:K) fF[,i]<-testPred[,i]>=threshold
	
			flF<-rep(0,K)
			for (i in 1:K) flF[i]<-mean(fF[,i],na.rm=T)			

			fA<-matrix(0,nrow(cvPred),K)
			for (i in 1:K) fA[,i]<-cvPred[,i]>=threshold

			flA<-matrix(0,K,K)
			for (l in 1:K) 
				for (k in 1:K) 
					flA[l,k]<-mean(fA[cvPredY==k,l],na.rm=T)
			flA<-flA+matrix(rnorm(nrow(flA)*ncol(flA),sd=0.0001),nrow(flA),ncol(flA))
			res<-tryCatch(quad.constrain(flF,flA),error = function(e) cat("medianSweep  warning\n"))
		 	if (!is.null(res)) betas[,counter]<-res
		 	else betas[,counter]<-rep(NA,K)
		 	counter<-counter+1
		}
		betas
}

mixtureL1<-function(testPred,cvPred, cvPredY,numberOfBins=15){
		K<-ncol(testPred)
		binProbs=seq(0, 1, by=1/numberOfBins)
		X<-matrix(0,(numberOfBins)*K,K)
		Y<-matrix(0,(numberOfBins)*K,1)
		
		for (k in 1:K){
				binNumsTrain<-findInterval(cvPred[,k],binProbs,rightmost.closed=T)
				binNumsTest<-findInterval(testPred[,k], binProbs,rightmost.closed=T)
				for (ki in 1:numberOfBins){
					Y[(k-1)*numberOfBins+ki]<-mean(binNumsTest<=ki,na.rm=T)
					for (j in 1:K) X[(k-1)*numberOfBins+ki,j]<-mean(binNumsTrain[cvPredY==j]<=ki,na.rm=T)
				}
		}
		X<-X+matrix(rnorm(numberOfBins*K*K,sd=0.0001),numberOfBins*K,K)
		X[X<0]<-0
		X[X>1]<-1
		Y[Y<0]<-0
		Y[Y>1]<-1
		tryCatch(LAD(Y,X),error = function(e)  cat("mixtureL1 LAD warning\n"))
		
}

naiveMethod<-function(trainMatrix,testMatrix,baseClassifierModelFunction=liblinearModelFunction,baseClassifierPredFunction=liblinearPredFunction){
	model<-baseClassifierModelFunction(trainMatrix)
	probs<-baseClassifierPredFunction(model,testMatrix[,-1])
	probs<-probs[,order(as.numeric(colnames(probs)))]
	colMeans(probs)
}


getNaiveProbs<-function(trainMatrix,testMatrix,baseClassifierModelFunction=liblinearModelFunction,baseClassifierPredFunction=liblinearPredFunction){
	model<-baseClassifierModelFunction(trainMatrix)
	probs<-baseClassifierPredFunction(model,testMatrix[,-1])
	probs[,order(as.numeric(colnames(probs)))]
}

#support functions

repmat<-function(x, n) matrix(rep(x, n), nrow = n, byrow = T)

cvRandomHalfSplit<-function(trainMatrix){
	trainMatrix<-trainMatrix[sample(1:nrow(trainMatrix),nrow(trainMatrix)),] 
	trainIndices<-tapply(1:nrow(trainMatrix), list(trainMatrix[,1]), function(s) s)
	cvTrainIndices<-unlist(lapply(trainIndices,function(x) x[1:(length(x)%/%2)]))
	cvTestIndices<-setdiff(unlist(trainIndices),cvTrainIndices)
	list(cvTrain=trainMatrix[cvTrainIndices,],cvTest=trainMatrix[cvTestIndices,])
}

getMean<-function(classMap,betas){
	if (length(betas)<3) return(betas)
	colnames(betas)<-classMap	
	colMeans(betas,na.rm=T)
}

getMedian<-function(classMap,betas){
	if (length(betas)<3) return(betas)
	colnames(betas)<-classMap	
	apply(betas,2,median,na.rm=T)
}

quantifyAll<-function(trainMatrix, testMatrix, baseClassifierModelFunction=liblinearModelFunction,baseClassifierPredFunction=liblinearPredFunction, trials=50){
	classMap<-levels(as.factor(trainMatrix[,1]))
	trainMatrix[,1]<-as.integer(as.factor(trainMatrix[,1]))
	testMatrix[,1]<-as.integer(factor(testMatrix[,1], levels=classMap))
	K<-length(classMap)
	friedmanBetas<-probAveragingBetas<-mixtureL2Betas<-mixtureL1Betas<-mixtureHPMFBetas<-adjCountBetas<-NULL
	medianSweepBetas<-array(0,c(trials,K,100))
	actual<-prop.table(table(factor(testMatrix[,1],levels=classMap)))
	for (trial in 1:trials){
		cvSplits<-cvRandomHalfSplit(trainMatrix)
		model<-baseClassifierModelFunction(cvSplits$cvTrain)
		cvPred<-baseClassifierPredFunction(model,cvSplits$cvTest[,-1])
		testPred<-baseClassifierPredFunction(model,testMatrix[,-1])
		trainingFractions<-prop.table(table(cvSplits$cvTrain[,1]))
		cvPredY<-cvSplits$cvTest[,1]
		
		res<-friedman(testPred,trainingFractions,cvPred, cvPredY)
		if (!is.null(res)) friedmanBetas<-rbind(friedmanBetas, res)
		
		res<-probAveraging(testPred,cvPred, cvPredY)
		if (!is.null(res)) probAveragingBetas<-rbind(probAveragingBetas, res)
		
		res<-mixtureL2(testPred,cvPred, cvPredY)
		if (!is.null(res)) mixtureL2Betas<-rbind(mixtureL2Betas, res)
		
		res<-mixtureL1(testPred,cvPred, cvPredY)
		if (!is.null(res)) mixtureL1Betas<-rbind(mixtureL1Betas, res)
		
		res<-hellingerPMF(testPred,cvPred, cvPredY,actual)
		if (!is.null(res)) mixtureHPMFBetas<-rbind(mixtureHPMFBetas, res)
		
		res<-medianSweep(testPred,cvPred, cvPredY)
		if (!is.null(res)) medianSweepBetas[trial,,]<-res
		
		res<-adjCount(testPred,cvPred, cvPredY)
		if (!is.null(res)) adjCountBetas<-rbind(adjCountBetas, res)
	}
	betaPerThreshold<-apply(medianSweepBetas,3,colMeans,na.rm=T)
	medianSweep<-prop.table(matrix(apply(betaPerThreshold,1,median,na.rm=T),ncol=K))
	
	dt<-data.frame(matrix(0,11,K+1))
	dt[,1]<-c("friedman","prob","mixtureL2","mixtureL1","mixtureHPMF","adjCount","medianSweep","hdx","va1","actual","naive")

	dt[1,-1]<-getMean(classMap,friedmanBetas)
	dt[2,-1]<-getMean(classMap,probAveragingBetas)
	dt[3,-1]<-getMean(classMap,mixtureL2Betas)
	dt[4,-1]<-getMean(classMap,mixtureL1Betas)
	dt[5,-1]<-getMean(classMap,mixtureHPMFBetas)
	dt[6,-1]<-getMean(classMap,adjCountBetas)
	dt[7,-1]<-medianSweep
	dt[8,-1]<-hdx(trainMatrix,testMatrix,start=actual)
	dt[9,-1]<-va1(trainMatrix,testMatrix,n=2)
	dt[10,-1]<-actual
	dt[11,-1]<-naiveMethod(trainMatrix,testMatrix)
	colnames(dt)<-c("alg",classMap)
	dt
}
getPrunedSize<-function (coef, size, minN = 0) {
    max.index <- which.max(coef)
    max.value <- max(coef)
    prunedSize <- round(coef * (size[max.index]/max.value))
    removalCounter <- 1
    while (any(((size - prunedSize)/max.value) < -0.01)) {
        prunedSize <- round(coef * ((size[max.index] - removalCounter)/max.value))
        removalCounter <- removalCounter + 1
    }
    prunedSize[prunedSize < min(size, minN)] <- min(size, minN)
    prunedSize
}

getErrorRate<-function(truth, probs, corrections=NULL, trainingFractions=NULL){
	if (is.null(corrections)){
		return(1-mean(truth==apply(probs,1,which.max)))
	} else {
		return(1-mean(apply(apply(repmat(as.numeric(corrections)/as.numeric(trainingFractions),nrow(probs))*probs,1,prop.table),2,which.max)==truth))
	}
}


stanfordSentimentTests<-function(){
	
	baseClassifierModelFunction=liblinearModelFunction
	baseClassifierPredFunction=liblinearPredFunction
	load("sentTermMatrix.rda") #termMatrix
	train<-termMatrix[1:8544,]
	trainMatrix<-train
	test<-termMatrix[8545:nrow(termMatrix),]
	load("sentTreeTest.probs") #probs
	load("sentTreeTest.truth") #truth
	indices<-tapply(1:nrow(test), list(test[,1]), function(s) s)
	size<-table(test[,1])
	x<-expand.grid(seq(0,1,0.1),seq(0,1,0.1),seq(0,1,0.1),seq(0,1,0.1),seq(0,1,0.1))
	pos<-x[which(round(rowSums(x),1)==1.0),]
	quantifications<-NULL
	errorRates<-NULL
	trainingFractions<-prop.table(table(trainMatrix[,1]))
	for (i in 1:1001){
		print(i)
		testSizes<-getPrunedSize(pos[i,], size)
		testIndices<-NULL
		for (j in 1:5) {
			if (testSizes[j]==0) next
			testIndices<-c(testIndices, sample(indices[[j]],as.numeric(testSizes[j])))
		}
		testMatrix<-test[testIndices,]
		d<-quantifyAll(trainMatrix,testMatrix)
		quantifications<-rbind(quantifications,d)
		ER<-NULL
		for (q in 1:nrow(d)){
			if (q==nrow(d)) ER<-c(ER,getErrorRate(truth[testIndices],probs[testIndices,]))
			else ER<-c(ER,getErrorRate(truth[testIndices],probs[testIndices,],d[q,-1],trainingFractions))
		}
		errorRates<-c(errorRates,ER)
		print(d)
		print(ER)
		
	}
	list(quantifications=quantifications, errorRates=errorRates)
}

marketingDataTests<-function(){
	load("nd.rda")
	getDataSets<-function(futureProps){
		subData<-nd
		subData<-subData[sample(1:nrow(subData), nrow(subData)),]
	
		trainingIndices<-tapply(1:nrow(subData), list(subData$Occupation), function(s) s)
	
		trainingIndicesModel<-unlist(lapply(trainingIndices,function(x) x[1:(2*length(x)/3)]))
	
		futureDataIndices<-lapply(trainingIndices,function(x) x[(1+2*length(x)/3):(length(x))])
	
		pRatio<-futureProps/prop.table(table(subData[unlist(futureDataIndices),]$Occupation))
	
		maxInd<-which.max(pRatio)
		futureCategoryCardinality<-floor((table(subData[unlist(futureDataIndices),]$Occupation)[maxInd]/futureProps[maxInd])*futureProps)
		testIndices<-NULL
		for (i in 1:length(futureProps))
			testIndices<-c(testIndices,sample(futureDataIndices[[i]], futureCategoryCardinality[i]))

		trainingIndicesModel<-sample(trainingIndicesModel,length(trainingIndicesModel))
		testIndices<-sample(testIndices,length(testIndices))
		list(trainMatrix=subData[trainingIndicesModel,],testMatrix=subData[testIndices,])
	}
	
	baseClassifierModelFunction=liblinearModelFunction
	baseClassifierPredFunction=liblinearPredFunction
	
	x<-expand.grid(seq(0,1,0.1),seq(0,1,0.1),seq(0,1,0.1),seq(0,1,0.1),seq(0,1,0.1))
	pos<-x[which(round(rowSums(x),1)==1.0),]
	quantifications<-NULL
	errorRates<-NULL
	
	for (i in 1:nrow(pos)){
		print(i)
		d<-getDataSets(as.numeric(pos[i,]))
		trainMatrix<-d$trainMatrix
		testMatrix<-d$testMatrix
		classMap<-levels(as.factor(trainMatrix[,1]))
		trainMatrix[,1]<-as.integer(as.factor(trainMatrix[,1]))
		testMatrix[,1]<-as.integer(factor(testMatrix[,1], levels=classMap))
		d<-quantifyAll(trainMatrix,testMatrix)
		quantifications<-rbind(quantifications,d)
		ER<-ER2<-NULL
		probs<-getNaiveProbs(trainMatrix,testMatrix)
		truth<-testMatrix[,1]
		trainingFractions<-prop.table(table(trainMatrix[,1]))
		for (q in 1:nrow(d)){
			if (q==nrow(d)) ER<-c(ER,getErrorRate(truth,probs))
			else ER<-c(ER,getErrorRate(truth,probs,d[q,-1], trainingFractions))
		}
		errorRates<-c(errorRates,ER)
		print(d)
		print(ER)
		
	}
	list(quantifications=quantifications, errorRates=errorRates)
}


#uncomment the following lines to run the tests
#results<-marketingDataTests()
#results<-stanfordSentimentTests()



