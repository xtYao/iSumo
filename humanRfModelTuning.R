h2o.init()
## TODO: split training set
## convention: label in last col, IDs in first 2 cols
dt = h2o.importFile(paste("data/",taxId,".iSumo.txt",sep=""))
dt.sp = h2o.splitFrame(dt, c(0.65,0.25), seed = 123)
train = dt.sp[[1]]; valid = dt.sp[[2]]; test = dt.sp[[3]]

## balance by split, not over-sampling
balanceBySplit = function(Dat){
	Labels = Dat[, ncol(Dat), with=F]
	Nlabels = table(Labels)
	Ratio = as.integer(Nlabels['FALSE']/Nlabels['TRUE'])
	Group = sample(1:Ratio, Nlabels['FALSE'], replace=TRUE)
	Negatives = cbind(Group, Dat[!Labels,])
	Positives = Dat[Labels,]
	BalancedDat = list()
	for (i in 1:Ratio){
		BalancedDat[[i]] = rbind(Negatives[Group==i, -1], Positives)
	}
	return(BalancedDat)
}

nc = ncol(dt)
## TODO: tune par
## default settings:
rf0 = h2o.randomForest(x = 3:(nc-1), y = nc,
					   training_frame = train, validation_frame = valid,
					   balance_classes = F,
					   nfolds = 10,
					   keep_cross_validation_predictions = T,
					   keep_cross_validation_fold_assignment = T,
					   ntrees = 200,
					   stopping_rounds = 2,
					   score_each_iteration = T,
					   seed = 123)

## default plus balancing
rf0.b = h2o.randomForest(x = 3:(nc-1), y = nc,
						 training_frame = train, validation_frame = valid,
						 balance_classes = T,
						 nfolds = 10,
						 keep_cross_validation_predictions = T,
						 keep_cross_validation_fold_assignment = T,
						 ntrees = 200,
						 stopping_rounds = 2,
						 score_each_iteration = T,
						 seed = 123)

## parameter tuning: larger max-depth per tree
rf0.b.md30 = h2o.randomForest(x = 3:(nc-1), y = nc,
							  training_frame = train, validation_frame = valid,
							  balance_classes = T,
							  nfolds = 10,
							  keep_cross_validation_predictions = T,
							  keep_cross_validation_fold_assignment = T,
							  ntrees = 200,
							  stopping_rounds = 2,
							  score_each_iteration = T,
							  seed = 123,
							  max_depth = 30)

## parameter tuning: stopping rounds > 2
rf1.b = h2o.randomForest(x = 3:(nc-1), y = nc,
						 training_frame = train, validation_frame = valid,
						 balance_classes = T,
						 nfolds = 10,
						 keep_cross_validation_predictions = T,
						 keep_cross_validation_fold_assignment = T,
						 ntrees = 200,
						 stopping_rounds = 3,
						 score_each_iteration = T,
						 seed = 123,
						 max_depth = 30)

## 
rf2.b.md30 = h2o.randomForest(x = 3:(nc-1), y = nc,
							  training_frame = train, validation_frame = valid,
							  balance_classes = T,
							  nfolds = 10,
							  keep_cross_validation_predictions = T,
							  keep_cross_validation_fold_assignment = T,
							  ntrees = 50,
							  #stopping_rounds = 3,
							  score_each_iteration = T,
							  seed = 123,
							  max_depth = 30)

rf2.b.md30.nt100 = h2o.randomForest(x = 3:(nc-1), y = nc,
							  training_frame = train, validation_frame = valid,
							  balance_classes = T,
							  nfolds = 10,
							  keep_cross_validation_predictions = T,
							  keep_cross_validation_fold_assignment = T,
							  ntrees = 100,
							  #stopping_rounds = 3,
							  score_each_iteration = T,
							  seed = 123,
							  max_depth = 30)

##
rf0.b.md30.mt10 = h2o.randomForest(x = 3:(nc-1), y = nc,
								   training_frame = train, validation_frame = valid,
								   balance_classes = T,
								   nfolds = 10,
								   keep_cross_validation_predictions = T,
								   keep_cross_validation_fold_assignment = T,
								   ntrees = 200,
								   mtries = 10,
								   stopping_rounds = 2,
								   score_each_iteration = T,
								   seed = 123,
								   max_depth = 30)

## TODO: compare the validation performance
rfs = list(rf0, rf0.b, rf0.b.md30, rf0.b.md30.mt10, ## early converging
		   rf1.b, ## late converging
		   rf2.b.md30, rf2.b.md30.nt100) ## hard-set nTree
metrics = do.call(rbind,
				  lapply(rfs,
				  	   function(x){
				  	   	tAuc = h2o.auc(x@model$training_metrics)
				  	   	vAuc = h2o.auc(x@model$validation_metrics)
				  	   	xAuc = h2o.auc(x@model$cross_validation_metrics)
				  	   	return(c(tAuc,vAuc,xAuc))
				  	   }))

## conclusion: the last one is the best, 