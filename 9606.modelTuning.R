h2o.init()
## TODO: split training set
## convention: label in last col, IDs in first 2 cols
dt = h2o.importFile(paste("data/",taxId,".iSumo.txt",sep=""))
dt.sp = h2o.splitFrame(dt, c(0.65,0.25), seed = 123)
train = dt.sp[[1]]; valid = dt.sp[[2]]; test = dt.sp[[3]]

nc = ncol(dt)
## TODO: tune par
## imbalanced, stopping round 2, default max depth 20
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

## imbalanced, stopping with 200 trees, default max depth 20
rf0.nt200 = h2o.randomForest(x = 3:(nc-1), y = nc,
                       training_frame = train, validation_frame = valid,
                       balance_classes = F,
                       nfolds = 10,
                       keep_cross_validation_predictions = T,
                       keep_cross_validation_fold_assignment = T,
                       ntrees = 200,
                       score_each_iteration = T,
                       seed = 123)

## balancing, stopping round 2, default max depth 20
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

## balancing, stopping rounds 2, max depth 30
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

## balancing, stopping round 2, max depth 30, mtries 10
rf0.b.md30.mt10 = h2o.randomForest(x = 3:(nc-1), y = nc,
								   training_frame = train,
								   validation_frame = valid,
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

## balancing, stopping rounds 3, max depth 30
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

## balancing, 50 trees, max depth 30, default mtries
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
## balancing, 100 trees, max depth 30, default mtries
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

## balancing, 200 trees, max depth 30, default mtries
rf2.b.md30.nt200 = h2o.randomForest(x = 3:(nc-1), y = nc,
									training_frame = train,
									validation_frame = valid,
									balance_classes = T,
									nfolds = 10,
									keep_cross_validation_predictions = T,
									keep_cross_validation_fold_assignment = T,
									ntrees = 200,
									#stopping_rounds = 3,
									score_each_iteration = T,
									seed = 123,
									max_depth = 30)

## looks like if we keep increasing the number of trees, we can acheive better
## validation and cross-validation AUC, but at the cost of long training time

## TODO: compare the validation performance
rfs = list(rf0,
           rf0.nt200,
           rf0.b,
           rf0.b.md30,
           rf0.b.md30.mt10, ## early converging
		   rf1.b, ## late converging
		   rf2.b.md30,
		   rf2.b.md30.nt100,
		   rf2.b.md30.nt200) ## hard-set nTree
names(rfs) = sapply(rfs, function(x) x@model_id)

## manifest of all these models
manifest =
    do.call('rbind',
            lapply(rfs, function(x) as.data.table(x@model$model_summary)))
manifest[, id := names(rfs)]
setkey(manifest, "id")

## gather metrics to do model selection
metrics = do.call(rbind,
                  lapply(rfs,
                         function(x){
                             ## AUCs
                             trAuc = h2o.auc(x, train = T)
                             vAuc = h2o.auc(x, valid = T)
                             xAuc = h2o.auc(x, xval = T)
                             tsAuc = h2o.auc(h2o.performance(x, test))
                             ## F1 and F2
                             perf = h2o.performance(x, valid = T)
                             maxPerf = as.data.table(
                                 perf@metrics$max_criteria_and_metric_scores)
                             vF1 = maxPerf[metric=="max f1", value]
                             vF2 = maxPerf[metric=="max f2", value]
                             return(c(trAuc,vAuc,xAuc,tsAuc, vF1, vF2))
                         }))
colnames(metrics) = c("trainingAuc", "validAuc", "xValidationAuc", "testAuc",
                      "validMaxF1", "validMaxF2")
metrics = data.table(metrics, keep.rownames = T, key = "rn")
manifest = merge(manifest, metrics, by.x = "id", by.y = "rn")

## save models
system(paste("mkdir -p models/", taxId, sep=""))
for (m in names(rfs)){
    fn = paste("models/", taxId, "/", m, ".h2o", sep="")
    h2o.saveModel(rfs[[m]], fn)
}
## save performance metrics
write.table(manifest,
            paste("models/",taxId,"/manifest.txt",sep=""),
            sep = "\t", quote = F, row.names = names(rfs))



## conclusion: the non-balanced 200 tree 20 max depth is the best
rf = rfs[['DRF_model_R_1481751746690_762']]
h2o.saveModel(rf, paste("models/",taxId,".rf.h2o",sep=""))

## fit null model the same way
rf.null = h2o.randomForest(x = names(gpsCount)[-1], y = "isSumo",
                           training_frame = train, validation_frame = valid,
                           balance_classes = F,
                           nfolds = 10,
                           keep_cross_validation_predictions = T,
                           keep_cross_validation_fold_assignment = T,
                           ntrees = 200,
                           score_each_iteration = T,
                           seed = 123)
h2o.saveModel(rf.null, paste("models/",taxId,".rf.null.h2o",sep=""))
