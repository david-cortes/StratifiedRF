#' @import C50
#' @import dplyr
#' @import parallel
#' @importFrom stats as.formula
#' @importFrom stats na.pass
`%>%` <- magrittr::`%>%`


#' Stratified Random Forest
#'
#' Random Forest that works with groups of predictor variables. When building a tree, a number of variables is taken from each group separately. Useful when rows contain information about different things (e.g. user information and product information) and it's not sensible to make a prediction with information from only one group of variables, or when there are far more variables from one group than the other and it's desired to have groups appear evenly on trees.
#'
#' Note that while this algorithm forces each tree to consider possible splits with variables from all groups, it doesn't guarantee that they will end up having splits with variables from different groups.
#'
#' The original Random Forest algorithm recommends a total number of sqrt(n_features), but this might not work so well when there are unequal groups of variables.
#'
#' Implementation of everything outside the tree-building is in native R code, thus might be slow. Trees are grown using the C5.0 algorithm from the 'C50' library, thus it can be used for classification only (not for regression). Refer to the 'C50' library for any documentation about the tree-building algorithm.

#' @param df Data to build the model (data.frame only).
#' @param targetvar String indicating the name of the target or outcome variable in the data. Character types will be coerced to factors.
#' @param groups Unnamed list, containing at each entry a group of variables (as a string vector with their names).
#' @param mtry A numeric vector indicating how many variables to take from each group when building each tree. If set to "auto" then, for each group, mtry=round(sqrt(m_total)*len(m_group)/len(m_total)) (with a minimum of 1 for each group).
#' @param ntrees Number of trees to grow. When setting multicore=TRUE, the number of trees should be a multiple of the number of cores, otherwise it will get rounded downwards to the nearest multiple.
#' @param multicore Whether to use multiple CPU cores to parallelize the construction of trees. Parallelization is done with the 'parallel' library's default settings.
#' @param class_quotas How many rows from each class to use in each tree (useful when there is a class imbalance). Must be a numeric vector or a named list with the number of desired rows to sample for each level of the target variable. Ignored when sample_weights is passed. Note that using more rows than the data originally had might result in incorrect out-of-bag error estimates.
#' @param sample_weights Probability of sampling each row when building a tree. Must be a numeric vector. If not defined, then all rows have the same probability. Note that, depending on the structure of the data, setting this might result in incorret out-of-bag error estimates.
#' @param fulldepth Whether to grow the trees to full depth. Ignored when passing c50_control.
#' @param replacement Whether to sample rows with replacement.
#' @param c50_control Custom parameters for growing trees. Must be a C5.0Control object compatible with the 'C50' package.
#' @param na.action A function indicating how to handle NAs. Default is to include missing values when building a tree (see 'C50' documentation).
#' @param drop_threshold Drop a tree whenever its resulting out-of-bag classification accuracy falls below a certain threshold specified here. Must be a number between 0 and 1.
#' @keywords stratified_rf
#' @export
#' @examples
#' data(iris)
#' groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
#' mtry <- c(1,1)
#' m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
#' summary(m)
#' @seealso 'C50' library: \url{https://cran.r-project.org/package=C50}
stratified_rf=function(df,targetvar,groups,mtry='auto',ntrees=500,multicore=TRUE,class_quotas=NULL,sample_weights=NULL,fulldepth=TRUE,replacement=TRUE,c50_control=NULL,na.action=na.pass,drop_threshold=NULL){
  # checking input parameters
  if (class(df)!="data.frame"){stop("Data to predict must be a data.frame")}
  if (class(targetvar)!="character"){stop('Target variable must be specified as a string')}
  if (class(groups)!="list"){stop('Predictor variables must be specified as an unnamed list of character vectors')}
  if (class(groups[[1]])!="character"){stop('Predictor variables must be specified as an unnamed list of character vectors')}
  if (class(df[,targetvar])=="character"){df[,targetvar]=as.factor(df[,targetvar])}
  if (class(df[,targetvar])!='factor'){stop('Target variable must be a factor')}

  if (class(mtry)=="character" && mtry=="auto"){
    m_total=sum(sapply(groups,length))
    frac_groups=sapply(groups,function(x) length(x)/m_total)
    mtry=round(sqrt(m_total)*frac_groups,0)
    mtry[mtry<1]=1}
  if (class(mtry)!="numeric"){stop('mtry must be a vector or list of numbers')}

  if ((class(ntrees)!="numeric")|(length(ntrees)!=1)){stop('ntrees must be a positive integer')}
  if (!is.null(sample_weights)){
    class_quotas=NULL
    if (class(sample_weights)!="numeric"){stop('sample_weights must be a numeric vector')}
    if (sum(sample_weights<0)>0){stop('sample_weights must contain only non-negative entries')}
    if (sum(sample_weights)<=0){stop('sample_weights must some non-zero entries')}
    }
  if (!is.null(class_quotas)){
    if (class(class_quotas)!='numeric' & class(class_quotas)!='list'){stop('class_quotas must be a numeric vector or list with one entry per level of the target variable')}
    levs_target=levels(eval(parse(text=paste0("df$",targetvar))))
    len_target=length(levs_target)
    if (length(class_quotas)!=len_target){stop('class_quotas must be a numeric vector or list with one entry per level of the target variable')}}
    if (class(class_quotas)=='numeric'){
        class_quotas=as.list(class_quotas)
        names(class_quotas)=levs_target
    }
  if (!is.null(drop_threshold)){
    if (class(drop_threshold)!='numeric'){stop('drop_threshold must be a number between 0 and 1')}
    if (drop_threshold<=0 | drop_threshold>=1){stop('drop_threshold must be a number between 0 and 1')}
  }

  if (multicore){return(
    stratified_rf_multi(df,targetvar,groups,mtry,ntrees,class_quotas,sample_weights,fulldepth,replacement,c50_control,na.action,drop_threshold)
  )}else{return(
    stratified_rf_single(df,targetvar,groups,mtry,ntrees,class_quotas,sample_weights,fulldepth,replacement,c50_control,na.action,drop_threshold)
  )}

}


#' Make predictions on new data
#'
#' Make predictions from a stratified_rf model on new data.
#'
#' Note that by default, for classification models the predictions are made quite differently from the original Random Forest algorithm.

#' @param object A stratified_rf model.
#' @param data New data on which to make predictions (data.frame only). Must have the same names as the data used to build the model.
#' @param type Prediction type. Either "class" to get the predicted class or "prob" to get the voting scores for each class.
#' @param agg_type How to combine the predictions from individual trees. Either "prob" to average the probabilities output from each tree or "class" to count the final predictions from each.
#' @param vote_type How to weight the outputs from each tree. Either "simple" to average them, or "weighted" for a weighted average according to their OOB classification accuracy.
#' @param  na.action Function indicating how to handle missing values (see the 'C50' documentation for details).
#' @param threshold Count only votes from trees whose out-of-bag classification accuracy is above this threshold. Must be a number between 0 and 1.
#' @param ... other options (not currently used)
#' @keywords predict.stratified_rf
#' @export
#' @examples
#' data(iris)
#' groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
#' mtry <- c(1,1)
#' m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
#' predict(m,iris)
#' @seealso 'C50' library: \url{https://cran.r-project.org/package=C50}
predict.stratified_rf=function(object,data,type='class',agg_type='prob',vote_type='simple',na.action=na.pass,threshold=NULL,...){

  if (!(type %in% c('class','prob','raw'))){stop("prediction type must be 'prob' or 'class'")}
  if (!(agg_type %in% c('class','prob'))){stop("aggregation type must be 'prob' or 'class'")}
  if (!(vote_type %in% c('simple','weighted'))){stop("vote type must be either simple or weighted")}
  if (class(data)!="data.frame"){stop("data to predict must be a data.frame")}
  if (sum(!(unlist(object$vargroups) %in% names(data)))>0){stop("data to predict doesn't have the same variables")}
  if (!is.null(threshold)){
    if ((class(threshold)!='numeric')|(length(threshold)>1)){stop('Threshold must be a number between 0 and 1')}
    if ((threshold>1)|(threshold<0)){stop('Threshold must be a number between 0 and 1')}
    ignored.trees=sapply(object$acc.trees, function(x) x<threshold)
    if (sum(ignored.trees)>0){
      object$trees=object$trees[-(which(ignored.trees,arr.ind=TRUE))]
      object$acc.trees=object$acc.trees[-(which(ignored.trees,arr.ind=TRUE))]
    }}


  if (agg_type=='prob') {
    if (vote_type=='simple'){preds=Reduce('+',lapply(object$trees,predict.C5.0,data,type='prob',na.action=na.action))}
    if (vote_type=='weighted'){
      preds=lapply(object$trees,predict.C5.0,data,type='prob',na.action=na.action)
      preds2=lapply(seq_along(preds),function(x) object$acc.trees[x]*preds[[x]])
      preds=Reduce('+',preds2)
    }
    if (type=='class'){return(object$levels_target[max.col(preds)])}
    if ((type=='prob')|(type=='raw')){return(preds/apply(preds,1,sum))}
  }

  if (agg_type=='class') {
    votes=Reduce(cbind,lapply(object$trees,predict.C5.0,data,type='class',na.action=na.action))

    if (vote_type=='simple'){
      if (type=='class'){return(object$levels_target[as.numeric(apply(votes,1,function(x) names(sort(-table(x)))[1]))])}
      if ((type=='prob')|(type=='raw')){
        vote_counts=data.frame(row.names=1:dim(data)[1])
        for (i in object$levels_target){vote_counts[,i]=apply(votes,1,function(x) sum(object$levels_target[x]==i))}
        return(vote_counts/apply(vote_counts,1,sum))
      }}

    if (vote_type=='weighted'){
      vote_counts=data.frame(row.names=1:dim(data)[1])
      for (i in object$levels_target){
        agg_func=function(x){
          cnt=(object$levels_target[x]==i)*object$acc.trees[seq_along(x)]
          return(sum(cnt))
        }
        vote_counts[,i]=apply(votes,1,agg_func)}

      if (type=='class'){return(object$levels_target[max.col(vote_counts)])}
      if ((type=='prob')|(type=='raw')){vote_counts/apply(vote_counts,1,sum)}
    }
  }

}


#' Heuristic on variable importance
#'
#' Heuristic on variable importance, taken as averages from the variable importances calculated for each tree.
#'
#' Methods are taken directly from the C5.0 trees. Currently doesn't support permutation tests.

#' @param model A stratified_rf model.
#' @param metric How to calculate the variable importance from each tree. Either "usage" or "splits".
#' @param agg_type How to aggregate the variable importances obtained from each tree. Either "simple" for a simple average, or "weighted" for an average weighted by each tree's accuracy.
#' @export
#' @examples
#' data(iris)
#' groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
#' mtry <- c(1,1)
#' m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
#' varimp_stratified_rf(m)
#' @return A named data frame with the importance score of each variable, sorted from largest to smallest.
varimp_stratified_rf=function(model,metric='usage',agg_type='simple'){
  if (class(model)!='stratified_rf'){stop("Model is not a stratified_rf")}
  if (agg_type=='simple'){
  tot=length(model$trees)
  lst.varimps=lapply(model$trees,C5imp,metric=metric,pct=FALSE)}
  if (agg_type=='weighted'){
  tot=sum(model$acc.trees)
  lst.varimps=mapply(function(tree,w) C5imp(tree,metric=metric,pct=FALSE)*w,model$trees,model$acc.trees,SIMPLIFY=FALSE,USE.NAMES = FALSE)
  }

  reduce.func=function(df1,df2){
    df1$key=row.names(df1)
    df2$key=row.names(df2)

    sum1=sum2=Overall=NULL
    res=df1 %>% rename(sum1=Overall) %>% full_join(df2 %>% rename(sum2=Overall),by=c('key'='key')) %>% rowwise() %>% mutate(Overall=sum(sum1,sum2,na.rm=TRUE))
    return(data.frame(row.names=res$key,Overall=res$Overall))
  }
  imps=Reduce(reduce.func,lst.varimps)/tot
  return(imps[order(-imps$Overall),,drop=FALSE])
}

#' Print summary statistics from a model
#' @param x A stratified_rf model.
#' @param ... other options (not currently used)
#' @export
#' @examples
#' data(iris)
#' groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
#' mtry <- c(1,1)
#' m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
#' print(m)
print.stratified_rf=function(x,...){
  summary(x)
}

#' Summary statistics from a model
#'
#'Calculates error statistics for out-of-bag samples from a stratified_rf model.
#'
#' Predictions for a class are made by averaging class probabilities across trees rather than by a majority vote. All trees are weighted equally.
#' @param object A stratified_rf model.
#' @param ... other options (not currently used)
#' @export
#' @examples
#' data(iris)
#' groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
#' mtry <- c(1,1)
#' m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
#' summary(m)
summary.stratified_rf=function(object,...){

  asperc=function(x){paste0(substr(round(100*x,2),1,4),'%')}
  cat('Stratified Random Forest object\n\n')
  cat('Out-of-bag prediction error: ',asperc(1-object$acc))
  cat('\n\nConfusion Matrix\n')
  tbl_errs=object$conf_matrix
  print(tbl_errs)
  cat('\n')
  class.errs=vector()
  for (i in 1:length(object$levels_target)){class.errs=c(class.errs,1-tbl_errs[i,i]/sum(tbl_errs[i,]))}
  for (i in 1:length(object$levels_target)){
    prec=tbl_errs[i,i]/sum(tbl_errs[,i])
    rec=tbl_errs[i,i]/sum(tbl_errs[i,])
    cat('Class ',object$levels_target[i],'- Precision:',asperc(prec),'Recall:',asperc(rec),'\n')
  }
  cat('\nPredictor Variables:\n')
  for (i in 1:length(object$vargroups)){
    cat('Group',i,': ',paste(object$vargroups[[i]],collapse=', '),'\n')
  }
  cat('\nTarget Variable: ',object$targetvar)
  cat('\n\nBuilt with',length(object$trees),'trees')

}



stratified_rf_single=function(df,targetvar,groups,mtry,ntrees,class_quotas=NULL,sample_weights=NULL,fulldepth=TRUE,replacement=TRUE,c50_control=NULL,na.action=na.pass,drop_threshold=NULL){

  # starting with the function
  tree.models=vector(mode="list", length=ntrees)
  tree.models.acc=vector(mode="list", length=ntrees)
  num.groups=length(groups)
  rowset=1:dim(df)[1]

  levels.target=levels(eval(parse(text=paste0("df$",targetvar))))
  oob.sumprob=data.frame(row.names=rowset)
  for (l in levels.target){
    oob.sumprob[l]=rep(0,dim(df)[1])
  }
  oob.times=rep(0,dim(df)[1])
  if (is.null(c50_control)){treecontrol=C5.0Control(noGlobalPruning = fulldepth)}
  else {treecontrol=c50_control}

  if (!is.null(class_quotas)){
    rowset.classes=list()
    for (i in levels.target){
      rowset.classes[[i]]=which(df[,targetvar]==i)
    }
  }

  # building all the trees
  for (i in 1:ntrees){

    # sampling columns
    vars=vector()
    for (j in 1:length(groups)){
      vars=c(vars,sample(groups[[j]],size=mtry[j],replace=FALSE))
    }

    # sampling rows
    if (is.null(class_quotas) & is.null(sample_weights)){
      trainrows=sample(rowset,replace=replacement)
    } else if (is.null(class_quotas) & !is.null(sample_weights)){
      trainrows=sample(rowset,replace=replacement,prob=sample_weights)
    } else {
      trainrows=list()
      for (l in levels.target){trainrows[[l]]=sample(rowset.classes[[l]],size=class_quotas[[l]],replace=TRUE)}
      trainrows=unlist(trainrows)
    }

    oob.rows=!(rowset %in% trainrows)

    traindata=df[trainrows,c(vars,targetvar)]
    form=as.formula(paste(targetvar, paste(vars, collapse=" + "), sep=" ~ "))

    # building and saving a tree
    tree=C5.0(form,data=traindata,control=treecontrol,rules=FALSE,trials=1,na.action=na.action)
    tree.models[[i]]=tree

    if (is.null(drop_threshold)){
      oob.sumprob[oob.rows,]=oob.sumprob[oob.rows,]+predict.C5.0(tree.models[[i]],df[oob.rows,],type='prob')
      oob.times[oob.rows]=oob.times[oob.rows]+1

      tree.models.acc[[i]]=mean(df[oob.rows,targetvar]==predict.C5.0(tree.models[[i]],df[oob.rows,],type='class'))}

    else {
      oob.acc=mean(df[oob.rows,targetvar]==predict.C5.0(tree.models[[i]],df[oob.rows,],type='class'))
      if (oob.acc>drop_threshold){
        oob.sumprob[oob.rows,]=oob.sumprob[oob.rows,]+predict.C5.0(tree.models[[i]],df[oob.rows,],type='prob')
        oob.times[oob.rows]=oob.times[oob.rows]+1
        tree.models.acc[[i]]=oob.acc
      }}
  }

  probs=oob.sumprob/oob.times
  pred=levels.target[max.col(probs)]
  real=eval(parse(text=paste0("df$",targetvar)))

  acc=mean(pred==real,na.rm = TRUE)
  conf_matrix=table(real,pred)

  if (is.null(drop_threshold)){
    return(structure(list(numtrees=length(tree.models),drop_threshold=0,trees=tree.models,acc.trees=unlist(tree.models.acc),oob.preds=probs,acc=acc,conf_matrix=conf_matrix,vargroups=groups,mtry=mtry,levels_target=levels.target,targetvar=targetvar,task='classification'),class='stratified_rf'))}
  else{
    if (sum(sapply(tree.models.acc,is.null))>0){
      tree.models=tree.models[-(which(sapply(tree.models.acc,is.null),arr.ind=TRUE))]
      tree.models.acc=tree.models.acc[-(which(sapply(tree.models.acc,is.null),arr.ind=TRUE))]
    }
    trees_after=length(tree.models)
    if (trees_after==0){stop('No tree had an accuracy above drop_threshold')}
    return(structure(list(numtrees=trees_after,drop_threshold=drop_threshold,trees=tree.models,acc.trees=unlist(tree.models.acc),oob.preds=probs,acc=acc,conf_matrix=conf_matrix,vargroups=groups,mtry=mtry,levels_target=levels.target,targetvar=targetvar,task='classification'),class='stratified_rf'))
  }
}


stratified_rf_multi=function(df,targetvar,groups,mtry,ntrees,class_quotas=NULL,sample_weights=NULL,fulldepth=TRUE,replacement=TRUE,c50_control=NULL,na.action=na.pass,drop_threshold=NULL){

  ncores=detectCores()
  ntrees=as.integer(ntrees/ncores)

  tree.models=vector(mode="list", length=ntrees)
  tree.models.acc=vector(mode="list", length=ntrees)
  num.groups=length(groups)
  rowset=1:dim(df)[1]

  levels.target=levels(eval(parse(text=paste0("df$",targetvar))))
  oob.sumprob=data.frame(row.names=rowset)
  for (l in levels.target){
    oob.sumprob[l]=rep(0,dim(df)[1])
  }
  oob.times=rep(0,dim(df)[1])
  if (is.null(c50_control)){treecontrol=C5.0Control(noGlobalPruning = fulldepth)}
  else {treecontrol=c50_control}

  if (!is.null(class_quotas)){
    rowset.classes=list()
    for (i in levels.target){
      rowset.classes[[i]]=which(df[,targetvar]==i)
    }
  }

  cl<-makeCluster(ncores)
  clusterExport(cl, c("df","targetvar","groups","mtry","ntrees","class_quotas","sample_weights","fulldepth","replacement","c50_control","na.action","drop_threshold"),envir=environment())
  clusterCall(cl, function() requireNamespace("C50"))
  clusterCall(cl, function() requireNamespace("stats"))
  res_list=parLapply(cl,seq(1,ncores),function(noarg){

    # building all the trees
    for (i in 1:ntrees){

      # sampling columns
      vars=vector()
      for (j in 1:length(groups)){
        vars=c(vars,sample(groups[[j]],size=mtry[j],replace=FALSE))
      }

      # sampling rows
      if (is.null(class_quotas) & is.null(sample_weights)){
        trainrows=sample(rowset,replace=replacement)
      } else if (is.null(class_quotas) & !is.null(sample_weights)){
        trainrows=sample(rowset,replace=replacement,prob=sample_weights)
      } else {
        trainrows=list()
        for (l in levels.target){trainrows[[l]]=sample(rowset.classes[[l]],size=class_quotas[[l]],replace=TRUE)}
        trainrows=unlist(trainrows)
      }

      oob.rows=!(rowset %in% trainrows)

      traindata=df[trainrows,c(vars,targetvar)]
      form=as.formula(paste(targetvar, paste(vars, collapse=" + "), sep=" ~ "))

      # building and saving a tree
      tree=C5.0(form,data=traindata,control=treecontrol,rules=FALSE,trials=1,na.action=na.action)
      tree.models[[i]]=tree

      if (is.null(drop_threshold)){
        oob.sumprob[oob.rows,]=oob.sumprob[oob.rows,]+predict.C5.0(tree.models[[i]],df[oob.rows,],type='prob')
        oob.times[oob.rows]=oob.times[oob.rows]+1

        tree.models.acc[[i]]=mean(df[oob.rows,targetvar]==predict.C5.0(tree.models[[i]],df[oob.rows,],type='class'))}

      else {
        oob.acc=mean(df[oob.rows,targetvar]==predict.C5.0(tree.models[[i]],df[oob.rows,],type='class'))
        if (oob.acc>drop_threshold){
          oob.sumprob[oob.rows,]=oob.sumprob[oob.rows,]+predict.C5.0(tree.models[[i]],df[oob.rows,],type='prob')
          oob.times[oob.rows]=oob.times[oob.rows]+1
          tree.models.acc[[i]]=oob.acc
        }}
    }
    list(tree.models.list=tree.models,tree.models.acc.list=tree.models.acc,oob.sumprob=oob.sumprob,oob.times=oob.times)
  })
  stopCluster(cl)
  tree.models=Reduce('c',lapply(res_list,function(x) x$tree.models.list))
  tree.models.acc=Reduce('c',lapply(res_list,function(x) x$tree.models.acc.list))
  oob.sumprob=Reduce('+',lapply(res_list,function(x) x$oob.sumprob))
  oob.times=Reduce('+',lapply(res_list,function(x) x$oob.times))

  probs=oob.sumprob/oob.times
  pred=levels.target[max.col(probs)]
  real=eval(parse(text=paste0("df$",targetvar)))

  acc=mean(pred==real,na.rm = TRUE)
  conf_matrix=table(real,pred)

  if (is.null(drop_threshold)){
    return(structure(list(numtrees=length(tree.models),drop_threshold=0,trees=tree.models,acc.trees=unlist(tree.models.acc),oob.preds=probs,acc=acc,conf_matrix=conf_matrix,vargroups=groups,mtry=mtry,levels_target=levels.target,targetvar=targetvar,task='classification'),class='stratified_rf'))}
  else{
    if (sum(sapply(tree.models.acc,is.null))>0){
      tree.models=tree.models[-(which(sapply(tree.models.acc,is.null),arr.ind=TRUE))]
      tree.models.acc=tree.models.acc[-(which(sapply(tree.models.acc,is.null),arr.ind=TRUE))]
    }
    trees_after=length(tree.models)
    if (trees_after==0){stop('No tree had an accuracy above drop_threshold')}
    return(structure(list(numtrees=trees_after,drop_threshold=drop_threshold,trees=tree.models,acc.trees=unlist(tree.models.acc),oob.preds=probs,acc=acc,conf_matrix=conf_matrix,vargroups=groups,mtry=mtry,levels_target=levels.target,targetvar=targetvar,task='classification'),class='stratified_rf'))
  }
}

