#written by Roger Mundry
#version from May 2014
#sole argument is 'model' which takes a model formula (in character mode) as handed over to, e.g., lm, glm, etc.
#which must not comprise the response
library(gtools)
all.models<-function(model){
  #browser()
  model=attr(terms(as.formula(paste(c("r", model), collapse="~"))), "term.labels")
  count.char<-Vectorize(function(xtext, char){
    sum(unlist(strsplit(xtext, split=""))==char)
  })
  #browser()
  model=model[order(count.char(model, "^"), count.char(model, ":"))]
  #build models with main effects only:
  to.treat=model[!grepl(model, pattern=":", fixed=T) & !grepl(model, pattern="^", fixed=T)]
  all.models=cbind(T, permutations(n=2, r=length(to.treat), v=c(F, T), repeats.allowed=T))
  to.treat=c(1, to.treat)
  all.models=lapply(1:nrow(all.models), function(x){to.treat[all.models[x,]]})
  #include squared terms (not involved in interactions):
  to.treat=model[!grepl(model, pattern=":", fixed=T) & grepl(model, pattern="^", fixed=T)]
  term.names=to.treat
  term.names=gsub(term.names, pattern="I(", replacement="", fixed=T)
  term.names=gsub(term.names, pattern="^2)", replacement="", fixed=T)
  if(length(to.treat)>0){
    for(i in 1:length(term.names)){
      all.models=c(all.models, lapply(all.models[unlist(lapply(all.models, function(x){length(intersect(term.names[i], x))==1}))], function(y){c(y, to.treat[i])}))
    }
  }
  #include interactions:
  to.treat=model[grepl(model, pattern=":", fixed=T)]
  #browser()
  if(length(to.treat)>0){
    for(i in 1:length(to.treat)){
      xx=unlist(strsplit(to.treat[i], split=":", fixed=T))
      yy=permutations(n=2, r=length(xx), v=c(F, T), repeats.allowed=T)
      yy=yy[-c(1, nrow(yy)),]
      yy=apply(yy, 1, function(x){paste(xx[x], collapse=":")})
      xx=yy[grepl(yy, pattern="^", fixed=T)]
      if(length(xx)>0){
        xx=c(xx, to.treat[i])
        xx=gsub(xx, pattern="I(", replacement="", fixed=T)
        xx=gsub(xx, pattern="^2)", replacement="", fixed=T)
      }
      yy=c(yy, xx)
      yy=unlist(lapply(all.models, function(x){length(intersect(x, yy))==length(yy)}))
      all.models=c(all.models, lapply(all.models[yy], function(x){c(x, to.treat[i])}))
    }
  }
  all.models=unlist(lapply(all.models, paste, collapse="+"))
  return(all.models)
}
