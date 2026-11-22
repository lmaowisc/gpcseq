# Exact cumulative-history multiset guard. Retains typed columns per channel.
# Joint radix order gives O(m log m) ordering with O(m) temporary storage;
# adjacent equality preserves duplicates and exact numeric values.
.gpc_history_records <- function(x,method) {
  normalize<-function(tab,fields) {
    if(is.null(tab) || !nrow(tab)) return(data.frame())
    for(f in setdiff(fields,names(tab))) tab[[f]]<-NA_real_
    tab[,fields,drop=FALSE]
  }
  if(method=='aj') return(list(mode='aj',tables=list(transitions=normalize(x$tr,c('id','time','from','to','available')))))
  if(is.null(x$channels)) {
    e<-x$e
    x$channels<-list(events=list(hospitalization=e[e$type=='hospitalization',c('id','time',intersect('available',names(e))),drop=FALSE]),
      measurements=list(biomarker=e[e$type=='biomarker',c('id','time','value',intersect('available',names(e))),drop=FALSE]))
  }
  tabs<-list()
  for(group in c('events','measurements')) for(name in names(x$channels[[group]])) {
    key<-paste0(group,':',nchar(enc2utf8(name),type='bytes'),':',enc2utf8(name))
    fields<-c('id','time',if(group=='measurements') 'value','available')
    tabs[[key]]<-normalize(x$channels[[group]][[name]],fields)
  }
  list(mode='named',tables=tabs)
}
# Canonical typed columns; no decimal conversion or row-string registry.
.gpc_column_table <- function(tab) {
  tab<-tab[,sort(names(tab)),drop=FALSE]
  for(n in names(tab)) {
    v<-tab[[n]]
    if(is.factor(v))v<-as.character(v)
    if(is.character(v))v<-enc2utf8(v)
    if(is.numeric(v)) {v<-as.double(v);v[v==0]<-0}
    if(!is.numeric(v) && !is.character(v) && !is.logical(v))stop('Unsupported field')
    tab[[n]]<-v
  }
  rownames(tab)<-NULL
  tab
}
.gpc_column_contains <- function(previous,current) {
  if(!nrow(previous))return(TRUE)
  if(!identical(names(previous),names(current)))return(FALSE)
  # Reject incompatible nonmissing types rather than coercing numeric to text.
  for(n in names(previous)) {
    a<-previous[[n]];b<-current[[n]]
    if(typeof(a)!=typeof(b) && !all(is.na(a)) && !all(is.na(b)))return(FALSE)
  }
  joined<-rbind(previous,current)
  ix<-do.call(order,c(unname(joined),list(method='radix',na.last=TRUE)))
  different<-rep(FALSE,length(ix)-1L)
  for(v in joined) {
    a<-v[ix[-length(ix)]];b<-v[ix[-1L]]
    equal<-(is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&a==b)
    different<-different | !equal
  }
  group<-cumsum(c(TRUE,different))
  n<-max(group)
  all(tabulate(group[ix<=nrow(previous)],n)<=tabulate(group[ix>nrow(previous)],n))
}
.gpc_history_check <- function(old,now) {
  counts<-lapply(now$tables,.gpc_column_table)
  if(!is.null(old)) {
    if(!identical(old$mode,now$mode))stop('History representation changed')
    for(name in names(old$counts)) {
      previous<-old$counts[[name]]
      if(!nrow(previous))next
      # Saved 0.3.0.9008 keys may contain the byte count of a Latin-1 name.
      # Recompute that prefix after UTF-8 conversion without changing the name.
      group<-sub(":.*$", "", name)
      channel<-enc2utf8(sub("^[^:]+:[0-9]+:", "", name))
      key<-if(old$mode=='named') paste0(group, ':', nchar(channel,type='bytes'), ':', channel) else name
      current<-counts[[key]]
      if(is.null(current) || !.gpc_column_contains(previous,current))
        stop('Previously observed history was removed or changed in ',name)
    }
  }
  list(mode=now$mode,counts=counts)
}
