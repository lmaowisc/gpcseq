test_that('cumulative histories are retained across methods and spending rules', {
 for(method in c('shared','ipcw','aj')) for(spend in list(c(.01,.05),'obrien-fleming','pocock','power','hsd')) {
  a<-gpc_example(1);b<-gpc_example(2)
  args<-list(looks=c(1,2),spending=spend,method=method)
  if(method!='shared')args$horizons<-c(1,2)
  if(method=='aj') {args$state_order<-c('well','hospitalized','dead');args$absorbing<-'dead';args$death_state<-'dead'}
  d<-do.call(gpc_design,args)
  update<-function(o,x)if(method=='aj')gpc_update(o,x$patients,transitions=x$transitions) else gpc_update(o,x$patients,events=x$events)
  g<-update(d,a);before<-serialize(g,NULL)
  bad<-b;bad$events<-b$events[0,];bad$transitions<-b$transitions[b$transitions$id!='C2',]
  expect_error(update(g,bad),'history was removed')
  expect_identical(serialize(g,NULL),before)
  expect_s3_class(update(unserialize(before),b),'gpcseq')
  legacy<-g;legacy$history_registry<-NULL
  expect_warning(res<-update(legacy,b),'Earlier cuts cannot be verified')
  ordinary<-update(g,b)
  expect_identical(res,ordinary)
 }
})
test_that('native and named reserved channels preserve histories', {
 a<-gpc_example(1);b<-gpc_example(2)
 prep<-function(x)gpc_trial_data(x$patients,events=list(hospitalization=x$events[,c('id','time')]),event='hospitalization')
 d<-gpc_design(c(1,2),c(.01,.05),'ipcw',horizons=c(1,2))
 n<-gpc_update(d,a$patients,a$events);p<-gpc_update(d,data=prep(a))
 expect_equal(gpc_results(gpc_update(n,data=prep(b))),gpc_results(gpc_update(p,b$patients,b$events)))
})
test_that('mixed dated channels retain dates and exact multiplicities', {
 a<-gpc_example(1);b<-gpc_example(2)
 prep<-function(x) {
  e<-x$events[,c('id','time')];e$available<-e$time
  gpc_trial_data(x$patients,events=list(hospitalization=e),measurements=list(marker=data.frame(id=x$patients$id,time=0,value=1)),event='hospitalization')
 }
 d<-gpc_design(c(1,2),c(.01,.05),'ipcw',horizons=c(1,2))
 g<-gpc_update(d,data=prep(a));good<-prep(b)
 expect_s3_class(gpc_update(g,data=good),'gpcseq')
 bad<-good;bad$channels$events$hospitalization$available[1]<-.75
 expect_error(gpc_update(g,data=bad),'history was removed')
 bad<-good;bad$channels$measurements$marker$value[1]<-2
 expect_error(gpc_update(g,data=bad),'history was removed')
 tab<-data.frame(id='a',time=1)
 expect_false(gpcseq:::.gpc_column_contains(tab[c(1,1),],tab))
 expect_true(gpcseq:::.gpc_column_contains(tab,tab[c(1,1),]))
})

test_that('channel encoding does not change cumulative history identity', {
 utf<-'caf\u00e9:1';latin<-iconv(utf,'UTF-8','latin1')
 expect_identical(utf,latin)
 prep<-function(k,name) {
  x<-gpc_example(k)
  gpc_trial_data(x$patients,events=setNames(list(x$events[,c('id','time')]),name),event=name,
   measurements=setNames(list(data.frame(id=x$patients$id,time=0,value=1)),name))
 }
 for(nm in list(c(latin,utf),c(utf,latin))) {
  d<-gpc_design(c(1,2),c(.01,.05),'ipcw',horizons=c(1,2))
  g<-gpc_update(d,data=prep(1,nm[1]))
  reference<-gpc_update(g,data=prep(2,nm[1]))
  result<-gpc_update(unserialize(serialize(g,NULL)),data=prep(2,nm[2]))
  expect_identical(gpc_results(result),gpc_results(reference))
  expect_identical(gpc_covariance(result),gpc_covariance(reference))
  # Reconstruct exactly the pre-normalization keys saved by 0.3.0.9008.
  names(g$history_registry$counts)<-paste0(c('events','measurements'),':',nchar(latin,type='bytes'),':',latin)
  expect_identical(gpc_results(gpc_update(g,data=prep(2,nm[2]))),gpc_results(reference))
  bad<-prep(2,nm[2]);bad$channels$measurements[[1]]$value[1]<-2
  expect_error(gpc_update(g,data=bad),'history was removed')
  bad<-prep(2,nm[2]);bad$channels$events[[1]]<-bad$channels$events[[1]][0,]
  expect_error(gpc_update(g,data=bad),'history was removed')
 }
})
