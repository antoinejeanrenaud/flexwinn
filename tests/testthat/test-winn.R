test_that("winn input", {
  expect_error(winn(NULL,100))
  met1<-c(rnorm(1000),NA)
  met1<-as.data.frame(met1)
  expect_error(winn(met1,c(100,300,600,800)))
  met1<-NA
  met1<-as.data.frame(met1)
  expect_error(winn(met1))
})

test_that("winn works", {
  set.seed(333)
  met1<-c(rnorm(200),rnorm(200,sd=2)+5+seq(0,3,length.out=200))
  met1<-as.data.frame(met1)
  result<-winn(met1,end.plates=200,graph = TRUE)
  expect_type(result,"list")

  met1<-c(rnorm(200),rnorm(200)+5+seq(0,3,length.out=200))
  met1<-as.data.frame(met1)
  result<-winn(met1,end.plates=200,graph = FALSE)
  expect_type(result,"list")

  met1<-c(rep(0,200),rep(5,200))
  met1<-as.data.frame(met1)
  result<-winn(met1,end.plates=200,graph = FALSE)
  expect_type(result,"list")

  met1<-c(rep(0,200),rnorm(200,mean=5))
  met1<-as.data.frame(met1)
  result<-winn(met1,end.plates=200,graph = FALSE)
  expect_type(result,"list")
})
