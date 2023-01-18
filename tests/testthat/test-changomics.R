test_that("fkPELT works", {
  expect_identical(fkPELT(rep(10,1000),NULL), numeric(0))
  expect_equal(fkPELT(c(rep(10,1000),rep(20,1000)),NULL), 1000)
  expect_identical(fkPELT(1:1000,NULL),numeric(0))
  expect_error(fkPELT(NULL,NULL))
})
test_that("changomics input", {
  expect_error(changomics(NULL))
  met1<-c(rnorm(1000),NA)
  met1<-as.data.frame(met1)
  expect_warning(changomics(met1))
  met1<-NA
  met1<-as.data.frame(met1)
  expect_error(changomics(met1))
})

  test_that("changomics works", {
  met1<-rep(15,1000)
  met1<-as.data.frame(met1)
  expect_identical(changomics(met1),met1-15)

  set.seed(333)
  met1<-c(rnorm(200),rnorm(200,sd=2)+5+seq(0,3,length.out=200))
  met1<-as.data.frame(met1)
  result<-changomics(met1,graph = TRUE)
  expect_type(result,"list")

  met1<-c(rnorm(200),rnorm(200)+5+seq(0,3,length.out=200))
  met1<-as.data.frame(met1)
  result<-changomics(met1,graph = FALSE)
  expect_type(result,"list")

  met1<-c(rep(0,200),rep(5,200))
  met1<-as.data.frame(met1)
  result<-changomics(met1,graph = FALSE)
  expect_type(result,"list")

  met1<-c(rep(0,200),rnorm(200,mean=5))
  met1<-as.data.frame(met1)
  result<-changomics(met1,graph = FALSE)
  expect_type(result,"list")
})


