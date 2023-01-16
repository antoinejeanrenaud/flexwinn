test_that("changomics.2 works", {
  expect_error(changomics.2(NULL))
  met1<-c(rnorm(1000),NA)
  met1<-as.data.frame(met1)
  expect_warning(changomics.2(met1))
  met1<-NULL
  met1<-as.data.frame(met1)
  expect_error(changomics.2(met1))
  set.seed(333)
  met1<-c(rnorm(200),rnorm(200,sd=2)+5+seq(0,3,length.out=200))
  met1<-as.data.frame(met1)
  result<-changomics.2(met1)
  expect_type(result,"list")
})

