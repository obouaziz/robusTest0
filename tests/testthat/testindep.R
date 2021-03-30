set.seed(1)
n<-40
x<-rnorm(n)
y<-x^2+0.3*rnorm(n)
robust_res=indeptest(x,y)

test_that("indeptest on Simu1",{
  expect_equal(robust_res$statistic,stat_indeptest(x,y))
  expect_equal(round(robust_res$statistic,4),0.8696)
  expect_equal(round(robust_res$p.value,2),0.01)
  expect_equal(round(robust_res$p.value,2),round(indeptest(x,y,N=50000)$p.value,2))
})

x1<-c(0.2, 0.3, 0.1, 0.4)
y1<-c(0.5, 0.4, 0.05, 0.2)

x2<-c(0.1, 0.3, 0.5, 0.3)
y2<-c(0.5, 0.45, 0.5, 0.25)

test_that("warning for indeptest when there are ties",{
  expect_warning(indeptest(x1,y1))
  expect_warning(indeptest(x2,y1))
  expect_warning(indeptest(x1,y2))
  expect_warning(indeptest(x2,y2))
  expect_silent(indeptest(x2,y2,ties.break = "random"))
})
