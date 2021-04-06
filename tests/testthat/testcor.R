test_that("Pearson correlation value, statistic and pvalue are correct",{
  expect_equal(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$estimate,as.numeric(with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1]))$estimate))
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$p.value,4),0.0174)
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$statistic,4),2.4126)
})

test_that("Kendall correlation value, statistic and pvalue are correct",{
  expect_warning(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="kendall")))
})

set.seed(1)
result=with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="spearman",ties.break="random"))

test_that("Spearman returns the same as the R version",{
  expect_warning(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="spearman")))
  expect_equal(round(result$statistic,2),2.57)
  expect_equal(round(result$p.value,2),0.01)
  expect_equal(round(result$estimate,2),0.26)
})

set.seed(1)
n<-4000
x<-rnorm(n)
y<-0.01*x^2+rnorm(n)
#result=cortest(x,y)
#result$statistic #-0.200384
#result$p.value#0.8411905
test_that("Pearson pvalue is correct on simulation 1",{
  expect_equal(round(cortest(x,y)$p.value,2),0.84)
})



set.seed(1)
n<-40
x<-rnorm(n)
y<-0.1*x^2+rnorm(n)
# result=cortest(x,y)
# result$statistic # 1.289596
# result$p.value# 0.22207

test_that("Pearson pvalue is correct on simulation 2",{
  expect_equal(round(cortest(x,y)$p.value,2),0.22)
})

set.seed(1)
n<-80
x<-rnorm(n)
y<-0.1*x^2+rnorm(n)
# result=cortest(x,y)
# result$statistic #-2.100994
# result$p.value# 0.03901

test_that("Pearson pvalue is correct on simulation 3",{
  expect_equal(round(cortest(x,y)$p.value,2),0.04)
})

set.seed(1)
n<-130
x<-rnorm(n)
y<-0.1*x^2+rnorm(n)
# result=cortest(x,y)
# result$statistic #1.161425
# result$p.value# 0.2476303

test_that("Pearson pvalue is correct on simulation 3",{
  expect_equal(round(cortest(x,y)$p.value,2),0.25)
})

set.seed(1)
n<-140
x<-rnorm(n)
y<-0.5*x^2+rnorm(n)
# result=cortest(x,y)
# result$statistic #0.3470794
# result$p.value# 0.7290604

test_that("Pearson pvalue is correct on simulation 3",{
  expect_equal(round(cortest(x,y)$p.value,2),0.73)
})

