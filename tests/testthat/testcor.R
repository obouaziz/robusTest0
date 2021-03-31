test_that("Pearson correlation value, statistic and pvalue are correct",{
  expect_equal(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$estimate,as.numeric(with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1]))$estimate))
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$p.value,4),0.0185)
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
