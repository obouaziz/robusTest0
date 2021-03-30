test_that("Pearson correlation value, statistic and pvalue are correct",{
  expect_equal(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$estimate,as.numeric(with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1]))$estimate))
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$p.value,4),0.0185)
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$statistic,4),2.4126)
})

test_that("Kendall correlation value, statistic and pvalue are correct",{
  expect_warning(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="kendall")))
})

