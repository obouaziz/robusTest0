data_rec=read.table("Recensement.txt",header=TRUE)
dataM=data_rec[(data_rec$SEXE=="M"),]

data3=data_rec[data_rec$CATEGORIE %in% c(3,4,5,10),]

test_that("vartest on Recensement data",{
  expect_equal(round(vartest(dataM$SAL_HOR~dataM$SYNDICAT)$p.value,4),0.0107)
  expect_equal(round(vartest(data3$SAL_HOR~data3$CATEGORIE)$p.value,3),0.041)
})
