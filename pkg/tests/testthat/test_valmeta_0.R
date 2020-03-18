context("valmeta 0. calculation of c-statistic")
skip_on_cran()

data(EuroSCORE)

test_that("Transformation of the c-statistic", {
  
  # Logit transformation
  logitc1 <- log(EuroSCORE$c.index/(1-EuroSCORE$c.index))
  logitc2 <- logit(EuroSCORE$c.index) #built-in function
  logitc3 <- (ccalc(cstat=c.index, cstat.se=se.c.index, data=EuroSCORE, g="log(cstat/(1-cstat))"))$theta
  expect_equal(logitc1, logitc2, logitc3)
  
  # Back-transformation
  ilogitc1 <- 1/(1+exp(-logitc1))
  ilogitc2 <- inv.logit(logitc2)
  expect_equal(ilogitc1, ilogitc2)

  # Arcsine square root transformed proportion
  ft1 <- asin(sqrt(EuroSCORE$c.index))
  ft2 <- (ccalc(cstat=c.index, cstat.se=se.c.index, data=EuroSCORE, g="asin(sqrt(cstat))"))$theta
  expect_equal(ft1, ft2)
  
})

test_that("Standard error of transformed c-statistic", {
  
  # Logit transformation
  logitcse1 <- EuroSCORE$se.c.index/(EuroSCORE$c.index*(1-EuroSCORE$c.index))
  logitcse2 <- (ccalc(cstat=c.index, cstat.se=se.c.index, data=EuroSCORE, g="log(cstat/(1-cstat))"))$theta.se
  expect_equal(logitcse1, logitcse2)
  
})


test_that("Anneke", {
  
  # Create test dataset provided by Anneke Damen
  testdat <- data.frame(Authoryear = c("Andersson 2015", "Chia 2014", "Goff 2014", "Goff 2014", "Jung 2015"),
                   Samplesize = c(3396, 307, 5041, 735, 114622),
                   N_events = c(284, 22, 539, 107, 7669),
                   Cstat = c(0.72, 0.55, 0.6843, 7669, 0.727),
                   Cstat95LB = c(NA, NA, NA, NA, 0.721),
                   Cstat95UB = c(NA, NA, NA, NA, 0.734),
                   CstatSE = c(0.014, 0.050, NA, NA, NA))
  
  ad <- ccalc(cstat = Cstat, cstat.cilb = Cstat95LB, cstat.ciub = Cstat95UB, 
              cstat.se = CstatSE, N = Samplesize, o = N_events, data = testdat, 
              slab = Authoryear)
  
  expect_equal(as.character(ad["Jung 2015","theta.se.source"]), "Confidence Interval")
  expect_equal(ad["Jung 2015","theta.cilb"], subset(testdat, Authoryear == "Jung 2015")$Cstat95LB)
  expect_equal(ad["Jung 2015","theta.ciub"], subset(testdat, Authoryear == "Jung 2015")$Cstat95UB)
  
})


