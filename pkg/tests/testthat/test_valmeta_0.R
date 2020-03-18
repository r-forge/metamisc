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


test_that("Preparation of untransformed aggregate data", {
  
  # Create test dataset provided by Anneke Damen
  testdat <- data.frame(Authoryear = c("Andersson 2015", "Chia 2014", "Goff 2014", "Goff 2014", "Jung 2015", 
                                       "Kavousi 2014", "Khalili 2015", "Lee 2015", "Muntner 2014"),
                   Samplesize = c(3396, 307, 5041, 735, 114622, 1513, 2353, 679, NA),
                   N_events = c(284, 22, 539, 107, 7669, 192, 200, 80, 376),
                   Cstat = c(0.72, 0.55, 0.6843, 7669, 0.727, 0.67, 0.74, 0.714, 0.65),
                   Cstat95LB = c(NA, NA, NA, NA, 0.721, 0.630, 0.710, 0.567, 0.620),
                   Cstat95UB = c(NA, NA, NA, NA, 0.734, 0.710, 0.780, 0.770, 0.680),
                   CstatSE = c(0.014, 0.050, NA, NA, NA, NA, NA, NA, NA))
  
  ad <- ccalc(cstat = Cstat, cstat.cilb = Cstat95LB, cstat.ciub = Cstat95UB, 
              cstat.se = CstatSE, N = Samplesize, o = N_events, data = testdat, 
              slab = Authoryear)
  
  # Check whether the prepared aggregate data for Jung 2015 is correct
  expect_equal(as.character(ad["Jung 2015","theta.se.source"]), "Confidence Interval")
  expect_equal(ad["Jung 2015","theta.cilb"], subset(testdat, Authoryear == "Jung 2015")$Cstat95LB)
  expect_equal(ad["Jung 2015","theta.ciub"], subset(testdat, Authoryear == "Jung 2015")$Cstat95UB)
  
  # Check whether the prepared aggregate data for Kavousi 2014 is correct
  expect_equal(as.character(ad["Kavousi 2014","theta.se.source"]), "Confidence Interval")
  expect_equal(ad["Kavousi 2014","theta.cilb"], subset(testdat, Authoryear == "Kavousi 2014")$Cstat95LB)
  expect_equal(ad["Kavousi 2014","theta.ciub"], subset(testdat, Authoryear == "Kavousi 2014")$Cstat95UB)
  
  # Check whether the prepared aggregate data for Nuntner 2014 is correct
  expect_equal(as.character(ad["Muntner 2014","theta.se.source"]), "Confidence Interval")
  expect_equal(ad["Muntner 2014","theta.cilb"], subset(testdat, Authoryear == "Muntner 2014")$Cstat95LB)
  expect_equal(ad["Muntner 2014","theta.ciub"], subset(testdat, Authoryear == "Muntner 2014")$Cstat95UB)
  
})


