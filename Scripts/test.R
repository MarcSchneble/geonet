document()
load_all()
library(spatstat)

X <- runifgn(50, small_gn)
delta <- "0"
h <- 0.2
r <- 2
formula <- ~1
model <- intensity_pspline(X, delta = delta, h = h, r = r)

X <- runifgn(50, small_gn)
delta <- 0.2
h <- 0.1
r <- 1
model <- intensity_pspline(X, delta = delta, h = h, r = 2, verbose = TRUE,
                           control = list(rho_start = 20,
                                          eps_theta = 1e-4,
                                          rho_max = 100))
summary(model)
plot(model)

X <- as_gnpp(chicago)
model <- intensity_pspline(X, formula = ~internal(y), scale = list(y = 1/1000),
                           delta = 10, h = 2, r = 2)
plot(model)
summary(model)


X <- as_gnpp(chicago)
delta <- 10
h <- 2
r <- 2
formula <- X ~ marks + x + y

start <- Sys.time()
model <- intensity_pspline(X, formula, delta = delta, h = h, r = r,
                           scale = list(x = 1/1000, y = 1/1000),
                           verbose = TRUE)
print(Sys.time() - start)

summary(model)
plot(model)


X <- montgomery
delta <- 0.1
h <- 0.05
r <- 2
formula <- ~ type + direction
formula <- ~1
model <- intensity_pspline(X, formula = formula, delta = delta, h = h, r = r, verbose = TRUE)
summary(model)
plot(model)


X <- runifgn(20, small_gn)
fit1 <- intensity_pspline(X, delta = 0.2, h = 0.1)
fit2 <- intensity_pspline(X, delta = 0.5, h = 0.1)
network_intensity(fit1, fit2)

X <- as_gnpp(chicago)
fit1 <- intensity_pspline(X)
X2 <- rgnpp(1000, fit)
fit2 <- intensity_pspline(X2, delta = 10, h = 2)
network_ISE(fit1, fit2)
