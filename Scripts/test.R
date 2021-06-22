document()
load_all()
library(spatstat)

X <- runifgn(50, small_gn)
delta <- 0.2
h <- 0.1
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

