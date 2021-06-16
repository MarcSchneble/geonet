library(spatstat)
document()
load_all()

X <- runifgn(50, small_gn)
delta <- 0.2
h <- 0.1
r <- 1
model <- intensity_pspline(X, delta = delta, h = h, r = 2)
summary(model)
plot(model)

X <- as_gnpp(chicago)
model <- intensity_pspline(X, delta = 10, h = 2, r = 2)
plot(model)


X <- as_gnpp(chicago)
delta <- 10
h <- 2
r <- 2
formula <- X ~ marks + internal(x) + internal(y)

start <- Sys.time()
model <- intensity_pspline(X, formula, delta = delta, h = h, r = r,
                           scale = list(x = 1/1000, y = 1/1000))
print(Sys.time() - start)

summary(model)
plot(model)


X <- montgomery
delta <- 0.05
h <- 0.025
r <- 2
formula <- ~ s(hour, k = 15) + internal(type) + internal(direction)
model <- intensity_pspline(X, formula, delta = delta, h = h, r = r)
sum(model$edf[model$ind[["hour"]]])
summary(model)
plot(model)


