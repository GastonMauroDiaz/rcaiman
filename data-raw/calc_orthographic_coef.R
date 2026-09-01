theta <- seq(0, pi/2, pi/180)
R <- seq(0,1, length.out = length(theta))
plot(theta, R, type = "n")

equi <- calc_relative_radius(theta, equidistant_coef, units = "radian")
lines(theta, equi)

R <- sin(theta)
lines(theta, R, col = 2)


fit <- lm(R ~ poly(theta, 4, raw = TRUE) - 1)
summary(fit)

validate_projection_function(signif(coef(fit), 4))
coef <- fix_projection_coef(signif(coef(fit), 4), tol = 0.0012)

or <- calc_relative_radius(theta, coef, units = "radian")
lines(theta, or, col = "blue")
