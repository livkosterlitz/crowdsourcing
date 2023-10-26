# Data point
x_data <- 2
y_data <- 3

# Identity line
x_identity <- seq(0, 5, by = 0.1)
y_identity <- x_identity

# Calculate perpendicular distance
distance <- abs(y_data - x_data) / sqrt(2)

# Calculate coordinates for the intersection point
x_intersection <- (x_data + y_data) / 2
y_intersection <- x_intersection

# Coordinates for the perpendicular line
x_perpendicular <- c(x_data, x_intersection)
y_perpendicular <- c(y_data, y_intersection)

Intersection point: (x_intersection, y_intersection) = ( (x + y) / 2, (x + y) / 2 )

# Plotting
plot(x_identity, y_identity, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "y", xlim = c(0, 5), ylim = c(0, 5))
abline(h = 0, v = 0, col = "black")
lines(x_perpendicular, y_perpendicular, col = "red", lwd = 2)
points(x_data, y_data, col = "green", pch = 19)
points(2.5, 2.5, col = "black", pch = 19)
text(x_data, y_data, " (2, 3)", pos = 2, col = "green")
segments(x_data, y_data, x_intersection, y_intersection, col = "purple", lwd = 2, lty = 2)
legend("topright", legend = c("Identity Line", "Perpendicular Line", "Data Point (2, 3)", "Perpendicular Distance"), col = c("blue", "red", "green", "purple"), lwd = c(2, 2, NA, 2), lty = c(1, 1, NA, 2))

# Coordinates of two points
x1 <- 2
y1 <- 3

x2 <- 2.5
y2 <- 2.5

# Calculate the differences
delta_x <- x2 - x1
delta_y <- y2 - y1

# Calculate the squared differences
delta_x_squared <- delta_x^2
delta_y_squared <- delta_y^2

# Sum of squared differences
sum_of_squares <- delta_x_squared + delta_y_squared

# Calculate the Euclidean distance
euclidean_distance <- sqrt(sum_of_squares)

# Print the Euclidean distance
print(euclidean_distance)

