## Playing in R to understand basic R functions and objects####

x <- 5
students <- data.frame(
  name = c("A","B","C"),
  height = c(62,68,72)
)

head(students)
class(students)
str(students)

students$name[2]

students[1,2]

mean(students$height)