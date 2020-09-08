# Right now, you are reading an R script (a doc with R code)

# Hashtags are used before text comments like this (green)
# This keeps them separate from your actual code 
### You can put more than one if you want!
####################  GO WILD! #######################

# Notice there are 4 panes in the R Studio window

# This pane (NW) is the script pane, where you can edit

# Notice that this script is named "GettingStarted.R" (see tab above)

# If the text in the tab is RED with a *, that means you haven't saved it
# If you don't see it in your Files (SE pane), you haven't saved a copy of your own

# Think of R files (.R) as text files that contain R code

#######################

# Put your cursor ANYWHERE on the line below and hit ctrl+enter
print("hello world")

## Notice how your output popped out on the bottom of your screen!
  
# Try a few more (just hit ctrl+enter or cmd+enter repeatedly)

9*9
9^2
1:100
1:100+1

# Notice that order of operations applies here!
6/2*(1+2)
6/(2*(1+2))

# Base R comes with many built-in calculator functions
log10(1000) # base 10 log
exp(1) # e^x
log(exp(3)) # natural log
sin(pi)

sum(1:100)
prod(1:100) # for ref, num of atoms in universe is around 6e+79
factorial(100)

choose(100,2) 
# Choose takes two arguments: if you have n things (e.g., 100 students), 
# it tells you how many unique ways there are to choose a group of size k (e.g., 2)
# For every student to meet every other student, it would take nearly 5000 pairings!
# If you forget what args a function takes, type the function name, click inside of it, and hit tab

choose()

# If you totally forget what a function does, put a question mark in front of the function name

?log

# Highlight JUST the "2+2" part below and hit ctrl+enter (or press "run")
is.numeric(2+2)

# Now try running the whole line (including the is.numeric thing)
# How was your output different?

# is.numeric() is a function that tells you if your data object is numeric or not
# Try this one

is.numeric("A")
is.character("A")

# Your results have been popping out in the console pane below (SW)
# Here you can directly run R commands!
# Type (or paste) the expression below into your console

100*101/2 #look at what this evaluates to; where have you seen it before?

# The pane in the NE is the environment pane. 
# This is where data, values, and functions will go!
# There's nothing there yet because we have created any objects!

## This thing is an assignment operator: <-
hi<- "hello world"

# It saves whatever's on the right ("hello world") as an object whose name is on the left (hi)
hi 

# See, it's still there!
# Take a look at the upper right pane (Environment tab). What do you see?

# You can also use the equals sign for assignment, but it is not considered good practice to do so
# Let's do it anyway! Let's use rep(), which is a function that takes two arguments. 
# The first argument is the thing to be repeated; the second is how many times to repeat it

hiTwice=rep(hi,2)
hiTwice

#############################################################################

# This creates a sequence of integers from 1 to 10 and stores it in the vector object "data"
data <- 1:10

# Does the same, but manually puts them in a vector with c(...)

data <- c(1,2,3,4,5,6,7,8,9,10)

# Does the same, but lets you make whatever regular sequences you want
data <- seq(from = 1, to = 10, by = 1)

# Notice that now "data" is a object in the environment pane! 
# It tells you the type, the range of elements [1:10], and the first 10 things

# Notice also that "data" never appeared in the console, it just got stored.
# To see it, just type the word data and hit enter (or run the line below)
data

# One reason why the equal sign shouldn't be used for assignment 
# is because the double equal sign (==) is used for logical tests of equivalance

2==2
2==3

# How many individual values does our "data" vector contain (i.e., how long is this vector)?
length(data)

# If we want ascending order, we can reverse it and then overwrite the old object. Let's do that.
data_rev <- rev(data)
data_rev

# Now let's remove it momentarily with the rm() function
rm(data_rev)

# To get it back, just put your cursor back on the line(s) where you assigned it (line 120)
# and hit ctrl+enter to rerun them

# Now that you've got it back, check to see what the first and last elements are
# We use brackets to index a vector (which is a list of elements of the same type)
data_rev[1]
data_rev[10]

# What's going on here?
data_rev[c(1,2,3)]

# How about here?
data_rev[data_rev>5] #what does data_rev>5 evaluate to?

# How about here?
data_rev[data>5] #what does data>5 evaluate to?

# How about here?
data_rev[c(F,F,F,F,F,T,T,T,T,T)]

# How about here?
data_rev[11]

#####################################################################################

# Let's go back to the original 1, 2, 3, ..., 10 data
# Add up all of them using the sum() function
sum(data)

# How about this one
cumsum(data)

# What about this?
cumsum(data)/sum(data)

# We can add 1 to every datapoint (many operations are vectorized)
dataPlusOne <- data+1
dataPlusOne
dataPlusOne-1

# Multiply every original datapoint by 2
doubleData <- data*2
doubleData

# OK, what does this do?
diff(doubleData)
diff(data)

# Bonus: what's going on here?!
# Note that you can nest functions (eval innermost to outermost)
diff(cumsum(data))

# Sum the vector containing the doubled data
sum(doubleData) # or sum(data*2)

# Sum of first n even numbers is n*(n+1)
length(doubleData)*(length(doubleData)+1)

# What's the average value? 
mean(doubleData)

# equivalently, sum and divide by number of elements
sum(doubleData)/length(doubleData)

# Is it the same as twice the mean of the original data?
mean(doubleData)==2*mean(data)


# How can we get all the first 10 odd numbers (from 1 to 19)? 
oddNums<-doubleData-1
oddNums

sum(oddNums)

# Interestingly, the sum of first n odd numbers is n^2
length(oddNums)^2

data=doubleData

#############################################################################

# The last pane (SE) shows your files, plots, packages, help etc.
# Let's make a quick plot
boxplot(data_rev, col="orange")

# how do the lines in this boxplot relate to the following quantities
fivenum(data)

# You can also flip it over if you want
boxplot(data1,horizontal = T)

# Or put two side-by-side
boxplot(data,doubleData)

# Or combine them together in a single vector using c()
# This concatenate function let's you combine values into a vector (if all the same type), or list (if different types)
allData<-c(data,doubleData)
allData
boxplot(allData)

## Let's make a prettier boxplot; we need to load a package first

library(ggplot2)

# If you are on your own computer, run install.packages("ggplot2") first!
# Wait for it to get done installing...It might take a second and look kinda crazy

# Now run this:
ggplot(data.frame(data=data1),aes(y=data))+geom_boxplot()+scale_x_discrete()


## Installing a package and :: instead of library()
# Many packages are already installed and ready to load for you; let's install a special one

devtools::install_github("hadley/emo")

emo::ji("sign of the horns") 

## the :: syntax allows you to use functions without loading the package that contains them
## that is, package::function() is equivalent to the following two-step process
## library(package)
## function()
