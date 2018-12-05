---
title: "Final_R_Simulations_DivyaP"
author: "Divya Purohith"
date: "December 5, 2018"
output: 
  html_document: 
    keep_md: yes
---



1 - Please do the DataCamp course (all chapters) on "Introduction to Bioconductor".

**DONE**

![Datacamp](Datacamp.png)


2- In class we developed the code to perform the deterministic simulation for a simple population genetics
model with one locus with 2 alleles in a haploid. Use what you learned in class, and the following information
to write a general simulator for diploids. Your code should be flexible enough (i.e. arguments in a function
call for instance) to try different sets of parameters including initial allele frequencies, fitness of each genotype,
number of generations. Your script should do the following at a minimum:
. Report whether the beneficial allele fixes by the end of your simulation.
. Plot the trajectory of the beneficial allele across generations.
Some things you will need:
mean population fitness in a given generation (t) is
. You can also use (if you want, but do not need to) the fact that p + q = 1 for allele frequencies, and
p2 + 2pq + q2 = 1 for genotypic frequencies.
. While mean fitness of the population changes with allele frequencies (according to the equation above),
individual genotypic fitnesses do not (i.e. no frequency dependent selection)

W = p2WAA + 2pqWAa + q2Waa
p(t+1) = p2 (t) WAA/W +P(t)Q(t) WAa/W
p(t+1) is the frequency of the A allele in the population in generation t + 1.
p(t) and q(t) are the allele frequencies of A and a in generation t.
p2 (t) is the genotypic frequency for the AA diploid genotype in generation t.
WAA, WAa & Waa are the fitnesses for each genotype AA, Aa and aa respectively.



```r
p2_t1 <- function(p,W_AA, W_Aa, W_aa) {
        W_bar <- p^2*W_AA + 2*p*(1-p)*W_Aa + (1-p)^2*W_aa # mean pop fitness
		p2_t1 <- p^2*W_AA/W_bar + p*(1-p)*W_Aa/W_bar
		return(p2_t1)
}
```


```r
rep <- replicate(n = 100,p2_t1(p=0.2,W_AA=0.25,W_Aa=0.5,W_aa=0.25))
head(rep)
```

```
## [1] 0.2727273 0.2727273 0.2727273 0.2727273 0.2727273 0.2727273
```


```r
diploid_selection <- function(p0=0.2,W_AA=0.25,W_Aa=0.5,W_aa=0.25 ,n = 100) {
    
    # Initialize vectors to store allele frequencies and mean pop fitness
    p <- rep(NA,n)  # a vector to store allele frequencies
   # q <- rep(NA,n)
    
    W_bar <- rep(NA, n)
     
    # starting conditions
   	p[1] <- p0 # starting allele frequencies
   #	q[1] <- q0


   	W_bar[1] <- p[1]^2*W_AA + 2*p[1]*(1-p)[1]*W_Aa + (1-p)[1]^2*W_aa 
   	
	# now we need to loop from generation to generation
	for ( i in 2:n) {
	  
	  W_bar[i-1] <- p[i-1]^2*W_AA + 2*p[i-1]*(1-p)[i-1]*W_Aa + (1-p)[i-1]^2*W_aa 
	  
	 # mean population fitness
	   p[i] <- (p[i - 1])^2*W_AA/W_bar[i - 1] + p[i-1]*(1-p)[i-1]*W_Aa/W_bar[i-1]
	}
    return(p)
	}
   	
diploid_selection()
```

```
##   [1] 0.2000000 0.2727273 0.3372781 0.3875487 0.4237468 0.4487673 0.4657249
##   [8] 0.4771141 0.4847321 0.4898182 0.4932112 0.4954739 0.4969825 0.4979883
##  [15] 0.4986589 0.4991059 0.4994039 0.4996026 0.4997351 0.4998234 0.4998823
##  [22] 0.4999215 0.4999477 0.4999651 0.4999767 0.4999845 0.4999897 0.4999931
##  [29] 0.4999954 0.4999969 0.4999980 0.4999986 0.4999991 0.4999994 0.4999996
##  [36] 0.4999997 0.4999998 0.4999999 0.4999999 0.4999999 0.5000000 0.5000000
##  [43] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [50] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [57] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [64] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [71] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [78] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [85] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [92] 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
##  [99] 0.5000000 0.5000000
```


```r
diploid.selection <- function(p0=0.2,W_AA=0.25,W_Aa=0.5,W_aa=0.25 ,n = 100) {
# Initialize vectors to store p, delta p and mean pop fitness
    
p <- rep(NA,n) 
    
delta_p <- rep(NA, n) 
	
w_bar <- rep(NA, n)
     
# starting conditions
p[1] <- p0 # starting allele frequencies

delta_p[1] <- 0 #change in allele frequency
	
W_bar[1] <- p[1]^2*W_AA + 2*p[1]*(1-p)[1]*W_Aa + (1-p)[1]^2*W_aa 
   	
	# now we need to loop from generation to generation
	for ( i in 2:n) {
	  
	  W_bar[i-1] <- p[i-1]^2*W_AA + 2*p[i-1]*(1-p[i-1])*W_Aa + (1-p[i-1])^2*W_aa 
	  
	 # mean population fitness
	   p[i] <- (p[i - 1])^2*W_AA/W_bar[i - 1] + p[i-1]*(1-p[i-1])*W_Aa/W_bar[i-1]
	   delta_p[i] <- p[i] - p[i-1] 
	}
    
 if (any(p > 0.9999)) {
    fixation <- min(which.max(p > 0.9999))
    cat("fixation for A1 occurs approximately at generation:", fixation )	
    } else {
        maxAlleleFreq <- max(p)
    	cat("fixation of A1 does not occur, max. allele frequency is:", print(maxAlleleFreq, digits = 3) )
    }
}
```


```r
diploid_fixation <- diploid_selection(p0 = 0.0001, W_AA = 1, W_Aa = 0.987, W_aa = 0.25,n = 1000) #The allele is not fixed
head(diploid_fixation)
```

```
## [1] 0.0001000000 0.0003945679 0.0015541473 0.0060802271 0.0231782239
## [6] 0.0806407565
```



```r
plt <- diploid_selection()
gen <- 1:length(plt)
plot(plt ~ gen, pch=20,
   ylab = "allele frequency", 
     xlab = "generation",
   xlim = c(1,100))
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


3. Genetic drift simulator. 

. Number of alleles in population (i.e. for diploids 2 times the number of individuals). #DONE
. Starting allele frequency #DONE
. Number of generations. #DONE
. plot allele frequency changes over time #DONE

#sample(x,size,replace)


```r
allele_freq <- function(A, alleles, gen) {
  allele <- sample(c("a", "A"), size = 40, replace = TRUE, prob = c(A, 1-A))
  
plot(1,0, type = "n",  xlim=c(1,40), ylim=c(0,20), xlab="number of generations", ylab="number of A alleles")

for(i in 1:gen) {
  alleles <- sample(allele, 40, replace = TRUE)
  
points(i, length(alleles[alleles=="A"]), pch = 18, col = "red")
}
}

allele_freq(0.7,14,40)
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

#OR


```r
allele_freq <- function(A, alle, gen) {
  allele_freq <- rep(NA,gen) #vector to store allele frequencies
  allele_freq[1] <- A
  
  allele_counts <- sample(c("a", "A"), size = alle, replace = TRUE, prob = c(A, 1-A))
plot(1,0, type = "n",  xlim=c(1,40), ylim=c(0,20), xlab="number of generations", ylab="number of A alleles")

for(i in 2:gen) {
  allele_counts <- sample(c("a","A"), size =  alle, replace = TRUE, prob = c(allele_freq[i-1],(1-allele_freq[i-1])))
  
  allele_freq[i] <- length(allele_counts[allele_counts=="A"])/length(allele_counts)
  
points(i, length(allele_counts[allele_counts=="A"]), pch = 18, col = "red")
}
}
allele_freq(0.7,14,40)
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

4 - Repeating your stochastic simulation. Now that you have a working simulator of genetic drift, you want
to use it to assess how likely it is for the allele to be lost (p = 0) after a certain number of generations of
drift (we will use 100 generations, but your function should be flexible). Using your function (you can modify
it if you need to), perform a simulation 1000 times each with starting allele frequencies of p = f(A) of 0.5,
0.25 and 0.1, with 200 diploid individuals in the population each generation. Have your function record the
proportion (out of 1000) of simulated runs that the A allele is lost from the population (p = 0).


```r
set.seed(720)
allele_freq <- function(A, alle, gen) {
  allele_freq <- rep(NA,gen) #Vecotr to store allele frequencies
  allele_freq[1] <- A
  
  allele_counts <- sample(c("a", "A"), size = alle, replace = TRUE, prob = c(A, 1-A))
#plot(1,0, type = "n",  xlim=c(1,40), ylim=c(0,20), xlab="number of generations", ylab="number of A alleles")

for(i in 2:gen) {
  allele_counts <- sample(c("a","A"), size =  alle, replace = TRUE, prob = c(allele_freq[i-1],(1-allele_freq[i-1])))
  
  allele_freq[i] <- length(allele_counts[allele_counts=="A"])/length(allele_counts)
  
#points(i, length(allele_counts[allele_counts=="A"]), pch = 18, col = "red")
}
  return(allele_freq[length(allele_freq)])
}

#Run for 100 generations & 200 individuals with different allele frequencies, one easy way to know is using replicate function.

mean(replicate(1000, allele_freq(0.5,200,100)) == 0)
```

```
## [1] 0.075
```

```r
mean(replicate(1000, allele_freq(0.25,200,100)) == 0)
```

```
## [1] 0.007
```

```r
mean(replicate(1000, allele_freq(0.1,200,100)) == 0)
```

```
## [1] 0
```


5 - Write some code that allows you to plot the allele trajectories for drift for 100 of the simulations starting
at p = 0.5. hint: I showed you an example of how to draw multiple lines with a single function last week. . .
It could look something like this.



```r
allele_freq <- function(A, alle, gen) {
  allele_freq <- rep(NA,gen)
  allele_freq[1] <- A
  
  allele_counts <- sample(c("a", "A"), size = alle, replace = TRUE, prob = c(A, 1-A))
#plot(1,0, type = "n",  xlim=c(1,40), ylim=c(0,20), xlab="number of generations", ylab="number of A alleles")

for(i in 2:gen) {
  allele_counts <- sample(c("a","A"), size =  alle, replace = TRUE, prob = c(allele_freq[i-1],(1-allele_freq[i-1]))) #freq of allels change
  
  allele_freq[i] <- length(allele_counts[allele_counts=="A"])/length(allele_counts)
  
  }
return(allele_freq)
}
allele_freq(0.5,200,100)
```

```
##   [1] 0.500 0.480 0.550 0.450 0.580 0.370 0.620 0.370 0.580 0.395 0.585
##  [12] 0.380 0.700 0.315 0.690 0.295 0.670 0.280 0.745 0.235 0.815 0.190
##  [23] 0.795 0.225 0.810 0.165 0.825 0.170 0.810 0.200 0.790 0.185 0.840
##  [34] 0.165 0.795 0.205 0.810 0.135 0.850 0.160 0.865 0.130 0.870 0.160
##  [45] 0.820 0.170 0.840 0.145 0.835 0.130 0.830 0.180 0.805 0.205 0.780
##  [56] 0.215 0.745 0.285 0.725 0.335 0.665 0.270 0.710 0.260 0.770 0.205
##  [67] 0.780 0.230 0.725 0.280 0.705 0.265 0.745 0.225 0.760 0.220 0.775
##  [78] 0.205 0.810 0.200 0.760 0.260 0.720 0.285 0.695 0.315 0.685 0.305
##  [89] 0.640 0.320 0.705 0.310 0.715 0.295 0.730 0.260 0.805 0.280 0.750
## [100] 0.270
```

```r
ran_colors <- sample(colors(), size = 100)  #random colors

p <- allele_freq(A = 0.5, alle = 200,gen = 100)
generations <- 1:length(p) #same as did in 2nd QA
plot(p ~ generations, lty = 1, lwd = 3, pch = 20,
     ylab = "Frequency of alleles",
     xlab = "Number of generations",
     ylim= c(0,1), xlim= c(0,100))

for(i in 1:100){
  lines(allele_freq(A = 0.5, alle = 200,gen = 100), type = "l", col= sample(ran_colors, rep = F, size = 1)) #replac is False,so one color is not replaced with other
}
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


6. A simple power analysis. One common use for stochastic simulations is to evaluate the statistical power of
a particular method or model given some known parameters such as effect size (like the slope of the regression
line), sample size, and residual standard error (the uncertainty in your estimator). While I kind of detest
teaching anything with p-values (it is not a great way of evaluating statistical models) we will use them for
this exercise.
We saw in class how to simulate data for the relationship between y and x and plot the regression fits (best fit
lines) for the simulated data. You are going to use that same approach to determine how often you observe a
"significant" p-value (p < 0.05) under some different scenarios.
In class we did something approximately like this (with a few changes)


#Given function


```r
x <- seq(from =1, to = 10, length.out = 20) # length.out is how many observations we will have
a <- 0.5 # intercept
b <- 0.1 # slope
y_deterministic <- a + b*x
y_simulated <- rnorm(length(x), mean = y_deterministic, sd = 2)
mod_sim <- lm(y_simulated ~ x)
p_val_slope <- summary(mod_sim)$coef[2,4] # extracts the p-value
p_val_slope
```

```
## [1] 0.8249865
```


#reset the formula


```r
diff_values <- function(a = int, b = slope, samplesize, rse) {
  x <- seq(from =1, to = 10, length.out = samplesize)
y_deterministic <- a + b*x
y_simulated <- rnorm(length(x), mean = y_deterministic, sd = rse)
mod_sim <- lm(y_simulated ~ x)
p_val_slope <- summary(mod_sim)$coef[2,4] # extracts the p-value
return(p_val_slope)
}
diff_values(0.5,0.1,20,2)
```

```
## [1] 0.3659154
```

#Run 1000 times and generate histogram 


```r
set.seed(720)
diff_values_1 <- replicate(1000, diff_values(0.5, 0.1, 20, 2))

hist(diff_values_1)
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

#Check the p value


```r
p_value <- (table(diff_values_1 <0.05)/length(diff_values_1)) [["TRUE"]]
p_value #The p-value is more than 0.05
```

```
## [1] 0.107
```

#Redo this by giving slope 0


```r
set.seed(720)
diff_values_2 <- replicate(1000, diff_values(0.5, 0, 20, 2))
hist(diff_values_2)
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

#Check the p value

```r
p_value_1 <- (table(diff_values_2<0.05)/length(diff_values_2)) [["TRUE"]]
p_value_1 
```

```
## [1] 0.061
```

Based on this (y_deterministic <- a + b*x)
If the the slope(b) is set zero and the intercept(a) 0.5, the output will be a straight line, Thats because when the slope is zero and if it interacts with x that also becomes zero. then y_deterministic value becomes constant showing no significance.


. Finally, using either a for loop or an Rish method (i.e. one of the apply family of functions) try a grid of
sample sizes from 10 to 100 by 5 (i.e. 10, 15, 20. . . ., 95, 100) with the slope being fixed at 0.1, intercept
= 0.5 and residual standard error as 1.5. At each different sample size, run the simulation 100 times
and report the frequency at which the p value is less than 0.05. What pattern do you observe. Plotting
the proportion of p-values less than 0.05 on the Y axis VS sample size (X axis) may be very helpful
3

#Used apply function (sapply)


```r
sample_grid <- seq(10,100, 5) #by 5 for sample size
fixed_sample_grid <- function(size, times = 100) {
  final_values <- replicate(100, diff_values(0.5, 0.1, size, 1.5))
  
  #for (i in 1:100) {
    fin_values <- diff_values(0.5,0.1,size,1.5)
  #}
  less_values <- (table(final_values <0.05)/length(final_values)) [["TRUE"]]
}
plot_x <- sapply(sample_grid, fixed_sample_grid, times = 100) #p value less than 0.05
head(plot_x)
```

```
## [1] 0.03 0.09 0.12 0.13 0.12 0.23
```


```r
plot(x = sample_grid, y = plot_x, main = "P-value of the different sample size", pch = 16, col = "red")
```

![](Final_R_Simulations_Assignment_DP_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

