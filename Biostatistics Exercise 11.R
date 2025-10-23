#joint PMF Gene expression  P(X=x, Y=y)
P <- matrix(
  c(0.05, 0.08, 0.12, 
    0.14, 0.19, 0.09, 
    0.22, 0.08, 0.03), 
  3, 3, byrow= TRUE
)
#add row/column names
  dimnames(P) <- list(X= c("1", "2", "3"), Y= c("1", "2", "3"))
)

#use P
P

#Matrix -> Dataframe
df <- as.data.frame(as.table(P))
install.packages("crayon")

#plot PMF
library(ggplot2)

ggplot(df, aes(x = X, y=Y, fill = Freq)) +
  geom_tile(color = "green") +
  geom_text(aes(label=round(Freq, 2)), color = "black", size=4) +
  scale_fill_gradient(low= "white", high = "steelblue") +
  labs(
    title= "Joint PMF Heatmap of Expression Levels of 
    two Chromosomes 'X and Y'", 
    x="X levels",
    y= "Y levels",
    fill = "P(X,Y)"
  ) +
  theme_minimal()
  
#marginal distribution of X and Y: 
marg_X <- rowSums(P)
marg_Y<- colSums(P)
print(marg_X)
print(marg_Y)

# Exercise C: we know from research that GENE X regulates Expression of GENE Y .
## Is it a up or down regulation?

# X up-regulates Y: high Expression of X= high Expression of Y
## -> Heatmap points from P(1,1) towards P(3.3)

#X down-regulates Y: high expression of X= low expression of Y
## -> Heatmap points from P(1,3) towards P(3,1)

#Create a plot of absolute Frequencies of Y and X Gene expression

##calculate absolute Frequencies for each Y level
Y1<- sum(P[1],P[2],P[3])
Y2<- sum(P[4],P[5],P[6])
Y3<- sum(P[7],P[8],P[9])

##calculate absolute Frequencies for each X level
X1<- sum(P[1],P[4],P[7])
X2<- sum(P[2],P[5],P[8])
X3<- sum(P[3],P[4],P[9])

#store X and Y values in a numeric vector
sum_Y <- c(Y1, Y2, Y3)
sum_X <- c(X1, X2, X3)

expression_levels <- data.frame(sum_X, sum_Y)
expression_levels$Level <- as.numeric(rownames(expression_levels))


##plot expression levels:
       
ggplot(data=expression_levels, aes(x=Level)) +
  geom_line(aes(y = sum_X, colour = "X chromosome")) +
  geom_line(aes(y= sum_Y, colour = "Y chromosome")) +
  ylim(0,1) +
  labs(
    title = "Absolute Frequencies of Chromosome Expressions", 
    x = "categorical expression levels", 
    y = "Expression frequencies") +
  theme_minimal()

