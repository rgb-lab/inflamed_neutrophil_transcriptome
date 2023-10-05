library(lme4)

expression <- read.csv("expression.csv", row.names = 1)
modeling <- function(df, parameter){
    lm0 <- lmer(get(parameter) ~ 1 + condition + (1|batch), data = df)
    lm1 = lmer(get(parameter) ~ 1 + (1|batch), data = df)
    a <- anova(lm0, lm1)
    return(a[[8]][[2]])
}

p_values = c()

for(gene in colnames(expression)){
    append(p_values, modeling(expression, gene))
}

names(p_values) <- colnames(expression)

saveRDS(p_values, "p_values.rds")