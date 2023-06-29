library(ggplot2)
x=read.table('conTr5.AU.genewisell.wpl.partlh',header=T)

x$LLsum[1]=x$deltaLL[1]
for (i in (2:length(x$deltaLL))){
x$LLsum[i]=x$LLsum[i-1]+x$deltaLL[i]
}

##accumulated likelihood
df <- data.frame(order = 1:length(x$Tree1), numbers = x$LLsum, color = ifelse(x$LLsum > 0, "red", "blue"))
pdf('conTr1.accumulatedLL.pdf',width=5,height=4)
ggplot(df, aes(x = order, y = numbers, fill = color)) +
	geom_bar(stat = "identity") +
	scale_fill_identity() +
	theme_bw() +
	labs(x = "Gene", y = "accumulated LL1-LL2")

dev.off()


##likelihood per tree
df <- data.frame(order = 1:length(x$Tree1), numbers = x$deltaLL, color = ifelse(x$deltaLL > 0, "red", "blue"))
pdf('conTr1.LLcomparison.pdf',width=5,height=4)
ggplot(df, aes(x = order, y = numbers, fill = color)) +
	geom_bar(stat = "identity") +
	scale_fill_identity() +
	theme_bw() +
	labs(x = "Gene", y = "LL1-LL2")

dev.off()
