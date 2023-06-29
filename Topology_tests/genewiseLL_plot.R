library(ggplot2)
df <- data.frame(order = 1:length(x$L.T1.T2.), numbers = x$L.T1.T2., color = ifelse(x$L.T1.T2. > 0, "red", "blue"))

pdf('conTr1.LLcomparison.pdf',width=5,height=4)
ggplot(df, aes(x = order, y = numbers, fill = color)) +
	geom_bar(stat = "identity") +
	scale_fill_identity() +
	theme_bw() +
	labs(x = "Gene", y = "LL1-LL2")

dev.off()
