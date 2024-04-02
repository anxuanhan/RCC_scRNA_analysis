df<- read.csv("cellcluster_count.csv")



library(ggplot2)
library(ggrepel)
library(tidyverse)


# Get the positions
my_colors <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#542788", "#807dba", "#c6dbef", "#efedf5", "#7CAEFA")
df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(freq))), 
         pos = freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), freq/2, pos))

p <- ggplot(df, aes(x = "" , y = freq, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = my_colors) +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(freq, "%")),
                   size = 3, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "group")) +
  theme_void()

ggsave(p, file = "pie_plot.pdf", dpi = 100, width = 7, height = 8)

dev.off()


# 加载gplots包
if(!require(gplots)) {
  install.packages("gplots")
}
library(gplots)

# 定义数据
overlap <- list(
  CCL2 = 1:165,
  CD74 = c(1:50, 1:1157),
  CA9 = c(1:65, 1:514)
)

# 绘制韦恩图
v <- venn(overlap)

# 定制韦恩图颜色
venn(overlap, small = 0.7, col = c("blue", "green", "red"), fill = c("skyblue", "pink", "mediumorchid"), cex = 1.5, cat.col = c("blue", "green", "red"),)
