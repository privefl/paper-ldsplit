corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))
corrT <- as(corr, "dgTMatrix")
THR_R2 <- 0.05
ind <- 81:103
library(magrittr)
df2 <- corrT[ind, ind] %>%
  { dplyr::filter(data.frame(i = .@i + ind[1], j = .@j + ind[1], r2 = .@x^2),
                  r2 >= THR_R2) }


library(ggplot2)
ggplot(df2, aes(i, j)) +
  theme_minimal(16) +
  geom_raster(aes(fill = r2), alpha = 0.6) +
  scale_fill_viridis_c(direction = -1, breaks = c(0.02, 0.1, 0.3, 0.5, 0.8)) +
  geom_text(aes(label = label),
            data = data.frame(i = c(80, 81), j = c(81, 80), label = "i")) +
  geom_text(aes(label = label), color = "blue",
            data = data.frame(i = c(80, 92), j = c(92, 80), label = "j")) +
  geom_text(aes(label = label), color = c("red", "blue", "blue"),
            data = data.frame(i = c(92.5, 100, 98), j = c(78.5, 86, 105),
                              label = c("C(i, k)", "E(i, j)", "C(j+1, k-1)"))) +
  coord_equal() +
  scale_x_continuous(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  scale_y_reverse(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  labs(x = NULL, y = NULL, fill = expression(r^2)) +
  geom_segment(aes(xend = c(80.5, 103.5), yend = c(103.5, 80.5)),
               data = data.frame(i = c(80.5, 80.5), j = c(80.5, 80.5)),
               linetype = 2) +
  geom_segment(aes(xend = c(103.5, 92.5), yend = c(92.5, 103.5)),
               data = data.frame(i = c(80.5, 92.5), j = c(92.5, 80.5)),
               linetype = 2, color = "blue") +
  geom_segment(aes(xend = c(103.5, 95.5), yend = c(95.5, 103.5)),
               data = data.frame(i = c(80.5, 95.5), j = c(95.5, 80.5)),
               linetype = 2, color = "red", alpha = 0.6) +
  geom_segment(aes(xend = c(104, 104), yend = c(104.5, 79.5)),
               data = data.frame(i = c(92.5, 80.5), j = c(104.5, 79.5)),
               linetype = 3, color = c("blue", "red"), arrow = arrow()) +
  guides(fill = guide_colorbar(barheight = 15, ticks.linewidth = 2))

# ggsave("illu.pdf", width = 6, height = 5.5)
ggsave("illu.png", width = 7, height = 6)
