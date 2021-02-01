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
  geom_raster(aes(fill = r2), alpha = 0.6) +
  scale_fill_viridis_b(direction = -1, breaks = c(0.02, 0.1, 0.3, 0.8)) +
  geom_text(aes(label = label),
            data = data.frame(i = c(80, 80, 81, 91), j = c(81, 91, 80, 80),
                              label = c("i", "j", "i", "j"))) +
  geom_text(aes(label = label), color = c("red", "blue", "blue"),
            data = data.frame(i = c(91.5, 99, 98), j = c(78.5, 86, 105),
                              label = c("C(i, k)", "E(i, j)", "C(j+1, k-1)"))) +
  coord_equal() +
  scale_x_continuous(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  scale_y_reverse(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  labs(x = NULL, y = NULL, fill = expression(r^2)) +
  geom_segment(aes(xend = c(80.5, 103.5, 103.5, 91.5), yend = c(103.5, 80.5, 91.5, 103.5)),
               data = data.frame(i = c(80.5, 80.5, 80.5, 91.5), j = c(80.5, 80.5, 91.5, 80.5)),
               linetype = 2) +
  geom_segment(aes(xend = c(104, 104), yend = c(104.5, 79.5)),
               data = data.frame(i = c(91.5, 80.5), j = c(104.5, 79.5)),
               linetype = 3, color = c("blue", "red"), arrow = arrow())

ggsave("illu.pdf", width = 6, height = 5.5)
