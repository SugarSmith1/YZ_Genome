# polar_concentric_with_outer_mean_labels.R
# 同心极环：Stage1→4（彩色） + 外层灰环（四时期平均），外层显示百分比数字

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(tibble)  # for tribble
})

# ---------- 1) 嵌入数据 ----------
df <- tribble(
  ~Stage, ~Type,   ~Total,  ~Expressed, ~Filtered,
  "stage1","scaf", 12107,   4602,       7505,
  "stage1","SO",   243853,  85592,      158261,
  "stage1","SS",   126589,  46273,      80316,
  "stage2","scaf", 12107,   4965,       7142,
  "stage2","SO",   243853,  93097,      150756,
  "stage2","SS",   126589,  50455,      76134,
  "stage3","scaf", 12107,   4986,       7121,
  "stage3","SO",   243853,  94278,      149575,
  "stage3","SS",   126589,  50698,      75891,
  "stage4","scaf", 12107,   5019,       7088,
  "stage4","SO",   243853,  94430,      149423,
  "stage4","SS",   126589,  50422,      76167
)

# ---------- 2) 预处理：同心环映射 ----------
df <- df %>%
  mutate(
    Stage = factor(Stage, levels = paste0("stage", 1:4)),
    ring  = as.numeric(Stage),                                # 内→外 = 1..4
    Type  = factor(Type, levels = c("SO","SS","scaf"),
                   labels = c("So","Ss","Sca"))
  )

# 四个阶段：先“表达(So,Ss,Sca)”，再“不表达(So,Ss,Sca)”
core_long <- df %>%
  pivot_longer(c(Expressed, Filtered), names_to = "Status", values_to = "Count") %>%
  mutate(
    Status = recode(Status,
                    "Expressed" = "Expression observed",
                    "Filtered"  = "No expression"),
    fill_key = paste0(Status, "|", Type)
  )

stack_levels_core <- c(
  "Expression observed|So","Expression observed|Ss","Expression observed|Sca",
  "No expression|So","No expression|Ss","No expression|Sca"
)
core_long$fill_key <- factor(core_long$fill_key, levels = stack_levels_core)

# ---------- 3) 外层灰环：四时期“平均比例”并放大到内层量级 ----------
stage_sum <- df %>%
  group_by(Stage) %>%
  summarise(Expr = sum(Expressed), Filt = sum(Filtered), .groups = "drop") %>%
  mutate(p_expr = Expr/(Expr+Filt), p_filt = 1 - p_expr)

p_expr_mean <- mean(stage_sum$p_expr)            # 算术平均比例；如需加权见注释
p_filt_mean <- 1 - p_expr_mean
mean_total_inner <- mean(stage_sum$Expr + stage_sum$Filt)  # 放大量级确保可见

outer_long <- tibble(
  ring = 5,
  Status = factor(c("Expression observed","No expression"),
                  levels = c("Expression observed","No expression")),
  Count = c(p_expr_mean, p_filt_mean) * mean_total_inner,
  fill_key = c("Expression observed|Overall","No expression|Overall")
) %>%
  mutate(
    perc   = Count / sum(Count),
    label  = scales::percent(perc, accuracy = 0.1),
    txt_col = ifelse(Status == "Expression observed", "#FFFFFF", "#333333")
  )

# # 若改为“加权平均”（按各期总量）：
# tot_expr <- sum(stage_sum$Expr); tot_filt <- sum(stage_sum$Filt)
# outer_long$Count <- c(tot_expr/(tot_expr+tot_filt), tot_filt/(tot_expr+tot_filt)) * mean_total_inner
# outer_long$perc  <- outer_long$Count / sum(outer_long$Count)
# outer_long$label <- scales::percent(outer_long$perc, accuracy = 0.1)

# ---------- 4) 颜色与图例 ----------
pal <- c(
  "Expression observed|So"     = "#E7962A",
  "No expression|So"           = "#F2D29C",
  "Expression observed|Ss"     = "#4B74B9",
  "No expression|Ss"           = "#B8CAE7",
  "Expression observed|Sca"    = "#2E6B3A",
  "No expression|Sca"          = "#AFCBB0",
  "Expression observed|Overall"= "#8C8C8C",
  "No expression|Overall"      = "#D9D9D9"
)

legend_labels <- c(
  "Expression observed|So"     = "So – Expression observed",
  "Expression observed|Ss"     = "Ss – Expression observed",
  "Expression observed|Sca"    = "Sca – Expression observed",
  "No expression|So"           = "So – No expression",
  "No expression|Ss"           = "Ss – No expression",
  "No expression|Sca"          = "Sca – No expression",
  "Expression observed|Overall"= "Overall mean – Expression observed",
  "No expression|Overall"      = "Overall mean – No expression"
)

# ---------- 5) 作图 ----------
ring_width_core  <- 0.80
ring_width_outer <- 0.70

p <- ggplot() +
  # 四个阶段（1..4）
  geom_bar(
    data = core_long,
    aes(x = ring, y = Count, fill = fill_key),
    stat = "identity", width = ring_width_core, color = "white", size = 0.5
  ) +
  # 最外层灰环（5）
  geom_bar(
    data = outer_long,
    aes(x = ring, y = Count, fill = fill_key),
    stat = "identity", width = ring_width_outer, color = "white", size = 0.6
  ) +
  # 外层百分比标签（居中）
  geom_text(
    data = outer_long,
    aes(x = ring, y = Count, label = label, color = txt_col),
    position = position_stack(vjust = 0.5),
    inherit.aes = FALSE,
    size = 4.2, fontface = "bold"
  ) +
  scale_color_identity(guide = "none") +
  coord_polar(theta = "y", clip = "off") +
  scale_x_continuous(limits = c(0.0, 6.2), breaks = 1:5) +  # 充足外缘防裁切
  scale_fill_manual(values = pal, breaks = names(legend_labels),
                    labels = unname(legend_labels),
                    name = "Gene class and status") +
  theme_void(base_size = 12) +
  labs(title = "Concentric polar rings: Stage1→4 (inner→outer) + outer gray = mean proportion") +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    plot.margin     = margin(8, 24, 8, 8)
  )

# 圈层标签（0°方向）
stage_labels <- tibble(
  ring = 1:5, y = 0,
  txt = c("Stage 1","Stage 2","Stage 3","Stage 4","Overall (mean)")
)
p <- p + geom_text(data = stage_labels, aes(x = ring, y = y, label = txt),
                   inherit.aes = FALSE, size = 4, fontface = "bold")

# ---------- 6) 导出 ----------
# ggsave("polar_concentric_with_outer_mean_labels.png", p, width = 9, height = 8, dpi = 300)
p
