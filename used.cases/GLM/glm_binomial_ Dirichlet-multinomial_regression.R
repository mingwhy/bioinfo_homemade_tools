
####################################################
#https://julianfaraway.github.io/faraway/
#Chap 2 Binomial Data from Julian Faraway's stat book: `Extending the Linear Model with R` 
#BiocManager::install('faraway')
library(faraway)
library(tidyverse)
library(emmeans)
data(orings) #from faraway package
orings #23 x 2 data.frame
with(orings,plot(damage/6~temp,pch=16,xlab='Temperature',ylab='Prob of damage'))

# bionimal regression model: glm( (success,fail) ~ predictors)
logitmod<-glm(cbind(damage,6-damage) ~ temp,
              family=binomial(link='logit'),orings)
summary(logitmod)

########################################################################
# https://www.nxn.se/valent/2020/11/28/s9jjv32ogiplagwx8xrkjk532p7k28
# Cell type abundance analysis in scRNA-seq
# github link: https://github.com/vals/Blog/tree/master/201127-cell-count-glm
library(tidyverse)
library(ggrepel)
library(emmeans)

cell_type_counts = read_csv('cell_type_counts.csv') %>% mutate(batch = as.character(batch)) 
cell_type_counts

cell_type_counts %>% group_by(Health, batch, Location, Subject)
cell_type_counts %>% filter(Total > 0) -> df_a

df_a 
unique(df_a$Cluster) #51 cell clusters

model0 <- glm(
  formula = cbind(Count, Other) ~ Cluster, #only use Cluster as predictor
  family = binomial(link = 'logit'),
  data = df_a
)
model0

emm0 <- emmeans(model0, specs = ~ Cluster) #grand mean per cell cluster
cell_type_probs=emm0 %>%
  summary(infer = TRUE, type = 'response') %>%
  arrange(prob)

cell_type_probs %>% head
cell_type_probs %>% pull(prob) %>% sum #equal to 1

df_a %>% pull(Count) %>% sum -> N
N #365492, total cell count
sum(df_a$Count) #the same

cell_type_fractions=df_a %>% 
  group_by(Cluster) %>% 
  summarise(fraction = sum(Count) / N) %>% 
  arrange(fraction) %>% 
  as.data.frame 

head(cell_type_fractions)
head(cell_type_probs)

cell_type_probs = cell_type_fractions %>% left_join(cell_type_probs) 
(
  ggplot(aes(x = fraction, y = prob), data = cell_type_probs)
  + geom_abline(color = 'red')
  + geom_point()
  + geom_segment(aes(xend = fraction, y = asymp.LCL, yend = asymp.UCL))
  + scale_x_log10()
  + scale_y_log10()
  + theme_minimal()
  + labs(x = 'Fraction of all cells', y = 'Probability from GLM', subtitle = 'Cell type fraction vs probability')
)

## with GLM, use healty, location, batch as predictors
unique(df_a$Health) #"Healthy"      "Inflamed"     "Non-inflamed"

df=df_a %>% filter(Health %in% c('Healthy', 'Inflamed')) 
formula = cbind(Count, Other) ~ Cluster * Health + Cluster * Location + Cluster * batch
model1 <- glm(formula = formula, family = binomial(link='logit'), data = df)
summary(model1)

# emeans usage: https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html
# within each cell cluster, show result of constrasts within Health variable: Inflamed - Healthy
emm1 <- emmeans(model1, specs = revpairwise ~ Health | Cluster)
emm1
#Results are averaged over the levels of: Location, batch 
#Results are given on the log odds ratio (not the response) scale. 

c_results = emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() 
dim(c_results) #51 cell clusters

c_results %>% arrange(desc(odds.ratio))
(
  ggplot(aes(x = odds.ratio, y = -log10(p.value)), data = c_results)
  + geom_point()
  + geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.value)))
  + geom_text_repel(aes(label = Cluster), data = c_results %>% filter(p.value < 1e-50))
  + scale_x_log10()
  + theme_minimal()
  + labs(x = 'Inflamed / Healthy', title = 'Cell type proportion change')
)


# show three cell clusters as examples: highly enriched in Inflamed, non-enriched, hily depleted cell type
c_results[c_results$Cluster %in% c('Inflammatory Fibroblasts','NKs','MT-hi'),]

(
  ggplot(aes(x = Health, y = Count / Total), data = df %>% filter(Cluster == 'Inflammatory Fibroblasts'))
  + geom_jitter(height = 0)
  + scale_y_log10()
  + theme_minimal()
  + labs(subtitle = 'Inflammatory Fibroblasts')
)

(
  ggplot(aes(x = Health, y = Count / Total), data = df %>% filter(Cluster == 'NKs'))
  + geom_jitter(height = 0)
  + scale_y_log10()
  + theme_minimal()
  + labs(subtitle = 'NKs')
)

(
  ggplot(aes(x = Health, y = Count / Total), data = df %>% filter(Cluster == 'MT-hi'))
  + geom_jitter(height = 0)
  + scale_y_log10()
  + theme_minimal()
  + labs(subtitle = 'MT-hi')
)

##
emm2 <- emmeans(model1, specs = ~ Cluster)
mean_probs=emm2 %>%
  summary(type = 'response') %>%
  select(Cluster, prob) 
m_results = c_results %>% left_join(mean_probs) 

(
  ggplot(aes(x = prob, y = odds.ratio, color = p.value < 0.001), data = m_results)
  + geom_point()
  + geom_text_repel(aes(label = Cluster), color = 'black', data = m_results %>% filter(abs(log(odds.ratio)) > log(2)))
  + scale_x_log10()
  + scale_y_log10()
  + theme_minimal()
  + labs(y = 'Inflamed / Healthy (odds ratio)', title = 'Cell type proportion change', x = 'Average abundance (probability)')
)

m_results %>% arrange(desc(odds.ratio))

########################################################################
# 用 Dirichlet-multinomial regression在单细胞数据中检测细胞成分变化
# https://zhuanlan.zhihu.com/p/341941329

########################################################################
# Anova – Type I/II/III SS explained
# https://md.psych.bio.uni-goettingen.de/mv/unit/lm_cat/lm_cat_unbal_ss_explained.html


