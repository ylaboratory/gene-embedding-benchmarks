# *******************************
# plots for the gene embedding
# benchmarking paper
# *******************************
library('ggplot2')
library('data.table')
library('RColorBrewer')
library('tidytext')
library('ggnewscale')
library('viridis')

base_dir = 'results/tsvs/'

meta <- fread('data/embed_meta.csv')
names(meta) <- c('method', 'algorithm', 'input', 'dims')
meta[, algorithm := tolower(algorithm)]
meta[, algorithm := ifelse(algorithm == 'skip-gram/cbow', 'skip-gram', algorithm)]
meta[, algorithm := ifelse(algorithm == 'matrix factorization', 'matrix_fact', algorithm)]
meta[, input := ifelse(input == 'gene expression (single cell)', 'exp_sc', input)]
meta[, input := ifelse(input == 'gene expression (bulk)', 'exp_bulk', input)]
meta[, input := ifelse(input == 'biomedical literature', 'literature', input)]
meta[, input := ifelse(input == 'amino acid sequence', 'amino_acid', input)]
meta[, input := ifelse(input == 'mutation profile, biomedical literature, PPI', 'multi', input)]

cols <- data.table(input=sort(unique(meta$input)), color=c('#f6d374', '#6d1615', '#d2526b', '#4287af', '#1a3a4b', '#61d2a3'))

# ------------------
# Figure 2: single gene performance
# ------------------

# AUROC plots with jitter
s_go <- fread(paste0(base_dir, 'single_gene_go_AUROC.tsv'), sep="\t")
s_go <- melt(s_go)
names(s_go) <- c('embedding', 'go.term', 'auroc')
s_go <- merge(s_go, meta, by.x='embedding', by.y='method')
s_go <- merge(s_go, cols, by='input')

go_mapping = fread('data/gmt/hsa_EXP_ALL_BP_direct.gmt', 
                   select = 1:2, sep='\t', header=F, fill=T)
go_mapping[, V2 := gsub("\\s*\\(\\d+\\)$", "", V2)]
names(go_mapping) = c('id', 'name')

ggplot(s_go, aes(x=reorder(input, -auroc, FUN = median), y=auroc, color = color)) + 
  geom_jitter(width = 0.2, height=0, size = 1, alpha = 0.65) + 
  theme_classic() + xlab('') + ylab('AUROC') + 
  scale_color_identity() + 
  stat_summary(fun = median, geom = "crossbar", color = "black", width=0.5) + ylim(0,1)

ggsave('results/plots/GO_dot_input.pdf', width=3.25, height=3.2)


# supplementary plot
s_go = merge(s_go, go_mapping, by.x='go.term', by.y = 'id')

ggplot(s_go, aes(y = reorder(embedding, auroc, FUN = median), 
                 x = reorder(name, -auroc, FUN = median), fill = auroc)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "AUROC", option = "viridis", direction = 1) +
  theme_minimal() + coord_equal() + 
  theme(axis.text = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title = element_blank())

ggsave('results/plots/GO_heatmap.pdf', width=8.5, height=6)


# disease gene analysis
s_omim <- fread(paste0(base_dir, 'single_gene_omim_AUROC.tsv'), sep="\t")
s_omim <- melt(s_omim)
names(s_omim) <- c('embedding', 'term', 'auroc')
s_omim <- merge(s_omim, meta, by.x='embedding', by.y='method')
s_omim <- merge(s_omim, cols, by='input')

ggplot(s_omim, aes(x=reorder(input, -auroc, FUN = median), y=auroc, color = color)) + 
  geom_jitter(width = 0.2, height=0, size = 1, alpha = 0.6) + 
  theme_classic() + xlab('') + ylab('AUROC') + 
  scale_color_identity() + 
  stat_summary(fun = median, geom = "crossbar", color = "black", width=0.5) + ylim(0,1)

ggsave('results/plots/omim_dot_input.pdf', width=3.25, height=3.2)

# show all omim - supplementary
omim_mapping = fread('data/gmt/omim_entrez.gmt', select = 1:2, sep='\t', header=F, fill=2854)
omim_mapping[, V2 := gsub("\\s*\\(\\d+\\)$", "", V2)]
names(omim_mapping) = c('id', 'name')
omim_mapping = rbind(omim_mapping,  data.table(id=c('DOID:0081337'), name=c('congenital_myopathy')))

s_omim = merge(s_omim, omim_mapping, by.x='term', by.y='id')#, all.x=T)

ggplot(s_omim, aes(y = reorder(embedding, auroc, FUN = median), 
                 x = reorder(name, -auroc, FUN = median), fill = auroc)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "AUROC", option = "viridis", direction = 1) +
  theme_minimal() + coord_equal() + 
  theme(axis.text = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title = element_blank())

ggsave('results/plots/omim_heatmap.pdf', width=8.5, height=6)

# stats
agg_omim = s_omim[, .(median_auroc = median(auroc, na.rm = TRUE)), by = .(term, input)]
wilcox.test(agg_omim[input == 'amino_acid']$median_auroc, agg_omim[input == 'exp_sc']$median_auroc, paired=T)
wilcox.test(agg_omim[input == 'amino_acid']$median_auroc, agg_omim[input == 'exp_bulk']$median_auroc, paired=T)
wilcox.test(agg_omim[input == 'exp_sc']$median_auroc, agg_omim[input == 'PPI']$median_auroc, paired=T)
wilcox.test(agg_omim[input == 'exp_sc']$median_auroc, agg_omim[input == 'exp_bulk']$median_auroc, paired=T)

median(s_omim[input == 'exp_sc']$auroc)
median(s_omim[input == 'exp_bulk']$auroc)

mean(s_omim[input == 'exp_sc']$auroc)
mean(s_omim[input == 'exp_bulk']$auroc)

# data leakage evaluation
s_go_pre_2024 <- fread(paste0(base_dir, 'single_gene_go_pre2024_AUROC.tsv'), sep="\t")
s_go_pre_2024 <- melt(s_go_pre_2024)
names(s_go_pre_2024) <- c('method', 'term', 'pre.auroc')
s_go_post_2024 <- fread(paste0(base_dir, 'single_gene_go_post2024_AUROC.tsv'), sep="\t")
s_go_post_2024 <- melt(s_go_post_2024)
names(s_go_post_2024) <- c('method', 'term', 'post.auroc')
go_dl = merge(s_go_pre_2024, s_go_post_2024, by=c('method', 'term'))
go_dl = merge(go_dl, meta, by='method')
go_dl <- merge(go_dl, cols, by='input')

go_dl[, delta := post.auroc - pre.auroc]

# supplementary 
ggplot(go_dl, aes(x = reorder(method, -delta, FUN = median), y = term, fill = delta)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato", midpoint = 0, name = "delta AUROC") +
  theme_minimal() + coord_equal() + 
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())

ggsave('results/plots/single_go_pre_and_post_heatmap.pdf', width=8.5, height=6)


# summary heatmap
res_sum = s_go[, .(go.avg.auroc = mean(auroc), color = unique(color), input=unique(input)), by = embedding]
o_sum = s_omim[, .(omim.avg.auroc = mean(auroc)), by = embedding]
res_sum = merge(res_sum, o_sum, by='embedding')
g_sum = go_dl[, .(go.temporal = mean(post.auroc)), by = method]
res_sum = merge(res_sum, g_sum, by.x='embedding', by.y='method')


# get the full embed versions and merge into one table
s_go_full <- fread(paste0(base_dir, 'single_gene_go_AUROC_full.tsv'), sep="\t")
s_go_full <- melt(s_go_full)
names(s_go_full) <- c('embedding', 'go.term', 'auroc')
s_go_full = s_go_full[, .(go.avg.full.auroc = mean(auroc)), by = embedding]
res_sum = merge(res_sum, s_go_full, by='embedding')

s_omim_full <- fread(paste0(base_dir, 'single_gene_omim_AUROC_full.tsv'), sep="\t")
s_omim_full <- melt(s_omim_full)
names(s_omim_full) <- c('embedding', 'go.term', 'auroc')
s_omim_full = s_omim_full[, .(omim.avg.full.auroc = mean(auroc)), by = embedding]
res_sum = merge(res_sum, s_omim_full, by='embedding')


times <- fread(paste0(base_dir, 'single_gene_times_intersectvfull.tsv'), sep="\t")
names(times) <- c('embedding', 'go.full.time', 'go.time', 'omim.full.time', 'omim.time')
res_sum = merge(res_sum, times, by='embedding')
res_sum[, time.full := ((go.full.time + omim.full.time) / 2)]
res_sum[, time := ((go.time + omim.time) / 2)]
res_sum[, overall := ((go.avg.auroc + go.avg.full.auroc + omim.avg.auroc + omim.avg.full.auroc + go.temporal) / 5)]

res_sum = res_sum[order(-overall)]
res_sum[, rank := seq(1,38)]
res_sum = res_sum[order(-overall)]

# plot
dt <- melt(res_sum, id.vars = c("embedding", 'input', 'color'),
           measure.vars = c("go.avg.auroc","go.avg.full.auroc", 'go.temporal', "omim.avg.auroc","omim.avg.full.auroc", "time.full", 
                            'time', 'rank'),
           variable.name = "metric",
           value.name = "value")

dt[, embedding := factor(embedding, levels = rev(res_sum$embedding))]

ggplot() +
  geom_tile(data=dt[metric %in% c('rank')], aes(x = metric, y = embedding, fill = color)) +
  scale_fill_identity() +
  geom_text(data=dt[metric=='rank'], aes(x=metric,y=embedding, label=as.integer(value)), color='white', size=2) +
  new_scale_fill() +
  geom_tile(data=dt[metric %in% c("go.avg.auroc","go.avg.full.auroc", 'go.temporal', "omim.avg.auroc","omim.avg.full.auroc")],
            aes(x = metric, y = embedding, fill = value)) +
  #scale_fill_gradient(name = "AUROC", low = "white", high = "steelblue", limits = c(0,1)) +
  scale_fill_viridis_c(name = "AUROC", begin = 0.2, end = 1.0, option = "magma", 
                       limits = c(0.5,1.0), oob = scales::squish) +
  new_scale_fill() +
  geom_tile(data=dt[metric %in% c('time.full')], aes(x = metric, y = embedding, fill = value)) +
  scale_fill_gradient(name = "Time (s)", low = "white", high = "#2166ac") +
  theme_classic() + scale_x_discrete(position = "top") + coord_equal() + 
  theme(axis.title.x.top = element_text(vjust = -1), axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid = element_blank(), axis.ticks = element_blank(), axis.line= element_blank()) +
  labs(x = "", y = "")

ggsave('results/plots/single_summary_heatmap.pdf', width=5, height=6)

# ------------------
# Figure 3: paired gene performance
# ------------------

sl_cat <- fread(paste0(base_dir, 'paired_gene_sl_concat_outerloop_avgs.tsv'), sep="\t")
sl_sum <- fread(paste0(base_dir, 'paired_gene_sl_sum_outerloop_avgs.tsv'), sep="\t")
sl_prod <- fread(paste0(base_dir, 'paired_gene_sl_product_outerloop_avgs.tsv'), sep="\t")
sl_all <- merge(sl_cat[,c('embedding', 'outer_AUC')], sl_sum[,c('embedding', 'outer_AUC')], by='embedding')
names(sl_all) <- c('method', 'cat', 'sum')
sl_all <- merge(sl_all, sl_prod[,c('embedding', 'outer_AUC')], by.x='method', by.y='embedding')
names(sl_all) <- c('method', 'cat', 'sum', 'prod')
sl_all <- merge(sl_all, meta, by='method')
sl_all <- merge(sl_all, cols, by='input')

sl_all <- melt(sl_all, id.vars = c("method", 'input','algorithm', 'color'),
               measure.vars = c('cat', 'sum', 'prod'),
               variable.name = "operation",
               value.name = "auroc")

ggplot(sl_all, aes(x=reorder(input, -auroc, FUN = median), y=auroc, color = color, shape=operation)) + 
  geom_jitter(width = 0.2, height=0, size = 1.5, alpha = 1) + ylim(0,1) + 
  theme_classic() + xlab('') + ylab('AUROC') + 
  scale_color_identity() 
ggsave('results/plots/SL_dot_input.pdf', width=4.75, height=4)

ng_cat <- fread(paste0(base_dir, 'paired_gene_ng_concat_outerloop_avgs.tsv'), sep="\t")
ng_sum <- fread(paste0(base_dir, 'paired_gene_ng_sum_outerloop_avgs.tsv'), sep="\t")
ng_prod <- fread(paste0(base_dir, 'paired_gene_ng_product_outerloop_avgs.tsv'), sep="\t")
ng_all <- merge(ng_cat[,c('embedding', 'outer_AUC')], ng_sum[,c('embedding', 'outer_AUC')], by='embedding')
names(ng_all) <- c('method', 'cat', 'sum')
ng_all <- merge(ng_all, ng_prod[,c('embedding', 'outer_AUC')], by.x='method', by.y='embedding')
names(ng_all) <- c('method', 'cat', 'sum', 'prod')
ng_all <- merge(ng_all, meta, by='method')
ng_all <- merge(ng_all, cols, by='input')

ng_all <- melt(ng_all, id.vars = c("method", 'input','algorithm', 'color'),
               measure.vars = c('cat', 'sum', 'prod'),
               variable.name = "operation",
               value.name = "auroc")

ggplot(ng_all, aes(x=reorder(input, -auroc, FUN = median), y=auroc, color = color, shape=operation)) + 
  geom_jitter(width = 0.2, height=0, size = 1.5, alpha = 1) + ylim(0,1) + 
  theme_classic() + xlab('') + ylab('AUROC') + 
  scale_color_identity()  
ggsave('results/plots/NG_dot_input.pdf',  width=4.75, height=4)


tf_cat <- fread(paste0(base_dir, 'paired_gene_tf_concat_outerloop_avgs.tsv'), sep="\t")
tf_sum <- fread(paste0(base_dir, 'paired_gene_tf_sum_outerloop_avgs.tsv'), sep="\t")
tf_prod <- fread(paste0(base_dir, 'paired_gene_tf_product_outerloop_avgs.tsv'), sep="\t")
tf_all <- merge(tf_cat[,c('embedding', 'outer_AUC')], tf_sum[,c('embedding', 'outer_AUC')], by='embedding')
names(tf_all) <- c('method', 'cat', 'sum')
tf_all <- merge(tf_all, tf_prod[,c('embedding', 'outer_AUC')], by.x='method', by.y='embedding')
names(tf_all) <- c('method', 'cat', 'sum', 'prod')
tf_all <- merge(tf_all, meta, by='method')
tf_all <- merge(tf_all, cols, by='input')

tf_all <- melt(tf_all, id.vars = c("method", 'input','algorithm', 'color'),
               measure.vars = c('cat', 'sum', 'prod'),
               variable.name = "operation",
               value.name = "auroc")

ggplot(tf_all, aes(x=reorder(input, -auroc, FUN = median), y=auroc, color = color, shape=operation)) + 
  geom_jitter(width = 0.2, height=0, size = 1.5, alpha = 1) + ylim(0,1) + 
  theme_classic() + xlab('') + ylab('AUROC') + 
  scale_color_identity()
ggsave('results/plots/TF_dot_input.pdf',  width=4.75, height=4)

# data leakage evaluation
pombe = fread(paste0(base_dir, 'paired_gene_pombe_ng_sum_outerloop_avgs.tsv'), sep='\t')
pombe_p = fread(paste0(base_dir, 'paired_gene_pombe_ng_product_outerloop_avgs.tsv'), sep='\t')
pombe_c = fread(paste0(base_dir, 'paired_gene_pombe_ng_concat_outerloop_avgs.tsv'), sep='\t')

pombe = pombe[,c('embedding', 'outer_AUC', 'outer_time')]
pombe[, type := 'sum']
pombe_p[, type := 'prod']
pombe_c[, type := 'cat']
pombe = rbind(pombe, pombe_p[,c('embedding', 'outer_AUC', 'outer_time', 'type')])
pombe = rbind(pombe, pombe_c[,c('embedding', 'outer_AUC', 'outer_time', 'type')])
names(pombe) = c('method', 'AUROC', 'time', 'type')

pombe = merge(pombe, meta, by='method')
pombe = merge(pombe, cols, by='input')

# add the full embedding results
ng_cat_f <- fread(paste0(base_dir, 'paired_gene_ng_concat_outerloop_avgs_full.tsv'), sep="\t")
ng_sum_f <- fread(paste0(base_dir, 'paired_gene_ng_sum_outerloop_avgs_full.tsv'), sep="\t")
ng_prod_f <- fread(paste0(base_dir, 'paired_gene_ng_product_outerloop_avgs_full.tsv'), sep="\t")
ng_all_f <- merge(ng_cat_f[,c('embedding', 'outer_AUC')], ng_sum_f[,c('embedding', 'outer_AUC')], by='embedding')
names(ng_all_f) <- c('method', 'cat', 'sum')
ng_all_f <- merge(ng_all_f, ng_prod_f[,c('embedding', 'outer_AUC')], by.x='method', by.y='embedding')
names(ng_all_f) <- c('method', 'cat', 'sum', 'prod')
ng_all_f <- merge(ng_all_f, meta, by='method')
ng_all_f <- merge(ng_all_f, cols, by='input')


sl_cat_f <- fread(paste0(base_dir, 'paired_gene_sl_concat_outerloop_avgs_full.tsv'), sep="\t")
sl_sum_f <- fread(paste0(base_dir, 'paired_gene_sl_sum_outerloop_avgs_full.tsv'), sep="\t")
sl_prod_f <- fread(paste0(base_dir, 'paired_gene_sl_product_outerloop_avgs_full.tsv'), sep="\t")
sl_all_f <- merge(sl_cat_f[,c('embedding', 'outer_AUC')], sl_sum_f[,c('embedding', 'outer_AUC')], by='embedding')
names(sl_all_f) <- c('method', 'cat', 'sum')
sl_all_f <- merge(sl_all_f, sl_prod_f[,c('embedding', 'outer_AUC')], by.x='method', by.y='embedding')
names(sl_all_f) <- c('method', 'cat', 'sum', 'prod')
sl_all_f <- merge(sl_all_f, meta, by='method')
sl_all_f <- merge(sl_all_f, cols, by='input')


tf_cat_f <- fread(paste0(base_dir, 'paired_gene_tf_concat_outerloop_avgs_full.tsv'), sep="\t")
tf_sum_f <- fread(paste0(base_dir, 'paired_gene_tf_sum_outerloop_avgs_full.tsv'), sep="\t")
tf_prod_f <- fread(paste0(base_dir, 'paired_gene_tf_product_outerloop_avgs_full.tsv'), sep="\t")
tf_all_f <- merge(tf_cat_f[,c('embedding', 'outer_AUC')], tf_sum_f[,c('embedding', 'outer_AUC')], by='embedding')
names(tf_all_f) <- c('method', 'cat', 'sum')
tf_all_f <- merge(tf_all_f, tf_prod_f[,c('embedding', 'outer_AUC')], by.x='method', by.y='embedding')
names(tf_all_f) <- c('method', 'cat', 'sum', 'prod')
tf_all_f <- merge(tf_all_f, meta, by='method')
tf_all_f <- merge(tf_all_f, cols, by='input')

ng_all_f <- melt(ng_all_f, id.vars = c("method", 'input','algorithm', 'color'),
                 measure.vars = c('cat', 'sum', 'prod'),
                 variable.name = "operation",
                 value.name = "auroc")

sl_all_f <- melt(sl_all_f, id.vars = c("method", 'input','algorithm', 'color'),
                 measure.vars = c('cat', 'sum', 'prod'),
                 variable.name = "operation",
                 value.name = "auroc")

tf_all_f <- melt(tf_all_f, id.vars = c("method", 'input','algorithm', 'color'),
                 measure.vars = c('cat', 'sum', 'prod'),
                 variable.name = "operation",
                 value.name = "auroc")

sl_all_f[, type := 'SL']
ng_all_f[, type := 'NG']
tf_all_f[, type := 'TF']

res_f = rbind(sl_all_f, ng_all_f)
res_f = rbind(res_f, tf_all_f)
res_f = dcast(res_f, method + operation ~ type, value.var = 'auroc')
names(res_f) <- c('method', 'operation', 'NG.f', 'SL.f', 'TF.f')

# summary heatmap
sl_all[, type := 'SL']
ng_all[, type := 'NG']
tf_all[, type := 'TF']
res = rbind(sl_all, ng_all)
res = rbind(res, tf_all)

res.df = dcast(res, method + operation ~ type, value.var = 'auroc')
res.df = merge(res.df, res_f, by=c('method', 'operation'))
res.df = merge(res.df, pombe[,c('method', 'AUROC', 'type')], by.x=c('method', 'operation'), by.y=c('method', 'type'))
names(res.df) <- c('method', 'operation', 'NG', 'SL', 'TF', 'NG.f', 'SL.f', 'TF.f', 'pombe')

res.df[, avg := (NG + SL + TF + NG.f + SL.f + TF.f + pombe) / 7]

res.df = res.df[order(-avg)]
res.df[, rank := seq(1,nrow(res.df))]
res.df = res.df[order(-avg)]
top.res = res.df[ , .SD[which.max(avg)], by = method]

# get the corresponding times
times = sl_sum[,c('embedding', 'outer_time')]
times = merge(times, sl_cat[, c('embedding', 'outer_time')], by='embedding')
names(times) = c('embedding', 'sl.sum', 'sl.cat') 
times = merge(times, sl_prod[, c('embedding', 'outer_time')], by='embedding')
names(times) = c('embedding', 'sl.sum', 'sl.cat', 'sl.prod')
times = merge(times, ng_sum[, c('embedding', 'outer_time')], by='embedding')
times = merge(times, ng_cat[, c('embedding', 'outer_time')], by='embedding')
times = merge(times, ng_prod[, c('embedding', 'outer_time')], by='embedding')
names(times) = c('embedding', 'sl.sum', 'sl.cat', 'sl.prod', 'ng.sum', 'ng.cat', 'ng.prod')
times = merge(times, tf_sum[, c('embedding', 'outer_time')], by='embedding')
times = merge(times, tf_cat[, c('embedding', 'outer_time')], by='embedding')
times = merge(times, tf_prod[, c('embedding', 'outer_time')], by='embedding')
names(times) = c('method', 'SL.sum', 'SL.cat', 'SL.prod', 'NG.sum', 'NG.cat', 'NG.prod',
                 'TF.sum', 'TF.cat', 'TF.prod')
times = melt(times)
names(times) = c('method', 'var', 'time')
times[ , c("type","operation") := tstrsplit(as.character(var), "\\.", fixed=FALSE) ]
times[, var := NULL]
res = merge(res, times, by=c('method', 'operation', 'type'))
top.res = merge(top.res, times[type == 'NG'], by=c('method', 'operation'))
top.res[, type := NULL]
names(top.res) = c('method', 'operation', 'NG', 'SL', 'TF', 'NG.f', 'SL.f', 'TF.f', 'pombe', 'avg', 'rank', 'ng.time')
top.res = merge(top.res, times[type == 'SL'], by=c('method', 'operation'))
top.res[, type := NULL]
names(top.res) = c('method', 'operation', 'NG', 'SL', 'TF', 'NG.f', 'SL.f', 'TF.f', 'pombe', 'avg', 'rank', 'ng.time', 'sl.time')
top.res = merge(top.res, times[type == 'TF'], by=c('method', 'operation'))
top.res[, type := NULL]
names(top.res) = c('method', 'operation', 'NG', 'SL', 'TF', 'NG.f', 'SL.f', 'TF.f', 'pombe', 'avg', 'rank', 'ng.time', 'sl.time', 'tf.time')
top.res = merge(top.res, meta, by='method')
top.res = merge(top.res, cols, by='input')
top.res[,op := rank]
top.res[, ng.time := log10(ng.time)]
top.res[, sl.time := log10(sl.time)]
top.res[, tf.time := log10(tf.time)]
top.res = top.res[order(rank)]

# plot
dt <- melt(top.res, id.vars = c("method", 'operation', 'input', 'color'),
           measure.vars = c("NG", 'NG.f',"SL", 'SL.f', "TF", 'TF.f', 'pombe', "ng.time", "sl.time", 'tf.time', 'rank', 'op'),
           variable.name = "metric",
           value.name = "value")

dt[, method := factor(method, levels = rev(top.res$method))]

ggplot() +
  geom_tile(data=dt[metric %in% c('rank')], aes(x = metric, y = method, fill = color)) +
  scale_fill_identity() +
  geom_text(data=dt[metric=='rank'], aes(x=metric,y=method, label=as.integer(value)), color='white', size=2) +
  new_scale_fill() +
  geom_tile(data=dt[metric %in% c('op')], aes(x = metric, y = method, fill=operation)) +
  scale_fill_brewer(palette = 'Dark2') +
  geom_text(data = dt[metric == 'op'],
            aes(x = metric, y = method,
                label = ifelse(operation == 'cat', '&',
                               ifelse(operation == 'sum', '+',
                                      ifelse(operation == 'prod', '*', '')))),
            color = 'white', size = 3) +
  new_scale_fill() +
  geom_tile(data=dt[metric %in% c("SL", 'SL.f', "NG", 'NG.f', "TF", 'TF.f', 'pombe')],
            aes(x = metric, y = method, fill = value)) +
 scale_fill_viridis_c(name = "AUROC",begin = 0.2,
    end = 1.0, option = "magma", limits = c(0.45, 1.0), oob = scales::squish) +
  new_scale_fill() +
  geom_tile(data=dt[metric %in% c('sl.time', 'ng.time', 'tf.time')], aes(x = metric, y = method, fill = value)) +
  scale_fill_gradient(name = "Time log(s)", low = "white", high = "#2166ac") +
  theme_classic() + scale_x_discrete(position = "top") + coord_equal() + 
  theme(axis.title.x.top = element_text(vjust = -1), axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid = element_blank(), axis.ticks = element_blank(), axis.line= element_blank()) +
  labs(x = "", y = "")

ggsave('results/plots/pair_summary_heatmap.pdf', width=5.5, height=6)

# supplementary
# compare the algorithm with operation and see if there are any
# correlations between overall performance
ggplot(res, aes(x=algorithm, y=auroc, fill = operation)) + geom_boxplot() + theme_classic() +
  xlab("") + scale_fill_brewer(palette = 'Dark2')

ggsave('results/plots/pair_algorithm_operation.pdf', width=5.5, height=4)

# stats
am = aov(auroc~type*operation*algorithm,data=res)
summary(am)

library(emmeans)
emmeans(am, pairwise ~ type, adjust="tukey")
emmeans(am, pairwise ~ operation, adjust="tukey")
emmeans(am, pairwise ~ algorithm, adjust="tukey")
emmeans(am, pairwise ~ operation | algorithm, adjust="tukey")

# -------------------
# Figure 4
# -------------------
kegg_go <- fread(paste0(base_dir, 'gene_set_andes_kegg_go_rank_no_overlap.tsv'), sep="\t")
kegg_go = melt(kegg_go)
names(kegg_go) <- c('method', 'sets', 'rank')
kegg_go[, type := 'no_overlap']

kg <- fread(paste0(base_dir, 'gene_set_andes_kegg_go_rank_overlap.tsv'), sep="\t")
kg = melt(kg)
names(kg) <- c('method', 'sets', 'rank')
kg[, type := 'overlap']
kegg_go = rbind(kegg_go, kg)
rm(kg)
kg <- fread(paste0(base_dir, 'gene_set_andes_kegg_go_rank_overlap_full.tsv'), sep="\t")
kg = melt(kg)
names(kg) <- c('method', 'sets', 'rank')
kg[, type := 'foverlap']
kegg_go = rbind(kegg_go, kg)
rm(kg)
kg <- fread(paste0(base_dir, 'gene_set_andes_kegg_go_rank_no_overlap_full.tsv'), sep="\t")
kg = melt(kg)
names(kg) <- c('method', 'sets', 'rank')
kg[, type := 'fno_overlap']
kegg_go = rbind(kegg_go, kg)
rm(kg)

kegg_go[ , c("s1","s2") := tstrsplit(as.character(sets), "\\_", fixed=FALSE)]

dis_tis <- fread(paste0(base_dir, 'gene_set_andes_disease_tissue_btoAUROC_no_overlap.tsv'), sep="\t")
dis_tis = melt(dis_tis)
names(dis_tis) <- c('method', 'sets', 'AUROC')
dis_tis[, type := 'no_overlap']

dts <- fread(paste0(base_dir, 'gene_set_andes_disease_tissue_btoAUROC_overlap.tsv'), sep="\t")
dts = melt(dts)
names(dts) <- c('method', 'sets', 'AUROC')
dts[, type := 'overlap']

dis_tis = rbind(dis_tis, dts)
rm(dts)

dts <- fread(paste0(base_dir, 'gene_set_andes_disease_tissue_btoAUROC_overlap_full.tsv'), sep="\t")
dts = melt(dts)
names(dts) <- c('method', 'sets', 'AUROC')
dts[, type := 'foverlap']

dis_tis = rbind(dis_tis, dts)
rm(dts)

dts <- fread(paste0(base_dir, 'gene_set_andes_disease_tissue_btoAUROC_no_overlap_full.tsv'), sep="\t")
dts = melt(dts)
names(dts) <- c('method', 'sets', 'AUROC')
dts[, type := 'fno_overlap']

dis_tis = rbind(dis_tis, dts)
rm(dts)

# calc normalized rank
kegg_go[, nrank := rank / length(unique(s1)), by = type]
kegg_go[, mtype := 'kegg.go']

kegg_go = merge(kegg_go, meta, by='method')
kegg_go = merge(kegg_go, cols, by='input')

dis_tis = merge(dis_tis, meta, by='method')
dis_tis = merge(dis_tis, cols, by='input')

ggplot(kegg_go[type=='no_overlap'], aes(x=reorder(input, nrank, FUN=median), y=1-nrank, color=color)) +
  geom_jitter(width = 0.2, height = 0.01, size = 1, alpha=0.7) +
  ylim(0,1) + theme_classic() +
  scale_color_identity() +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.5) + 
  xlab('') + ylab('normalized rank')

ggsave('results/plots/geneset_kegg-go_dot_input.pdf', width=3.25, height=3.2)

ggplot(dis_tis[type=='no_overlap'], aes(x=reorder(input, -AUROC, FUN=median), y=AUROC, color=color)) +
  geom_jitter(width = 0.2, height = 0.01, size = 1, alpha=0.7) +
  ylim(0,1) + theme_classic() +
  scale_color_identity() +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.5) + 
  xlab('') + ylab('AUROC') 

ggsave('results/plots/geneset_dis-tis_dot_input.pdf', width=3.25, height=3.2)

# data leakage
go_go <- fread(paste0(base_dir, 'gene_set_andes_go_go_rank_no_overlap.tsv'))
go_go = melt(go_go)
names(go_go) <- c('method', 'sets', 'rank')
go_go[, type := 'intersect']

gg2 <- fread(paste0(base_dir, 'gene_set_andes_go_go_rank_no_overlap_full.tsv'))
gg2 = melt(gg2)
names(gg2) <- c('method', 'sets', 'rank')
gg2[, type := 'full']

go_go = rbind(go_go, gg2)
rm(gg2)

go_go[ , c("s1","s2") := tstrsplit(as.character(sets), "\\_", fixed=FALSE)]
go_go[, nrank := rank / length(unique(s1)), by = type]

go_go = merge(go_go, meta, by='method')
go_go = merge(go_go, cols, by='input')
go_go = merge(go_go, go_mapping, by.y = 'id', by.x = 's1')

go_go[, go_avg_rank := mean(rank), by = c('s1', 'type')] 
go_go[, go_med_rank := median(rank), by = c('s1', 'type')] 

go_tmp = go_go[order(rank), .SD[1], by = c('s1', 'type')]

# supplementary 
ggplot(go_go[type=='intersect'], aes(x = reorder(method, rank, FUN = median), 
                 y = reorder(name, rank, FUN = median), fill = rank)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "rank", option = "viridis", direction = -1) +
  theme_minimal() + coord_equal() + 
  theme(axis.text = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title = element_blank())

ggsave('results/plots/geneset_generalization_heatmap.pdf', width=8.5, height=6)


# summmary plot
tmpk = kegg_go[type == 'no_overlap']
tmpd = dis_tis[type == 'no_overlap']
res_gs = tmpk[, .(kegg_go = mean(nrank), color = unique(color), input=unique(input)), by = method]
res_gs = merge(res_gs, tmpd[, .(dis_tis = mean(AUROC)), by=method], by='method')

tmpk = kegg_go[type == 'fno_overlap']
tmpd = dis_tis[type == 'fno_overlap']
res_gs = merge(res_gs, tmpd[, .(dis_tis_fno = mean(AUROC)), by=method], by='method')
res_gs = merge(res_gs, tmpk[, .(kegg_go_fno = mean(nrank)), by=method], by='method')

rm(tmpk, tmpd)

res_gs = merge(res_gs, go_go[type == 'intersect', .(temporal = mean(nrank)), by = method], by='method')

res_gs[, avg := (kegg_go + kegg_go_fno + temporal)/ 3]
res_gs[, avg2 := (dis_tis + dis_tis_fno)/ 2]
res_gs = res_gs[order(avg)]
res_gs[, rank.kg := seq(1,nrow(res_gs))]


res_gs = res_gs[order(-avg2)]
res_gs[, rank.dt := seq(1,nrow(res_gs))]

res_gs[, avg.rank := (rank.kg + rank.dt) / 2]
res_gs = res_gs[order(avg.rank)]
res_gs[, rank := seq(1,nrow(res_gs))]

# plot
dt <- melt(res_gs, id.vars = c("method", 'input', 'color'),
           measure.vars = c("kegg_go", 'kegg_go_fno', 'temporal', "dis_tis", 'dis_tis_fno', 'rank'),
           variable.name = "metric",
           value.name = "value")

dt[, method := factor(method, levels = rev(res_gs$method))]

ggplot() +
  geom_tile(data=dt[metric %in% c('rank')], aes(x = metric, y = method, fill = color)) +
  scale_fill_identity() +
  geom_text(data=dt[metric=='rank'], aes(x=metric, y=method, label=as.integer(value)), color='white', size=2) +
  new_scale_fill() +
  geom_tile(data=dt[metric %in% c("kegg_go", 'kegg_go_fno', 'temporal')],
              aes(x = metric, y = method, fill = 1-value)) +
  scale_fill_viridis_c(name = "normalized rank", option = "viridis", direction = 1) +
  new_scale_fill() + 
  geom_tile(data=dt[metric %in% c("dis_tis",'dis_tis_fno')], aes(x = metric, y = method, fill = value)) +
  scale_fill_viridis_c(name = "AUROC",
                       begin = 0.2, end = 1.0, option = "magma",
                       limits = c(0.49,1.0), oob = scales::squish) +
  theme_classic() + scale_x_discrete(position = "top") + coord_equal() + 
  theme(axis.title.x.top = element_text(vjust = -1), axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid = element_blank(), axis.ticks = element_blank(), axis.line= element_blank()) +
  labs(x = "", y = "")

ggsave('results/plots/geneset_summary_heatmap.pdf', width=5, height=6)

# **************************
# Figure 5
# **************************
# cca 
cca <- fread(paste0(base_dir, 'cca_pca10.tsv'))
mat <- cca[,-1]
row_clust <- hclust(dist(mat))
col_clust <- hclust(dist(t(mat)))
mat_ordered <- as.matrix(mat)
row.names(mat_ordered) <- cca$V1
mat_ordered <- mat_ordered[row_clust$order, col_clust$order]
cca <- melt(mat_ordered)
cca <- as.data.table(cca)
names(cca) <- c('method1', 'method2', 'cca')

# plot CCA matrix
library(patchwork)

baseheat <- ggplot(cca, aes(x = method2, y = method1, fill = cca)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = "CCA") +
  theme_minimal() + ylab("") + xlab("") + coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

cl = merge(cca, merge(meta, cols, by='input'), by.x='method1', by.y='method')
cl$method1 = factor(cl$method1, levels=levels(cca$method1))
col_labs = ggplot(cl,
                  aes(method1, y=1, fill=color)) + geom_tile() +
  scale_fill_identity() + coord_equal() + ylab("") + xlab("") +
  theme_void() 

col_labs2 = ggplot(cl,
                  aes(method1, y=1, fill=algorithm)) + geom_tile() +
  coord_equal() + ylab("") + xlab("") + scale_fill_brewer(palette = 'Dark2') + 
  theme_void() 

(col_labs / col_labs2) / baseheat

ggsave('results/plots/cca_heatmap.pdf', width=8, height=8)

# jaccard gene set coverage stats
jac <- fread(paste0(base_dir, 'jaccard_genes.tsv'))
jmat <- as.matrix(jac[,-1])
mean(jmat[upper.tri(jmat)], na.rm = TRUE)

mjac <- melt(jac)
names(mjac) <- c('method1', 'method2', 'jac')
mjac = mjac[method1 != method2]

# supplementary jaccard plot
mat2 <- jac[,-1]
row_clust <- hclust(dist(mat2))
col_clust <- hclust(dist(t(mat2)))
mat2_ordered <- as.matrix(mat2)
row.names(mat2_ordered) <- jac$V1
mat2_ordered <- mat2_ordered[row_clust$order, col_clust$order]
mjac2 <- melt(mat2_ordered)
mjac2 <- as.data.table(mjac2)
names(mjac2) <- c('method1', 'method2', 'jac')

basejac <- ggplot(mjac2, aes(x = method2, y = method1, fill = jac)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "C", direction=-1, limits = c(0, 1), name = "Jac") +
  theme_minimal() + ylab("") + xlab("") + coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


sizes = meta[,c('method', 'dims')]
sizes[, n_genes := as.integer(sub("^\\((\\d+),.*", "\\1", dims))]

sizes$method = factor(sizes$method, levels = levels(mjac2$method1))

upper_counts = ggplot(sizes, aes(x = method, y=1, fill = n_genes)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "G", direction=-1, limits = c(15000, 40000), name = "N") +
  theme_void() + coord_equal() + ylab('') + xlab('') 

upper_counts / basejac

ggsave('results/plots/jac_heatmap.pdf', width=8, height=8)


# get pairwise correlations of performance
tmp = dcast(s_go, go.term~embedding, value.var = 'auroc')
tmp = cor(tmp[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 's.go.cor')

perf_cor = merge(cca, tmp, by=c('method1', 'method2'))

tmp = dcast(s_omim, term~embedding, value.var = 'auroc')
tmp = cor(tmp[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 's.omim.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

tmp = dcast(go_dl, term~method, value.var = 'post.auroc')
tmp = cor(tmp[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 's.temporal.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

# paired
sl_cat <- fread(paste0(base_dir, 'paired_gene_sl_concat_outerloop_AUC.tsv'), sep="\t")
sl_sum <- fread(paste0(base_dir, 'paired_gene_sl_sum_outerloop_AUC.tsv'), sep="\t")
sl_prod <- fread(paste0(base_dir, 'paired_gene_sl_product_outerloop_AUC.tsv'), sep="\t")

sl_cat[, outer_AUC_average := NULL]
sl_sum[, outer_AUC_average := NULL]
sl_prod[, outer_AUC_average := NULL]

sl_cat = melt(sl_cat)
names(sl_cat) <- c('embedding', 'fold', 'cat')
sl_cat = dcast(sl_cat, fold~embedding)

sl_sum = melt(sl_sum)
names(sl_sum) <- c('embedding', 'fold', 'sum')
sl_sum = dcast(sl_sum, fold~embedding)

sl_prod = melt(sl_prod)
names(sl_prod) <- c('embedding', 'fold', 'prod')
sl_prod = dcast(sl_prod, fold~embedding)

ng_cat <- fread(paste0(base_dir, 'paired_gene_ng_concat_outerloop_AUC.tsv'), sep="\t")
ng_sum <- fread(paste0(base_dir, 'paired_gene_ng_sum_outerloop_AUC.tsv'), sep="\t")
ng_prod <- fread(paste0(base_dir, 'paired_gene_ng_product_outerloop_AUC.tsv'), sep="\t")

ng_cat[, outer_AUC_average := NULL]
ng_sum[, outer_AUC_average := NULL]
ng_prod[, outer_AUC_average := NULL]

ng_cat = melt(ng_cat)
names(ng_cat) <- c('embedding', 'fold', 'cat')
ng_cat = dcast(ng_cat, fold~embedding)

ng_sum = melt(ng_sum)
names(ng_sum) <- c('embedding', 'fold', 'sum')
ng_sum = dcast(ng_sum, fold~embedding)

ng_prod = melt(ng_prod)
names(ng_prod) <- c('embedding', 'fold', 'prod')
ng_prod = dcast(ng_prod, fold~embedding)

tf_cat <- fread(paste0(base_dir, 'paired_gene_tf_concat_outerloop_AUC.tsv'), sep="\t")
tf_sum <- fread(paste0(base_dir, 'paired_gene_tf_sum_outerloop_AUC.tsv'), sep="\t")
tf_prod <- fread(paste0(base_dir, 'paired_gene_tf_product_outerloop_AUC.tsv'), sep="\t")

tf_cat[, outer_AUC_average := NULL]
tf_sum[, outer_AUC_average := NULL]
tf_prod[, outer_AUC_average := NULL]

tf_cat = melt(tf_cat)
names(tf_cat) <- c('embedding', 'fold', 'cat')
tf_cat = dcast(tf_cat, fold~embedding)

tf_sum = melt(tf_sum)
names(tf_sum) <- c('embedding', 'fold', 'sum')
tf_sum = dcast(tf_sum, fold~embedding)

tf_prod = melt(tf_prod)
names(tf_prod) <- c('embedding', 'fold', 'prod')
tf_prod = dcast(tf_prod, fold~embedding)

pair_cat = rbind(sl_cat, ng_cat)
pair_cat = rbind(pair_cat, tf_cat)

pair_sum = rbind(sl_sum, ng_sum)
pair_sum = rbind(pair_sum, tf_sum)

pair_prod = rbind(sl_prod, ng_prod)
pair_prod = rbind(pair_prod, tf_prod)


p_cat <- fread(paste0(base_dir, 'paired_gene_pombe_ng_concat_outerloop_AUC.tsv'), sep="\t")
p_sum <- fread(paste0(base_dir, 'paired_gene_pombe_ng_sum_outerloop_AUC.tsv'), sep="\t")
p_prod <- fread(paste0(base_dir, 'paired_gene_pombe_ng_product_outerloop_AUC.tsv'), sep="\t")

p_cat[, outer_AUC_average := NULL]
p_sum[, outer_AUC_average := NULL]
p_prod[, outer_AUC_average := NULL]

p_cat = melt(p_cat)
names(p_cat) <- c('embedding', 'fold', 'cat')
p_cat = dcast(p_cat, fold~embedding)

p_sum = melt(p_sum)
names(p_sum) <- c('embedding', 'fold', 'sum')
p_sum = dcast(p_sum, fold~embedding)

p_prod = melt(p_prod)
names(p_prod) <- c('embedding', 'fold', 'prod')
p_prod = dcast(p_prod, fold~embedding)

pair_cat = rbind(pair_cat, p_cat)
pair_sum = rbind(pair_sum, p_sum)
pair_prod = rbind(pair_prod, p_prod)

tmp = cor(pair_cat[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 'p.cat.cor')
 
perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

tmp = cor(pair_sum[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 'p.sum.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

tmp = cor(pair_prod[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 'p.prod.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))


# sets
tmp = dcast(dis_tis[type == 'no_overlap'], sets~method, value.var='AUROC')
tmp = cor(tmp[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 'dis.tis.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

tmp = dcast(kegg_go[type == 'no_overlap'], sets~method, value.var='rank')
tmp = cor(tmp[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 'kegg.go.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

tmp = dcast(go_go[type == 'intersect'], sets~method, value.var='rank')
tmp = cor(tmp[, -1, with = FALSE], method='spearman')
tmp = melt(tmp)
names(tmp) <- c('method1', 'method2', 'go.go.temporal.cor')

perf_cor = merge(perf_cor, tmp, by=c('method1', 'method2'))

# keep only one direction of the pair since
# directionality isn't important
perf_cor = perf_cor[method1 != method2]
perf_cor[, method1 := as.character(method1)]
perf_cor[, method2 := as.character(method2)]
perf_cor[, pair_id := paste0(pmin(method1, method2),'-', pmax(method1, method2))]
perf_cor_uni = perf_cor[!duplicated(pair_id)]
perf_cor_uni = merge(perf_cor_uni, meta[,c('method','input')], by.x='method1', by.y = 'method')
perf_cor_uni = merge(perf_cor_uni, meta[,c('method','input')], by.x='method2', by.y = 'method')
perf_cor_uni[, pair_type := paste0(pmin(input.x, input.y),'-', pmax(input.x, input.y))]
perf_cor_uni[, pair_type_bin := ifelse(input.x == input.y, 'same', 'cross')]

# averages for reducing the number of values
cor_res_pt = character()
cor_res_task = character()
cor_res_value = numeric()
tasks = c('s.go.cor', 's.omim.cor', 's.temporal.cor', 'p.cat.cor', 'p.sum.cor', 'p.prod.cor', 'dis.tis.cor', 'kegg.go.cor', 'go.go.temporal.cor', 'cca')

for(pt in unique(perf_cor_uni$pair_type))
{
  for(task in tasks)
  {
    tmp = perf_cor_uni[pair_type == pt]
    cor_res_pt = append(cor_res_pt, pt)
    cor_res_task = append(cor_res_task, task)
    cor_res_value = append(cor_res_value, median(as.numeric(unlist(tmp[,task,with=F]))))
  }
}
for(pt in unique(perf_cor_uni$pair_type_bin))
{
  for(task in tasks)
  {
    tmp = perf_cor_uni[pair_type_bin == pt]
    cor_res_pt = append(cor_res_pt, pt)
    cor_res_task = append(cor_res_task, task)
    cor_res_value = append(cor_res_value, median(as.numeric(unlist(tmp[,task,with=F]))))
  }
}


cor_res_avg = data.table(pair=cor_res_pt, task=cor_res_task, cor=cor_res_value)

pair_order_man = c('same', 'cross', 'amino_acid-amino_acid', 
                   'literature-literature', 'multi-multi', 'PPI-PPI',
                   'exp_bulk-exp_bulk',
                   'exp_sc-exp_sc', 'exp_bulk-exp_sc',
                   'amino_acid-exp_bulk', 'amino_acid-exp_sc', 'amino_acid-literature',
                   'amino_acid-multi', 'amino_acid-PPI', 'exp_bulk-literature',
                   'exp_sc-literature', 'exp_bulk-multi', 'exp_sc-multi', 'exp_bulk-PPI',
                   'exp_sc-PPI', 'literature-multi', 'literature-PPI', 'multi-PPI')
cor_res_avg$pair = factor(cor_res_avg$pair, levels = pair_order_man)

cor_res_avg$task = factor(cor_res_avg$task, levels = 
                            rev(c('s.go.cor', 's.omim.cor', 's.temporal.cor',
                                  'p.cat.cor', 'p.sum.cor', 'p.prod.cor',
                                  'kegg.go.cor','dis.tis.cor','go.go.temporal.cor',
                                  'cca')))


lower = ggplot(cor_res_avg[task != 'cca'], aes(x = pair, y = task, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = "cor") +
  theme_classic() + coord_equal() + ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

upper = ggplot(cor_res_avg[task == 'cca'], aes(x = pair, y = task, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = rev(colorRampPalette(brewer.pal(5, "PuBu"))(100))) +
  theme_void() + coord_equal() + ylab('') + xlab('') 

upper / lower

ggsave('results/plots/cor_sum_heatmap.pdf', width=5, height=6)


