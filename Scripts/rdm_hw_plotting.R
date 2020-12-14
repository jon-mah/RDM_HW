library(ggplot2)
setwd("C:/Users/jonat/Desktop/GitHub/RDM_HW")

exact_vs_approx = read.csv('./Data/exact_vs_approximate_exp_homo_b.csv')
max_diff_data = read.csv('./Data/exact_vs_approximate_max_diff.csv')

p_exact_vs_approx_diff <- ggplot(data=exact_vs_approx, aes(x=n_ton, y=diff)) +
  geom_point(color='black') +
  ggtitle('Absolute difference between exact and approximate expected homozygotes') +
  xlab('n-ton count for a sample of 100 individuals') +
  ylab('Absolute difference in expected homozygotes')
  
p_exact_vs_approx_diff

p_exact_vs_approx_fraction <- ggplot(data=exact_vs_approx, aes(x=n_ton, y=fraction)) +
  geom_point(color='black') +
  ggtitle('Fractional difference between exact and approximate homozygotes') +
  xlab('n-ton count for a sample of 100 individuals') +
  ylab('Fractional difference in expected homozygotes') +
  scale_y_continuous(trans='log10')
p_exact_vs_approx_fraction

p_exact_vs_approx_percent_diff <- ggplot(data=exact_vs_approx, aes(x=n_ton, y=percent_diff)) +
  geom_point(color='black') +
  ggtitle('Difference / n_ton for exact and approximate homozygotes') +
  xlab('n-ton count for a sample of 100 individuals') +
  ylab('Difference / n_ton count')
p_exact_vs_approx_percent_diff

p_max_diff <- ggplot(data=max_diff_data, aes(x=num_ind, y=max_diff)) +
  geom_point(color='black') +
  ggtitle('Maximum difference vs. number of individuals') +
  xlab('Number of individuals') +
  ylab('Maximum difference between exact and proximate expectation')
p_max_diff

p_fractional_max_diff <- ggplot(data=max_diff_data, aes(x=num_ind, fractional_max_diff)) +
  geom_point(color='black') +
  ggtitle('Maximum difference / number of indvidiuals') +
  xlab('Number of individuals') +
  ylab('Maximum difference vs. number of individauls')
p_fractional_max_diff

s_0_rec_10 = read.csv('./Data/s_0_rec/sample_10_summary.csv')

s_0_001_add_1000 = read.csv('./Data/s_0_001_add/sample_1000_summary.csv') # 0_001_add_1000
s_0_001_rec_1000 = read.csv('./Data/s_0_001_rec/sample_1000_summary.csv') # 0_001_rec_1000
s_0_01_add_1000 = read.csv('./Data/s_0_01_add/sample_1000_summary.csv') # 0_01_add_1000
s_0_01_rec_1000 = read.csv('./Data/s_0_01_rec/sample_1000_summary.csv') # 0_01_rec_1000
s_0_1_add_1000 = read.csv('./Data/s_0_1_add/sample_1000_summary.csv') # 0_1_add_1000
s_0_1_rec_1000 = read.csv('./Data/s_0_1_rec/sample_1000_summary.csv') # 0_1_rec_1000
s_0_add_1000 = read.csv('./Data/s_0_add/sample_1000_summary.csv') # 0_1000
s_0_rec_1000 = read.csv('./Data/s_0_rec/sample_1000_summary.csv') # 0_1000
# chi-squared p-values
chi_p_s_0_add_1000 = s_0_add_1000$chi_pr_values
chi_p_s_0_001_add_1000 = s_0_001_add_1000$chi_pr_values
chi_p_s_0_01_add_1000 = s_0_01_add_1000$chi_pr_values
chi_p_s_0_1_add_1000 = s_0_1_add_1000$chi_pr_values
chi_p_s_0_rec_1000 = s_0_rec_1000$chi_pr_values
chi_p_s_0_001_rec_1000 = s_0_001_rec_1000$chi_pr_values
chi_p_s_0_01_rec_1000 = s_0_01_rec_1000$chi_pr_values
chi_p_s_0_1_rec_1000 = s_0_1_rec_1000$chi_pr_values

# fisher p-values (leq)
fisher_p_s_0_add_1000 = s_0_add_1000$fisher_p_n_homo_leq
fisher_p_s_0_001_add_1000 = s_0_001_add_1000$fisher_p_n_homo_leq
fisher_p_s_0_01_add_1000 = s_0_01_add_1000$fisher_p_n_homo_leq
fisher_p_s_0_1_add_1000 = s_0_1_add_1000$fisher_p_n_homo_leq
fisher_p_s_0_rec_1000 = s_0_rec_1000$fisher_p_n_homo_leq
fisher_p_s_0_001_rec_1000 = s_0_001_rec_1000$fisher_p_n_homo_leq
fisher_p_s_0_01_rec_1000 = s_0_01_rec_1000$fisher_p_n_homo_leq
fisher_p_s_0_1_rec_1000 = s_0_1_rec_1000$fisher_p_n_homo_leq

log_chi_p_s_0_add_1000 = -log10(sort(chi_p_s_0_add_1000))
log_chi_p_s_0_rec_1000 = -log10(sort(chi_p_s_0_rec_1000))
log_chi_p_s_0_001_add_1000 = -log10(sort(chi_p_s_0_001_add_1000))
log_chi_p_s_0_001_rec_1000 = -log10(sort(chi_p_s_0_001_rec_1000))
log_chi_p_s_0_01_add_1000 = -log10(sort(chi_p_s_0_01_add_1000))
log_chi_p_s_0_01_rec_1000 = -log10(sort(chi_p_s_0_01_rec_1000))
log_chi_p_s_0_1_add_1000 = -log10(sort(chi_p_s_0_1_add_1000))
log_chi_p_s_0_1_rec_1000 = -log10(sort(chi_p_s_0_1_rec_1000))

log_fisher_p_s_0_add_1000 = -log10(sort(fisher_p_s_0_add_1000))
log_fisher_p_s_0_rec_1000 = -log10(sort(fisher_p_s_0_rec_1000))
log_fisher_p_s_0_001_add_1000 = -log10(sort(fisher_p_s_0_001_add_1000))
log_fisher_p_s_0_001_rec_1000 = -log10(sort(fisher_p_s_0_001_rec_1000))
log_fisher_p_s_0_01_add_1000 = -log10(sort(fisher_p_s_0_01_add_1000))
log_fisher_p_s_0_01_rec_1000 = -log10(sort(fisher_p_s_0_01_rec_1000))
log_fisher_p_s_0_1_add_1000 = -log10(sort(fisher_p_s_0_1_add_1000))
log_fisher_p_s_0_1_rec_1000 = -log10(sort(fisher_p_s_0_1_rec_1000))

y_axis_0_001_add_1000 = -log10(1:nrow(s_0_001_add_1000) / nrow(s_0_001_add_1000))
y_axis_0_01_add_1000 = -log10(1:nrow(s_0_01_add_1000) / nrow(s_0_01_add_1000))
y_axis_0_1_add_1000 = -log10(1:nrow(s_0_1_add_1000) / nrow(s_0_1_add_1000))
y_axis_0_1000 = -log10(1:nrow(s_0_add_1000) / nrow(s_0_add_1000))
y_axis_0_001_rec_1000 = -log10(1:nrow(s_0_001_rec_1000) / nrow(s_0_001_rec_1000))
y_axis_0_01_rec_1000 = -log10(1:nrow(s_0_01_rec_1000) / nrow(s_0_01_rec_1000))
y_axis_0_1_rec_1000 = -log10(1:nrow(s_0_1_rec_1000) / nrow(s_0_1_rec_1000))

data_0_001_add_1000 = data.frame(y_axis_0_001_add_1000, log_chi_p_s_0_001_add_1000, log_fisher_p_s_0_001_add_1000)
data_0_1000 = data.frame(y_axis_0_1000, log_chi_p_s_0_add_1000, log_chi_p_s_0_rec_1000,log_fisher_p_s_0_add_1000, log_fisher_p_s_0_rec_1000)
data_0_001_rec_1000 = data.frame(y_axis_0_001_rec_1000, log_chi_p_s_0_001_rec_1000, log_fisher_p_s_0_001_rec_1000)
data_0_01_add_1000 = data.frame(y_axis_0_01_add_1000, log_chi_p_s_0_01_add_1000, log_fisher_p_s_0_01_add_1000)
data_0_1_add_1000 = data.frame(y_axis_0_1_add_1000, log_chi_p_s_0_1_add_1000, log_fisher_p_s_0_1_add_1000)
data_0_01_rec_1000 = data.frame(y_axis_0_01_rec_1000, log_chi_p_s_0_01_rec_1000, log_fisher_p_s_0_01_rec_1000)
data_0_1_rec_1000 = data.frame(y_axis_0_1_rec_1000, log_chi_p_s_0_1_rec_1000, log_fisher_p_s_0_1_rec_1000)

p_0_001_add_1000 = ggplot(data=data_0_001_add_1000, aes(y=log_chi_p_s_0_001_add_1000, x=y_axis_0_001_add_1000, color='chi_additive')) + 
  geom_point() +
  geom_point(aes(y=log_fisher_p_s_0_001_add_1000, x=y_axis_0_001_add_1000, color='fisher_additive')) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green'),
                     name='P-value calculation',
                     breaks=c('chi_additive', 'fisher_additive'),
                     labels=c('chi_additive', 'fisher_additive')) +
  ggtitle('s = 0.001') +
  geom_abline(slope=1) +
  xlim(0, 4) +
  ylim(0, 4) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_001_add_1000

p_0_001_rec_1000 = ggplot(data=data_0_001_rec_1000, aes(y=log_chi_p_s_0_001_rec_1000, x=y_axis_0_001_rec_1000, color='chi_recessive')) + 
  geom_point() +
  geom_point(aes(y=log_fisher_p_s_0_001_rec_1000, x=y_axis_0_001_rec_1000, color='fisher_recessive')) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green'),
                     name='P-value calculation',
                     breaks=c('chi_recessive', 'fisher_recessive'),
                     labels=c('chi_recessive', 'fisher_recessive')) +
  ggtitle('s = 0.001') +
  geom_abline(slope=1) +
  xlim(0, 3) +
  ylim(0, 3) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_001_rec_1000

p_0_1000 = ggplot(data=data_0_1000, aes(y=log_chi_p_s_0_add_1000, x=y_axis_0_1000, color='chi_additive')) + 
  geom_point(shape=3) +
  geom_point(aes(y=log_fisher_p_s_0_add_1000, x=y_axis_0_1000, color='fisher_additive'), shape=3) +
  geom_point(aes(y=log_chi_p_s_0_rec_1000, x=y_axis_0_1000, color='chi_recessive'), shape=4) +
  geom_point(aes(y=log_fisher_p_s_0_rec_1000, x=y_axis_0_1000, color='fisher_recessive'), shape=4) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green', 'purple', 'orange'),
                     name='P-value calculation',
                     breaks=c('chi_additive', 'fisher_additive', 'chi_recessive', 'fisher_recessive'),
                     labels=c('chi_additive', 'fisher_additive', 'chi_recessive', 'fisher_recessive')) +
  ggtitle('s = 0') +
  geom_abline(slope=1) +
  xlim(0, 4) +
  ylim(0, 4) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_1000

p_0_01_add_1000 = ggplot(data=data_0_01_add_1000, aes(y=log_chi_p_s_0_01_add_1000, x=y_axis_0_01_add_1000, color='chi_additive')) + 
  geom_point() +
  geom_point(aes(y=log_fisher_p_s_0_01_add_1000, x=y_axis_0_01_add_1000, color='fisher_additive')) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green'),
                     name='P-value calculation',
                     breaks=c('chi_additive', 'fisher_additive'),
                     labels=c('chi_additive', 'fisher_additive')) +
  ggtitle('s = 0.01') +
  geom_abline(slope=1) +
  xlim(0, 4) +
  ylim(0, 4) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_01_add_1000

p_0_1_add_1000 = ggplot(data=data_0_1_add_1000, aes(y=log_chi_p_s_0_1_add_1000, x=y_axis_0_1_add_1000, color='chi_additive')) + 
  geom_point() +
  geom_point(aes(y=log_fisher_p_s_0_1_add_1000, x=y_axis_0_1_add_1000, color='fisher_additive')) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green'),
                     name='P-value calculation',
                     breaks=c('chi_additive', 'fisher_additive'),
                     labels=c('chi_additive', 'fisher_additive')) +
  ggtitle('s = 0.1') +
  geom_abline(slope=1) +
  xlim(0, 2) + 
  ylim(0, 2) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_1_add_1000

p_0_1_rec_1000 = ggplot(data=data_0_1_rec_1000, aes(y=log_chi_p_s_0_1_rec_1000, x=y_axis_0_1_rec_1000, color='chi_recessive')) + 
  geom_point() +
  geom_point(aes(y=log_fisher_p_s_0_1_rec_1000, x=y_axis_0_1_rec_1000, color='fisher_recessive')) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green'),
                     name='P-value calculation',
                     breaks=c('chi_recessive', 'fisher_recessive'),
                     labels=c('chi_recessive', 'fisher_recessive')) +
  ggtitle('s = 0.1') +
  geom_abline(slope=1) +
  xlim(0, 3)+ 
  ylim(0, 3) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_1_rec_1000

p_0_01_rec_1000 = ggplot(data=data_0_01_rec_1000, aes(y=log_chi_p_s_0_01_rec_1000, x=y_axis_0_01_rec_1000, color='chi_recessive')) + 
  geom_point() +
  geom_point(aes(y=log_fisher_p_s_0_01_rec_1000, x=y_axis_0_01_rec_1000, color='fisher_recessive')) +
  xlab('-log10 of expected p-values') + 
  ylab('-log10 of observed p-values') + 
  scale_color_manual(values=c('blue', 'green'),
                     name='P-value calculation',
                     breaks=c('chi_recessive', 'fisher_recessive'),
                     labels=c('chi_recessive', 'fisher_recessive')) +
  ggtitle('s = 0.01') +
  geom_abline(slope=1) +
  xlim(0, 3) +
  ylim(0, 3) +
  geom_hline(color = 'red', yintercept = -log10(0.05))
p_0_01_rec_1000

s_0_rec_10$allele_count = as.factor(s_0_rec_10$allele_count)
p_fisher_boxplot = ggplot(data=s_0_rec_10, aes(x=allele_count, y=fisher_p_n_homo_leq)) + 
  geom_boxplot() +
  ggtitle("s=0, Fisher's Exact Test") +
  xlab('allele count') +
  ylab('p-values')
p_fisher_boxplot

p_chi_boxplot = ggplot(data=s_0_rec_10, aes(x=allele_count, y=chi_pr_values)) +
  geom_boxplot() +
  ggtitle("s=0, chi-squared Test") +
  xlab('allele count') +
  ylab('p-values')
p_chi_boxplot
