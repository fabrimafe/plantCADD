import pandas as pd
import sys
import matplotlib.pyplot as plt

df_name = sys.argv[1]  #scores.txt

# Inits
df = pd.read_csv(df_name, sep='\t', header=None, names=['score', 'AF_pop'])
common_af_threshold = 0.01
rare_af_threshold = 0.001

# Initialize lists for metrics
quantiles = []
precision_metrics = []
sensitivity_metrics = []
specificity_metrics = []
f1_metrics = []

# Iterate over quantile values from 0 to 1 with a step of 0.1
for quantile_value in [0.1,0.5,0.9,0.99,0.995,0.999,0.9999]:
    print(f'Processing quantile: {quantile_value}')
    
    # Calculate rates
    quantile_threshold = df['score'].quantile(quantile_value)
    
    true_pos  = len(df[(df['score'] >= quantile_threshold) & (df['AF_pop'] <= rare_af_threshold)])
    false_pos = len(df[(df['score'] >= quantile_threshold) & (df['AF_pop'] >= common_af_threshold)])
    false_neg = len(df[(df['score'] <  quantile_threshold) & (df['AF_pop'] <= rare_af_threshold)])
    true_neg  = len(df[(df['score'] <  quantile_threshold) & (df['AF_pop'] >= common_af_threshold)])
    
    # Calculate metrics
    percision=true_pos/(true_pos+false_pos)
    sensitivity=true_pos/(true_pos+false_neg)
    specificity=true_neg/(true_neg+false_pos)
    f1=2*(percision*sensitivity)/(percision+sensitivity)
    
    # Append metrics to lists
    quantiles.append(quantile_value)
    precision_metrics.append(percision)
    sensitivity_metrics.append(sensitivity)
    specificity_metrics.append(specificity)
    f1_metrics.append(f1)
    
    
# Plot
plt.figure(figsize=(8, 8))

bar_width = 0.1
index = range(len(quantiles))

plt.bar(index, precision_metrics, width=bar_width, label='Precision')
plt.bar([i + bar_width for i in index], sensitivity_metrics, width=bar_width, label='Sensitivity')
plt.bar([i + 2 * bar_width for i in index], specificity_metrics, width=bar_width, label='Specificity')
plt.bar([i + 3 * bar_width for i in index], f1_metrics, width=bar_width, label='F1 Score')

plt.xlabel('Quantiles')
plt.ylabel('Metrics')
plt.xticks([i + 1.5 * bar_width for i in index], quantiles)
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig('odds_ratio_plot.svg', dpi=300)  
