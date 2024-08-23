import pandas as pd
import sys


df_name = sys.argv[1]  #scores.txt

# Inits
df = pd.read_csv(df_name, sep='\t', header=None, names=['score', 'AF_pop'])
common_af_threshold = 0.01
rare_af_threshold = 0.001

# Initialize lists for metrics
quantiles = []
accuracy_metrics=[]
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
    accuracy=(true_pos+true_neg)/(true_pos+false_pos+true_neg+false_neg)
    percision=true_pos/(true_pos+false_pos)
    sensitivity=true_pos/(true_pos+false_neg)
    specificity=true_neg/(true_neg+false_pos)
    f1=2*(percision*sensitivity)/(percision+sensitivity)
    
    # Append metrics to lists
    quantiles.append(quantile_value)
    accuracy_metrics.append(accuracy)
    precision_metrics.append(percision)
    sensitivity_metrics.append(sensitivity)
    specificity_metrics.append(specificity)
    f1_metrics.append(f1)
    
metrics_df = pd.DataFrame({
    'Quantile': quantiles,
    'Accuracy': accuracy_metrics, 
    'Precision': precision_metrics,
    'Sensitivity': sensitivity_metrics,
    'Specificity': specificity_metrics,
    'F1 Score': f1_metrics
})
metrics_df.to_csv('metrics.txt', sep='\t', index=False, header=True)