import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("libsyn_schedule_dist.csv")
actual_makespan = 2070.74  # exact operation times
df['delta_percentage'] = ((df['makespan'] - actual_makespan) / actual_makespan) * 100


plt.figure(figsize=(6, 5))
sns.boxplot(data=df, x='dist_cv', y='delta_percentage', hue='dist_type')
# plt.axhline(0, color='red', linestyle='--', label='Actual Makespan (0% delta)')
# plt.title(r'Distribution of Percentage \Delta Makespan by dist_cv and dist_type')
plt.xlabel('Coefficient of variation in operation time')
plt.ylabel(r'$\Delta$ makespan (%)')
plt.legend(title='Distribution type')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("libsyn_schedule_dist_plot.png", dpi=600)
