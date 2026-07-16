import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("LJ_energy.csv")
df = df.dropna()
plt.figure(figsize=(10, 5))
plt.scatter(df['time'], df['Kinetic energy'], color='#004f77', s=10, alpha=0.6, label="Kinetic energy")
plt.scatter(df['time'], df['Potential energy'], color="#ff0000", s=10, alpha=1.0, label="Potential energy")
plt.scatter(df['time'], df['Total energy'], color="#1a7700", s = 10, label="Total energy")

plt.title("Energy Stores of the LJ System", fontsize=14)
plt.xlabel("Time (seconds)", fontsize=12)
plt.ylabel("Energy (Joules)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.ylim(bottom=0)
plt.legend()
plt.ylim(np.min(df['Potential energy']) + 0.1, np.max(df['Kinetic energy']) + 0.1)

# Save and show the plot
plt.savefig('energy_plot.png', dpi=300, bbox_inches='tight')
plt.show()