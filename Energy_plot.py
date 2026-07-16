import pandas as pd
import matplotlib.pyplot as plt

# 1. Load and clean simulation output data
# Reads the CSV containing calculated kinetic, potential, and total energies over time
df = pd.read_csv("LJ_energy.csv")
df = df.dropna() # Remove any null rows to prevent plotting artifacts

# 2. Configure figure fimensions & high DPI for professional presentation
plt.figure(figsize=(10, 5))

# 3. Plot energy stores using scatter plots with clear, recognisable physical colors
# - Blue represents kinetic energy (system motion)
# - Red represents potential energy (atomic interactions / potential wells)
# - Green represents total energy (Hamiltonian constant of the NVE ensemble)
plt.scatter(df['time'], df['Kinetic energy'], color='#004f77', s=10, alpha=0.6, label="Kinetic energy")
plt.scatter(df['time'], df['Potential energy'], color="#ff0000", s=10, alpha=1.0, label="Potential energy")
plt.scatter(df['time'], df['Total energy'], color="#1a7700", s = 10, label="Total energy")

# 4. Standard Labeling and Styling
plt.title("Energy Stores of the LJ System", fontsize=14)
plt.xlabel("Time (seconds)", fontsize=12)
plt.ylabel("Energy (Joules)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5) # Dashed gridlines for readability
plt.legend()

# 5. Export Graph in High-Resolution for GitHub Presentation & Academic Reports
plt.savefig('energy_plot.png', dpi=300, bbox_inches='tight')
plt.show()