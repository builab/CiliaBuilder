import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv('cilia_coordinates.csv')

# Extract X and Y coordinates and convert to numpy arrays
x = df['X'].values
y = df['Y'].values

# Create the plot
plt.figure(figsize=(10, 12))
plt.plot(x, y, 'b-', linewidth=2, label='Cilia Curve')
plt.scatter(x[0], y[0], color='green', s=100, zorder=5, label='Start (0,0)')
plt.scatter(x[-1], y[-1], color='red', s=100, zorder=5, label='End')

# Customize the plot
plt.xlabel('X', fontsize=12)
plt.ylabel('Y', fontsize=12)
plt.title('Cilia Curve Plot', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.legend()
plt.axis('equal')
plt.tight_layout()

# Display the plot
plt.show()

# Optional: Save the plot to a file
# plt.savefig('cilia_curve.png', dpi=300, bbox_inches='tight')