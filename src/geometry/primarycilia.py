# Read CSV

# CILIA_LENTGTH=10000

# SCALE CSV TO THE RANGE

# INITIAL_LENGTH = 200

# EXTEND THE CSV with INITIAL_LENGTH

# BTUBULE= 25%
# A & BTUBULE3 = 30-35% of Z

# At least 4 go to 90%
# Other from 55-85 random (55, 65, 75, 85)

"""
Zero column Idx_A and Idx_B
Find the max Z value from everything.
Read the csv file, for DoubletNumber=3, pick a value between 0.3 - 0.35 of Z and set  Idx_A & Idx_B to 1 from Z=0 to that Z value.

For all other DoubletNumber, set the Idx_B from Z=20 to Z = 0.25*maxZ  of value 1.
At least 4 of the DoubletNumber (except DoubletNumber 3), set Idx_A randomly of 0.90, 0.925, 0.95 and 0.975*maxZ, no duplication
The rest of the 4 DoubletNumber, set Idx_A randomonly of .55, .65, .75, .85*maxZ

then remove all line with both Idx_A=0 and Idx_B=0.

Write the new csv as primarycilia.csv

All the way, 1, 2, 7, 8
Mid way, 4, 5, 6, 9

[0.95, 0.9, 0.4, 0.65, 0.80, 0.55, 0.975, 0.925, 0.75]
[0.25, 0.25, 0.4, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
"""
"""
import pandas as pd
import numpy as np

# Read the CSV file
df = pd.read_csv('primarycilia_template.csv')

# Set Angle values for each DoubletNumber
# Doublet1 = 90, Doublet2 = 130, Doublet3 = 170, etc. (increment by 40)
for doublet_num in range(1, 10):
    angle_value = 90 + (doublet_num - 1) * 40
    df.loc[df['DoubletNumber'] == doublet_num, 'Angle'] = angle_value
    print(f"DoubletNumber {doublet_num}: Angle = {angle_value}")

print()

# Zero out Idx_A and Idx_B columns
df['Idx_A'] = 0
df['Idx_B'] = 0

# Find maximum Z value
max_z = df['Z'].max()
print(f"Maximum Z value: {max_z}")

# Process DoubletNumber = 3
# Pick a random value between 0.3 and 0.35 of max_z
z_threshold_3 = np.random.uniform(0.35, 0.45) * max_z
print(f"DoubletNumber 3 threshold: {z_threshold_3}")

# Set Idx_A and Idx_B to 1 for DoubletNumber=3 where Z is between 0 and threshold
mask_doublet3 = (df['DoubletNumber'] == 3) & (df['Z'] >= 0) & (df['Z'] <= z_threshold_3)
df.loc[mask_doublet3, ['Idx_A', 'Idx_B']] = 1

# Process all other DoubletNumbers (1, 2, 4, 5, 6, 7, 8, 9)
# Set Idx_B to 1 for Z between 20 and 0.25*max_z
other_doublets = [1, 2, 4, 5, 6, 7, 8, 9]
z_threshold_b = 0.25 * max_z

for doublet in other_doublets:
    mask_b = (df['DoubletNumber'] == doublet) & (df['Z'] >= 20) & (df['Z'] <= z_threshold_b)
    df.loc[mask_b, 'Idx_B'] = 1

# Randomly assign Idx_A values to the 8 other doublets
# First 4 get high values: 0.90, 0.925, 0.95, 0.975
# Next 4 get lower values: 0.55, 0.65, 0.75, 0.85
np.random.shuffle(other_doublets)

high_values = [0.90, 0.925, 0.95, 0.975]
low_values = [0.55, 0.65, 0.72, 0.80]

# Assign high values to first 4 doublets
# Set Idx_A to 1 for Z from 0 to the threshold
for i, doublet in enumerate(other_doublets[:4]):
    z_threshold_a = high_values[i] * max_z
    mask_a = (df['DoubletNumber'] == doublet) & (df['Z'] >= 0) & (df['Z'] <= z_threshold_a)
    df.loc[mask_a, 'Idx_A'] = 1
    print(f"DoubletNumber {doublet}: Idx_A set for Z from 0 to {z_threshold_a:.2f} (high: {high_values[i]})")

# Assign low values to next 4 doublets
# Set Idx_A to 1 for Z from 0 to the threshold
for i, doublet in enumerate(other_doublets[4:]):
    z_threshold_a = low_values[i] * max_z
    mask_a = (df['DoubletNumber'] == doublet) & (df['Z'] >= 0) & (df['Z'] <= z_threshold_a)
    df.loc[mask_a, 'Idx_A'] = 1
    print(f"DoubletNumber {doublet}: Idx_A set for Z from 0 to {z_threshold_a:.2f} (low: {low_values[i]})")

# Remove all rows where both Idx_A=0 and Idx_B=0
df_filtered = df[~((df['Idx_A'] == 0) & (df['Idx_B'] == 0))]

print(f"\nOriginal rows: {len(df)}")
print(f"Filtered rows: {len(df_filtered)}")
print(f"Rows removed: {len(df) - len(df_filtered)}")

# Write to new CSV file
df_filtered.to_csv('primarycilia.csv', index=False)
print("\nProcessed file saved as 'primarycilia.csv'")
"""

import csv
import pandas as pd
import numpy as np

def process_cilia_csv(df_in, idx_a_thresholds, idx_b_thresholds):
    """
    Process primary cilia CSV file by setting Angle, Idx_A, and Idx_B values.
    
    Parameters:
    input_file (str): Path to input CSV file
    output_file (str): Path to output CSV file
    """
    
    
    # Set Angle values for each DoubletNumber
    # Doublet1 = 90, Doublet2 = 130, Doublet3 = 170, etc. (increment by 40)
    print("Setting Angle values:")
    for doublet_num in range(1, 10):
        angle_value = 90 + (doublet_num - 1) * 40
        df.loc[df['DoubletNumber'] == doublet_num, 'Angle'] = angle_value
        print(f"  DoubletNumber {doublet_num}: Angle = {angle_value}°")
    
    print()
    
    # Zero out Idx_A and Idx_B columns
    df['Idx_A'] = 0
    df['Idx_B'] = 0
    
    # Find maximum Z value
    max_z = df['Z'].max()
    print(f"Maximum Z value: {max_z}")
    print()
    
    # Process each DoubletNumber
    print("Setting Idx_A and Idx_B values:")
    for doublet_num in range(1, 10):
        idx = doublet_num - 1  # Array index (0-8)
        
        # Calculate Z thresholds
        z_threshold_a = idx_a_thresholds[idx] * max_z
        z_threshold_b = idx_b_thresholds[idx] * max_z
        
        # Set Idx_A to 1 for Z from 0 to z_threshold_a
        mask_a = (df['DoubletNumber'] == doublet_num) & (df['Z'] >= 0) & (df['Z'] <= z_threshold_a)
        df.loc[mask_a, 'Idx_A'] = 1
        
        # Set Idx_B to 1 for Z from 20 to z_threshold_b
        mask_b = (df['DoubletNumber'] == doublet_num) & (df['Z'] >= 20) & (df['Z'] <= z_threshold_b)
        df.loc[mask_b, 'Idx_B'] = 1
        
        print(f"  DoubletNumber {doublet_num}:")
        print(f"    Idx_A = 1 for Z: 0 to {z_threshold_a:.2f} ({idx_a_thresholds[idx]} * max_z)")
        print(f"    Idx_B = 1 for Z: 20 to {z_threshold_b:.2f} ({idx_b_thresholds[idx]} * max_z)")
        
    # Remove all rows where both Idx_A=0 and Idx_B=0
    df_filtered = df[~((df['Idx_A'] == 0) & (df['Idx_B'] == 0))]
    
    print(f"Original rows: {len(df)}")
    print(f"Filtered rows: {len(df_filtered)}")
    print(f"Rows removed: {len(df) - len(df_filtered)}")
    

    
    return df_filtered


# Run the function
if __name__ == "__main__":
    # Define thresholds for Idx_A and Idx_B for each doublet (1-9)
    idx_a_thresholds = [0.95, 0.92, 0.4, 0.94, 1.0, 0.55, 0.80, 0.65, 0.75]
    idx_b_thresholds = [0.20, 0.20, 0.4, 0.20, 0.20, 0.22, 0.23, 0.45, 0.20]
    
    # Read the CSV file
    df = pd.read_csv('primarycilia_template.csv')
    df_filtered = process_cilia_csv(df, idx_a_thresholds, idx_b_thresholds)

    # Write to new CSV file
    df_filtered.to_csv('primarycilia_fixed.csv', index=False)
