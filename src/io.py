import pandas as pd
import csv
from . import default_config

# --- MODIFIED save_curves_to_csv FUNCTION ---
def save_curves_to_csv(curves, filename='cilia_doublets.csv'):
    """
    Save curve coordinates to CSV file with the specific doublet format:
    'DoubletNumber', 'X', 'Y', 'Z', 'Idx_A', 'Idx_B','Angle','A_Shift', 'B_Shift'
    """
    all_data = []
    
    # Define the maximum Z-coordinate for the '1' region
    Z_MAX_IDX_B = INITIAL_LENGTH + TRANSITION_LENGTH # 300 + 2000 = 2300
    
    for i, curve_data in enumerate(curves):
        curve = curve_data['curve']
        doublet_number = i + 1 
        
        # Calculate Angle
        angle = 90 + 40 * (doublet_number - 1)
        
        # Define constant values
        idx_a = 1
        a_shift = -CILIA_DOUBLET_SHIFT
        b_shift = CILIA_DOUBLET_SHIFT
        
        # --- MODIFICATION 3: Randomized Idx_B ---
        # 1. Randomly select the transition Z-coordinate for this specific curve
        idx_b_transition_z = random.uniform(0, Z_MAX_IDX_B)
        
        # 2. Iterate over all points and assign Idx_B conditionally
        for point in curve:
            x, y, z = point
            
            # Idx_B is 1 if Z is within the randomized region, 0 otherwise
            idx_b = 1 if z <= idx_b_transition_z else 0
            
            # Append data in the required column order
            all_data.append([
                doublet_number, x, y, z, 
                idx_a, idx_b, angle, 
                a_shift, b_shift
            ])
    
    # Create the DataFrame
    columns = [
        'DoubletNumber', 'X', 'Y', 'Z', 
        'Idx_A', 'Idx_B', 'Angle', 
        'A_Shift', 'B_Shift'
    ]
    df = pd.DataFrame(all_data, columns=columns)
    
    # Save to CSV
    df.to_csv(filename, index=False)
    print(f"Curves saved to {filename} in the requested Doublet format.")
    return filename
    
# --- Create mixed CSV directly from in-memory data ---
def create_mixed_csv_from_memory(all_curves_data, output_filename='tip_cilia_doublets.csv'):
    """
    Create mixed CSV file by randomly selecting curves from different tip lengths.
    Works directly with in-memory curve data without creating intermediate files.
    """
    
    mixed_data = []
    num_doublets = 9
    
    # Define the maximum Z-coordinate for the '1' region
    Z_MAX_IDX_B = INITIAL_LENGTH + TRANSITION_LENGTH
    
    for line_num in range(1, num_doublets + 1):
        # Randomly select one of the tip length variations
        selected_data = random.choice(all_curves_data)
        selected_tip_length = selected_data['tip_length']
        selected_curves = selected_data['curves']
        
        print(f"Line {line_num}: Selected from tip_length={selected_tip_length}")
        
        # Get the specific curve for this line number
        curve_data = selected_curves[line_num - 1]
        curve = curve_data['curve']
        doublet_number = line_num
        
        # Calculate Angle
        angle = 90 + 40 * (doublet_number - 1)
        
        # Define constant values
        idx_a = 1
        a_shift = -70
        b_shift = 70
        
        # Randomized Idx_B transition
        idx_b_transition_z = random.uniform(0, Z_MAX_IDX_B)
        
        # Process all points for this curve
        for point in curve:
            x, y, z = point
            idx_b = 1 if z <= idx_b_transition_z else 0
            
            mixed_data.append([
                doublet_number, x, y, z, 
                idx_a, idx_b, angle, 
                a_shift, b_shift
            ])
    
    # Create DataFrame and save
    columns = [
        'DoubletNumber', 'X', 'Y', 'Z', 
        'Idx_A', 'Idx_B', 'Angle', 
        'A_Shift', 'B_Shift'
    ]
    mixed_df = pd.DataFrame(mixed_data, columns=columns)
    mixed_df.to_csv(output_filename, index=False)
    
    print(f"\nMixed CSV saved to {output_filename}")
    print(f"Total points: {len(mixed_df)}")
    print("=" * 60)
    
    return output_filename


def generate_tip_csv(cilia_radius=default_config.CILIA_RADIUS, tip_length_end=default_config.TIP_LENGTH, transition_radius=default_config.TIP_TRANSITION_RADIUS, final_radius=default_config.TIP_FINAL_RADIUS):
    
    tip_length_start = INITIAL_LENGTH + TRANSITION_LENGTH
    all_curves_data = generate_multiple_tip_lengths_in_memory(
        cilia_radius, transition_radius, final_radius,
        tip_length_start, tip_length_end, number_of_steps=10
    )
            
    # Create mixed CSV directly from in-memory data
    mixed_filename = create_mixed_csv_from_memory(all_curves_data)
            

# Example usage
# generate_tip_csv(cilia_radius=875, tip_length_end=5000, transition_radius=656, final_radius=460)