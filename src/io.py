"""
Input/output utilities for reading and writing microtubule geometry data 
using Pandas DataFrames and NumPy arrays.

Handles loading of 2D centerline templates and 3D microtubule geometry 
CSVs, as well as formatting and writing the final 3D geometry output.
"""
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

# Assuming default_config contains necessary constants like MAX_MT_LENGTH, etc.
from . import default_config
# Assuming this module exists for tip length generation
from .geometry.tip import generate_multiple_tip_lengths_in_memory 


# --- Global Constants for Microtubule Data Structure ---

REQUIRED_COLUMNS = [
    "DoubletNumber", "X", "Y", "Z", 
    "Idx_A", "Idx_B", "Angle", 
    "A_Shift", "B_Shift"
]
"""Minimum set of columns required for a base DMT/TMT geometry file (Doublet/Triplet)."""

EXTENDED_COLUMNS = [
    'DoubletNumber', 'X', 'Y', 'Z', 
    'Idx_A', 'Idx_B', 'Idx_C', 
    'Angle', 
    'A_Shift', 'B_Shift', 'C_Shift'
]
"""The full set of columns expected for output/storage, including the 'C' tubule data."""

MIN_2D_POINTS = 100
"""Minimum number of points required for a 2D centerline template."""


# --- Core I/O Functions ---

def load_template_data(template_filename: str) -> pd.DataFrame:
    """
    Loads microtubule geometry data from a template CSV file located 
    relative to the script's directory.

    The function assumes the template file is located in a 'geometry'
    subdirectory relative to the current script's parent directory.
    It then calls `read_3d_csv` to validate and load the data.

    Args:
        template_filename: The name of the template CSV file (e.g., 'dmt_template.csv').

    Returns:
        pandas.DataFrame: The DataFrame containing the validated 3D geometry data
                          (limited to REQUIRED_COLUMNS).

    Raises:
        FileNotFoundError: If the template file cannot be found at the calculated path.
    """
    # 1. Get the path to the current Python file (where this function is defined)
    script_path = Path(__file__).resolve()

    # 2. Navigate to the 'geometry' directory to find the template file.
    template_file_path = script_path.parent / 'geometry' / template_filename
    
    if not template_file_path.exists():
        # Raise an exception if the file is not found
        raise FileNotFoundError(f"Template file not found at: {template_file_path}")
    
    # 3. Load the data using the absolute path
    df = read_3d_csv(template_file_path)
    
    return df

def read_3d_csv(file_path: str) -> pd.DataFrame | None:
    """
    Reads a 3D microtubule geometry CSV file and validates the presence 
    of all columns listed in `REQUIRED_COLUMNS`.

    On validation failure, an error message is printed to stderr and None is returned.

    Args:
        file_path: The path to the input CSV file.

    Returns:
        pandas.DataFrame | None: The DataFrame containing the validated data, 
                                 or None if the file is missing, empty, unparsable, or 
                                 lacks required columns.
    """
    # --- 1. Check if file exists ---
    if not os.path.exists(file_path):
        print(f"Error: File not found at path '{file_path}'", file=sys.stderr)
        return None

    # --- 2. Read the CSV ---
    try:
        df = pd.read_csv(file_path)
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{file_path}' is empty.", file=sys.stderr)
        return None
    except pd.errors.ParserError:
        print(f"Error: Could not parse '{file_path}'. Check if it is a valid CSV format.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}", file=sys.stderr)
        return None

    # --- 3. Validate Columns ---
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if missing_columns:
        print("\n" + "="*50, file=sys.stderr)
        print("Error: The CSV file is missing critical columns.", file=sys.stderr)
        print(f"Missing columns: {missing_columns}", file=sys.stderr)
        print("="*50 + "\n", file=sys.stderr)
        return None
    
    # Return only the validated columns
    return df[REQUIRED_COLUMNS] 
    
def read_2d_csv(file_path: str) -> np.ndarray:
    """
    Reads a 2D centerline CSV file (e.g., for flagellar curvature), 
    validates it contains 'X' and 'Y' columns, and ensures it has a 
    minimum number of points (defined by `MIN_2D_POINTS`).

    Args:
        file_path: The path to the input CSV file, expected to have 'X' and 'Y' headers.

    Returns:
        np.ndarray: A NumPy array (N, 2) containing only the 'X' and 'Y' coordinates.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the file is corrupted, missing required headers, or has too few points.
    """
    REQUIRED_2D_COLUMNS = ['X', 'Y']

    # --- 1. Load and Validate Template Data ---
    try:
        # Load the CSV file into a Pandas DataFrame
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        # Re-raise the FileNotFoundError with context
        raise FileNotFoundError(f"Template file not found: '{file_path}'")
    except Exception as e:
        # Catch other potential loading errors (e.g., corrupted file)
        raise ValueError(f"Could not load template file '{file_path}'. Check format: {e}")

    # --- 2. Header and Column Validation ---
    
    # Check for required columns
    missing_columns = [col for col in REQUIRED_2D_COLUMNS if col not in df.columns]
    if missing_columns:
        raise ValueError(
            f"Template file '{file_path}' is missing required header columns: "
            f"{missing_columns}. Expected columns: {REQUIRED_2D_COLUMNS}."
        )

    # Check for point count
    if len(df) < MIN_2D_POINTS:
        raise ValueError(
            f"Template file must have at least {MIN_2D_POINTS} points, got {len(df)}."
        )
    
    # --- 3. Extract Data as NumPy Array ---
    # Select only the 'X' and 'Y' columns and convert to a NumPy array
    template_data = df[REQUIRED_2D_COLUMNS].values

    return template_data

def write_3d_csv(df: pd.DataFrame, file_path: str):
    """
    Writes a DataFrame to a CSV file, ensuring all columns from 
    `EXTENDED_COLUMNS` are present, correctly ordered, and formatted.
    
    If 'Idx_C' or 'C_Shift' are missing from the input DataFrame, they 
    are automatically added with a default value of 0.

    Position columns ('X', 'Y', 'Z') are formatted to two decimal places.

    Args:
        df: The DataFrame containing the microtubule geometry data.
        file_path: The full path where the CSV file should be saved.
    """
    
    # 1. Create a working copy of the DataFrame to modify
    df_output = df.copy()

    # 2. Ensure all required columns for the final output are present.
    # Add 'Idx_C' and 'C_Shift' with default 0 if they don't exist (Triplets)
    if 'Idx_C' not in df_output.columns:
        df_output['Idx_C'] = 0
        
    if 'C_Shift' not in df_output.columns:
        df_output['C_Shift'] = 0
        
    # Check for missing core geometry columns ('X', 'Y', 'Z' etc.)
    missing_core_columns = [col for col in EXTENDED_COLUMNS if col not in df_output.columns]
    if missing_core_columns:
        error_msg = (
            f"Cannot write CSV: Input DataFrame is missing core geometry columns: "
            f"{missing_core_columns}"
        )
        print(error_msg, file=sys.stderr)
        return

    # 3. Check for extra columns and print a warning
    extra_columns = [col for col in df_output.columns if col not in EXTENDED_COLUMNS]
    if extra_columns:
        print("="*50, file=sys.stderr)
        print("WARNING: Extra columns found in the input DataFrame.", file=sys.stderr)
        print(f"Ignored columns: {extra_columns}", file=sys.stderr)
        print(f"Only columns: {EXTENDED_COLUMNS} will be written.", file=sys.stderr)
        print("="*50, file=sys.stderr)

    # 4. Select only the required columns, format, and write to CSV
    try:
        # Select columns in the required order from EXTENDED_COLUMNS
        df_output = df_output[EXTENDED_COLUMNS]
        
        # Apply formatting to position columns
        df_output['X'] = df_output['X'].apply(lambda x: '{:.2f}'.format(x))
        df_output['Y'] = df_output['Y'].apply(lambda x: '{:.2f}'.format(x))
        df_output['Z'] = df_output['Z'].apply(lambda x: '{:.2f}'.format(x))

        df_output.to_csv(file_path, index=False)
        print(f"\nSuccessfully wrote data to '{file_path}'.")
        print(f"File includes only the required columns.")
        
    except Exception as e:
        print(f"Error occurred while writing to CSV: {e}", file=sys.stderr)