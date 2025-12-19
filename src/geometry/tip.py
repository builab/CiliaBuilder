import numpy as np
from scipy.interpolate import CubicSpline

from .base import calculate_doublet_angle, calculate_radial_position

import pandas as pd
import sys
import os
import random

# --- CONSTANTS (Global variables for easy modification) ---
        # Default interval for dynamic point calculation
MAX_INTERVAL = 20 # Inveral
TIP_LENGTH=5000
CILIA_RADIUS=875
TRANSITION_RADIUS=656
FINAL_RADIUS=460
INITIAL_LENGTH = 300     # Z-coordinate for control point p2
TRANSITION_LENGTH = 2000 # Z-coordinate for control point p3
CILIA_DOUBLET_SHIFT = 70

def generate_tip_curves(tip_length=TIP_LENGTH, cilia_radius=CILIA_RADIUS, transition_radius=TRANSITION_RADIUS, final_radius=FINAL_RADIUS, initial_length=INITIAL_LENGTH, transition_length=TRANSITION_LENGTH,interval=MAX_INTERVAL):
    """
    Generate 9 cilia curves centered around the Z axis.
    
    Parameters:
    -----------
    tip_length : float
        Total length of the tip
    cilia_radius : float
        Radius of the circle at Z=0 (starting points)
    transition_radius : float
        Radius of the circle at Z=1/5*tip_length (third points)
    final_radius : float
        Radius of the circle at Z=tip_length (final points)
    
    Returns:
    --------
    curves : list of dicts
        List of 9 curves, each containing 'curve' (array of shape (n_points, 3)) 
        and 'control_points'.
    """
    num_doublets = 9
    curves = []
    
    # --- MODIFICATION 1: Dynamic Point Calculation ---
    
    # Segment 1: p1 to p2 (Length = INITIAL_LENGTH)
    n_points_linear = int(np.ceil(initial_length/ interval)) + 1
    n_points_linear = max(n_points_linear, 5) # Ensure a minimum of 5 points

    # Segment 2: p2 to p4 (Length = tip_length - INITIAL_LENGTH)
    spline_sampling_length = tip_length - initial_length
    n_points_spline = int(np.ceil(spline_sampling_length / interval)) + 1
    n_points_spline = max(n_points_spline, 5) # Ensure a minimum of 5 points
    
    # Z-coordinates for sampling the spline (for the smooth curve fix)
    z_spline_fine = np.linspace(initial_length, tip_length, n_points_spline)

    for i in range(num_doublets):
        angle = (i / num_doublets) * 2 * np.pi + np.pi # Temporary fix to match curve.py
        
        # Define 4 control points
        p1 = np.array([cilia_radius * np.cos(angle), cilia_radius * np.sin(angle), 0])
        p2 = np.array([p1[0], p1[1], initial_length])
        p3 = np.array([transition_radius * np.cos(angle), transition_radius * np.sin(angle), transition_length])
        p4 = np.array([final_radius * np.cos(angle), final_radius * np.sin(angle), tip_length])
        control_points = np.array([p1, p2, p3, p4])
        
        # --- Linear segment from p1 to p2 ---
        t_linear = np.linspace(0, 1, n_points_linear)
        linear_segment = np.array([
            p1 + t * (p2 - p1) for t in t_linear
        ])
        
        # --- MODIFICATION 2: Smooth Spline Segment (Z-Parametrization Fix) ---
        
        # 1. Define control points for X(Z) and Y(Z) splines
        z_control = np.array([p2[2], p3[2], p4[2]])
        x_control = np.array([p2[0], p3[0], p4[0]])
        y_control = np.array([p2[1], p3[1], p4[1]])
        
        # 2. Create X(Z) and Y(Z) Cubic Splines
        cs_x = CubicSpline(z_control, x_control)
        cs_y = CubicSpline(z_control, y_control)
        
        # 3. Sample the splines using the fine, increasing Z-coordinates
        x_new = cs_x(z_spline_fine)
        y_new = cs_y(z_spline_fine)
        
        # 4. Combine X, Y, Z to form the smooth segment
        smooth_segment = np.column_stack([x_new, y_new, z_spline_fine])
        
        # Combine segments (skip first point of smooth segment to avoid duplication)
        curve = np.vstack([linear_segment, smooth_segment[1:]])
        
        curves.append({
            'curve': curve,
            'control_points': control_points
        })
    
    return curves



# ---  Generate multiple tip lengths and return data in memory ---
def generate_multiple_tip_lengths_in_memory(cilia_radius=CILIA_RADIUS, transition_radius=TRANSITION_LENGTH, final_radius=FINAL_RADIUS, 
                                            tip_length_start=4000, tip_length_end=8000, number_of_steps=10):
    """
    Generate multiple tip length variations and return all curve data in memory.
    Does not create intermediate CSV files.
    """
    all_curves_data = []
    tip_lengths = np.linspace(tip_length_start, tip_length_end, number_of_steps)
    
    print("GENERATING MULTIPLE TIP LENGTH VARIATIONS (IN MEMORY)")
    
    for tip_length in tip_lengths:
        print(f"\nGenerating curves for tip_length = {tip_length}...")
        curves = generate_tip_curves(tip_length, cilia_radius, transition_radius, final_radius)
        all_curves_data.append({
            'tip_length': tip_length,
            'curves': curves
        })
    
    print(f"Generated {len(all_curves_data)} tip length variations")
    
    return all_curves_data


