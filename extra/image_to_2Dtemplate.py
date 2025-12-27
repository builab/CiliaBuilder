import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize
from scipy import ndimage
from scipy.interpolate import UnivariateSpline

# Load the image
image = cv2.imread('cilia_image.png', cv2.IMREAD_GRAYSCALE)

# Convert to binary (black and white)
# Invert if needed so the curve is white on black background
_, binary = cv2.threshold(image, 127, 255, cv2.THRESH_BINARY_INV)

# Skeletonize the binary image
skeleton = skeletonize(binary // 255).astype(np.uint8) * 255

# Extract skeleton coordinates
y_coords, x_coords = np.where(skeleton > 0)

# Sort points by y coordinate (from bottom to top)
sorted_indices = np.argsort(y_coords)[::-1]  # Reverse to go from bottom to top
x_sorted = x_coords[sorted_indices]
y_sorted = y_coords[sorted_indices]

# Flip y-coordinates so (0,0) is at bottom
y_sorted = image.shape[0] - y_sorted

# Normalize to start at (0,0)
x_sorted = x_sorted - x_sorted[0]
y_sorted = y_sorted - y_sorted[0]

# Remove duplicate y values and ensure strictly increasing
unique_indices = []
prev_y = -1
for i, y in enumerate(y_sorted):
    if y > prev_y:
        unique_indices.append(i)
        prev_y = y

x_sorted = x_sorted[unique_indices]
y_sorted = y_sorted[unique_indices]

# Fit smoothing spline (adjust 's' parameter for more/less smoothing)
# Higher s = more smoothing, lower s = closer to original points
# s=0 would be like CubicSpline (passes through all points)
# Try values like 100, 500, 1000, 5000 depending on desired smoothness
smoothing_factor = 1000
spline = UnivariateSpline(y_sorted, x_sorted, s=smoothing_factor)

# Generate smooth curve points (at least 200 points)
num_points = 200
y_smooth = np.linspace(y_sorted.min(), y_sorted.max(), num_points)
x_fit = spline(y_smooth)

# Save to CSV
csv_data = np.column_stack((x_fit, y_smooth))
np.savetxt('cilia_coordinates.csv', csv_data, delimiter=',', 
           header='X,Y', comments='', fmt='%.4f')

# Plotting
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Original binary image
axes[0, 0].imshow(binary, cmap='gray')
axes[0, 0].set_title('Binary Image')
axes[0, 0].axis('off')

# Skeletonized image
axes[0, 1].imshow(skeleton, cmap='gray')
axes[0, 1].set_title('Skeletonized')
axes[0, 1].axis('off')

# Extracted points
axes[1, 0].scatter(x_sorted, y_sorted, s=1, c='blue', alpha=0.5)
axes[1, 0].set_title('Extracted Points')
axes[1, 0].set_xlabel('X')
axes[1, 0].set_ylabel('Y')
axes[1, 0].grid(True, alpha=0.3)

# Smoothing spline fit
axes[1, 1].scatter(x_sorted, y_sorted, s=1, c='gray', alpha=0.3, label='Original points')
axes[1, 1].plot(x_fit, y_smooth, 'b-', linewidth=2, label=f'Smoothing Spline (s={smoothing_factor})')
axes[1, 1].set_title('Smoothing Spline Fit')
axes[1, 1].set_xlabel('X')
axes[1, 1].set_ylabel('Y')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('curve_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Smoothing spline fit complete (smoothing factor: {smoothing_factor})")
print(f"CSV file 'cilia_coordinates.csv' saved with {num_points} points along the smoothed spline")
print(f"Tip: Increase smoothing_factor for more smoothing, decrease for less smoothing")