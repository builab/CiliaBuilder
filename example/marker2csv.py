import xml.etree.ElementTree as ET
import numpy as np
import csv
from scipy.interpolate import CubicSpline

input_cmm = "primarycilia.cmm"
output_csv = "primarycilia_template.csv"

Z_OFFSET = 1620.0
Z_TARGET_MAX = 10000.0
Z_STEP = 20.0

# --------------------------------------------------
# Parse ChimeraX marker set
# --------------------------------------------------
tree = ET.parse(input_cmm)
root = tree.getroot()

lines = {}
all_z = []

for marker in root.findall("marker"):
    mid = int(marker.attrib["id"])
    x = float(marker.attrib["x"])
    y = float(marker.attrib["y"])
    z = float(marker.attrib["z"]) - Z_OFFSET  # subtract first

    line_number = mid // 100  # 100–199 → 1, etc.

    lines.setdefault(line_number, []).append((x, y, z))
    all_z.append(z)

# --------------------------------------------------
# Global Z rescaling
# --------------------------------------------------
all_z = np.array(all_z)
z_max = all_z.max()

scale = Z_TARGET_MAX / z_max

# Z sampling grid (shared by all lines)
z_sample = np.arange(0, Z_TARGET_MAX + Z_STEP, Z_STEP)

# --------------------------------------------------
# Process each line
# --------------------------------------------------
rows_out = []

for line_number, pts in sorted(lines.items()):
    pts = np.array(pts)

    # Sort by Z
    pts = pts[np.argsort(pts[:, 2])]

    x = pts[:, 0]
    y = pts[:, 1]
    z = pts[:, 2] * scale

    # Use Z as spline parameter (monotonic)
    t = z

    # Fit smooth splines
    sx = CubicSpline(t, x)
    sy = CubicSpline(t, y)
    sz = CubicSpline(t, z)

    # Resample every 20 units
    x_new = sx(z_sample)
    y_new = sy(z_sample)
    z_new = sz(z_sample)

    for xi, yi, zi in zip(x_new, y_new, z_new):
        rows_out.append([
            line_number,
            xi,
            yi,
            zi,
            1, 0, 0, -70, 70
        ])

# --------------------------------------------------
# Write CSV
# --------------------------------------------------
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["DoubletNumber", "X", "Y", "Z", "Idx_A", "Idx_B", "Angle", "A_Shift", "B_Shift"])
    writer.writerows(rows_out)

print(f"Wrote {output_csv}")

