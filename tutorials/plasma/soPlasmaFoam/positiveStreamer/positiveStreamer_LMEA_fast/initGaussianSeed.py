import numpy as np
import fluidfoam
import os
import re

# --- Case 1 Parameters from Paper ---
N0 = 5e18        # Peak seed density 
sigma = 0.0004   # 0.4 mm width 
x0 = 0.0   
y0 = 0.01      # 1.0 cm height 
bg = 7.9e12      # N2+ share of the 1e13 background (0.79); O2+ carries 2.1e12

sol_path = './'
# POSITIVE IONS ONLY. The Gaussian is a seed of positive ions on a uniform
# n_e = n_i = 1e13 background; the resulting net space charge is what enhances
# the field locally and initiates the streamer. That is the benchmark
# definition (Bagheri et al., Plasma Sources Sci. Technol. 27 (2018) 095002,
# eq. for n_i; and Pasolari & Kourtzanidis, arXiv:2607.05137 sec. 6.1) -- a
# neutral seed produces no space-charge enhancement and a different, later
# initiation.
field_paths = [os.path.join(sol_path, '0', 'n_N2p')]

# 1. Read Mesh for coordinates
# Returns x, y, z arrays for cell centers
x, y, z = fluidfoam.readmesh(sol_path)

# 2. Calculate Gaussian values
r2 = (x-x0)**2 + z**2
y_dist2 = (y - y0)**2
n_values = bg + N0 * np.exp(-(r2 + y_dist2) / (sigma**2))

# 3. Read original file
for field_path in field_paths:
  with open(field_path, 'r') as f:
    lines = f.readlines()

# 4. Filter and Reconstruct the file
  new_lines = []
  skip = False

  for line in lines:
      # Identify the internalField line (uniform or nonuniform)
      if line.strip().startswith("internalField"):
          new_lines.append(f"internalField   nonuniform List<scalar> {len(n_values)}\n(\n")
          # Add all calculated values
          for val in n_values:
              new_lines.append(f"{val}\n")
          new_lines.append(");\n")
          skip = True # Start skipping the old internalField data
          continue
    
      # Identify the start of the boundaryField to stop skipping
      if "boundaryField" in line:
          skip = False
    
      if not skip:
          # Avoid double-adding the internalField if the original was multi-line
          if line.strip() != ");" or "boundaryField" not in "".join(new_lines[-2:]):
              new_lines.append(line)

  new_lines = [l.replace('$internalField', f'uniform {bg}')
               for l in new_lines]

  # 5. Write back
  with open(field_path, 'w') as f:
      f.writelines(new_lines)

  print(f"Applied Gaussian seed + {bg} background to {field_path}")
