import os
import cairosvg
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# Parameters
svg_filename = "triple_dot_layout.svg"
png_filename = "triple_dot_layout.png"
Nx, Ny = 300, 150  # Resolution of the domain

# 1. Define SVG layout string
svg_content = f'''
<svg width="{Nx}" height="{Ny}" xmlns="http://www.w3.org/2000/svg">
  <!-- Reservoirs -->
  <rect id="Reservoir_Left" x="0" y="50" width="30" height="50" fill="#cccccc"/>
  <rect id="Reservoir_Right" x="{Nx-30}" y="50" width="30" height="50" fill="#cccccc"/>

  <!-- Barrier Gates -->
  <rect id="Barrier_Left" x="60" y="40" width="10" height="70" fill="#ff0000"/>
  <rect id="Barrier_Center" x="145" y="40" width="10" height="70" fill="#ff0000"/>
  <rect id="Barrier_Right" x="230" y="40" width="10" height="70" fill="#ff0000"/>

  <!-- Plunger Gates -->
  <rect id="Plunger_Left" x="80" y="30" width="30" height="90" fill="#0000ff"/>
  <rect id="Plunger_Center" x="160" y="30" width="30" height="90" fill="#0000ff"/>
  <rect id="Plunger_Right" x="240" y="30" width="30" height="90" fill="#0000ff"/>
</svg>
'''

# 2. Save the SVG to disk
with open(svg_filename, "w") as f:
    f.write(svg_content)

# 3. Convert SVG to PNG using cairosvg
cairosvg.svg2png(url=svg_filename, write_to=png_filename, output_width=Nx, output_height=Ny)

# 4. Load PNG as RGB array
img = Image.open(png_filename).convert("RGB")
img_array = np.array(img)

# 5. Define color-to-voltage mapping
color_voltage_map = {
    (204, 204, 204): 0.0,   # Reservoirs
    (255, 0, 0): 0.2,       # Barriers
    (0, 0, 255): 0.1        # Plungers
}

# 6. Generate potential mask V_gate
V_gate = np.zeros((Ny, Nx), dtype=float)
for color, voltage in color_voltage_map.items():
    mask = np.all(img_array == color, axis=-1)
    V_gate[mask] = voltage

# 7. Visualize the resulting potential
plt.imshow(V_gate, origin='lower', cmap='viridis', extent=[0, Nx, 0, Ny])
plt.colorbar(label="Gate Voltage (V)")
plt.title("Gate Voltage Map from SVG Layout")
plt.xlabel("x [a.u.]")
plt.ylabel("y [a.u.]")
plt.tight_layout()
plt.show()

# Optional cleanup
# os.remove(svg_filename)
# os.remove(png_filename)
