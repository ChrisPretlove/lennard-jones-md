import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

df = pd.read_csv("LJ_data.csv")

# Defining the exact box dimensions from your C code
BOX_WIDTH = 5.0e-9  
BOX_HEIGHT = 5.0e-9

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.5)
ax.set_title("Lennard-Jones Potential Simulation")

# --- THE FIX: Explicitly set the limits to match your nanometer box ---
ax.set_xlim(0, BOX_WIDTH)
ax.set_ylim(0, BOX_HEIGHT)

# Optional quality of life: make the axes ticks show in "nm" instead of "0.000000001"
ax.set_xticklabels([f"{x*1e9:.1f}" for x in ax.get_xticks()])
ax.set_yticklabels([f"{y*1e9:.1f}" for y in ax.get_yticks()])
ax.set_xlabel("X (nm)")
ax.set_ylabel("Y (nm)")
# ----------------------------------------------------------------------

scat = ax.scatter([], [], s=50, color="black")
grouped = {time: group[['x', 'y']].values for time, group in df.groupby('time')}

timestamps = sorted(list(grouped.keys()))[::10]

def init():
    scat.set_offsets(np.empty((0, 2)))
    return [scat]

def update(frame):
    current_time = timestamps[frame]
    coordinates = grouped[current_time]
    scat.set_offsets(coordinates)
    return [scat]

ani = animation.FuncAnimation(fig, update, frames=len(timestamps), init_func=init, blit=True, interval=10)

print("Saving animation... (this might take a minute)")
ani.save("lj_simulation.gif", writer="Pillow", fps=30)
print("Animation successfully saved")
plt.show()