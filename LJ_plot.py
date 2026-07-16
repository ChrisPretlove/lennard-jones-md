import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

# 1. Load particle trajectories
# Expects a CSV containing spatial coordinates of all particles over time
df = pd.read_csv("LJ_data.csv")

# Defining the exact enviroment dimensions
BOX_WIDTH = 5.0e-9  
BOX_HEIGHT = 5.0e-9

# 2. Setup figure canvas
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal') # Maintain 1:1 spatial aspect ratio
ax.set_xlim(0, BOX_WIDTH) # Box boundaries as definined in the LJ_Potential.c script
ax.set_ylim(0, BOX_HEIGHT)
ax.set_xlabel("X Position (meters)")
ax.set_ylabel("Y Position (meters)")
ax.set_title("2D Lennard-Jones potential simulation of Argon atoms")

# Make the axes ticks show in nanometres instead of metres
ax.set_xticklabels([f"{x*1e9:.1f}" for x in ax.get_xticks()])
ax.set_yticklabels([f"{y*1e9:.1f}" for y in ax.get_yticks()])
ax.set_xlabel("X (nm)")
ax.set_ylabel("Y (nm)")

# 3. Create the particle instance
# This placeholder scatter plot will be modified dynamically frame-by-frame
scat = ax.scatter([], [], s=50, color="black")
grouped = {time: group[['x', 'y']].values for time, group in df.groupby('time')}
timestamps = sorted(list(grouped.keys()))[::10]

# 4. Frame initialization function (Required by Matplotlib FuncAnimation)
def init():
    scat.set_offsets(np.empty((0, 2)))
    return [scat]

# 5. Update function called on every frame of the simulation
def update(frame):
    # Filter the trajectory dataframe to only grab coordinates for the current timestep
    current_time = timestamps[frame]

    # Extract coordinate positions into an (N, 2) array for scatter updating
    coordinates = grouped[current_time]
    scat.set_offsets(coordinates)
    return [scat]

# 6. Generate the physics animation 
# Determines unique timeframes dynamically, pacing the rendering with a 20ms interval
ani = animation.FuncAnimation(fig, update, frames=len(timestamps), init_func=init, blit=True, interval=10)

# 7. Save as a high-quality GIF to embed directly in the GitHub README
print("Saving animation... (this might take a minute)")
ani.save("lj_simulation.gif", writer="pillow", fps=30)
print("Animation successfully saved")
plt.show()