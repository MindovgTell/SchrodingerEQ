import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the initial CSV file
z_data_list = []
t_points = np.arange(0, 300)

# Load multiple CSV files into z_data_list
for i in range(0, 300):  # Adjust the range as needed
    filename = f"../Matrice/matrix{i}.csv"  # Use f-string to create filenames like matrix0.csv, matrix1.csv, etc.
    data = np.genfromtxt(filename, delimiter=',', dtype=float)  # Use dtype=float for consistency

    # Pad data to ensure it is 125,125
    padded_data = np.pad(data, ((0, max(0, 125 - data.shape[0])), (0, max(0, 200 - data.shape[1]))), 'constant')
    z_data_list.append(padded_data[:125, :125])  # Ensure it's exactly 100x100

# Load potential V1
V1 = np.array(np.genfromtxt('../Matrice/potential.csv', delimiter=',', dtype=float))  # Ensure dtype=float

# Ensure V1 is reshaped or cropped to (100, 100) as well
V1 = V1[:125, :125]

# Define x and y points based on the 100x100 dimensions
x_points = np.arange(125)  # 100 columns
y_points = np.arange(125)  # 100 row

fontsize = 12
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create the figure and axis
fig, ax = plt.subplots()

# Mask V1 where values are greater than 10
V1 = np.ma.masked_where(V1 > 10.00, V1)

# Set up the colormap normalization based on the first frame
norm = matplotlib.cm.colors.Normalize(vmin=np.min(z_data_list[0]), vmax=np.max(z_data_list[0]))

# Plot the static V1 background
pb = ax.imshow(V1, extent=[x_min, x_max, y_min, y_max], cmap='Greys_r', norm=norm)

# Plot the dynamic z_data_list[0]
img = ax.imshow(z_data_list[0], extent=[x_min, x_max, y_min, y_max], cmap='viridis', norm=norm)

# Add labels and colorbar
plt.xlabel("x [Units of distance /1]", fontsize=fontsize)
plt.ylabel("y [Units of distance /1]", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colorbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("Data Value", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Animation function
def animation(i):
    # Normalize the color scale to the current frame
    norm = matplotlib.cm.colors.Normalize(vmin=np.min(z_data_list[i]), vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data (Sums with V1 to overlay boundary conditions)
    img.set_data(np.add(z_data_list[i], V1))

    return img,

# Create the animation object
anim = FuncAnimation(fig, animation, interval=10, frames=len(z_data_list), repeat=False, blit=False)

# Save the animation as a GIF (make sure ffmpeg is installed and configured)
anim.save('./animation.gif', writer="ffmpeg", fps=30)

# Show the animation (optional)
plt.show()
