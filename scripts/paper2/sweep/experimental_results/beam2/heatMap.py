import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

fileName = str(sys.argv[1])

img_array = np.load(fileName)

vmin_ = np.min(img_array[:])
vmax_ = np.max(img_array[:])

# Print for debugging
print("Maximum value in img_array:", vmax_)
print("Minimum value in img_array:", vmin_)
# Filter out values above 30% of the maximum
filtered_img_array = np.where(img_array > vmax_ * 0.3, vmax_ * 0.3, img_array)

# Flip the image vertically
flipped_img_array = np.rot90(filtered_img_array, 2)


# Calculate max and min y values
max_y = np.max(len(flipped_img_array[:]))
min_y = 0


xStrech=4
max_x = np.max(len(flipped_img_array[0,:]))*xStrech
min_x = 0

print("Maximum Y in img_array:", max_y)
print("Minimum Y in img_array:", min_y)

# Desired tick values and number of ticks
desired_maxY_tick = 50
desired_maxX_tick = 1.8*2+0.65

num_ticks = 5  # Number of ticks you want to generate

# Calculate tick positions based on max and min y values
tick_yPositions = np.linspace(max_y-1, min_y, num_ticks)
yTicks=np.linspace(desired_maxY_tick, 0, num_ticks)
xTicks=[0,1.8,1.8+0.65,1.8*2+0.65]
#tick_xPositions = [0,1.8,1.8+0.65,1.8*2+0.65]*(max_x/(1.8*2+0.65))
tick_xPositions=[(max_x/(1.8*2+0.65))*n for n in xTicks]


# Define extent based on max and min y values
extent = [min_x, max_x, min_y, max_y]


# Desired number of intervals for the colorbar and image
num_intervals = 5
boundaries = np.linspace(0, 50, num_intervals + 1)  # Define boundaries for intervals

# Create a BoundaryNorm instance for both image and colorbar
norm = colors.BoundaryNorm(boundaries, plt.cm.viridis.N)



# Plotting
plt.imshow(flipped_img_array, cmap='plasma', extent=extent,norm=norm,interpolation='none')

# Add colorbar with discrete intervals
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='plasma'), boundaries=boundaries, ticks=boundaries)
cbar.set_label('Volumefraction [%]')  # Set label for colorbar

# Custom axis labels
plt.xlabel('Width[mm]')
plt.ylabel('Length[mm]')

plt.yticks(tick_yPositions, yTicks)
plt.xticks(tick_xPositions, xTicks)
plt.axvline(x = 1.8*(max_x/(1.8*2+0.65)), color = 'b', label = 'membrane interface')
plt.axvline(x = (1.8+0.65)*(max_x/(1.8*2+0.65)), color = 'b', label = 'membrane interface')




#plt.colorbar()  # Optionally, add a colorbar

plt.show()

