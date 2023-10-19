import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

colormap = LinearSegmentedColormap.from_list('white_gist_earth', [
    (0,     (1,        1,        1       )),
    (1e-20, (0.965882, 0.915975, 0.913378)),
    (0.2,   (0.772885, 0.646409, 0.444171)),
    (0.4,   (0.568932, 0.677541, 0.340330)),
    (0.6,   (0.249216, 0.576471, 0.342046)),
    (0.8,   (0.143740, 0.396564, 0.488306)),
    (1,     (0.013067, 0.000000, 0.348089)),
    ], N=256)



# Define the range of values for the colorbar
min_value = 0
max_value = 1

# Create a figure and axes
fig, ax = plt.subplots()

# Create a ScalarMappable object that maps the colormap to the range of values
sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=min_value, vmax=max_value))
sm.set_array([])  # You need to set an empty array-like object

# Create the colorbar using the ScalarMappable object
cbar = plt.colorbar(sm, ax=ax, orientation='vertical')

# Add a label to the colorbar
cbar.set_label('Colorbar Label')

# Show the plot
plt.show()

