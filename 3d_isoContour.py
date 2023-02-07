import yt
from yt.units import kpc
from yt.units import dimensions
from yt import YTArray
import sys
from PIL import Image
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
from matplotlib import rc_context
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt



redshift = 0

#smoothing = 0.5 Mpc
filepath_z_0_s_0 = 'bpData/BP_0214_densities_1024_0.bin'
filepath_z_5_s_0 = 'bpData/BP_0170_densities_1024_0.bin'
filepath_z_1_s_0 = 'bpData/BP_0136_densities_1024_0.bin'
filepath = ''

#Bolshio-P Parameters
box_size = 1024
voxel = 0.3603515625 #Mpc
hubble_param = 0.6775067751
sim_dim = 369 # Mpc
#Constants
lambda_0 = 0.693
matter_0 = 0.307
density_conv = 10e-29

#Setup parameters
density_ratios = [0.25,0.5,1,2,4]
smoothing = 0.5 #mpc

scale_factor = YTArray([(1 / (1 + redshift ))])

lambda_a = (lambda_0) / ( (matter_0/(scale_factor**3)) + lambda_0)

if (smoothing == 0.5):
    if (redshift == 0):
        filepath = filepath_z_0_s_0
    elif (redshift == 0.5):
        filepath = filepath_z_5_s_0
    elif (redshift == 1):
        filepath = filepath_z_1_s_0



#grab binary data file, convert to 3d array
simulation_data = np.fromfile(filepath,np.float64).reshape((box_size,box_size,box_size))
##First Question, can I convert the raw density values into g/cm**3 by multiplying by density_conv variable?
data = dict(density = (simulation_data*density_conv, "g/cm**3"),particle_position_x = (simulation_data, 'Mpc'),particle_position_y = (simulation_data, 'Mpc'), particle_position_z = (simulation_data, 'Mpc'))
bbox = np.array([[0, sim_dim], [0, sim_dim], [0, sim_dim]])

#load array into yt dataset
ds = yt.load_uniform_grid(data, simulation_data.shape, length_unit="Mpc", bbox=bbox, nprocs=box_size)

def matter_density_field(field, data):
    matter_a = (matter_0 / (scale_factor**3)) / ( (matter_0 / (scale_factor**3)) + lambda_0)
    matter_a = data['density'] * (matter_0 / lambda_0 )
    return ( matter_a    / lambda_a).in_units("g/cm**3")



def dark_energy_density_field(field, data):
    yt_lambda_a = yt.YTArray(lambda_a,"g/cm**3")
    return ((yt_lambda_a[0]).in_units("g/cm**3"))



ds.add_field(('gas','matter_density_field'), function=matter_density_field, units="g/cm**3")
ds.add_field(('gas','dark_energy_density_field'), function=dark_energy_density_field, units="g/cm**3")


# Create a sphere object centered on the highest density point in the simulation
# with radius 1 Mpc
sphere = ds.sphere("max", (5.0, "Mpc"))

# Identify the isodensity surface in this sphere with density = 1e-24 g/cm^3
surface = ds.surface(sphere, "matter_density_field", 3.12626e-25)

# Color this isodensity surface according to the log of the temperature field
# colors = yt.apply_colormap((surface["matter_density_field"]), cmap_name="hot")

# Create a 3D matplotlib figure for visualizing the surface
fig = plt.figure()
ax = fig.gca(projection='3d')
p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

# Set the surface colors in the right scaling [0,1]
p3dc.set_facecolors([0,1,0])
ax.add_collection(p3dc)

# Let's keep the axis ratio fixed in all directions by taking the maximum
# extent in one dimension and make it the bounds in all dimensions
# max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
# centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
# bounds = np.zeros([3,2])
# bounds[:,0] = centers[:] - max_extent/2
# bounds[:,1] = centers[:] + max_extent/2
# ax.auto_scale_xyz(bounds[0,:], bounds[1,:], bounds[2,:])

# Save the figure
plt.savefig("3D/%s_Surface.png" % ds)
