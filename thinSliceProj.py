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
data = dict(matter_density = (simulation_data*density_conv, "g/cm**3"),particle_position_x = (simulation_data, 'Mpc'),particle_position_y = (simulation_data, 'Mpc'), particle_position_z = (simulation_data, 'Mpc'))
bbox = np.array([[0, sim_dim], [0, sim_dim], [0, sim_dim]])

#load array into yt dataset
ds = yt.load_uniform_grid(data, simulation_data.shape, length_unit="Mpc", bbox=bbox, nprocs=box_size)

def matter_density_field(field, data):
    matter_a = (matter_0 / (scale_factor**3)) / ( (matter_0 / (scale_factor**3)) + lambda_0)
    matter_a = data['matter_density'] * (matter_0 / lambda_0 )
    return ( matter_a    / lambda_a).in_units("g/cm**3")



def dark_energy_density_field(field, data):
    yt_lambda_a = yt.YTArray(lambda_a,"g/cm**3")
    return ((yt_lambda_a[0]).in_units("g/cm**3"))


ds.add_field(('gas','matter_density_field'), function=matter_density_field, units="g/cm**3")
ds.add_field(('gas','dark_energy_density_field'), function=dark_energy_density_field, units="g/cm**3")
slc = yt.SlicePlot(ds,'z',  ['matter_density_field'],origin="native",center=[0,0,0],width=(369,'Mpc'))


slc.set_zlim('matter_density_field', density_ratios[0]*density_conv, (density_ratios[4])*density_conv) #4.1*density_conv)
slc.set_cmap('matter_density_field', "Rainbow")

slc.hide_colorbar()
slc.annotate_grids(cmap=None)
slc.annotate_contour("matter_density_field")


# slc6.set_unit('z', 'Mpc/h')
if (smoothing == 0.5):
    if (redshift == 0):
        slc.annotate_title('DM/DE = '+'multicolor'+' z = 0')
        slc.save('z0_de/'+filepath+'multicolor'+"z=0"+'rainbow')
    elif (redshift == 0.5):
        slc.annotate_title('DM/DE = '+'multicolor'+' z = 0.5')
        slc.save('z5_de/'+filepath+'multicolor'+"z=05"+'rainbow')
    elif (redshift == 1):
        slc.annotate_title('DM/DE = '+str(density_ratios[index])+' z = 1')
        slc.save('z1_de/'+filepath+'multicolor'+"z=1"+'rainbow')
