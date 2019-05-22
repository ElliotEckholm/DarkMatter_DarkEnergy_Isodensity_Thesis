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
import matplotlib.colorbar as cb
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


redshift = 1.0

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

# print("----lambda_a------")
# print(lambda_a)
# print("----------")


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
# ds = yt.load_uniform_grid(data, simulation_data.shape, length_unit=(voxel, "Mpc"), mass_unit=(1.0,"Msun"),bbox=bbox, nprocs=box_size)
ds = yt.load_uniform_grid(data, simulation_data.shape, length_unit="Mpc", bbox=bbox, nprocs=box_size)
# ad = ds.all_data()

def matter_density_field(field, data):
    matter_a = (matter_0 / (scale_factor**3)) / ( (matter_0 / (scale_factor**3)) + lambda_0)
    # print("----matter_a------")
    # print(matter_a)
    # print(np.average(data['matter_density'] ))
    # print("----------")

    matter_a = data['matter_density'] * (matter_0 / lambda_0 )

    # return (data['matter_density']    / lambda_a).in_units("g/cm**3")
    return ( matter_a    / lambda_a).in_units("g/cm**3")



def dark_energy_density_field(field, data):
    yt_lambda_a = yt.YTArray(lambda_a,"g/cm**3")
    # print("-----yt lambda_a-----")
    # print(yt_lambda_a)
    # print("----------")
    # yt_lambda_a = YTArray([lambda_a]).in_units("g/cm**3")
    return ((yt_lambda_a[0]).in_units("g/cm**3"))


# print("----matter_0 / lambda_0 ------")
# print(matter_0 / lambda_0)
# print("----------")


ds.add_field(('gas','matter_density_field'), function=matter_density_field, units="g/cm**3")

ds.add_field(('gas','dark_energy_density_field'), function=dark_energy_density_field, units="g/cm**3")



L = [1,1,0] # vector normal to cutting plane



for index in range(len(density_ratios)):

    # fig = plt.figure()
    # grid = AxesGrid(fig, (0.075,0.075,0.85,0.85,0.95),
    #                 nrows_ncols = (2, 2),
    #                 axes_pad = 1.0,
    #                 label_mode = "1",
    #                 share_all = True,
    #                 cbar_location="right",
    #                 cbar_mode="each",
    #                 cbar_size="3%",
    #                 cbar_pad="10%")

    slc = yt.SlicePlot(ds,'z',  ['matter_density_field'],origin="native",center=[0,0,0],width=(369,'Mpc'))



    #
    #
    # plot = slc.plots["matter_density_field"]
    # # plot.figure = fig
    # # plot.axes = grid[index].axes
    # # plot.cax = grid.cbar_axes[index]
    # colorbar = plot.cb
    # #
    # # # colorbar.set_ticks([density_ratios[0]*density_conv,density_ratios[1]*density_conv,density_ratios[2]*density_conv,density_ratios[3]*density_conv,density_ratios[4]*density_conv])
    # # # colorbar.set_ticks([density_ratios[index]*density_conv])
    # # # colorbar.set_ticklabels(['$420$'])
    # #
    # colorbar.set_ticks([1e-28])
    # colorbar.set_ticklabels(['$10^{-28}$'])  #

    # slc.set_zlim('matter_density_field', density_ratios[index]*density_conv, (density_ratios[index]+0.1)*density_conv) #4.1*density_conv)
    # slc.set_cmap('matter_density_field', "Blues")
    # slc.hide_colorbar()

    # This is the plot object.




    #
    # cbar.set_ticks([5e-29])
    # cbar.set_ticklabels(['$420$'])
    # slc.annotate_grids(cmap=None)
    # slc6.set_unit('z', 'Mpc/h')


    if (smoothing == 0.5):
        if (redshift == 0):
            slc.annotate_title('No Cuts z = 0')
            # slc._setup_plots()
            slc.save('z0/'+filepath+str(density_ratios[index])+"z=0")
        elif (redshift == 0.5):
            slc.annotate_title('No Cuts z = 0.5')
            # slc._setup_plots()
            slc.save('z5/'+filepath+str(density_ratios[index])+"z=05")
        elif (redshift == 1):
            slc.annotate_title('No Cuts z = 1')
            # slc._setup_plots()
            slc.save('z1/'+filepath+str(density_ratios[index])+"z=1")


#
# slc = yt.SlicePlot(ds,'z',  ['matter_density_field'],origin="native",center=[0,0,0],width=(369,'Mpc'))
#
# cm = yt.visualization.color_maps.make_colormap([('red', 10), ('orange', 10), ('yellow', 10), ('green', 10), ('blue', 10)], name="steps",interpolate=False)
# slc.set_zlim('matter_density_field', density_ratios[0]*density_conv, (density_ratios[4])*density_conv) #4.1*density_conv)
# slc.set_cmap('matter_density_field', cm)
#
# slc.hide_colorbar()
# # slc.annotate_grids(cmap=None)
# # slc6.set_unit('z', 'Mpc/h')
# if (smoothing == 0.5):
#     if (redshift == 0):
#         slc.annotate_title('DM/DE = '+'multicolor'+' z = 0')
#         slc.save('z0_de/'+filepath+'multicolor'+"z=0"+'custom')
#     elif (redshift == 0.5):
#         slc.annotate_title('DM/DE = '+'multicolor'+' z = 0.5')
#         slc.save('z5_de/'+filepath+'multicolor'+"z=05"+'custom')
#     elif (redshift == 1):
#         slc.annotate_title('DM/DE = '+str(density_ratios[index])+' z = 1')
#         slc.save('z1_de/'+filepath+'multicolor'+"z=1"+'custom')
#






# slc = yt.SlicePlot(ds,'z',  ['matter_density_field'],origin="native",center=[0,0,0],width=(369,'Mpc'))
#
#
# slc.set_zlim('matter_density_field', density_ratios[0]*density_conv, (density_ratios[4])*density_conv) #4.1*density_conv)
# slc.set_cmap('matter_density_field', "Rainbow")
#
# slc.hide_colorbar()
# slc.annotate_grids(cmap=None)
# slc.annotate_contour("matter_density_field")
# # slc6.set_unit('z', 'Mpc/h')
# if (smoothing == 0.5):
#     if (redshift == 0):
#         slc.annotate_title('DM/DE = '+'multicolor'+' z = 0')
#         slc.save('z0_de/'+filepath+'multicolor'+"z=0"+'rainbow')
#     elif (redshift == 0.5):
#         slc.annotate_title('DM/DE = '+'multicolor'+' z = 0.5')
#         slc.save('z5_de/'+filepath+'multicolor'+"z=05"+'rainbow')
#     elif (redshift == 1):
#         slc.annotate_title('DM/DE = '+str(density_ratios[index])+' z = 1')
#         slc.save('z1_de/'+filepath+'multicolor'+"z=1"+'rainbow')
#
