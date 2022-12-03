import netCDF4 as nc4
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

f0 = "./ref_files/A201606npp_ref.nc"
f1 = "./output_files/A201606npp_vgpm.nc"

f0nc = nc4.Dataset(f0, 'r')
f1nc = nc4.Dataset(f1, 'r')

npp_ref = np.array(f0nc.variables['npp'][:])
npp = np.array(f1nc.variables['npp_vgpm'][:])

f0nc.close()
f1nc.close()

a = np.isnan(npp)
aref = np.isnan(npp_ref)
pixel_plot = np.logical_and(~a, ~aref)
npixel_plot = np.sum(pixel_plot)
nan_both = np.logical_and(a, aref)
n_nan_both = np.sum(nan_both)

n_diff = 2160*4320-npixel_plot-n_nan_both

diff = n_diff/npixel_plot
diff_percent = '{:.2%}'.format(diff)
print(diff_percent)

xx = np.arange(0, 50000, 10000)
plt.scatter(npp, npp_ref, marker='o', s=0.1)
plt.xlabel('npp')
plt.ylabel('npp_ref')
plt.savefig('./output_files/r_comp.jpg')
plt.show()

print("done")
