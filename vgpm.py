import netCDF4 as nc4
import numpy as np
from PIL import Image

def opp_befa(chl, irr, sst, dayL):
    if chl < 1.0:
        chl_tot = 38.0 * np.power(chl, 0.425)
    else:
        chl_tot = 40.2 * np.power(chl, 0.507)

    z_eu = 200.0 * np.power(chl_tot, (-0.293))

    if z_eu <= 102.0:
        z_eu = 568.2 * np.power(chl_tot, -0.746)

    if sst < -10.0:
        pb_opt = 0.0
    elif sst < -1.0:
        pb_opt = 1.13
    elif sst > 28.5:
        pb_opt = 4.0
    else:
        pb_opt = 1.2956 + 2.749e-1*sst + 6.17e-2*np.power(sst, 2) - \
            2.05e-2*np.power(sst, 3) + 2.462e-3*np.power(sst, 4) - \
            1.348e-4*np.power(sst, 5) + 3.4132e-6*np.power(sst, 6) - \
            3.27e-8*np.power(sst, 7)

    irrFunc = 0.66125 * irr / (irr + 4.1)

    npp = pb_opt * chl * dayL * irrFunc * z_eu

    return npp

def cal_dayL(lat, yDay):
    gamma = lat/180.0 * np.pi
    psi = yDay/365.0 * 2.0 * np.pi
    solarDec = (0.39637 - 22.9133*np.cos(psi) + 4.02543*np.sin(psi) - \
                0.38720*np.cos(2*psi) + 0.05200*np.sin(2*psi)) * np.pi/180.0
    r = -np.tan(gamma) * np.tan(solarDec)

    if r<=-1:
        return 24.0
    elif np.fabs(r)<1:
        return 24.0 * np.arccos(r) / np.pi
    else:
        return 0


if __name__=="__main__":

    filename_chl = "./input_files/A20161532016182.L3m_MO_CHL_chlor_a_9km.nc"
    filename_sst = "./input_files/AQUA_MODIS.20160601_20160630.L3m.MO.SST.sst.9km.nc"
    filename_par = "./input_files/A20161532016182.L3m_MO_PAR_par_9km.nc"

    nc_chl = nc4.Dataset(filename_chl, 'r')
    nc_sst = nc4.Dataset(filename_sst, 'r')
    nc_par = nc4.Dataset(filename_par, 'r')

    dim_dic = nc_chl.dimensions
    nlat = dim_dic['lat'].size
    nlon = dim_dic['lon'].size

    chl = np.array(nc_chl.variables['chlor_a'][:])
    par = np.array(nc_par.variables['par'][:])
    sst = np.array(nc_sst.variables['sst'][:])

    chl = np.where(chl == -32767.0, np.nan, chl)
    par = np.where(par == -32767.0, np.nan, par)
    sst = np.where(sst == -32767.0, np.nan, sst)

    nrows, ncols = chl.shape
    latarray = np.array(nc_chl.variables['lat'][:])
    lonarray = np.array(nc_chl.variables['lon'][:])
    npparray = np.full((nrows, ncols), np.nan)

    nc_chl.close()
    nc_sst.close()
    nc_par.close()

    ndays = 30
    dayL_ave_max = -1000.0
    for i in range(nrows):
        dayL_sum = 0.0
        lat = latarray[i]
        for yDay in range(153, 183):
            dayL = cal_dayL(lat, yDay)
            dayL_sum += dayL

        dayL_ave = dayL_sum / ndays
        if dayL_ave > dayL_ave_max:
            dayL_ave_max = dayL_ave

        for j in range(ncols):
            npparray[i, j] = opp_befa(chl[i, j], par[i, j], sst[i, j], dayL_ave)

        print("complete one row")

    print("calculation done")


    Image.fromarray(npparray).save('./output_files/A201606npp_vgpm.tif')
    f = nc4.Dataset('./output_files/A201606npp_vgpm.nc', 'w', format='NETCDF4')

    f.createDimension('lat', nlat)
    f.createDimension('lon', nlon)

    npp_nc = f.createVariable('npp_vgpm', 'f4', ('lat', 'lon'))
    lat_nc = f.createVariable('lat', 'f4', 'lat')
    lon_nc = f.createVariable('lon', 'f4', 'lon')

    lat_nc[:] = latarray
    lon_nc[:] = lonarray
    npp_nc[:, :] = npparray
    f.close()

    print("done")











