import xarray as xr
import numpy as np
import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt

PgCperppm = 2.123

#--------------------------
# Mixed-Layer Response Functions
#--------------------------

def HILDA_response(years, OceanMLDepth): 
    """
    Mixed layer response function of the HILDA box-diffusion model 
    (Shaffer and Sarmiento 1992, Sigenthaler and Joos 1996)
    
    Parameters
    ----------
    years : `numpy.array` or `xarray.DataArray` 
        Points in time to at which the user wants calculate ocean flux.
    OceanMLDepth : `int` or `float`
        Ocean mixed layer depth, used as a tuning parameter.
        
    Returns
    -------
    return_val : `numpy.array`
        An array of values of the response function at any point in time, 
        times a conversion factor
    """
    

    ocean_area = 3.62E14  # ocean area in square meters
    g_cper_mole = 12.0113  # molar mass of carbon.
    sea_water_dens = 1.0265E3  # sea water density in kg/m^3.
    # c = 1.722E17 umol m^3 kg^-1 ppm^-1; Joos et al. 1996 pg. 402
    c = (1E21 * PgCperppm / g_cper_mole) / (sea_water_dens) 

    years = years-np.min(years)
    return_val = np.zeros(len(years))
    
    for yr_ind,yr in enumerate(years):
        if yr < 2.0:
            value = 0.12935 + 0.21898 * np.exp(-yr / 0.034569) + 0.17003 * np.exp(-yr / 0.26936) + 0.24071 * np.exp(
                -yr / 0.96083) + 0.24093 * np.exp(-yr / 4.9792)
        else:
            value = 0.022936 + 0.24278 * np.exp(-yr / 1.2679) + 0.13963 * np.exp(-yr / 5.2528) + 0.089318 * np.exp(
                -yr / 18.601) + 0.037820 * np.exp(-yr / 68.736) + 0.035549 * np.exp(-yr / 232.30)

        # scale values to umol kg^-1 ppm^-1
        return_val[yr_ind] = value * c / (OceanMLDepth * ocean_area)


    return return_val


def Sarmiento92_response(years, OceanMLDepth):
    """
    Mixed layer response function (GFDL model, Sarmiento et al. 1992) as described in Joos et al., 1996.
    
    Parameters
    ----------
    years : `numpy.array` or `xarray.DataArray` 
        Points in time to at which the user wants calculate ocean flux.
    OceanMLDepth : `int` or `float`
        Ocean mixed layer depth, used as a tuning parameter.
        
    Returns
    -------
    return_val : `numpy.array`
        An array of values of the response function at any point in time, 
        times a conversion factor
    """

    ocean_area = 3.62E14  # ocean area in square meters
    g_cper_mole = 12.0113  # molar mass of carbon.
    sea_water_dens = 1.0265E3  # sea water density in kg/m^3.
    # c = 1.722E17 umol m^3 kg^-1 ppm^-1; Joos et al. 1996 pg. 402
    c = (1E21 * PgCperppm / g_cper_mole) / (sea_water_dens) 

    years = years-np.min(years)
    return_val = np.zeros(len(years))
    
    for yr_ind,yr in enumerate(years):
        if yr <= 1.0:
            value = 1
        else:
            value = (0.014819 + 0.70367*np.exp(-yr/0.70177)
                     + 0.24966*np.exp(-yr/2.3488) + 0.066485
                     *np.exp(-yr/15.281) + 0.038344
                     *np.exp(-yr/65.359) + 0.019439
                     *np.exp(-yr/347.55))
            
        

        # scale values to umol kg^-1 ppm^-1
        return_val[yr_ind] = value * c / (OceanMLDepth * ocean_area)

    return return_val

#--------------------------
# Chemistry
#--------------------------

def delta_pco2_ocean(surface_ocean_ddic):
    """
    Perturbation pCO2 calculated from parameterized carbonate chemistry
    as described in Joos et al., 1996.
    
    Parameters
    ----------
    surface_ocean_ddic : `float`
        Anthropogenic carbon in the surface ocean (umol/kg)
        
    Returns
    -------
    return_val : `float`
        perturbation pCO2 in units of ppm
    """
    
    TC = 18.1716  # Effective Ocean temperature for carbonate chemistry in deg C.
    A1 = (1.5568 - 1.3993E-2 * TC)
    A2 = (7.4706 - 0.20207 * TC) * 1E-3
    A3 = -(1.2748 - 0.12015 * TC) * 1E-5
    A4 = (2.4491 - 0.12639 * TC) * 1E-7
    A5 = -(1.5468 - 0.15326 * TC) * 1E-10
    # from Joos et al. 1996, pg. 402
    return_val = surface_ocean_ddic * (
            A1 + surface_ocean_ddic * (A2 + surface_ocean_ddic * (A3 + surface_ocean_ddic * (A4 + surface_ocean_ddic * A5))))

    return return_val

def delta_pco2_ocean_lin(surface_ocean_ddic):
    """
    Perturbation pCO2 based on constant buffer factor. The buffer factor is calculated
    as the slope of delta_pco2_ocean at an ocean perturbation pCO2 of 0.001 ppm and
    anthropogenic carbon concentration (dDIC) of 0.001
    
    Parameters
    ----------
    surface_ocean_ddic : `float`
        Anthropogenic carbon in the surface ocean (umol/kg)
        
    Returns
    -------
    return_val : `float`
        perturbation pCO2 in units of ppm
    """
    
    ddic = 0.001
    dpco2 = delta_pco2_ocean(0.001) - delta_pco2_ocean(0)
    
    R = dpco2/ddic
    
    return_val = (R*surface_ocean_ddic)
        

    return return_val

#--------------------------
# Flux Calculations (and idealized flux calculations)
#--------------------------

def flux_TC_CV(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
               surface_ocean_ddic,dpco2_oc): #ocean_flux
    
    """
    air-sea flux at constant temperature but variable carbonate chemistry
    
    Parameters
    ----------
    yr_id : `int`
        Array index corresponding to the time value of the current timestep
    datmos_co2 : `float`
        Atmospheric perturbation pCO2 (ppm)
    air_sea_gas_exchange_coeff: `float`
        Air-sea gas exchange coefficient (1/year)
    surface_ocean_ddic : `float`
        Anthropogenic carbon in the surface ocean (umol/kg)   
    dpco2_oc: `numpy.array`
        Ocean perturbation pCO2 (ppm)
        
    Returns
    -------
    `float`
        air-sea flux in units of ppm/timestep
    """

    if yr_ind > 0:
        dpco2_oc[yr_ind] = delta_pco2_ocean(surface_ocean_ddic[yr_ind])

    return air_sea_gas_exchange_coeff * (datmos_co2[yr_ind] - dpco2_oc[yr_ind])

#--------------------------

def flux_TV_CC(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
               surface_ocean_ddic,dpco2_oc,pco2_oc,pco2_oc_pi,
               atmos_co2,DT): #ocean_flux_dt_lin  
    """
    Air-sea flux at variable temperature but constant carbonate chemistry
    
    Parameters
    ----------
    yr_id : `int`
        Array index corresponding to the time value of the current timestep
    datmos_co2 : `float`
        Atmospheric perturbation pCO2 (ppm)
    air_sea_gas_exchange_coeff: `float`
        Air-sea gas exchange coefficient (1/year)
    surface_ocean_ddic : `float`
        Anthropogenic carbon in the surface ocean (umol/kg)   
    dpco2_oc : `numpy.array`
        Ocean perturbation pCO2 (ppm)
    pco2_oc : `numpy.array`
        Ocean pCO2 (ppm)
    pco2_oc_pi : `int` or `float`
        Preindustrial ocean pCO2 (ppm)
    atmos_co2 :  `numpy.array`
        Atmopsheric pCO2 (ppm)
    DT : `numpy.array` or `xarray.DataArray`
        Warming since preindustrail (degrees celsius)
        
    Returns
    -------
    `float`
        air-sea flux in units of ppm/timestep
    """

    dpco2_oc[yr_ind] = delta_pco2_ocean_lin(surface_ocean_ddic[yr_ind])

    pco2_oc[yr_ind] = (pco2_oc_pi + dpco2_oc[yr_ind])*np.exp(0.0423*DT[yr_ind])

    return air_sea_gas_exchange_coeff * (atmos_co2[yr_ind] - pco2_oc[yr_ind])

#--------------------------

def flux_TC_CC(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
               surface_ocean_ddic,dpco2_oc): #ocean_flux_lin   
    """
    Air-sea flux at constant temperature and constant carbonate chemistry
    
    Parameters
    ----------
    yr_id : `int`
        Array index corresponding to the time value of the current timestep
    datmos_co2 : `float`
        Atmospheric perturbation pCO2 (ppm)
    air_sea_gas_exchange_coeff: `float`
        Air-sea gas exchange coefficient (1/year)
    surface_ocean_ddic : `float`
        Anthropogenic carbon in the surface ocean (umol/kg)   
    dpco2_oc : `numpy.array`
        Ocean perturbation pCO2 (ppm)
        
    Returns
    -------
    `float`
        air-sea flux in units of ppm/timestep 
    """

    if yr_ind > 0:
        dpco2_oc[yr_ind] = delta_pco2_ocean_lin(surface_ocean_ddic[yr_ind])

    return air_sea_gas_exchange_coeff * (datmos_co2[yr_ind] - dpco2_oc[yr_ind])

#--------------------------

def flux_TV_CV(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
               surface_ocean_ddic,dpco2_oc,pco2_oc,pco2_oc_pi,
               atmos_co2,DT): #ocean_flux_dt
    """
    Actual air-sea flux (variable temperature, variable carbonate chemistry)
    
    Parameters
    ----------
    yr_id : `int`
        Array index corresponding to the time value of the current timestep
    datmos_co2 : `float`
        Atmospheric perturbation pCO2 (ppm)
    air_sea_gas_exchange_coeff: `float`
        Air-sea gas exchange coefficient (1/year)
    surface_ocean_ddic : `float`
        Anthropogenic carbon in the surface ocean (umol/kg)   
    dpco2_oc : `numpy.array`
        Ocean perturbation pCO2 (ppm)
    pco2_oc : `numpy.array`
        Ocean pCO2 (ppm)
    pco2_oc_pi : `int` or `float`
        Preindustrial ocean pCO2 (ppm)
    atmos_co2 :  `numpy.array`
        Atmopsheric pCO2 (ppm)
    DT : `numpy.array` or `xarray.DataArray`
        Warming since preindustrail (degrees celsius)
        
    Returns
    -------
    `float`
        air-sea flux in units of ppm/timestep
    """

    dpco2_oc[yr_ind] = delta_pco2_ocean(surface_ocean_ddic[yr_ind])

    pco2_oc[yr_ind] = (pco2_oc_pi + dpco2_oc[yr_ind])*np.exp(0.0423*DT[yr_ind])

    return air_sea_gas_exchange_coeff * (atmos_co2[yr_ind] - pco2_oc[yr_ind])

#--------------------------

def ocean_flux(atmos_co2,
                   OceanMLDepth=51, HILDA=True,
                   DT=None,
                   temperature='constant', chemistry='variable',
                  ):   
    """
    Calculate ocean carbon uptake as in Joos et al. 1996
    
    Parameters
    ----------
    atmos_co2 :  `xarray.DataArray`
        Atmopsheric pCO2 (ppm)
    OceanMLDepth : `int` or `float`
        Ocean mixed layer depth, used as a tuning parameter.
    HILDA : bool
        When set to True, use the HILDA response function. When False, use Sarmiento et al. 1992 response.
        Default value is True
    DT : `numpy.array` or `xarray.DataArray`
        Warming since preindustrail (degrees celsius)
    temperature : str
        'variable' or 'constant', determines whether surface ocean experiences temperature increase
    chemistry: str
        'variable' or 'constant', determines whether carbonate chemistry is linear or variable
        
    Returns
    -------
    ds: `xarray.Dataset`
        This Dataset contains time (year), ocean flux ((Pg C)/yr), surface ocean anthropogenic
        carbon (umol/kg), surface ocean perturbation pCO2 (ppm)
    """
    
    if temperature not in ['constant','variable']:
        raise ValueError('temperature flag must be either \'contstant\' or \'variable\'')
    if chemistry not in ['constant','variable']:
        raise ValueError('chemistry flag must be either \'contstant\' or \'variable\'')

#     air_sea_gas_exchange_coeff = 0.1042  # kg m^-2 year^-1
    ocean_area = 3.62E14
    air_sea_gas_exchange_coeff = 1/9.06 # year^-1
    surface_ocean_ddic = np.zeros(len(atmos_co2))
    dpco2_oc = np.zeros(len(atmos_co2))
    pco2_oc = np.zeros(len(atmos_co2))
    air_sea_flux = np.zeros(len(atmos_co2))
    
    pco2_oc_pi = 284.4
    datmos_co2 = atmos_co2-284.4
    
    
    if HILDA:

        ocean_response = HILDA_response(atmos_co2.year, OceanMLDepth)
    
    else:
        
        ocean_response = Sarmiento92_response(atmos_co2.year, OceanMLDepth)

    
    for yr_ind in range(len(atmos_co2) - 1):

        # if block for various experiments
        if (temperature == 'variable') and (chemistry == 'variable'):
            
            air_sea_flux[yr_ind] = flux_TV_CV(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
                                              surface_ocean_ddic,dpco2_oc,pco2_oc,pco2_oc_pi,atmos_co2,DT)
        
        if (temperature == 'constant') and (chemistry == 'variable'):

            air_sea_flux[yr_ind] = flux_TC_CV(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
                                              surface_ocean_ddic,dpco2_oc)
        
        if (temperature == 'constant') and (chemistry == 'constant'):
            
            air_sea_flux[yr_ind] = flux_TC_CC(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
                                              surface_ocean_ddic,dpco2_oc)
        
        if (temperature == 'variable') and (chemistry == 'constant'):
            
            air_sea_flux[yr_ind] = flux_TV_CC(yr_ind,datmos_co2,air_sea_gas_exchange_coeff,
                                              surface_ocean_ddic,dpco2_oc,pco2_oc,pco2_oc_pi,atmos_co2,DT)

        # Accumulate committments of these fluxes to all future times.
        for j in range(yr_ind + 1, len(surface_ocean_ddic)):
            surface_ocean_ddic[j] = surface_ocean_ddic[j] + air_sea_flux[yr_ind] * ocean_response[j - yr_ind]
            
    if (temperature == 'variable'):
        dpco2_oc = pco2_oc - pco2_oc_pi
        
    ds = xr.Dataset(data_vars={'F_as':('year',(air_sea_flux*PgCperppm)),
                               'dDIC':('year',surface_ocean_ddic),
                               'dpCO2':('year',dpco2_oc)
                               },
                coords={'year':atmos_co2.year}) 

    # annually average/sum output
    year_bins = np.unique(atmos_co2.year.astype(int))

    ds['F_as'] = ds.F_as.groupby_bins(atmos_co2.year, bins=year_bins, right=False).sum(dim='year')
    ds['dDIC'] = ds.dDIC.groupby_bins(atmos_co2.year, bins=year_bins, right=False).mean(dim='year')
    ds['dpCO2'] = ds.dpCO2.groupby_bins(atmos_co2.year, bins=year_bins, right=False).mean(dim='year')
    
    ds = ds.drop('year')
    ds = ds.rename({'year_bins':'year'})
    ds['year'] = year_bins[0:-1]


    return ds

def plot_experiments(atmos_co2, DT, OceanMLDepth=51):

    """
    Run multiple experiments and plot the results. Same as the example on readthedocs
    """

    plt.figure(dpi=300)

    # linear buffering and constant solubility
    ds = qoccm.ocean_flux(atmos_co2,
                          OceanMLDepth=OceanMLDepth, HILDA=True,
                          DT=None,
                          temperature='constant', chemistry='constant',
                         )
    flux = ds.F_as
    plt.plot(atmos_co2.year,flux,label = 'Fixed Temperature and PI Buffer Factor')

    # linear buffering
    ds = qoccm.ocean_flux(atmos_co2,
                          OceanMLDepth=OceanMLDepth, HILDA=True,
                          DT=DT,
                          temperature='variable', chemistry='constant',
                         )
    flux = ds.F_as
    plt.plot(atmos_co2.year,flux,label='Constant PI Buffer Capacity')

    # constant solubility
    ds = qoccm.ocean_flux(atmos_co2,
                          OceanMLDepth=OceanMLDepth, HILDA=True,
                          DT=None,
                          temperature='constant', chemistry='variable',
                         )
    flux = ds.F_as
    plt.plot(atmos_co2.year,flux,label = 'Fixed Temperature',color='tab:green')

    # control
    ds = qoccm.ocean_flux(atmos_co2,
                          OceanMLDepth=OceanMLDepth, HILDA=True,
                          DT=DT,
                          temperature='variable', chemistry='variable',
                         )
    flux = ds.F_as
    plt.plot(atmos_co2.year,flux,label='Control',color='k')

    plt.ylabel('Pg C yr$^{-1}')
    plt.grid()
    plt.xlim(1850.5,2080)
    plt.legend()

    ax = plt.gca()
    
    return(ax)
