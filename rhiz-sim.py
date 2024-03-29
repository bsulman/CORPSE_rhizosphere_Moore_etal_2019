import CORPSE
import pandas
from pylab import *

# 5% clay
params={
    'vmaxref':[1000,50,600], # Relative maximum enzymatic decomp rates. Multiply by 0.022 to get actual turnover rates (to take moisture equation into account)
    'Ea':[37e3,54e3,50e3],    # Activation energy
    'kC':[0.01,0.01,0.01],    # Michaelis-Menton parameter
    'gas_diffusion_exp':2.5,  # Determines suppression of decomp at high soil moisture
    'minMicrobeC':1e-3,       #Minimum microbial biomass (fraction of total C)
    'Tmic':0.25,       # Microbial lifetime
    'et':0.5,          # Fraction of turnover not converted to CO2
    'eup':[0.6,0.05,0.6], # Carbon uptake efficiency
    'tProtected':75.0,    # Protected C turnover time (years)
    'protection_rate':[0.25,0.00005,0.75], # Protected carbon formation rate (year-1)
    'Resp_uses_total_C':False
}

dt=1.0/365 # Daily time step
depth=5 #cm

# Read daily temperature data for upland site
## NOTE: Weather data downloaded from SPRUCE data repository: http://dx.doi.org/10.3334/CDIAC/spruce.001
weatherdata=pandas.read_csv('../SPRUCE-data/EM123_Combined_public_data_2010_2015.csv',index_col='datetime',parse_dates=True,na_values=[-9999.0])
soilT=weatherdata['EM1_Hummock10cm']+273.15
soilmoisture=weatherdata['EM1_VW1_TopofHummock']
soilmoisture[soilmoisture<0.0]=0.0
soilmoisture=soilmoisture/soilmoisture.max()

# Soil moisture from observations
theta=soilmoisture.fillna(method='backfill').resample('1D').mean()[365:]
T=soilT.fillna(method='backfill').resample('1D').mean()[365:]

# After spinup, run with no fresh inputs
inputs = array([0.3,0.7,0.0])*0.0e-3 # gC/g soil/year

# Calculate ranges of root density
rootdata=pandas.read_excel('All Site SEM data SULMAN with volume.xlsx')
bulkdensity=1.26 #g/cm3

rootdensity=rootdata['roots per volume']*1e-3 #g/cm3 of soil
rootdensity[rootdensity<0]=nan
rootmass=rootdata['roots']*1e-3 # g roots/g soil

srl=rootdata['SRL (m g-1)'].fillna(rootdata['SRL (m g-1)'].mean()) # m/g root
rootlength = rootmass*srl # m root/g soil
exudationrate=rootlength*100*1e-6*24*365/1000.0 # Phillips et al: 1 ugC/cm root length/hour

do_spinup=False
soc=rootdata['SOC g kg-1']
soc_kg_m2=soc/(bulkdensity*1e3)*100**2*depth*1e-3

# Set up initial cohorts
if do_spinup:
    nbins=1
else:
    nbins=10

rootlength_bins=rootlength.quantile(linspace(0.5/nbins,1.0-0.5/nbins,nbins)) # m/g soil
rootdensity_cm_gsoil=array(rootlength_bins)*100
exudate_ugC_cmroot_hour=0.25
exudate_gC_gsoil_hour=exudate_ugC_cmroot_hour*1e-6*rootdensity_cm_gsoil
exudate_gC_gsoil_year=exudate_gC_gsoil_hour*24*365
exudationrate=exudate_gC_gsoil_year # Phillips et al: 1 ugC/cm root length/hour

# c=CORPSE.soil_carbon_cohort(litterC=[6e-05,0.012,9.6e-05], protectedC=[0.0037,4.2e-05,0.0054], livingMicrobeC=8.4e-05, params=params)
c=CORPSE.soil_carbon_cohort(litterC=[0.0001,0.011,8.7e-05], protectedC=[0.0032,4.1e-05,0.0061], livingMicrobeC=0.00021, params=params)
c2=c.copy()
c2.params.update({'Resp_uses_total_C':True})

cohorts=[]
cohorts2=[]
for ii in range(nbins):
    cohorts.append(c.copy())
    cohorts2.append(c2.copy())

if do_spinup:
    nsteps=365*50
else:
    nsteps=365*4


outputs={'unprotectedC':zeros((nsteps,nbins,3)),
            'protectedC':zeros((nsteps,nbins,3)),
            'decomp':zeros((nsteps,nbins,3)),
            'microbeC':zeros((nsteps,nbins)),
            'CO2':zeros((nsteps,nbins))
            }
outputs2={'unprotectedC':zeros((nsteps,nbins,3)),
            'protectedC':zeros((nsteps,nbins,3)),
            'decomp':zeros((nsteps,nbins,3)),
            'microbeC':zeros((nsteps,nbins)),
            'CO2':zeros((nsteps,nbins))
            }

doy=T.index.dayofyear
exudation_ts=cos((doy-0.67*365)*2*pi/365)+1

for step in range(nsteps):
    if step%365==0:
        print ('Year %d of %d'%(floor(step/365),floor(nsteps/365)))
    for cc in range(len(cohorts)):
        out=cohorts[cc].update(T[step%len(T)],theta[step%len(T)],dt)
        cohorts[cc].check_validity()
        outputs['unprotectedC'][step,cc,:]=cohorts[cc].litterC
        outputs['decomp'][step,cc,:]=out['decomp']
        outputs['protectedC'][step,cc,:]=cohorts[cc].protectedC
        outputs['microbeC'][step,cc]=cohorts[cc].livingMicrobeC
        outputs['CO2'][step,cc]=cohorts[cc].CO2

        cohorts[cc].add_carbon((inputs+array([exudationrate[cc],0.0,0.0]))*dt*exudation_ts[step%len(T)])

        out2=cohorts2[cc].update(T[step%len(T)],theta[step%len(T)],dt)
        cohorts2[cc].check_validity()
        outputs2['unprotectedC'][step,cc,:]=cohorts2[cc].litterC
        outputs2['decomp'][step,cc,:]=out2['decomp']
        outputs2['protectedC'][step,cc,:]=cohorts2[cc].protectedC
        outputs2['microbeC'][step,cc]=cohorts2[cc].livingMicrobeC
        outputs2['CO2'][step,cc]=cohorts2[cc].CO2

        cohorts2[cc].add_carbon((inputs+array([exudationrate[cc],0.0,0.0]))*dt*exudation_ts[step%len(T)])

# Plot results


if do_spinup:
    t=arange(nsteps)/365.0
    figure(1);clf()
    subplot(211)
    plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1),'b-',label='Unprotected')
    plot(t,outputs['protectedC'][:,0,:].sum(axis=1),'r-',label='Protected')
    plot(t,outputs['unprotectedC'][:,0,0],':g',label='Fast')
    plot(t,outputs['unprotectedC'][:,0,1],'b:',label='Slow')
    plot(t,outputs['unprotectedC'][:,0,2],'r:',label='Microbe necro')
    plot(t,outputs['microbeC'][:,0],'m-',label='Microbe')
    plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1)+outputs['protectedC'][:,0,:].sum(axis=1),'k-',label='Total')
    legend(loc='best',fontsize='medium')
    draw()

    subplot(212)
    plot(t,outputs2['unprotectedC'][:,0,:].sum(axis=1),'b-',label='Unprotected')
    plot(t,outputs2['protectedC'][:,0,:].sum(axis=1),'r-',label='Protected')
    plot(t,outputs2['unprotectedC'][:,0,0],':g',label='Fast')
    plot(t,outputs2['unprotectedC'][:,0,1],'b:',label='Slow')
    plot(t,outputs2['unprotectedC'][:,0,2],'r:',label='Microbe necro')
    plot(t,outputs2['microbeC'][:,0],'m-',label='Microbe')
    plot(t,outputs2['unprotectedC'][:,0,:].sum(axis=1)+outputs2['protectedC'][:,0,:].sum(axis=1),'k-',label='Total')

else:
    t=T.index[:nsteps]
    figure(1);clf()
    # subplot(131)
    lowpoint=0
    plot(t,outputs['unprotectedC'][:,lowpoint,:].sum(axis=1)*1e3,'b-',label='Unprotected')
    # plot(t,outputs['protectedC'][:,0,:].sum(axis=1),'r-',label='Protected')
    plot(t,outputs['unprotectedC'][:,lowpoint,0]*1e3,'-g',label='Fast')
    # plot(t,outputs['unprotectedC'][:,0,1],'c:',label='Slow')
    plot(t,outputs['unprotectedC'][:,lowpoint,2]*1e3,'y-',label='Dead microbe')
    plot(t,outputs['microbeC'][:,lowpoint]*1e3,'m-',label='Live microbe')
    # plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1)+outputs['protectedC'][:,0,:].sum(axis=1),'k-',label='Total')
    leg=legend(loc='best',fontsize='medium');leg.get_frame().set_alpha(0.5)

    highpoint=4
    plot(t,outputs['unprotectedC'][:,highpoint,:].sum(axis=1)*1e3,'b--',label='Unprotected')
    # plot(t,outputs['protectedC'][:,-1,:].sum(axis=1),'r--',label='Protected')
    plot(t,outputs['unprotectedC'][:,highpoint,0]*1e3,'--g',label='Fast')
    # plot(t,outputs['unprotectedC'][:,-1,1],'c-.',label='Slow')
    plot(t,outputs['unprotectedC'][:,highpoint,2]*1e3,'y--',label='Microbe necro')
    plot(t,outputs['microbeC'][:,highpoint]*1e3,'m--',label='Microbe')
    # plot(t,outputs['unprotectedC'][:,0-1,:].sum(axis=1)+outputs['protectedC'][:,0-1,:].sum(axis=1),'k--',label='Total')
    xlabel('Time (years)')
    ylabel('Carbon content(mgC/g soil)')
    title('Labile carbon pools')

    plot(t,outputs2['unprotectedC'][:,highpoint,:].sum(axis=1)*1e3,'b:')
    # plot(t,outputs['protectedC'][:,0,:].sum(axis=1),'r-',label='Protected')
    plot(t,outputs2['unprotectedC'][:,highpoint,0]*1e3,':g')
    # plot(t,outputs['unprotectedC'][:,0,1],'c:',label='Slow')
    plot(t,outputs2['unprotectedC'][:,highpoint,2]*1e3,'y:')
    plot(t,outputs2['microbeC'][:,highpoint]*1e3,'m:')

    draw()

    figure(2,figsize=(8,5.3));clf()

    timepoint=nonzero((t.week==25)&(t.year==2015))[0]


    subplot(122)
    semilogx(rootlength_bins*1e3,outputs['microbeC'][timepoint,:].mean(axis=0)*1e6,'ko-',label='June')
    semilogx(rootlength_bins*1e3,outputs2['microbeC'][timepoint,:].mean(axis=0)*1e6,'k^--')
    xlabel('Root density (mm g soil$^{-1}$)')
    ylabel('Microbial biomass (mg C kg soil$^{-1}$)')
    title('Living microbial biomass')
    # gca().set_ylim(bottom=-0.1,top=200)
    # legend(loc='upper left',fontsize='medium')

    data_out=pandas.DataFrame({'H2 Microbial biomass (mgC/kgsoil)':outputs['microbeC'][timepoint,:].mean(axis=0)*1e6,
        'H1 Microbial biomass (mgC/kgsoil)':outputs2['microbeC'][timepoint,:].mean(axis=0)*1e6,
        'rootdensity (mm/gsoil)':rootlength_bins*1e3})

    subplot(121)
    zeropoint=outputs['decomp'][timepoint,0,0].mean(axis=0)/outputs['unprotectedC'][timepoint,0,0].mean(axis=0)
    zeropoint=1.0
    semilogx(rootlength_bins*1e3,outputs['decomp'][timepoint,:,0].mean(axis=0)/outputs['unprotectedC'][timepoint,:,0].mean(axis=0)/10,'go-',label='Labile*0.1')
    data_out['H2 Labile Decomp rate (year-1)']=outputs['decomp'][timepoint,:,0].mean(axis=0)/outputs['unprotectedC'][timepoint,:,0].mean(axis=0)

    zeropoint=outputs['decomp'][timepoint,0,1].mean(axis=0)/outputs['unprotectedC'][timepoint,0,1].mean(axis=0)
    zeropoint=1.0
    semilogx(rootlength_bins*1e3,outputs['decomp'][timepoint,:,1].mean(axis=0)/outputs['unprotectedC'][timepoint,:,1].mean(axis=0),'bo-',label='Resistant')
    data_out['H2 Resistant Decomp rate (year-1)']=outputs['decomp'][timepoint,:,1].mean(axis=0)/outputs['unprotectedC'][timepoint,:,1].mean(axis=0)

    zeropoint=outputs2['decomp'][timepoint,0,0].mean(axis=0)/outputs2['unprotectedC'][timepoint,0,0].mean(axis=0)
    zeropoint=1.0
    semilogx(rootlength_bins*1e3,outputs2['decomp'][timepoint,:,0].mean(axis=0)/outputs2['unprotectedC'][timepoint,:,0].mean(axis=0)/10,'g^--')
    data_out['H1 Labile Decomp rate (year-1)']=outputs2['decomp'][timepoint,:,0].mean(axis=0)/outputs2['unprotectedC'][timepoint,:,0].mean(axis=0)

    zeropoint=outputs2['decomp'][timepoint,0,1].mean(axis=0)/outputs2['unprotectedC'][timepoint,0,1].mean(axis=0)
    zeropoint=1.0
    semilogx(rootlength_bins*1e3,outputs2['decomp'][timepoint,:,1].mean(axis=0)/outputs2['unprotectedC'][timepoint,:,1].mean(axis=0),'b^--')
    data_out['H1 Resistant Decomp rate (year-1)']=outputs2['decomp'][timepoint,:,1].mean(axis=0)/outputs2['unprotectedC'][timepoint,:,1].mean(axis=0)


    xlabel('Root density (mm g soil$^{-1}$)')
    ylabel('Decomposition rate (year$^{-1}$)')
    title('SOC turnover rate')
    legend(fontsize='medium')
    # ylim(-.5,20.2)

    tight_layout()
    draw()

    data_out.to_csv('plotted_data.csv',index_label='Root density percentile')


    figure(3);clf()
    subplot(111)
    totalC=(outputs['unprotectedC'][timepoint,:,:].mean(axis=0).sum(axis=1)+outputs['protectedC'][timepoint,:,:].mean(axis=0).sum(axis=1))*1000
    plot(rootlength_bins,totalC-totalC[0],'go-')
    totalC2=(outputs2['unprotectedC'][timepoint,:,:].mean(axis=0).sum(axis=1)+outputs2['protectedC'][timepoint,:,:].mean(axis=0).sum(axis=1))*1000
    plot(rootlength_bins,totalC2-totalC2[0],'go--')
    xlabel('Root length (m/g soil)')
    ylabel('Difference in soil C (mgC/g soil)')
    title('4-year difference in soil C')

    subplots_adjust(left=0.07,right=0.95,wspace=0.25)

    draw()

show()
