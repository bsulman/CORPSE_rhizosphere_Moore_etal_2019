import CORPSE
import pandas
from pylab import *

# 5% clay
params={
    'vmaxref':[1000,50,600], #Relative maximum enzymatic decomp rates
    'Ea':[37e3,54e3,50e3],    # Activation energy
    'kC':[0.01,0.01,0.01],    # Michaelis-Menton parameter
    'gas_diffusion_exp':2.5,  # Determines suppression of decomp at high soil moisture
    'minMicrobeC':1e-3,       #Minimum microbial biomass (fraction of total C)
    'Tmic':0.25,       # Microbial lifetime
    'et':0.5,          # Fraction of turnover not converted to CO2
    'eup':[0.6,0.05,0.6], # Carbon uptake efficiency
    'tProtected':75.0,    # Protected C turnover time (years)
    'protection_rate':[0.75,0.00005,0.75], # Protected carbon formation rate (year-1)
}

dt=1.0/365 # Daily time step
depth=5 #cm

# Read daily temperature data for upland site
# datafile='../RDS-2005-0001_UTF8/UTF8/weather/dailyWeather.txt'
# weatherdata=pandas.read_csv(datafile,index_col=('date','watershed','position'),encoding='utf-8-sig',parse_dates=True)
weatherdata=pandas.read_csv('../SPRUCE-data/EM123_Combined_public_data_2010_2015.csv',index_col='datetime',parse_dates=True,na_values=[-9999.0])
soilT=weatherdata['EM1_Hummock10cm']+273.15
soilmoisture=weatherdata['EM1_VW1_TopofHummock']
soilmoisture[soilmoisture<0.0]=0.0
soilmoisture=soilmoisture/soilmoisture.max()

# Constant soil moisture for now
# theta=zeros(len(soilT))+0.7
theta=soilmoisture.fillna(method='backfill').resample('1D')[365:]
T=soilT.fillna(method='backfill').resample('1D')[365:]

inputs = array([0.1,0.9,0.0])*2.0e-3 # gC/g soil/year

# Calculate ranges of root density
rootdata=pandas.read_excel('../All Site SEM data SULMAN with volume.xlsx')
bulkdensity=1.26 #g/cm3
# rootdensity=0.1e-3 # g/cm3 soil
rootdensity=rootdata['roots per volume']*1e-3 #g/cm3 of soil
rootdensity[rootdensity<0]=nan
rootmass=rootdata['roots']*1e-3 # g roots/g soil
# srl=40.0      # m/g
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

rootlength_bins=rootlength.quantile(linspace(0.5/nbins,1.0-0.5/nbins,nbins))
exudationrate=array(rootlength_bins*100*1e-6*24*365/1000.0) # Phillips et al: 1 ugC/cm root length/hour

c=CORPSE.soil_carbon_cohort(litterC=[6e-05,0.012,9.6e-05], protectedC=[0.0037,4.2e-05,0.0054], livingMicrobeC=8.4e-05, params=params)

cohorts=[]
for ii in xrange(nbins):
    cohorts.append(c.copy())

if do_spinup:
    nsteps=365*50
else:
    nsteps=365*4


outputs={'unprotectedC':zeros((nsteps,nbins,3)),
            'protectedC':zeros((nsteps,nbins,3)),
            'microbeC':zeros((nsteps,nbins)),
            'CO2':zeros((nsteps,nbins))
            }

doy=T.index.dayofyear
exudation_ts=cos((doy-0.67*365)*2*pi/365)+1

for step in xrange(nsteps):
    if step%365==0:
        print 'Year %d of %d'%(floor(step/365),floor(nsteps/365))
    for cc in xrange(len(cohorts)):
        out=cohorts[cc].update(T[step%len(T)],theta[step%len(T)],dt)
        cohorts[cc].check_validity()
        outputs['unprotectedC'][step,cc,:]=cohorts[cc].litterC
        outputs['protectedC'][step,cc,:]=cohorts[cc].protectedC
        outputs['microbeC'][step,cc]=cohorts[cc].livingMicrobeC
        outputs['CO2'][step,cc]=cohorts[cc].CO2

        cohorts[cc].add_carbon((inputs+array([exudationrate[cc],0.0,0.0]))*dt*exudation_ts[step%len(T)])


# Plot results
t=arange(nsteps)/365.0

if do_spinup:
    figure(1);clf()
    plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1),'b-',label='Unprotected')
    plot(t,outputs['protectedC'][:,0,:].sum(axis=1),'r-',label='Protected')
    plot(t,outputs['unprotectedC'][:,0,0],':g',label='Fast')
    plot(t,outputs['unprotectedC'][:,0,1],'b:',label='Slow')
    plot(t,outputs['unprotectedC'][:,0,2],'r:',label='Microbe necro')
    plot(t,outputs['microbeC'][:,0],'m-',label='Microbe')
    plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1)+outputs['protectedC'][:,0,:].sum(axis=1),'k-',label='Total')
    legend(loc='best',fontsize='medium')
    draw()

else:
    figure(1,figsize=(13.6,5.3));clf()
    subplot(131)
    # plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1),'b-',label='Unprotected')
    # plot(t,outputs['protectedC'][:,0,:].sum(axis=1),'r-',label='Protected')
    plot(t,outputs['unprotectedC'][:,0,0]*1e6,'-g',label='Fast')
    # plot(t,outputs['unprotectedC'][:,0,1],'c:',label='Slow')
    plot(t,outputs['unprotectedC'][:,0,2]*1e6,'y-',label='Dead microbe')
    plot(t,outputs['microbeC'][:,0]*1e6,'m-',label='Live microbe')
    # plot(t,outputs['unprotectedC'][:,0,:].sum(axis=1)+outputs['protectedC'][:,0,:].sum(axis=1),'k-',label='Total')
    leg=legend(loc='best',fontsize='medium');leg.get_frame().set_alpha(0.5)

    # plot(t,outputs['unprotectedC'][:,-1,:].sum(axis=1),'b--',label='Unprotected')
    # plot(t,outputs['protectedC'][:,-1,:].sum(axis=1),'r--',label='Protected')
    plot(t,outputs['unprotectedC'][:,-1,0]*1e6,'--g',label='Fast')
    # plot(t,outputs['unprotectedC'][:,-1,1],'c-.',label='Slow')
    plot(t,outputs['unprotectedC'][:,-1,2]*1e6,'y--',label='Microbe necro')
    plot(t,outputs['microbeC'][:,-1]*1e6,'m--',label='Microbe')
    # plot(t,outputs['unprotectedC'][:,0-1,:].sum(axis=1)+outputs['protectedC'][:,0-1,:].sum(axis=1),'k--',label='Total')
    xlabel('Time (years)')
    ylabel('Carbon content(mg/kg soil)')
    title('Labile carbon pools')

    subplot(132)
    plot(log10(rootlength_bins/(srl.mean())*1e3+1),outputs['microbeC'][-140,:]*1e6,'bo-',label='August')
    plot(log10(rootlength_bins/(srl.mean())*1e3+1),outputs['microbeC'][-201,:]*1e6,'go-',label='June')
    xlabel('Log(root biomass [g/g soil] + 1)')
    ylabel('Microbial biomass (mg/kg soil)')
    title('Living microbial biomass')
    legend(loc='upper left',fontsize='medium')

    subplot(133)
    plot(log10(rootlength_bins/(srl.mean())*1e3+1),(outputs['unprotectedC'][-140,:,:].sum(axis=1)+outputs['protectedC'][-140,:,:].sum(axis=1))*1000,'bo-')
    plot(log10(rootlength_bins/(srl.mean())*1e3+1),(outputs['unprotectedC'][-201,:,:].sum(axis=1)+outputs['protectedC'][-201,:,:].sum(axis=1))*1000,'go-')
    xlabel('Log(root biomass [g/g soil] + 1)')
    ylabel('Soil C (gC/kg soil)')
    title('Exudation effect on soil C')

    subplots_adjust(left=0.07,right=0.95,wspace=0.25)

    draw()

show()
