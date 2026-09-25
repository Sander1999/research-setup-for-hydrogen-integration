"""Conservative H2/H2O/NaCl flash for an ideal isothermal source-rock case.

Adapted from the verified bachelor DARTS implementation; capillary and
relative-permeability conventions follow the minor Bourakébougou model.

H2 solubility: Kerkache et al. (2024), doi:10.1016/j.molliq.2024.124497.
The published true-species mole fraction is converted to a salt-free basis.
Water and NaCl are nonvolatile here; aqueous density/viscosity are declared
engineering inputs, not a brine EOS. No mineral equilibrium is implied.
"""
from pathlib import Path
import json
import numpy as np
from functools import lru_cache
from iapws import IAPWS97
from darts.physics.base.property_container import PropertyContainer
from darts.physics.properties.basic import ConstFunc, RockCompactionEvaluator

COMPONENTS=['H2','H2O','NaCl']
MW=np.array([2.01588,18.01528,58.4428])  # kg/kmol
FE_MW=55.845
G=9.80665e-5  # bar/m per kg/m3


def load_config(parameters=None):
    c=json.loads(Path(__file__).with_name('config.json').read_text())
    if isinstance(parameters,(str,Path)):parameters=json.loads(Path(parameters).read_text())
    if parameters:
        if set(parameters)-set(c):raise ValueError('Unknown config keys: '+str(set(parameters)-set(c)))
        c.update(parameters)
    for key in ['nx','ny','nz']:
        if not isinstance(c[key],int) or isinstance(c[key],bool) or c[key]<2:raise ValueError('Grid dimensions must be integers >=2')
    for k,v in c.items():
        if isinstance(v,(int,float)) and not isinstance(v,bool) and (not np.isfinite(v) or v<0):
            raise ValueError(k+' must be finite and nonnegative')
    for k in ['porosity','length_x_m','length_y_m','thickness_m','runtime_days','half_time_days',
              'rock_density_kg_m3','host_density_kg_m3','host_permeability_md','source_permeability_md',
              'damage_stress_scale_mpa','first_timestep_days','max_timestep_days','obl_pressure_step_bar',
              'obl_hydrogen_step','obl_water_step','brine_density_kg_m3','brine_viscosity_cp']:
        if c[k]<=0:raise ValueError(k+' must be positive')
    for k in ['fe_mass_fraction','ferrous_fraction','accessible_fraction','porosity']:
        if not 0<c[k]<=(1 if k!='porosity' else .99):raise ValueError('Invalid fraction '+k)
    if not 298.15<=c['temperature_k']<=453.15 or not 10<c['producer_bhp_bar']<1000:
        raise ValueError('Kerkache correlation requires 298.15–453.15K and pressures below1000bar')
    if not .1<=c['salinity_molal']<=6:raise ValueError('Brine scenario requires0.1–6mol/kg water; the correlation has a small pure-water discontinuity')
    if c['boundary_mode'] not in ['closed','natural_flow']:raise ValueError('Unknown boundary mode')
    if c['residual_water_saturation']+c['residual_gas_saturation']>=1:raise ValueError('Residual saturations sum to >=1')
    for key in ['source_centre_m','source_semiaxes_m']:
        arr=np.asarray(c[key],float)
        if arr.shape!=(3,) or not np.isfinite(arr).all():raise ValueError('Invalid '+key)
    if min(c['source_semiaxes_m'])<=0:raise ValueError('Positive source axes required')
    lo=np.array([0.,0.,c['top_depth_m']]);hi=lo+[c['length_x_m'],c['length_y_m'],c['thickness_m']]
    ctr=np.asarray(c['source_centre_m']);rad=np.asarray(c['source_semiaxes_m'])
    if np.any(ctr-rad<lo) or np.any(ctr+rad>hi):raise ValueError('Source ellipsoid must fit inside domain')
    return c


A_H2=np.array([.05888460,-.06136111,-.002650473,.002731125,.001802374,-.001150707,.9588528e-4,-.1109040e-6,.1264403e-9])
B_H2=np.array([1.325,1.87,2.5,2.8,2.938,3.14,3.37,3.75,4.])
D_H2=np.array([1.,1.,2.,2.,2.42,2.63,3.,4.,5.])

@lru_cache(maxsize=32)
def saturation_pressure_mpa(temperature_k):
    return float(IAPWS97(T=float(temperature_k),x=0).P)


def hydrogen_fugacity_coefficient(pressure_bar,temperature_k):
    """Thermodynamic integral ln(phi)=integral_0^p(Z-1)/p dp of Lemmon Z."""
    p=np.asarray(pressure_bar)
    return np.exp(np.sum(A_H2*(100/temperature_k)**B_H2*(p[...,None]/10)**D_H2/D_H2,axis=-1))


def solubility_x(pressure_bar,temperature_k,molality,particles_per_nacl=2.):
    """Salt-FREE aqueous xH2 from Kerkache2024 Eqs9–15, Table3.

    T298.15–453.15K, P<=1000bar, NaCl<=6mol/kg. Pressure numerator
    and Henry constant use MPa; the Poynting exponent uses Pa and SI R.
    Adopt a dissociated-species basis for published true mole fraction.
    This is an explicit convention inference; use particles_per_nacl=1
    to test the alternative formula-unit convention.
    Return nH2/(nH2+nH2O), converting that basis explicitly.
    Positive-brine Eq12 has ln(gamma)=m*(A0*T+A1/T+A3+A4*m*T)+A2*T.
    At exactly zero salt gamma=1 by pure-water definition (a ~1% fitted
    discontinuity); simulations therefore require m>=0.1 and audit m<=6.
    """
    p=np.asarray(pressure_bar)/10.;t=temperature_k;m=np.asarray(molality)
    h=np.exp(np.polyval([-4.40347e-10,8.33536e-7,-6.05506e-4,.190837,-12.7991],t))
    ln_gamma=np.where(m>0,m*(1.59054e-3*t+2.63936e2/t-1.14604-1.18419e-5*m*t)+3.79192e-5*t,0.)
    ps=saturation_pressure_mpa(float(t));poynting=np.exp(19e-6*(p-ps)*1e6/(8.314462618*t))
    # Water vapour correction only enters equilibrium fugacity; the native
    # flow flash neglects its small transported mass in the H2 gas phase.
    true_x=(p-ps)*hydrogen_fugacity_coefficient(pressure_bar,t)/(h*np.exp(ln_gamma)*poynting)
    ion_ratio=particles_per_nacl*m*MW[1]/1000
    return true_x*(1+ion_ratio)/(1+ion_ratio*true_x)


def hydrogen_density(p,t):
    """NIST pressure-explicit H2 density correlation, Lemmon et al.2008."""
    a=np.array([.05888460,-.06136111,-.002650473,.002731125,.001802374,-.001150707,.9588528e-4,-.1109040e-6,.1264403e-9])
    b=np.array([1.325,1.87,2.5,2.8,2.938,3.14,3.37,3.75,4.])
    d=np.array([1.,1.,2.,2.,2.42,2.63,3.,4.,5.])
    z=1+np.sum(a*(100/t)**b*(p/10)**d)
    return p*1e5*(MW[0]/1000)/(z*8.314462618*t)


class BrineFlash:
    def __init__(self,cfg):self.cfg=cfg
    def split(self,p,z):
        z=np.maximum(np.asarray(z,float),1e-14);z/=z.sum()
        h,w,s=z
        m=1000*s/(w*MW[1])
        xs=float(solubility_x(p,self.cfg['temperature_k'],m,self.cfg['salt_basis_particles_per_nacl']))
        if not 0<xs<1:raise ValueError('Invalid solubility outside property domain')
        h_diss=min(h,xs/(1-xs)*w)
        gas=max(h-h_diss,0.)
        nu=np.array([gas,1-gas]);x=np.zeros((2,3))
        x[0]=[1,0,0];x[1]=[h_diss,w,s];x[1]/=x[1].sum()
        if np.max(np.abs(nu@x-z))>1e-12:raise RuntimeError('Analytical flash did not conserve components')
        return nu,x
    def state(self,p,z):
        nu,x=self.split(p,z);t=self.cfg['temperature_k']
        rho=np.array([hydrogen_density(p,t),self.cfg['brine_density_kg_m3']*(1+self.cfg['brine_compressibility_per_bar']*(p-self.cfg['producer_bhp_bar']))])
        rm=rho/(x@MW);vol=nu/rm;sat=vol/vol.sum()
        return {'nu':nu,'x':x,'rho':rho,'rho_m':rm,'sat':sat,
                'molality':1000*x[1,2]/(x[1,1]*MW[1])}


class Density:
    def __init__(self,cfg,gas):self.cfg,self.gas=cfg,gas
    def evaluate(self,p,t,x):
        return hydrogen_density(p,t) if self.gas else self.cfg['brine_density_kg_m3']*(1+self.cfg['brine_compressibility_per_bar']*(p-self.cfg['producer_bhp_bar']))


class RelPerm:
    def evaluate(self,s):return float(np.clip(s,0,1)**2)


class Properties(PropertyContainer):
    def __init__(self,cfg):
        self.cfg=cfg;self.system=BrineFlash(cfg)
        super().__init__(['gas','brine'],COMPONENTS,MW,eps_z=1e-14,rock_comp=cfg['rock_compressibility_per_bar'],temperature=cfg['temperature_k'])
        self.rock_compr_ev=RockCompactionEvaluator(pref=cfg['producer_bhp_bar'],compres=cfg['rock_compressibility_per_bar'])
        self.density_ev={'gas':Density(cfg,True),'brine':Density(cfg,False)}
        self.viscosity_ev={'gas':ConstFunc(.0094),'brine':ConstFunc(cfg['brine_viscosity_cp'])}
        swr,sgr=cfg['residual_water_saturation'],cfg['residual_gas_saturation']
        self.rel_perm_ev={'gas':RelativePermeability(sgr,swr,2.),'brine':RelativePermeability(swr,sgr,3.)}
        self.capillary_pressure_ev=CapillaryPressure(cfg)
        self.diffusion_ev={'gas':ConstFunc(np.zeros(3)), 'brine':ConstFunc(np.ones(3)*cfg['aqueous_effective_diffusion_m2_s']*86400)}
        # Additional hydration is a declared congruent-alteration mass sink.
        redox_water_per_kg_rock=cfg['fe_mass_fraction']*cfg['ferrous_fraction']*MW[1]/(3*FE_MW)
        self.water_stoichiometry=1+cfg['hydration_water_kg_per_kg_altered_rock']/redox_water_per_kg_rock
    def run_flash(self,pressure,temperature,zc,evaluate_PT=False):
        self.nu,self.x=self.system.split(pressure,zc)
        return np.flatnonzero(self.nu>0)
    def evaluate_mass_source(self,pressure,temperature,zc):
        # Reference reaction:3FeO+H2O -> Fe3O4+H2. Positive native KIN is
        # a sink; mesh.kin_factor supplies kmol extent/(m3 bulk day) relative
        # to this unit reference. All solids remain in the separate Fe ledger.
        self.mass_source[:]=[-1.,self.water_stoichiometry,0.]
        return self.mass_source


class RelativePermeability:
    def __init__(self,own,other,exponent):self.own,self.other,self.exponent=own,other,exponent
    def evaluate(self,s):return float(np.clip((s-self.own)/(1-self.own-self.other),0,1)**self.exponent)

class CapillaryPressure:
    """Regularised modified Brooks–Corey; inherited minor-model convention Pg−Pw."""
    def __init__(self,cfg):self.cfg=cfg
    def evaluate(self,s):
        c=self.cfg;se=np.clip((s[1]-c['residual_water_saturation'])/(1-c['residual_water_saturation']),1e-8,1)
        pc=min(c['capillary_max_bar'],c['capillary_entry_bar']*(se**(-1/c['capillary_exponent'])-1))
        return np.array([0.,pc])
