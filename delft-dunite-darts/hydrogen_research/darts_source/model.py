"""Native 3-D reactive flow with finite Fe and explicitly reduced damage.

Damage and reaction rates are empirical scenarios. No displacement equation,
fracture network, mineral equilibrium or calibrated serpentinization rate is
claimed. Geometric pore volume is fixed; pressure compliance is conservative.
"""
from __future__ import annotations
import numpy as np
from darts.engines import value_vector,well_control_iface
from darts.models.darts_model import DartsModel
from darts.nonlinear_solvers import NewtonSolver
from darts.physics.base.physics import PhysicsBase
from darts.reservoirs.struct_reservoir import StructReservoir
from .properties import COMPONENTS,MW,FE_MW,G,Properties,load_config,solubility_x


def grid_geometry(c):
    """Layer-aligned z mesh; nx/ny exact, nz a requested minimum resolution."""
    dx=c['length_x_m']/c['nx'];dy=c['length_y_m']/c['ny']
    top=c['top_depth_m'];bottom=top+c['thickness_m'];target_dz=c['thickness_m']/c['nz']
    boundaries=sorted(set([top,bottom]+[float(z) for l in c['layers'] for z in (l['top_depth_m'],l['bottom_depth_m']) if top<z<bottom]))
    edges=[]
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        edges.extend(np.linspace(lo,hi,max(1,int(np.ceil((hi-lo)/target_dz)))+1)[:-1].tolist())
    edges=np.array(edges+[bottom]);zc=(edges[:-1]+edges[1:])/2;dz=np.diff(edges)
    z,y,x=np.meshgrid(zc,(np.arange(c['ny'])+.5)*dy,(np.arange(c['nx'])+.5)*dx,indexing='ij')
    xyz=np.column_stack((x.ravel(),y.ravel(),z.ravel()))
    dims=np.column_stack((np.full(len(xyz),dx),np.full(len(xyz),dy),np.repeat(dz,c['nx']*c['ny'])))
    return xyz,dims


def source_fractions(c,xyz,dims):
    """Subcell Gaussian x/y integration of analytic ellipsoid z intersections.

    Partial cells are corrected conservatively to the exact ellipsoid volume.
    All fractions remain in [0,1]; this avoids changing Fe mass on refinement.
    """
    ctr=np.asarray(c['source_centre_m']);r=np.asarray(c['source_semiaxes_m'])
    q,w=np.polynomial.legendre.leggauss(12);fractions=np.zeros(len(xyz))
    for i,wi in zip(q,w):
        for j,wj in zip(q,w):
            xy=xyz[:,:2]+np.column_stack((i*dims[:,0]/2,j*dims[:,1]/2))
            inside=1-np.sum(((xy-ctr[:2])/r[:2])**2,axis=1)
            dz=r[2]*np.sqrt(np.maximum(inside,0))
            lo=np.maximum(xyz[:,2]-dims[:,2]/2,ctr[2]-dz)
            hi=np.minimum(xyz[:,2]+dims[:,2]/2,ctr[2]+dz)
            fractions+=wi*wj/4*np.maximum(hi-lo,0)/dims[:,2]*(inside>0)
    target=4*np.pi*np.prod(r)/3;volumes=np.prod(dims,axis=1)
    for _ in range(10):
        delta=target-fractions@volumes
        if abs(delta)<1e-11*target:break
        partial=(fractions>1e-12)&(fractions<1-1e-12)
        weights=fractions*(1-fractions)*partial
        if weights.sum()==0:raise ValueError('Grid too coarse to represent the source')
        fractions=np.clip(fractions+delta*weights/(weights@volumes),0,1)
    if abs(fractions@volumes-target)>1e-9*target:raise RuntimeError('Source-volume quadrature adjustment failed')
    return fractions


def damage_closure(progress,source_fraction,c):
    """Irreversible scalar damage target from an assumed constrained eigenstrain.

    progress is fraction of the accessible reactive rock, NOT all source rock.
    This strain is a diagnostic: it is not imposed on geometric pore storage.
    """
    strain=c['solid_expansion_fraction']*(1-c['porosity'])*source_fraction*c['accessible_fraction']*progress
    stress=c['effective_confinement_modulus_mpa']*strain
    damage=-np.expm1(-np.maximum(stress-c['tensile_threshold_mpa'],0)/c['damage_stress_scale_mpa']) if c['damage_enabled'] else np.zeros_like(progress)
    return strain,stress,damage


class Model(DartsModel):
    def __init__(self,parameters=None):
        self.cfg=load_config(parameters);c=self.cfg
        super().__init__();self.timer.node['initialization'].start()
        self.xyz,self.dims=grid_geometry(c);self.nb=len(self.xyz);self.depths=self.xyz[:,2]
        self.cell_volume=np.prod(self.dims,axis=1);self.bulk_volume=float(self.cell_volume.sum())
        self.nz=self.nb//(c['nx']*c['ny']);self.z_centres=np.unique(self.depths)
        self.height=c['thickness_m'];self.source_fraction=source_fractions(c,self.xyz,self.dims)
        self.layer_ids=np.zeros(self.nb,dtype=int);self.porosity=np.full(self.nb,c['porosity'])
        self.host_permeability=np.full(self.nb,c['host_permeability_md'])
        self.vertical_ratio=np.full(self.nb,c['vertical_permeability_ratio'])
        self.host_density=np.full(self.nb,c['host_density_kg_m3'])
        for i,layer in enumerate(c['layers']):
            mask=(self.depths>=layer['top_depth_m'])&(self.depths<layer['bottom_depth_m'])
            self.layer_ids[mask]=i;self.porosity[mask]=layer['porosity']
            self.host_permeability[mask]=layer['permeability_md'];self.vertical_ratio[mask]=layer['vertical_permeability_ratio']
            self.host_density[mask]=layer['bulk_density_kg_m3']
        # Source porosity is a volume mixture where the body crosses a layer.
        self.porosity=self.porosity*(1-self.source_fraction)+c['porosity']*self.source_fraction
        self.source_rock_mass=self.source_fraction*self.cell_volume*(1-c['porosity'])*c['rock_density_kg_m3']
        self.fe_initial=self.source_rock_mass*c['fe_mass_fraction']*c['ferrous_fraction']*c['accessible_fraction']/FE_MW
        self.fe_initial_kmol=float(self.fe_initial.sum());self.fe_remaining=self.fe_initial.copy()
        self.generated_kmol=0.;self._pending_extent=np.zeros(self.nb);self.damage=np.zeros(self.nb)
        self.strain=np.zeros(self.nb);self.stress=np.zeros(self.nb)
        # Harmonic lithological mixture for partially intersected source cells.
        self.base_permeability=1/((1-self.source_fraction)/self.host_permeability+self.source_fraction/c['source_permeability_md'])
        self.corridor_fraction=np.zeros(self.nb)
        if c['corridor_enabled']:
            bounds=np.array([c['corridor_x_bounds_m'],c['corridor_y_bounds_m'],c['corridor_depth_bounds_m']])
            lo=np.maximum(self.xyz-self.dims/2,bounds[:,0]);hi=np.minimum(self.xyz+self.dims/2,bounds[:,1])
            self.corridor_fraction=np.prod(np.maximum(hi-lo,0)/self.dims,axis=1)
            vertical=self.base_permeability*self.vertical_ratio
            self.base_permeability=(1-self.corridor_fraction)*self.base_permeability+self.corridor_fraction*c['corridor_horizontal_permeability_md']
            vertical=(1-self.corridor_fraction)*vertical+self.corridor_fraction*c['corridor_vertical_permeability_md']
            self.vertical_ratio=vertical/self.base_permeability
        self.permeability=self.base_permeability.copy()
        self.reservoir=StructReservoir(self.timer,nx=c['nx'],ny=c['ny'],nz=self.nz,
            dx=self.dims[:,0],dy=self.dims[:,1],dz=self.dims[:,2],
            permx=self.base_permeability,permy=self.base_permeability,
            permz=self.base_permeability*self.vertical_ratio,
            poro=self.porosity,depth=self.depths,op_num=self.layer_ids,cache=False)
        self.property_container=Properties(c);self.system=self.property_container.system
        self.physics=PhysicsBase(COMPONENTS,['gas','brine'],self.timer,
            axes_step=[c['obl_pressure_step_bar'],c['obl_hydrogen_step'],c['obl_water_step']],
            axes_origin=[1.,1e-14,1e-14],epsilon_z=1e-14,extrapolation_flag=False,
            state_spec=PhysicsBase.StateSpecification.P,cache=False)
        self.regional_properties={}
        for i in range(max(1,len(c['layers']))):
            cfg_region=c.copy()
            if c['layers']:cfg_region['capillary_entry_bar']=c['layers'][i]['capillary_entry_bar']
            self.regional_properties[i]=Properties(cfg_region)
            self.physics.add_property_region(self.regional_properties[i],region=i)
        w=1/(1+c['salinity_molal']*MW[1]/1000)
        self.inj_composition=np.array([1e-12,w*(1-1e-12),(1-w)*(1-1e-12)])
        self.timer.node['initialization'].stop()
    def hydrostatic_pressure(self,z):
        c=self.cfg;d=np.asarray(z)-c['datum_depth_m'];beta=c['brine_compressibility_per_bar']
        if beta:return c['producer_bhp_bar']+np.expm1(beta*G*c['brine_density_kg_m3']*d)/beta
        return c['producer_bhp_bar']+G*c['brine_density_kg_m3']*d
    def set_wells(self):
        c=self.cfg
        if c['boundary_mode']=='closed':return
        # One virtual boundary segment per depth, perforated across y. These
        # represent distant aquifer heads; they are not installed pumping wells.
        for k in range(1,self.nz+1):
            for side,i in [('LEFT',1),('RIGHT',c['nx'])]:
                name=f'{side}_{k}';self.reservoir.add_well(name)
                for j in range(1,c['ny']+1):
                    cell=(k-1)*c['nx']*c['ny']+(j-1)*c['nx']+(i-1)
                    # Exact half-cell face conductance in DARTS field units,
                    # rather than a radius-dependent Peaceman pumping index.
                    conductance=.00852671467191601*self.base_permeability[cell]*self.dims[cell,1]*self.dims[cell,2]/(.5*self.dims[cell,0])
                    self.reservoir.add_perforation(name,res_cell_idx=(i,j,k),
                        well_diameter=2*c['boundary_well_radius_m'],well_index=conductance,well_indexD=0.)
    def set_well_controls(self):
        c=self.cfg
        delta=G*c['brine_density_kg_m3']*c['length_x_m']*c['hydraulic_gradient']/2
        for w in self.reservoir.wells:
            left=w.name.startswith('LEFT');k=int(w.name.split('_')[1])-1
            z=self.z_centres[k]
            target=float(self.hydrostatic_pressure(z))+(delta if left else -delta)
            self.physics.set_well_controls(wctrl=w.control,control_type=well_control_iface.BHP,is_inj=True,
                target=target,inj_composition=self.inj_composition.tolist())
    def set_initial_conditions(self):
        c=self.cfg;p=self.hydrostatic_pressure(self.depths)
        if c['boundary_mode']=='natural_flow':
            p=p+G*c['brine_density_kg_m3']*c['hydraulic_gradient']*(c['length_x_m']/2-self.xyz[:,0])
        self.physics.set_initial_conditions_from_array(self.reservoir.mesh,
            {'pressure':p,'H2':np.full(self.nb,1e-12),'H2O':np.full(self.nb,self.inj_composition[1])})
        self.reservoir.mesh.kin_factor=value_vector(np.zeros(self.nb))
    def set_solver(self):
        super().set_solver();c=self.cfg
        self.ts_control.dt_first=c['first_timestep_days'];self.ts_control.dt_min=1e-10
        self.ts_control.dt_max=c['max_timestep_days'];self.ts_control.dt_mult=1.5
        self.ts_control.runtime=c['runtime_days']
        self.nonlinear_solver=NewtonSolver(tolerance=c['nonlinear_tolerance'],max_iterations=20)
        self.linear_solver.spec.tolerance=1e-11
    def init(self,*args,**kwargs):
        result=super().init(*args,**kwargs)
        self._base_trans=np.array(self.reservoir.mesh.tran,copy=True)
        self._bm=np.array(self.reservoir.mesh.block_m,dtype=int);self._bp=np.array(self.reservoir.mesh.block_p,dtype=int)
        self._initial_pressure=self.reservoir_state()[:,0].copy()
        self._initial_fluid_mass=self.cell_fluid_mass()
        return result
    def progress(self):
        return np.divide(self.fe_initial-self.fe_remaining,self.fe_initial,out=np.zeros(self.nb),where=self.fe_initial>0)
    def reaction_factor(self,q):
        c=self.cfg
        if not c['feedback_enabled']:return np.ones(self.nb)
        xi=self.progress();contact=np.clip(q['sat'][:,1],0,1)**c['water_contact_exponent']
        dissolved=q['x'][:,1,0]/(q['x'][:,1,0]+q['x'][:,1,1])
        sat=solubility_x(self.reservoir_state()[:,0],c['temperature_k'],q['molality'],c['salt_basis_particles_per_nacl'])
        inhibition=1/(1+c['inhibition_strength']*dissolved/sat)
        surface=1+c['surface_area_gain']*self.damage
        passivation=np.exp(-c['passivation_strength']*xi)
        return contact*inhibition*surface*passivation
    def update_native_permeability(self):
        c=self.cfg
        self.permeability=self.base_permeability*(1+c['permeability_gain']*self.damage**c['damage_permeability_exponent'])
        a,b=self._bm,self._bp;inside=(a<self.nb)&(b<self.nb)
        multipliers=np.ones(len(a));ai,bi=a[inside],b[inside]
        # Uniform spacings and common directional anisotropy make harmonic
        # cell-k ratios exact for this TPFA grid. Reverse connections match.
        vertical=np.abs(self.xyz[bi,2]-self.xyz[ai,2])>1e-8
        ra=np.where(vertical,self.vertical_ratio[ai],1.);rb=np.where(vertical,self.vertical_ratio[bi],1.)
        axis=np.argmax(np.abs(self.xyz[bi]-self.xyz[ai]),axis=1)
        da=self.dims[ai,axis];db=self.dims[bi,axis]
        multipliers[inside]=(da/(self.base_permeability[ai]*ra)+db/(self.base_permeability[bi]*rb))/(da/(self.permeability[ai]*ra)+db/(self.permeability[bi]*rb))
        for i in np.flatnonzero(~inside):
            cells=[j for j in (a[i],b[i]) if j<self.nb]
            if cells:multipliers[i]=self.permeability[cells[0]]/self.base_permeability[cells[0]]
        np.asarray(self.reservoir.mesh.tran)[:]=self._base_trans*multipliers
    def run_timestep(self,dt,t,verbose=0):
        q=self.cell_properties();c=self.cfg
        fraction=-np.expm1(-np.log(2)*dt/c['half_time_days']*self.reaction_factor(q)) if c['reaction_enabled'] else np.zeros(self.nb)
        extent=self.fe_remaining/3*fraction
        pv=self.current_pore_volume()
        water=pv*np.sum(q['sat']*q['rho_m']*q['x'][:,:,1],axis=1)
        self._pending_extent=np.minimum(extent,.1*water/self.property_container.water_stoichiometry)
        self.reservoir.mesh.kin_factor=value_vector(self._pending_extent/(self.cell_volume*dt))
        # Immutable accepted damage: repeated rejected attempts do not increase k.
        self.update_native_permeability()
        return super().run_timestep(dt,t,verbose)
    def after_converged_timestep(self):
        self.fe_remaining-=3*self._pending_extent;self.generated_kmol+=float(self._pending_extent.sum())
        if self.fe_remaining.min()<-1e-14:raise RuntimeError('Finite Fe inventory became negative')
        self.strain,self.stress,target=damage_closure(self.progress(),self.source_fraction,self.cfg)
        self.damage=np.maximum(self.damage,target)
        # Do NOT update tran here: the accepted-step ledger needs the k that
        # actually entered this backward-Euler solve, not next step's damage.
        super().after_converged_timestep()
    def reservoir_state(self):return np.asarray(self.physics.engine.X)[:3*self.nb].reshape(self.nb,3).copy()
    def cell_properties(self,states=None):
        states=self.reservoir_state() if states is None else states
        if hasattr(self,'_property_state') and np.array_equal(states,self._property_state):return self._property_result
        q=[self.system.state(float(s[0]),np.r_[s[1:],1-s[1:].sum()]) for s in states]
        self._property_state=states.copy();self._property_result={k:np.array([v[k] for v in q]) for k in q[0]}
        return self._property_result
    def current_pore_volume(self):
        p=self.reservoir_state()[:,0]
        return self.cell_volume*self.porosity*(1+self.cfg['rock_compressibility_per_bar']*(p-self.cfg['producer_bhp_bar']))
    def cell_fluid_mass(self):
        q=self.cell_properties()
        return self.current_pore_volume()*np.sum(q['sat']*q['rho'],axis=1)
    def operators(self):
        states=np.asarray(self.physics.engine.X).reshape(-1,3)
        if hasattr(self,'_operator_state') and np.array_equal(states,self._operator_state):return self._operator_result
        ev=self.physics.reservoir_operators[0];ops=[]
        for s in states:
            a=value_vector(np.zeros(ev.n_ops))
            (self.physics.acc_flux_itor[int(self.layer_ids[len(ops)])] if len(ops)<self.nb else self.physics.acc_flux_w_itor).evaluate(value_vector(s),a)
            ops.append(np.array(a))
        self._operator_state=states.copy();self._operator_result=np.array(ops)
        return self._operator_result
    def inventory(self,interpolated=False):
        if interpolated:
            ev=self.physics.reservoir_operators[0];ops=self.operators()[:self.nb]
            kmol=((self.cell_volume*self.porosity)[:,None]*ops[:,ev.ACC_OP:ev.ACC_OP+3]).sum(axis=0)
            return {'component_kmol':kmol,'component_kg':kmol*MW}
        q=self.cell_properties();phase=self.current_pore_volume()[:,None,None]*(q['sat']*q['rho_m'])[:,:,None]*q['x']
        return {'component_kmol':phase.sum(axis=(0,1)),'component_kg':phase.sum(axis=(0,1))*MW,
                'free_H2_kg':float(phase[:,0,0].sum()*MW[0]),'dissolved_H2_kg':float(phase[:,1,0].sum()*MW[0])}
    def face_phase_flows(self):
        """Native Darcy flux m3/day for all directed faces, from m toward p."""
        mesh=self.reservoir.mesh;s=np.asarray(self.physics.engine.X).reshape(-1,3);op=self.operators()
        ev=self.physics.reservoir_operators[0];a,b=self._bm,self._bp
        flows=np.zeros((len(a),2));component=np.zeros((len(a),2,3))
        for ph in range(2):
            dp=s[b,0]-s[a,0]+.5*(op[a,ev.GRAV_OP+ph]+op[b,ev.GRAV_OP+ph])*np.asarray(mesh.grav_coef)-op[b,ev.PC_OP+ph]+op[a,ev.PC_OP+ph]
            upstream=np.where(dp<0,a,b);flows[:,ph]=-np.asarray(mesh.tran)*dp*op[upstream,ev.LAMBDA_OP+ph]
            component[:,ph,:]=flows[:,ph,None]*op[upstream,ev.FLUX_OP+3*ph:ev.FLUX_OP+3*(ph+1)]*MW
        return flows,component
    def perforation_rates(self):
        _,component=self.face_phase_flows();out={}
        for w in self.reservoir.wells:
            phase=np.zeros((2,3))
            for seg,cell,_wi,_wid in w.perforations:
                wc=w.well_body_idx+seg;idx=np.flatnonzero((self._bm==cell)&(self._bp==wc))
                if len(idx)!=1:raise RuntimeError('Unexpected well connection topology')
                phase+=component[idx[0]]
            out[w.name]={'component_kg_day':phase.sum(axis=0),'phase_component_kg_day':phase}
        return out
    def internal_velocity(self):
        flows,_=self.face_phase_flows();v=np.zeros((self.nb,3));n=np.zeros((self.nb,3))
        for i,(a,b) in enumerate(zip(self._bm,self._bp)):
            if a>=self.nb or b>=self.nb or a>b:continue
            dx=self.xyz[b]-self.xyz[a];axis=int(np.argmax(abs(dx)))
            value=flows[i,1]/(self.cell_volume[a]/self.dims[a,axis])*np.sign(dx[axis])
            for cell in [a,b]:v[cell,axis]+=value;n[cell,axis]+=1
        return np.divide(v,n,out=np.zeros_like(v),where=n>0)
