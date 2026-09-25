"""Reproducible native DARTS source-rock scenarios and conservative ledgers.

Run with the registered DARTS Python, from the research-project directory:
python -m hydrogen_research.darts_source.run --suite --output results/darts
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from pathlib import Path
import numpy as np
from darts.engines import redirect_darts_output,set_num_threads
from .model import Model
from .properties import COMPONENTS,MW,FE_MW


def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


class RecordedModel(Model):
    def __init__(self,parameters=None):super().__init__(parameters);self.records=[];self.progress_path=None
    def record(self):
        t=float(self.physics.engine.t);s=self.reservoir_state();q=self.cell_properties()
        if not np.isfinite(s).all() or s[:,0].min()<10 or s[:,0].max()>1000:
            raise RuntimeError('Accepted pressure outside 10–1000 bar property domain: '+str([s[:,0].min(),s[:,0].max()]))
        if q['molality'].min()<0 or q['molality'].max()>6:raise RuntimeError('Accepted salinity outside 0–6 molal domain')
        exact=self.inventory();native=self.inventory(True)
        rates=self.perforation_rates() if t else {}
        prev=self.records[-1] if self.records else None;dt=t-prev['time_days'] if prev else 0.
        total_in=np.zeros(3);total_out=np.zeros(3);gas_out=aqueous_out=0.
        for rate in rates.values():
            signed=rate['component_kg_day'];total_in+=np.maximum(-signed,0);total_out+=np.maximum(signed,0)
            gas_out+=max(rate['phase_component_kg_day'][0,0],0);aqueous_out+=max(rate['phase_component_kg_day'][1,0],0)
        prevval=lambda name:prev[name] if prev else 0.
        row={'time_days':t,'pressure_min_bar':float(s[:,0].min()),'pressure_max_bar':float(s[:,0].max()),
            'generated_H2_kg':self.generated_kmol*MW[0],'unreacted_accessible_Fe_kg':float(self.fe_remaining.sum()*FE_MW),
            'free_H2_kg':exact['free_H2_kg'],'dissolved_H2_kg':exact['dissolved_H2_kg'],
            'gas_H2_out_kg':prevval('gas_H2_out_kg')+dt*gas_out,
            'aqueous_H2_out_kg':prevval('aqueous_H2_out_kg')+dt*aqueous_out,
            'max_gas_saturation':float(q['sat'][:,0].max()),'salinity_min_molal':float(q['molality'].min()),'salinity_max_molal':float(q['molality'].max()),
            'max_damage':float(self.damage.max()),'source_volume_mean_damage':float((self.damage*self.source_fraction)@self.cell_volume/(self.source_fraction@self.cell_volume)),
            'max_bulk_expansion_strain':float(self.strain.max()),'max_trial_stress_mpa':float(self.stress.max()),
            'max_permeability_used_md':float(self.permeability.max()),
            'reaction_water_consumed_kg':self.generated_kmol*MW[1]*self.property_container.water_stoichiometry,
            'H2_out_rate_kg_day':float(total_out[0]),'H2O_in_rate_kg_day':float(total_in[1])}
        reaction=np.array([1.,-self.property_container.water_stoichiometry,0.])*self.generated_kmol*MW
        for j,name in enumerate(COMPONENTS):
            row[name+'_kg']=float(exact['component_kg'][j]);row['native_'+name+'_kg']=float(native['component_kg'][j])
            row[name+'_in_kg']=prevval(name+'_in_kg')+dt*total_in[j]
            row[name+'_out_kg']=prevval(name+'_out_kg')+dt*total_out[j]
            first=self.records[0] if self.records else row
            ledger=row[name+'_out_kg']-row[name+'_in_kg']-reaction[j]
            row[name+'_balance_kg']=row['native_'+name+'_kg']-first['native_'+name+'_kg']+ledger
            row[name+'_exact_balance_kg']=row[name+'_kg']-first[name+'_kg']+ledger
        self.records.append(row)
        if self.progress_path and len(self.records)%20==0:
            self.progress_path.write_text(json.dumps({'status':'running','time_days':t,'target_days':self.cfg['runtime_days'],'accepted_steps':len(self.records)-1})+'\n')
    def after_converged_timestep(self):super().after_converged_timestep();self.record()


def spatial_rows(m):
    q=m.cell_properties();s=m.reservoir_state();v=m.internal_velocity()
    h2conc=np.sum(q['sat']*q['rho_m']*q['x'][:,:,0],axis=1)*MW[0]
    fluid_mass=m.cell_fluid_mass()
    solid_gain=(m.fe_initial-m.fe_remaining)/3*(m.property_container.water_stoichiometry*MW[1]-MW[0])
    k_next=m.base_permeability*(1+m.cfg['permeability_gain']*m.damage**m.cfg['damage_permeability_exponent'])
    return [{'cell':i,'x_m':m.xyz[i,0],'y_m':m.xyz[i,1],'z_depth_m':m.xyz[i,2],
        'dx_m':m.dims[i,0],'dy_m':m.dims[i,1],'dz_m':m.dims[i,2],
        'source_fraction':m.source_fraction[i],'corridor_fraction':m.corridor_fraction[i],
         'density_contrast_kg_m3':m.source_fraction[i]*((1-m.cfg['porosity'])*m.cfg['rock_density_kg_m3']+m.cfg['porosity']*m.cfg['brine_density_kg_m3']-m.host_density[i]),
        'layer_id':int(m.layer_ids[i]),'porosity_reference':m.porosity[i],
        'host_bulk_density_kg_m3':m.host_density[i],
        'fluid_mass_initial_kg':m._initial_fluid_mass[i],
        'fluid_mass_final_kg':fluid_mass[i],
        'solid_mass_gain_kg':solid_gain[i],
        'pressure_bar':s[i,0],'H2_kg_m3_pore_fluid':h2conc[i],'gas_saturation':q['sat'][i,0],
        'aqueous_H2_salt_free_mole_fraction':q['x'][i,1,0]/(q['x'][i,1,0]+q['x'][i,1,1]),
        'salinity_molal':q['molality'][i],'conversion_fraction_accessible_rock':m.progress()[i],
        'damage':m.damage[i],'bulk_expansion_strain_diagnostic':m.strain[i],
        'permeability_used_md':m.permeability[i],'permeability_next_step_md':k_next[i],
        'darcy_x_m_day':v[i,0],'darcy_y_m_day':v[i,1],'darcy_z_m_day':v[i,2]} for i in range(m.nb)]


def code_fingerprint():
    h=hashlib.sha256()
    for name in ['model.py','properties.py','run.py','config.json']:
        h.update(name.encode());h.update(Path(__file__).with_name(name).read_bytes())
    return h.hexdigest()


def layer_inventories(m):
    q=m.cell_properties();phase=m.current_pore_volume()[:,None,None]*(q['sat']*q['rho_m'])[:,:,None]*q['x']*MW
    rows=[]
    for i,layer in enumerate(m.cfg['layers']):
        mask=m.layer_ids==i
        rows.append({'layer_id':i,'name':layer['name'],'retained_free_H2_kg':float(phase[mask,0,0].sum()),
                     'retained_dissolved_H2_kg':float(phase[mask,1,0].sum()),
                     'generated_H2_kg':float(((m.fe_initial-m.fe_remaining)/3*MW[0])[mask].sum())})
    return rows


def run_case(parameters=None,output=None,*,plots=True):
    output=Path(output or 'results/darts/natural_flow').resolve();output.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();status=output/'status.json';status.write_text('{"status":"running"}\n')
    try:
        redirect_darts_output(str(output/'native.log'));set_num_threads(2)
        m=RecordedModel(parameters);m.progress_path=status;m.init();m.set_output(output_folder=str(output));m.record()
        config_text=json.dumps(m.cfg,sort_keys=True,indent=2)+'\n';(output/'parameters.json').write_text(config_text)
        milestones=sorted(set([float(t) for t in m.cfg['snapshot_days'] if 0<=t<=m.cfg['runtime_days']]+[0.,float(m.cfg['runtime_days'])]))
        write_csv(output/'spatial_initial.csv',spatial_rows(m))
        for stop in milestones[1:]:
            m.run(stop-float(m.physics.engine.t),verbose=0)
            if abs(m.physics.engine.t-stop)>1e-7:raise RuntimeError('DARTS stopped early')
            write_csv(output/f'spatial_day_{stop:g}.csv',spatial_rows(m))
        rows=m.records;first,last=rows[0],rows[-1];balance={};exact_balance={}
        for name in COMPONENTS:
            scale=max(first[name+'_kg']+last[name+'_in_kg']+(last['generated_H2_kg'] if name=='H2' else 0),1e-8)
            balance[name]=max(abs(r[name+'_balance_kg']) for r in rows)/scale
            exact_balance[name]=max(abs(r[name+'_exact_balance_kg']) for r in rows)/scale
        capacity=m.fe_initial_kmol/3*MW[0]
        ideal=capacity*(-np.expm1(-np.log(2)*m.cfg['runtime_days']/m.cfg['half_time_days'])) if m.cfg['reaction_enabled'] else 0
        solid_oxygen=last['reaction_water_consumed_kg']-last['generated_H2_kg']
        summary={'scope':'Ideal 3D olivine-rich source body: native multiphase flow, finite Fe reaction, empirical damage/permeability. Scenario, not a calibrated field prediction.',
            'configuration_sha256':hashlib.sha256(config_text.encode()).hexdigest(),'code_sha256':code_fingerprint(),'parameters':m.cfg,
            'runtime_seconds':time.monotonic()-start,'days':last['time_days'],'cells':m.nb,'actual_nz':m.nz,
            'bulk_volume_m3':m.bulk_volume,'source_bulk_volume_m3':float(m.source_fraction@m.cell_volume),
            'source_rock_mass_kg':float(m.source_rock_mass.sum()),'reactive_source_rock_mass_kg':float(m.source_rock_mass.sum()*m.cfg['accessible_fraction']),
            'initial_reference_pore_volume_m3':float(m.cell_volume@m.porosity),
            'capacity_H2_kg':capacity,'generated_H2_kg':last['generated_H2_kg'],
            'retained_free_H2_kg':last['free_H2_kg'],'retained_dissolved_H2_kg':last['dissolved_H2_kg'],
            'exported_H2_kg':last['H2_out_kg'],'gas_H2_exported_kg':last['gas_H2_out_kg'],'aqueous_H2_exported_kg':last['aqueous_H2_out_kg'],
            'numerical_initial_H2_kg':first['H2_kg'],'numerical_boundary_H2_in_kg':last['H2_in_kg'],
            'reaction_water_consumed_kg':last['reaction_water_consumed_kg'],'solid_mass_gain_kg':solid_oxygen,
            'water_in_kg':last['H2O_in_kg'],'water_out_kg':last['H2O_out_kg'],
            'max_native_relative_balance_error':balance,'max_exact_relative_balance_error':exact_balance,
            'native_H2_balance_error_kg':last['H2_balance_kg'],
            'generated_H2_export_lower_bound_kg':max(0.,last['H2_out_kg']-first['H2_kg']-last['H2_in_kg']),
            'layer_inventories':layer_inventories(m),
            'intrinsic_analytical_H2_kg':ideal,'max_damage':last['max_damage'],
            'source_volume_mean_damage':last['source_volume_mean_damage'],
            'max_diagnostic_expansion_to_reference_pore_volume_ratio':float(np.max(m.strain/m.porosity)),
            'total_diagnostic_expansion_m3':float(m.strain@m.cell_volume),
            'final_pressure_range_bar':[last['pressure_min_bar'],last['pressure_max_bar']],
            'accepted_timesteps':len(rows)-1,
            'coupling':'Accepted-step Lie splitting. Fe and damage commit only after convergence; rate, hydration sink and native harmonic TPFA transmissibility are frozen during each attempted timestep. Ledger integrates the exact transmissibility used by that accepted solve.',
            'storage_limit':'Reference geometric porosity remains fixed. Pressure-dependent native accumulation includes reversible pore compliance. Reaction expansion is a damage diagnostic only; its displacement of pore water and a mechanical displacement/stress solution are omitted.',
            'energy_boundary':'Head boundaries represent an assumed naturally maintained aquifer gradient. No artificial pumping is prescribed. Geological head maintenance and actual water flux require field evidence.',
            'export_limit':'Export is across deep model boundaries, not arrival at the surface or at hydrogen-oxidising bacteria. No migration efficiency to the surface is inferred.',
            'rate_limit':'At 80 C the 10-year half-time is an assumed sensitivity, not a rate fitted to high-temperature experiments or the XRF composition.'}
        if max(balance.values())>2e-5:raise RuntimeError('Native component balance failed: '+str(balance))
        if not m.cfg['feedback_enabled'] and abs(summary['generated_H2_kg']-ideal)>max(capacity*1e-9,1e-8):
            raise RuntimeError('Intrinsic finite-source integration failed or water limiting became active')
        if summary['generated_H2_kg']>capacity*(1+1e-10):raise RuntimeError('Finite Fe capacity exceeded')
        write_csv(output/'history.csv',rows);fields=spatial_rows(m);write_csv(output/'spatial_final.csv',fields)
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        if plots:plot_case(m,rows,fields,output)
        status.write_text(json.dumps({'status':'complete','days':last['time_days'],'runtime_seconds':time.monotonic()-start})+'\n')
        return summary,rows
    except BaseException as e:
        status.write_text(json.dumps({'status':'failed','error':str(e)})+'\n');raise


def plot_case(m,rows,fields,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=np.array([r['time_days'] for r in rows])/365.25
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for key,label in [('generated_H2_kg','Generated'),('dissolved_H2_kg','Retained dissolved'),('free_H2_kg','Retained gas'),('H2_out_kg','Deep-boundary export')]:
        axes[0].plot(t,[r[key]/1000 for r in rows],label=label)
    axes[0].set_ylabel('Hydrogen (tonne)');axes[0].legend(fontsize=8)
    axes[1].plot(t,[r['max_damage'] for r in rows],label='Maximum');axes[1].plot(t,[r['source_volume_mean_damage'] for r in rows],label='Source-volume mean')
    axes[1].set_ylabel('Regularised damage proxy');axes[1].legend(fontsize=8)
    for key,label in [('pressure_min_bar','Minimum'),('pressure_max_bar','Maximum')]:axes[2].plot(t,[r[key] for r in rows],label=label)
    axes[2].set_ylabel('Pressure (bar absolute)');axes[2].legend(fontsize=8)
    for ax in axes:ax.set_xlabel('Elapsed time (year)');ax.grid(alpha=.2)
    fig.suptitle('Ideal source rock | '+m.cfg['boundary_mode'].replace('_',' ')+' | uncalibrated 80°C rate and damage assumptions')
    fig.savefig(output/'history.png',dpi=180);plt.close(fig)
    fig=plt.figure(figsize=(15,6),layout='constrained')
    xyz=m.xyz;source=m.source_fraction>0
    for panel,(key,label,cmap) in enumerate([('H2_kg_m3_pore_fluid','H₂ (kg/m³ pore fluid)','viridis'),('damage','Damage proxy (dimensionless)','magma'),('permeability_next_step_md','Permeability after accepted damage (mD)','cividis')],1):
        ax=fig.add_subplot(1,3,panel,projection='3d');values=np.array([f[key] for f in fields])
        mask=source if panel>1 else ((xyz[:,1]<m.cfg['length_y_m']/2)|source)
        dots=ax.scatter(*xyz[mask].T,c=values[mask],cmap=cmap,s=14,alpha=.8,depthshade=False)
        fig.colorbar(dots,ax=ax,shrink=.55,pad=.05,label=label)
        ax.set(xlabel='x (m)',ylabel='y (m)',zlabel='Depth (m)',zlim=(m.cfg['top_depth_m']+m.height,m.cfg['top_depth_m']))
        ax.set_box_aspect((m.cfg['length_x_m'],m.cfg['length_y_m'],m.height));ax.view_init(elev=21,azim=-62)
        ax.tick_params(labelsize=7);ax.set_title(['Fluid concentration (front cutaway)','Source cells','Source cells'][panel-1],fontsize=10)
    fig.suptitle(f'Native 3D fields at {m.cfg["runtime_days"]/365.25:.2f} years; {m.nb:,} cells\nThe damage field is a continuum proxy, not a resolved fracture network')
    fig.savefig(output/'source_rock_3d.png',dpi=190);plt.close(fig)
    # A physical x-z cross-section shows transport direction without depicting
    # unresolved cracks or drawing artificial surface migration paths.
    jy=m.cfg['ny']//2;indices=np.array([i for i in range(m.nb) if (i//m.cfg['nx'])%m.cfg['ny']==jy])
    xx=xyz[indices,0].reshape(m.nz,m.cfg['nx']);zz=xyz[indices,2].reshape(m.nz,m.cfg['nx'])
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,key,label in [(axes[0],'H2_kg_m3_pore_fluid','H₂ (kg/m³ pore fluid)'),(axes[1],'pressure_bar','Pressure (bar)')]:
        v=np.array([f[key] for f in fields])[indices].reshape(xx.shape)
        im=ax.pcolormesh(xx,zz,v,shading='nearest');fig.colorbar(im,ax=ax,label=label)
        ux=np.array([f['darcy_x_m_day'] for f in fields])[indices].reshape(xx.shape)
        uz=np.array([f['darcy_z_m_day'] for f in fields])[indices].reshape(xx.shape)
        ax.quiver(xx,zz,ux,uz,color='white',width=.003)
        ax.invert_yaxis();ax.set(xlabel='x (m)',ylabel='Depth (m)')
    fig.suptitle('Central cross-section; arrows show native brine Darcy velocity (common scale within each panel)')
    fig.savefig(output/'transport_section.png',dpi=180);plt.close(fig)


def plot_comparison(root,summaries):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=list(summaries);x=np.arange(len(names));fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for offset,key,label in [(-.3,'generated_H2_kg','Generated'),(-.1,'retained_dissolved_H2_kg','Dissolved'),(.1,'retained_free_H2_kg','Free gas'),(.3,'exported_H2_kg','Deep export')]:
        axes[0].bar(x+offset,[summaries[n][key]/1000 for n in names],width=.2,label=label)
    axes[0].set_ylabel('Hydrogen (tonne)');axes[0].legend(fontsize=8)
    axes[1].bar(x,[summaries[n]['source_volume_mean_damage'] for n in names],color='#496f82');axes[1].set_ylabel('Source-volume mean damage proxy')
    for ax in axes:ax.set_xticks(x,names,rotation=25,ha='right');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Identical ideal source-rock inventory | conditional rate feedback and deep transport')
    fig.savefig(root/'comparison.png',dpi=180);plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--suite',action='store_true');p.add_argument('--validate',action='store_true');p.add_argument('--config',type=Path);p.add_argument('--output',type=Path,default=Path('results/darts'));p.add_argument('--no-plots',action='store_true')
    p.add_argument('--preset',choices=['pilot','standard','detailed'],default='standard');args=p.parse_args()
    if args.validate:
        from .validation.convergence import run_validation
        print(json.dumps(run_validation(args.output),indent=2));return
    override=json.loads(args.config.read_text()) if args.config else {}
    if args.preset=='pilot':override={'nx':6,'ny':4,'nz':4,'runtime_days':365.,'max_timestep_days':10.,'snapshot_days':[0.,365.]}|override
    elif args.preset=='detailed':override={'nx':20,'ny':16,'nz':24,'max_timestep_days':5.}|override
    cases=[('closed',{'boundary_mode':'closed'}),('natural_flow',{'boundary_mode':'natural_flow'}),
           ('flow_no_damage',{'boundary_mode':'natural_flow','damage_enabled':False}),
           ('flow_blocked_path',{'boundary_mode':'natural_flow','corridor_enabled':False}),
           ('intrinsic_closed',{'boundary_mode':'closed','feedback_enabled':False,'damage_enabled':False}),
           ('intrinsic_flow',{'boundary_mode':'natural_flow','feedback_enabled':False,'damage_enabled':False}),
           ('zero_source',{'boundary_mode':'natural_flow','reaction_enabled':False}),
           ('flow_salt_basis1',{'boundary_mode':'natural_flow','salt_basis_particles_per_nacl':1.})] if args.suite else [(override.get('boundary_mode','natural_flow'),{})]
    summaries={}
    for name,extra in cases:
        result,_=run_case(override|extra,args.output/name,plots=not args.no_plots);summaries[name]=result
        print(json.dumps({'case':name,'generated_H2_kg':result['generated_H2_kg'],'retained_dissolved_H2_kg':result['retained_dissolved_H2_kg'],'exported_H2_kg':result['exported_H2_kg'],'seconds':result['runtime_seconds']}),flush=True)
    (args.output/'suite_summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    if args.suite and not args.no_plots:plot_comparison(args.output,summaries)

if __name__=='__main__':main()
