"""Reconstruct native snapshot fluxes to diagnose water renewal and rate factors.

No simulation is rerun and no core parameter is changed. Cumulative exchange
uses four saved snapshots and trapezoids, so is a diagnostic approximation,
not the accepted-timestep component conservation ledger.
"""
from pathlib import Path
import csv,json
import numpy as np
from darts.engines import value_vector,redirect_darts_output,set_num_threads
from ..model import Model
from ..properties import MW,hydrogen_density,solubility_x
from ..run import code_fingerprint


def reconstruct(case_dir,filename):
    case_dir=Path(case_dir);c=json.loads((case_dir/'parameters.json').read_text())
    rows=list(csv.DictReader((case_dir/filename).open()))
    val=lambda k:np.array([float(r[k]) for r in rows])
    m=Model(c);m.init();p=val('pressure_bar');sg=val('gas_saturation');xs=val('aqueous_H2_salt_free_mole_fraction');sal=val('salinity_molal')
    h_over_w=xs/(1-xs);s_over_w=sal*MW[1]/1000
    xw=np.column_stack((h_over_w,np.ones(len(rows)),s_over_w));xw/=xw.sum(axis=1)[:,None]
    rw=c['brine_density_kg_m3']*(1+c['brine_compressibility_per_bar']*(p-c['producer_bhp_bar']))/(xw@MW)
    rg=np.array([hydrogen_density(float(pi),c['temperature_k']) for pi in p])/MW[0]
    nwater=(1-sg)*rw;ngas=sg*rg
    amounts=nwater[:,None]*xw;amounts[:,0]+=ngas;z=amounts/amounts.sum(axis=1)[:,None]
    states=np.array(m.physics.engine.X).reshape(-1,3)
    states[:m.nb,0]=p;states[:m.nb,1:]=z[:,:2];m.physics.engine.X=value_vector(states.ravel())
    multiplier=val('permeability_used_md')/m.base_permeability
    m.damage=(np.maximum(multiplier-1,0)/c['permeability_gain'])**(1/c['damage_permeability_exponent'])
    m.update_native_permeability()
    velocity=m.internal_velocity();stored=np.column_stack([val('darcy_'+axis+'_m_day') for axis in 'xyz'])
    scale=max(float(np.max(abs(stored))),1e-10)
    reconstruction_error=float(np.max(abs(velocity-stored))/scale)
    if reconstruction_error>1e-4:raise RuntimeError('Snapshot-native velocity reconstruction failed: '+str(reconstruction_error))
    flow,_=m.face_phase_flows();a,b=m._bm,m._bp;mask=(a<m.nb)&(b<m.nb)&(a<b)
    weighted=flow[mask,1]*(m.source_fraction[b[mask]]-m.source_fraction[a[mask]])
    qin=float(np.maximum(weighted,0).sum());qout=float(np.maximum(-weighted,0).sum())
    f=m.source_fraction;vol=m.cell_volume;weight=f*vol;normal=weight.sum();alpha=val('conversion_fraction_accessible_rock');damage=val('damage')
    contact=(1-sg)**c['water_contact_exponent']
    inhibition=1/(1+c['inhibition_strength']*xs/solubility_x(p,c['temperature_k'],sal,c['salt_basis_particles_per_nacl']))
    area=1+c['surface_area_gain']*damage;passivation=np.exp(-c['passivation_strength']*alpha)
    factor=contact*inhibition*area*passivation
    avg=lambda v:float(v@weight/normal)
    remaining_weight=weight*(1-alpha)
    return {'source_reference_pore_volume_m3':float((weight*m.porosity).sum()),
        'water_in_to_fractional_source_m3_day':qin,'water_out_of_fractional_source_m3_day':qout,
        'source_volume_mean_darcy_speed_m_day':avg(np.linalg.norm(velocity,axis=1)),
        'source_volume_mean_water_contact':avg(contact),'source_volume_mean_inhibition':avg(inhibition),
        'source_volume_mean_area_factor':avg(area),'source_volume_mean_passivation':avg(passivation),
        'source_volume_mean_rate_factor':avg(factor),
        'remaining_Fe_weighted_rate_factor':float(factor@remaining_weight/remaining_weight.sum()),
        'source_volume_mean_H2_saturation_ratio':avg(xs/solubility_x(p,c['temperature_k'],sal,c['salt_basis_particles_per_nacl'])),
        'velocity_reconstruction_relative_error':reconstruction_error}


def audit(root):
    root=Path(root);redirect_darts_output(str(root/'mechanism_audit_native.log'));set_num_threads(2)
    cases={}
    for name in ['closed','natural_flow']:
        folder=root/name;summary=json.loads((folder/'summary.json').read_text())
        if summary['code_sha256']!=code_fingerprint():raise RuntimeError('Stale case: '+name)
        snapshot_map={0.:'spatial_initial.csv'}
        snapshot_map.update({float(p.stem.removeprefix('spatial_day_')):p.name for p in folder.glob('spatial_day_*.csv')})
        snapshot_map[float(summary['days'])]='spatial_final.csv'
        snapshots=list(snapshot_map.items())
        records=[]
        for t,path in sorted(snapshots):records.append({'time_days':t}|reconstruct(folder,path))
        t=np.array([r['time_days'] for r in records]);incoming=np.array([r['water_in_to_fractional_source_m3_day'] for r in records]);outgoing=np.array([r['water_out_of_fractional_source_m3_day'] for r in records])
        vin=float(np.trapezoid(incoming,t));vout=float(np.trapezoid(outgoing,t));pv=records[0]['source_reference_pore_volume_m3']
        cases[name]={'snapshots':records,'approximate_cumulative_water_in_m3':vin,'approximate_cumulative_water_out_m3':vout,
            'approximate_incoming_source_pore_volumes':vin/pv,'approximate_outgoing_source_pore_volumes':vout/pv}
    result={'code_sha256':code_fingerprint(),'cases':cases,
        'method':'Native brine face fluxes reconstructed from accepted snapshot states and the actual permeability used for each saved state. For a face a→b, q*(f_source,b−f_source,a) enters the fractional-source control volume; positive/negative parts distinguish inward/outward exchange. Face pairs are counted once.',
        'integration_limit':'Cumulative water exchange uses only saved time0,365,1825,3650day snapshots and trapezoidal integration. It is a coarse water-renewal diagnostic, not the exact timestep conservation ledger and not net geogenic H2 export.',
        'interpretation':'Hydration can draw water toward the source even with closed exterior boundaries. A higher Darcy speed does not guarantee complete source-water replacement; fast flow in the upper sandstone can bypass the source. The rate law changes through contact, dissolved-H2 inhibition, passivation and damage-related area, not through a direct flow-rate multiplier.'}
    (root/'mechanism_audit.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path('results/darts'));args=parser.parse_args()
    print(json.dumps(audit(args.output),indent=2))
