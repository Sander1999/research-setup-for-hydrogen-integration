"""Test the water-renewal bottleneck without changing the chemical source law.

The 100 mD deeper host is a hypothetical fractured-host end member, not a
reinterpretation of Delft geology. Sedimentary layers above 2400 m are unchanged.
Run with the registered DARTS interpreter, after the reference suite.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .properties import solubility_x
from .run import run_case


def feedback_diagnostics(folder):
    folder=Path(folder)
    cfg=json.loads((folder/'parameters.json').read_text())
    f=pd.read_csv(folder/'spatial_final.csv')
    weights=(f.source_fraction*f.dx_m*f.dy_m*f.dz_m).to_numpy()
    mean=lambda values:float(np.average(values,weights=weights))
    xs=solubility_x(f.pressure_bar.to_numpy(),cfg['temperature_k'],f.salinity_molal.to_numpy(),cfg['salt_basis_particles_per_nacl'])
    contact=(1-f.gas_saturation)**cfg['water_contact_exponent']
    inhibition=1/(1+cfg['inhibition_strength']*f.aqueous_H2_salt_free_mole_fraction/xs)
    area=1+cfg['surface_area_gain']*f.damage
    passivation=np.exp(-cfg['passivation_strength']*f.conversion_fraction_accessible_rock)
    speed=np.linalg.norm(f[['darcy_x_m_day','darcy_y_m_day','darcy_z_m_day']],axis=1)
    return {'water_contact_mean':mean(contact),'inhibition_multiplier_mean':mean(inhibition),
            'area_multiplier_mean':mean(area),'passivation_multiplier_mean':mean(passivation),
            'combined_multiplier_mean':mean(contact*inhibition*area*passivation),
            'source_brine_Darcy_speed_m_day':mean(speed),
            'source_horizontal_Darcy_speed_m_day':mean(abs(f.darcy_x_m_day)),
            'source_reaction_fraction_mean':mean(f.conversion_fraction_accessible_rock),
            'diagnostic_scope':'Final source-volume-weighted cell means; Darcy speed is not an integrated fresh-water renewal measurement.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results/flow_sensitivity'))
    parser.add_argument('--reference',type=Path,default=Path('results/darts'))
    args=parser.parse_args()
    cfg=json.loads((args.reference/'natural_flow/parameters.json').read_text())
    original=cfg['layers'][-1]['permeability_md']
    cfg['layers'][-1]['permeability_md']=100.
    result={'scope':'Hypothetical connectivity end member; no calibrated fractured host or observed Delft groundwater gradient.',
            'changed_parameter':'Deepest hypothetical host horizontal permeability',
            'reference_md':original,'alternative_md':100.,
            'unchanged':'Source inventory, kinetics, caprock, sandstone, corridor, head gradient, temperature and salinity.',
            'diagnostic_limit':'Endpoint velocities and rate multipliers diagnose response; integrated fresh-water turnover is not inferred.',
            'cases':{},'feedback':{}}
    for name in ['closed','natural_flow']:
        result['feedback'][name]=feedback_diagnostics(args.reference/name)
        s,_=run_case(cfg|{'boundary_mode':name},args.output/name)
        result['cases'][name]=s
        result['feedback'][name+'_permeable_host']=feedback_diagnostics(args.output/name)
    a,b=result['cases']['closed'],result['cases']['natural_flow']
    result['generation_change_percent']=100*(b['generated_H2_kg']/a['generated_H2_kg']-1)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'sensitivity_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(result['feedback']).T.to_csv(args.output/'feedback_diagnostics.csv')
    print(json.dumps({'generation_change_percent':result['generation_change_percent'],
                      'generated_kg':{k:v['generated_H2_kg'] for k,v in result['cases'].items()}},indent=2))


if __name__=='__main__':main()
