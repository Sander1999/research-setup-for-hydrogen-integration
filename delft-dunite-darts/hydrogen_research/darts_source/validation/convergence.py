"""Native verification, including resolution of the FLOW MINUS CLOSED effect."""
from pathlib import Path
import csv,json,math,unittest
from ..run import run_case,code_fingerprint
from .properties_reference import benchmarks


def run_validation(root):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    suite=unittest.defaultTestLoader.loadTestsFromName('hydrogen_research.darts_source.validation.test_model')
    test_result=unittest.TextTestRunner(verbosity=1).run(suite)
    if not test_result.wasSuccessful():raise RuntimeError('Native unit or conservation tests failed')
    benchmarks(root/'fluid_benchmarks.json')
    base={}
    for case in ['natural_flow','closed']:
        p=root/case/'summary.json'
        if p.exists():
            base[case]=json.loads(p.read_text())
            if base[case]['code_sha256']!=code_fingerprint():raise RuntimeError('Baseline results have a stale code fingerprint; rerun the suite')
        else:base[case],_=run_case({'boundary_mode':case},root/case)
    results={};refined={};metrics=['generated_H2_kg','retained_free_H2_kg','retained_dissolved_H2_kg','exported_H2_kg']
    for name in ['timestep_refinement','grid_refinement']:
        refined[name]={}
        for case in ['natural_flow','closed']:
            c=base[case]['parameters']
            cfg=c|({'max_timestep_days':c['max_timestep_days']/2} if name=='timestep_refinement' else {'nx':max(20,math.ceil(1.25*c['nx'])),'ny':max(16,math.ceil(1.25*c['ny'])),'nz':max(24,math.ceil(1.25*c['nz']))})
            output_name=name if case=='natural_flow' else name+'_closed'
            cached=root/output_name/'summary.json';state=root/output_name/'status.json'
            result=json.loads(cached.read_text()) if cached.exists() else None
            complete=state.exists() and json.loads(state.read_text()).get('status')=='complete'
            if not (result and complete and result.get('code_sha256')==code_fingerprint() and result['parameters']==cfg):
                result,_=run_case(cfg,root/output_name)
            refined[name][case]=result;change={}
            for metric in metrics:
                reference=base[case][metric];absolute=abs(result[metric]-reference)
                relative=absolute/max(abs(result[metric]),abs(reference),1e-30)
                floor=max(base[case]['generated_H2_kg'],result['generated_H2_kg'])*1e-5
                limit=.03 if name=='timestep_refinement' else .10
                change[metric]={'baseline':reference,'refined':result[metric],
                    'relative_difference':relative,'absolute_difference_kg':absolute,
                    'absolute_resolution_floor_kg':floor,'relative_target':limit,
                    'within_target':bool(absolute<=max(limit*max(abs(result[metric]),abs(reference)),floor))}
            results[output_name]={'metrics':change,'all_within_target':all(v['within_target'] for v in change.values()),
                'cells':result['cells'],'max_native_relative_balance_error':result['max_native_relative_balance_error']}
    effects={}
    # Keep the paired difference explicit: checking each million-kg inventory
    # against a 10% target cannot resolve a much smaller flow effect.
    for name,pair in [('baseline',base)]+list(refined.items()):
        value=pair['natural_flow']['generated_H2_kg']-pair['closed']['generated_H2_kg']
        effects[name]={'flow_minus_closed_generated_H2_kg':value,
            'fraction_of_closed_generation':value/pair['closed']['generated_H2_kg']}
    b=effects['baseline']['flow_minus_closed_generated_H2_kg']
    deltas=[abs(effects[n]['flow_minus_closed_generated_H2_kg']-b) for n in refined]
    same_sign=all(e['flow_minus_closed_generated_H2_kg']*b>0 for e in effects.values()) if b else False
    magnitude_resolved=bool(same_sign and max(deltas)<.5*abs(b))
    report={'code_sha256':code_fingerprint(),'unit_tests_run':test_result.testsRun,'unit_tests_passed':True,
        'native_conservation_threshold':2e-5,'timestep_target':.03,'grid_target':.10,
        'tiny_output_absolute_floor_fraction_of_generation':1e-5,
        'fluid_benchmarks_passed':True,'comparisons':results,
        'flow_minus_closed':{'results':effects,'same_sign_in_all_tested_pairs':same_sign,
            'maximum_change_in_effect_kg':max(deltas),'effect_change_less_than_half_baseline_effect':magnitude_resolved,
            'interpretation':'A consistent tested numerical effect remains conditional on uncalibrated kinetics, damage, geology and boundary assumptions. If the paired effect changes by at least half its baseline magnitude, do not report a resolved quantitative flow advantage.'},
        'scope':'Numerical verification only. Grid sensitivity targets concern integrated outputs, not a resolved fracture geometry or field validation.',
        'all_sensitivity_targets_met':all(r['all_within_target'] for r in results.values()),
        'accepted_state_domains':accepted_state_domains(root)}
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def accepted_state_domains(root):
    result={}
    for path in sorted(Path(root).glob('*/history.csv')):
        rows=list(csv.DictReader(path.open()))
        if not rows:continue
        initial=float(rows[0]['salinity_min_molal'])
        low=min(float(r['salinity_min_molal']) for r in rows);high=max(float(r['salinity_max_molal']) for r in rows)
        result[path.parent.name]={'pressure_min_bar':min(float(r['pressure_min_bar']) for r in rows),
            'pressure_max_bar':max(float(r['pressure_max_bar']) for r in rows),
            'salinity_min_molal':low,'salinity_max_molal':high,
            'maximum_absolute_salinity_drift_molal':max(abs(low-initial),abs(high-initial)),
            'maximum_fractional_salinity_drift':max(abs(low-initial),abs(high-initial))/initial}
    return result
