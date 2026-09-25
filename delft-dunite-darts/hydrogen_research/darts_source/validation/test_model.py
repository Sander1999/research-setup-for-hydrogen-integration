"""Independent finite-source, geometry, native-transport and ledger checks."""
import tempfile,unittest
from pathlib import Path
import numpy as np
from darts.engines import value_vector
from ..model import Model,damage_closure,grid_geometry,source_fractions
from ..properties import load_config,MW,FE_MW,Properties
from ..run import run_case

SMALL={'nx':4,'ny':4,'nz':4,'runtime_days':30.,'max_timestep_days':2.,'snapshot_days':[0.,30.]}

class PhysicsTests(unittest.TestCase):
    def test_exact_source_volume_under_refinement(self):
        c=load_config();target=4*np.pi*np.prod(c['source_semiaxes_m'])/3
        for n in [4,8,12]:
            q=c|{'nx':n,'ny':n,'nz':n};xyz,d=grid_geometry(q);f=source_fractions(q,xyz,d)
            self.assertAlmostEqual((f@np.prod(d,axis=1))/target,1,places=10)
            self.assertTrue(np.all((f>=0)&(f<=1)))
    def test_reactive_fraction_scales_expansion(self):
        c=load_config();xi=np.ones(3);f=np.array([0.,.5,1.]);strain,stress,d=damage_closure(xi,f,c)
        np.testing.assert_allclose(strain,c['solid_expansion_fraction']*(1-c['porosity'])*c['accessible_fraction']*f)
        self.assertEqual(d[0],0);self.assertTrue(np.all((d>=0)&(d<=1)))
    def test_hydration_and_redox_stoichiometry(self):
        c=load_config();props=Properties(c)
        perkgH2=c['fe_mass_fraction']*c['ferrous_fraction']/(3*FE_MW)*MW[0]
        consumed=perkgH2/MW[0]*MW[1]*props.water_stoichiometry
        self.assertAlmostEqual(consumed,c['hydration_water_kg_per_kg_altered_rock']+perkgH2*MW[1]/MW[0],places=12)
        np.testing.assert_allclose(props.evaluate_mass_source(100,323.15,np.array([1e-5,.99,.00999])),[-1,props.water_stoichiometry,0])
    def test_harmonic_native_permeability_update(self):
        m=Model(SMALL|{'boundary_mode':'natural_flow'});m.init();before=np.asarray(m.reservoir.mesh.tran).copy();flows,_=m.face_phase_flows()
        m.damage[:]=(9/m.cfg['permeability_gain'])**(1/m.cfg['damage_permeability_exponent']);m.update_native_permeability()
        connected=(m._bm<m.nb)|(m._bp<m.nb)
        np.testing.assert_allclose(np.asarray(m.reservoir.mesh.tran)[connected],before[connected]*10,rtol=1e-12)
        after,_=m.face_phase_flows();np.testing.assert_allclose(after[connected],flows[connected]*10,rtol=1e-11,atol=1e-8)
        m.damage[:]=0;m.damage[0]=1;m.update_native_permeability()
        for i,(a,b) in enumerate(zip(m._bm,m._bp)):
            if a<m.nb and b<m.nb:
                axis=int(np.argmax(abs(m.xyz[b]-m.xyz[a])));da,db=m.dims[a,axis],m.dims[b,axis]
                ra=m.vertical_ratio[a] if axis==2 else 1.;rb=m.vertical_ratio[b] if axis==2 else 1.
                expected=before[i]*(da/(m.base_permeability[a]*ra)+db/(m.base_permeability[b]*rb))/(da/(m.permeability[a]*ra)+db/(m.permeability[b]*rb))
                self.assertAlmostEqual(np.asarray(m.reservoir.mesh.tran)[i],expected,places=12)
    def test_darcy_velocity_matches_imposed_hydraulic_gradient(self):
        from ..properties import G
        c=SMALL|{'layers':[],'host_permeability_md':10.,'source_permeability_md':10.,
                 'corridor_enabled':False,'reaction_enabled':False,'boundary_mode':'natural_flow'}
        m=Model(c);m.init();velocity=m.internal_velocity()
        expected=.00852671467191601*10.*G*m.cfg['brine_density_kg_m3']*m.cfg['hydraulic_gradient']/m.cfg['brine_viscosity_cp']
        np.testing.assert_allclose(velocity[:,0],expected,rtol=1e-7)
    def test_rejected_step_does_not_commit_fe_or_damage(self):
        m=Model(SMALL);m.init();old=m.fe_remaining.copy();d=m.damage.copy();oldgen=m.generated_kmol
        # Isolate rejection hook; source preparation still executes fully.
        m.nonlinear_solver.solve_timestep=lambda dt,t,verbose:False
        result=m.run_timestep(2.,0.);self.assertFalse(result)
        np.testing.assert_array_equal(m.fe_remaining,old);np.testing.assert_array_equal(m.damage,d);self.assertEqual(m.generated_kmol,oldgen)
    def test_native_closed_and_open_ledgers(self):
        with tempfile.TemporaryDirectory() as tmp:
            for mode in ['closed','natural_flow']:
                summary,rows=run_case(SMALL|{'boundary_mode':mode},Path(tmp)/mode,plots=False)
                self.assertLess(max(summary['max_native_relative_balance_error'].values()),2e-5)
                self.assertGreater(summary['generated_H2_kg'],0)
                self.assertLess(summary['generated_H2_kg'],summary['capacity_H2_kg'])
                if mode=='closed':self.assertEqual(summary['exported_H2_kg'],0)
    def test_intrinsic_flow_has_same_finite_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=[]
            for mode in ['closed','natural_flow']:
                summary,_=run_case(SMALL|{'boundary_mode':mode,'feedback_enabled':False,'damage_enabled':False},Path(tmp)/mode,plots=False)
                result.append(summary['generated_H2_kg'])
                self.assertAlmostEqual(summary['generated_H2_kg']/summary['intrinsic_analytical_H2_kg'],1,places=9)
            self.assertAlmostEqual(result[0]/result[1],1,places=9)
    def test_zero_source_cannot_generate_hydrogen(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary,_=run_case(SMALL|{'reaction_enabled':False},Path(tmp),plots=False)
            self.assertEqual(summary['generated_H2_kg'],0);self.assertEqual(summary['max_damage'],0)
    def test_independent_high_pressure_fluid_references(self):
        from .properties_reference import benchmarks
        result=benchmarks();self.assertEqual(len(result['rows']),21)
    def test_strict_domain_guards(self):
        for override in [{'nx':2.5},{'temperature_k':273.15},{'salinity_molal':6.1},{'boundary_mode':'pump'},{'source_semiaxes_m':[10000.,2.,3.]}]:
            with self.assertRaises(ValueError):load_config(override)

if __name__=='__main__':unittest.main()
