# An ideal hydrogen-generating body below a layered Delft sandstone analogue

This model uses native Open-DARTS to couple a finite Fe inventory, water consumption, a scalar reaction-damage law and three-dimensional H₂/brine transport. The Delft-inspired formation order is a literature-constrained analogue. The dunite body and the connecting corridor are inserted hypotheses. They are not an interpretation of an observed dunite body beneath Delft.

The model can test a conditional question: **with the same source rock, does a maintained natural water-head gradient change generation, retention or deep export?** The answer depends on the reaction feedback. With intrinsic first-order conversion alone, both boundary cases must generate the same amount. A flowing model is not allowed to create extra hydrogen merely because it transports more water.

Use a compatible, separate DARTS interpreter through `run_model.sh`; see the root README for setup. The module also runs directly from the project directory:

```
"$DARTS_PYTHON" -B -m hydrogen_research.darts_source.run --suite --output results/darts
```

The standard mesh uses 16 × 12 horizontal cells and at least 18 vertical cells, with each formation boundary honoured exactly. `--preset detailed` uses 20 × 16 and at least 24 vertical cells and halves the maximum timestep. `--preset pilot` is a quick 365-day calculation. `--config path.json` accepts documented overrides and rejects unknown parameter names. The ten-year standard duration and ten-year intrinsic reaction half-time are assumptions, not an experimentally calibrated forecast.

## Geometry and geology

The domain spans 4 × 3 km and depths 1900–2800 m. The ideal ellipsoid has centre (2000, 1500, 2580) m and semi-axes (800, 500, 180) m. It is an Fe-bearing olivine-rich endmember with 6 wt% elemental Fe in its solids, 1% accessible reactive rock, 3300 kg/m³ grain density and 5% source porosity. Subcell quadrature intersects the ellipsoid with the actual finite volumes; a bounded correction preserves the exact ellipsoid volume on every mesh.

The prescribed formation sequence comprises a Rodenrijs claystone cap, upper Delft sandstone, an internal fine-grained baffle, lower Delft sandstone, Alblasserdam underburden and a deeper hypothetical host. All numerical layer tops, thicknesses and effective flow properties in `config.json` are coarsened scenario choices. They are not measured log picks. In particular, the assumed lateral geometry is not a mapped Delft reservoir.

The optional 250 × 250 m vertical corridor exists only below the reservoir base, between depths 2300 and 2600 m. Its horizontal and vertical permeabilities are scenario inputs. Subcell corridor area is preserved using face-parallel effective conductances. `flow_blocked_path` disables it while retaining the same source and overlying formations. No artificial bypass is put through the overlying caprock.

## Source and water balance

For each cell, the initial available iron is

\[
N_{Fe,0}=V_b f_s(1-\phi_s)\rho_s w_{Fe} f_{II}f_a/M_{Fe}.
\]

Here \(f_s\) is the ellipsoid volume fraction and \(f_a\) is the fraction of source solids treated as accessible. The maximum H₂ inventory is \(N_{Fe,0}/3\), using the reference redox reaction \(3FeO+H_2O\rightarrow Fe_3O_4+H_2\). XRF elemental Fe alone does not determine either ferrous fraction or accessibility.

A timestep proposes

\[
\Delta \xi=\frac{N_{Fe}}{3}\left[1-\exp(-k\,F\,\Delta t)\right],\quad k=\ln 2/t_{1/2},
\]

with an empirical factor

\[
F=S_w^{n_c}\frac{(1+a_D D)\exp(-b\alpha)}{1+b_H x_{H_2,aq}/x_{H_2,sat}}.
\]

\(\alpha\) is conversion of the **accessible reactive rock**, not conversion of all rock. Water contact, passivation, fresh area and H₂ inhibition are bounded or regularised scenario closures. The H₂ inhibition law is not a thermodynamic affinity calculation. Its coefficient is not fitted to laboratory data. Turning `feedback_enabled` off sets \(F=1\).

The additional hydration demand is 0.15 kg water per kg of altered reactive rock. The native water sink per mole H₂ is

\[
\nu_w=1+\frac{h_w}{w_{Fe}f_{II}M_{H_2O}/(3M_{Fe})}.
\]

The first term is redox water and the second is additional hydration. Each attempted step is capped at 10% of the current cell water inventory to avoid exhausting water in the split source update. Iron consumption and damage are committed only after the nonlinear timestep converges. Rejected attempts leave both histories unchanged.

The exported solid mass gain is \(\Delta \xi(\nu_wM_{H_2O}-M_{H_2})\). This includes transferred fluid mass and is needed for time-lapse gravity: simply removing produced H₂ from density would give the wrong mass change.

## Expansion and cracking proxy

The bulk-volume expansion diagnostic is

\[
\epsilon_v=\epsilon_s(1-\phi_s)f_s f_a\alpha,
\quad \sigma_{trial}=K_{eff}\epsilon_v,
\]

and the irreversible damage target is

\[
D_{new}=\max\left[D_{old},1-\exp\left(-\frac{\max(\sigma_{trial}-\sigma_t,0)}{\sigma_D}\right)\right].
\]

The reactive fraction appears explicitly in the strain. Depleting the 1% accessible Fe inventory must not imply 40% expansion of every solid in the cell. The stiffness and tensile threshold are imposed effective parameters. The trial stress is not the result of a displacement or force-equilibrium solution. The model resolves neither fracture surfaces nor aperture, fracture energy, fault failure or measured crack propagation.

The permeability used in the next accepted flow solve is

\[
k=k_0(1+a_kD^{n_k}).
\]

It enters the **native** TPFA connection transmissibilities through the harmonic average of the two cell permeabilities and the actual half-cell lengths. Reverse connections receive the same multiplier. This is a permeability law for a continuum damage proxy, not a universal cubic fracture law.

The reference geometric pore volume is held fixed against reaction-induced displacement. A pressure-dependent native accumulation factor \(1+c_r(P-P_{ref})\) represents reversible pore compliance. The expansion-to-reference-pore-volume ratio is exported to show the size of the omitted displacement. Thus the implementation couples generation, damage, permeability and phase transport, but it is a reduced model rather than a complete chemo-poroelastic calculation. It should not be used to predict fracture aperture or subsidence.

## Flow and accounting

Native DARTS solves three component balances, H₂, H₂O and NaCl, with an H₂ gas phase and a brine phase. Water and salt are nonvolatile. Relative permeabilities, regularised capillary pressure and aqueous Fick transport are explicit scenario properties. The model is isothermal; heat released by reaction and geothermal temperature differences across the domain are not solved.

The closed case starts hydrostatic and has no external fluid connections. It is initially stationary, not constrained to remain motionless: reaction and buoyancy can produce internal flow. The open case adds a horizontal hydraulic gradient of 0.05. Each lateral aquifer connection uses an exact half-cell planar conductance. Its pressure is maintained by the assumed geological head; the model prescribes no artificial pumping. Groundwater flow is not automatically available at this magnitude in a real field.

All exports are across **deep model boundaries**. Neither surface arrival nor a flux available to hydrogen-oxidising bacteria is inferred. A separate migration/soil model must describe that missing pathway and its losses.

`history.csv` contains generated, retained, exported and imported component inventories. Native OBL accumulation and independent TPFA component fluxes supply the conservation ledger. Fluxes are integrated at the right endpoint of each accepted backward-Euler step, using the permeability that was actually used by that step. The exact-property ledger is also retained so interpolation error is visible.

`spatial_final.csv` contains absolute cell depths and dimensions, source and corridor fractions, density contrast, phase concentration, damage, permeability, Darcy velocity, and initial/final fluid mass plus reaction-related solid gain. These are the authoritative fields for the gravity calculation. Each result has a parameter snapshot, SHA-256 hash and completion status.

## Tests and comparison controls

The suite includes `closed`, `natural_flow`, `flow_no_damage`, `flow_blocked_path`, `intrinsic_closed`, `intrinsic_flow`, `zero_source` and `flow_salt_basis1`. The intrinsic pair must agree with the exact finite-source solution. The zero-source case must generate no H₂ or damage. Unit tests additionally check source volume under refinement, hydration stoichiometry, rejected-step rollback, native permeability scaling and component conservation. Timestep and grid checks are required before interpreting the comparative results.

Run numerical verification and the water-renewal audit from the project directory:

```
"$DARTS_PYTHON" -B -m hydrogen_research.darts_source.run --validate --output results/darts
"$DARTS_PYTHON" -B -m hydrogen_research.darts_source.validation.mechanisms --output results/darts
```

The delivered ten-year reference suite passed eleven native tests and its component-conservation checks. The 3,648-cell reference gives a flow-minus-closed generation difference of +0.3373%; halving the maximum timestep gives +0.3375%, and refining to 8,640 cells gives +0.2934%. The direction survives these tests, but this small conditional difference is not a field-calibrated prediction. The salt-basis alternative alone changes total generation by about −0.47%.

**Free-gas inventory is not grid converged.** Relative to the refined result, its discrepancy is 20.1% for natural flow and 15.7% for closed boundaries, exceeding the declared 10% target. Total generation differs by 0.27–0.31%. `verification.json` therefore records `all_sensitivity_targets_met: false`; a successful verification command means the checks completed, not that every sensitivity target passed. Quantitative gas accumulation and plume shape require further spatial refinement before design use.

The mechanism audit reconstructs the native brine face flux from saved states. It estimates incoming source-pore-volume renewal over ten years as 0.068 for natural flow and 0.040 for the closed case, where reaction can still draw water internally. These integrals use the saved snapshots and trapezoidal integration, so they are approximate diagnostics rather than the accepted-step conservation ledger. Water contact remains almost complete in both cases; only a small fraction of source water is renewed. This explains why an imposed regional gradient need not produce a large generation increase in this geometry. The supplementary, uncalibrated 100 mD deeper-host case is reproduced with `bash run_model.sh sensitivity` after the reference pair.

The maximum diagnostic expansion is about 5.39% of reference pore volume. Fixed reaction-related storage is therefore a material approximation even though the discrete component balances close. The accepted primary-suite states span about 195–305 bar and 2.500–2.683 mol/kg, a maximum salinity increase of 7.32%; the finer-grid checks extend the pressure range to about 194–306 bar.

Scientific provenance and the limits of the low-temperature rate evidence are recorded in `docs/source_rock_science.md` in the project root. Numerical input values without a cited measurement remain scenarios even when the governing mechanism is supported by the literature.

## High-pressure fluid reference

The deep case uses a provisional implementation of the gamma–phi correlation of Kerkache et al. (2024), Eqs.9–15 and Table3, with SI pressure conversion in the Poynting term. Fugacity follows the thermodynamic pressure integral of the Lemmon hydrogen compressibility correlation; IAPWS supplies water vapour pressure. The salt denominator convention is an explicit inference: the default counts dissociated ions, and `flow_salt_basis1` tests NaCl formula-unit counting. The two conventions are not fitted separately; this is not a claim of independently verified exact reproduction of the published model. `fluid_benchmarks.json` reports both against all retained primary observations.

The default 2.5 mol/kg brine is a scenario, not a Delft chemical analysis. Molality evolves through conservative water and salt balances. Independent recent freshwater points differ by about 2–3%; brine points near the used temperature/pressure domain differ by up to about 5% under the default convention. Larger discrepancies at 150°C are retained in the benchmark report. These checks do not calibrate geological permeability or hydrogen-generation rate.

The source model's lower bound on generated-H₂ deep export subtracts the entire initial numerical H₂ seed and boundary seed input. It therefore cannot mistake the tiny positivity floor for a physical geological source. The positivity floor is 1e-12 mole fraction, explicitly reported as an inventory.

Primary fluid references: [Kerkache et al. (2024)](https://doi.org/10.1016/j.molliq.2024.124497), [ELEGANCY brine measurements, Table 1](https://www.sintef.no/globalassets/project/elegancy/deliverables/elegancy_d2.1.6_h2-solubility-in-brine_v2.pdf), and [Wolff et al. (2026), Table 5](https://www.frontiersin.org/journals/energy-research/articles/10.3389/fenrg.2026.1919332/full).
