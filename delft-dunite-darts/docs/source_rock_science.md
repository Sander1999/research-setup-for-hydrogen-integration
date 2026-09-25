# Source-rock model: evidence and assumptions

The executable source-rock model is an **idealised scenario experiment**, not a calibrated prediction of hydrogen production at a particular site. A dunite body, its depth and density contrast are geological hypotheses. The native Delft analogue uses an isothermal temperature of 80 °C, an initial reference pressure of 265.78 bar at 2500 m depth and hydrostatic variation through a 1900–2800 m domain. These are explicit scenario choices, informed by the local temperature/depth scale; they do not reconstruct an actual source or measured natural flow. Gravity may constrain geometry and density contrast, but does not supply a reaction rate or a hydrogen yield.

## What is defensible at 80 °C

There is no verified rate calibration in the supplied literature for the proposed 80 °C, deep pressurised body. Neubeck et al. measured approximately stable headspace H2 at 30 and 50 °C, and increasing H2 at 70 °C in carbonate-rich batch experiments. The 70 °C amount rose from 2.65 ± 0.44 to 3.79 ± 0.35 **nmol per gram of olivine** between days 43 and 315. The initial autoclaving-related pulse must not be treated as sustained generation. These amounts describe the experiment and include possible loss processes; they are not an intrinsic mineral reaction constant. [Neubeck et al. (2014), section 3.3](https://liu.diva-portal.org/smash/get/diva2:732907/FULLTEXT01.pdf).

Mayhew et al. show that low-temperature H2 production depends on mineral assemblage and surface reactions, including spinel-associated processes, in experiments at 55 and 100 °C. That evidence supports testing mechanism sensitivity, but does not establish a universal dunite rate at 80 °C. [Mayhew et al. (2013)](https://www.nature.com/articles/ngeo1825).

Consequently, every kinetic timescale, reactive area, inhibition constant and cracking threshold used at 80 °C must retain `assumed_scenario` status. Vary these inputs independently; report response surfaces or ranges instead of a single field forecast. A high-temperature calibration must not be silently extrapolated.

## Chemistry and finite inventories

Let `r` denote moles of a declared reaction extent per **bulk m³ per second**, with signed molar coefficients `nu_i`. Conservation requires

```
source_i = nu_i * r                         [mol / bulk m³ / s]
delta_solid_volume = sum_solid(nu_i * Vbar_i) * delta_extent
```

Two atom-balanced end-member bookkeeping reactions are

```
2 Mg2SiO4 + 3 H2O -> Mg3Si2O5(OH)4 + Mg(OH)2
3 Fe2SiO4 + 2 H2O -> 2 Fe3O4 + 3 SiO2 + 2 H2
```

These are idealised stoichiometric limits, not an assertion that this product assemblage forms at 80 °C. Mg-only hydration produces no H2. In the magnetite-only iron oxidation limit, H2/Fe is 1/3 mol/mol; the absolute electron-balance ceiling is 1/2 mol H2 per mol initially ferrous iron if all of it becomes ferric. An actual yield must specify iron partition among residual olivine, serpentine, brucite and oxides. Do not add an H2 source independently of the consumed Fe and water. Stop reaction on exhaustion of either limiting inventory, without clipping material after it has already been generated.

Fauguerolles et al. demonstrate that changing H2 mobility changes product mineralogy and redox state. Their 300/350 °C, 50 MPa capsule experiments compare H2-permeable AgPd with effectively impermeable Au. Magnetite abundance ceases to quantify total progress when hematite forms; ferric serpentine also complicates any fixed magnetite-to-H2 proxy. These experiments do not calibrate a reservoir flow-rate or inhibition constant. [Fauguerolles et al. (2024), sections 5.2.1–5.2.2](https://ejm.copernicus.org/articles/36/555/2024/).

## Kinetics: useful forms and their limits

Evans et al. use a transparent continuum closure (their Eq. 53):

```
R_mass = M_product * k_surface * a0 * phi * phi_ol * A/(Rg*T)
```

`R_mass` is kg product/(bulk m³ s); `k_surface` is mol/(m² s), `a0` is m²/bulk m³, `phi` and `phi_ol` are volume fractions, and `A` is reaction affinity in J/mol. They hold affinity constant for their chosen model. The availability factor makes reaction vanish when water-filled porosity or olivine is exhausted. Their porosity mass balance is

```
dphi/dt = (1-phi)*div(v_s) - sum_solid(nu_mass_i/rho_i)*R_mass
```

Solid-volume increase can coexist with **net contraction of solid plus consumed fluid**. Their microcrack closure is `k = C*w²*phi` (Eq. 44); a particular constant-aspect-ratio geometry gives a phi-squared scaling. These choices require microstructure assumptions and are not universal calibrations. [Evans et al. (2018), Eqs. 9, 44, 53](https://www.ldeo.columbia.edu/sites/default/files/u988/2018%20Evans%20et%20al%202018%20JGR%20reaction%20driven%20volume%20change.pdf).

For a separate high-temperature benchmark, Malvoisin et al. give

```
f(T) = 808.3 * exp(-3640/T) * (1-exp(-8759*(1/T-1/623.6)))
z(t,T,R) = 1-exp(-k*f(T)*(t/R²)^n)
```

Here T is kelvin, t is **hours**, R is initial grain radius in **micrometres**; k and n are grain-size-regime fits, not SI constants. Appendix A restricts the temperature function to 513–623.6 K, derived from hydrothermal experiments at 500 bar. Fits exclude grains below 5 µm; the larger-grain case requires multiple regimes. Reaction progress was corrected for an initial experimental stage and rescaled for armouring-limited completion. This law is unsuitable for direct evaluation at 353.15 K. [Malvoisin et al. (2012), Appendix A](https://www.normalesup.org/~malvoisin/pdf/Malvoisin_et_al_2012b_jgr.pdf).

For the low-temperature implementation, an explicitly assumed effective rate with inventory and water limitation is acceptable for sensitivity experiments. If a decreasing H2 factor such as `1/(1+C_H2/C50)` is used, C50 is an **uncalibrated scenario parameter**, not a literature-derived half-inhibition constant. A thermodynamic alternative uses `A = Rg*T*ln(K/Q)` and a specified affinity function, but requires consistent activities, fugacities, reaction definition, standard states and equilibrium data. Concentration cannot silently replace fugacity in a high-pressure reaction quotient.

## Cracking, flow and feedback

A model that increases permeability from accumulated reaction or a scalar damage variable is a **damage proxy**. It cannot predict tensile stress, crack orientation, aperture or fracture energy unless it solves displacement/stress equilibrium and a failure law. A volume-based proxy should retain separate ledgers for precipitate volume (clogging) and assumed mechanical dilation (opening), and state the fixed-bulk-volume approximation. Do not create pore space without recording its constitutive origin.

Evans et al. couple poroelasticity, chemical mass balance and phase-field brittle failure. Their permeability includes matrix and oriented crack components; the latter enhances flow parallel to cracks. The exact permeability relation is acknowledged as uncertain, and their full reaction outcome depends on fluid access as well as failure. Their study supports a coupled feedback hypothesis, not an automatic rule that more reaction always produces more permeability or complete serpentinization. [Evans et al. (2020), sections 2.2 and 5](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019JB018614).

Use the same source chemistry and inventories for stagnant and flowing scenarios. Darcy advection, molecular diffusion, phase partition and boundary export then change H2 residence time. Water replenishment and H2 removal can alter reaction through the explicitly chosen availability/affinity functions. A lower local H2 concentration under stronger flow can accompany higher cumulative generation. Report separately: generated H2, remaining dissolved/gas H2, exported H2, consumed H2 and numerical balance error. Record water and reactive-rock consumption alongside those quantities.

## Connection to gravity, shallow microbiology and later LiDAR

The source body supplies a geological scenario to a gravity forward model. The transport result supplies a flux boundary or ensemble to the shallow model; it should not be inferred from gravity amplitude directly. Keep density, source capacity and reaction rate as different parameters.

H2-oxidising bacteria can consume atmospheric H2. Khdhiri et al. compare 0.5 and 10,000 ppmv exposures and find soil-dependent responses, so community composition alone is not a unique signature of a deep source. [Khdhiri et al. (2017)](https://journals.asm.org/doi/10.1128/AEM.00275-17). Soil H2 uptake, soil-gas profiles and independent geological evidence are needed to test the source-to-surface chain. Later LiDAR can characterize vegetation and terrain covariates; it is not a direct H2 sensor, and vegetation correlations require controls for moisture, soil, nutrients and land use.

## Verification targets

- Closed-cell stoichiometric tests conserve atoms and enforce nonnegative inventories.
- Zero kinetic rate and zero Fe cases generate no H2; no water disables hydration.
- With no source or sink, the transport ledger closes under both closed and open boundaries.
- Compare fixed permeability, clogging only, and explicit damage-proxy feedback using identical initial material inventories.
- Refine timestep and grid; compare generated and exported amounts, not only concentration pictures.
- Keep numerical verification separate from experimental calibration and field validation.
