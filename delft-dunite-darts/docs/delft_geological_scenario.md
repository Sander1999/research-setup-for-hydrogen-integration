# Delft sedimentary architecture with a hypothetical hydrogen source

The Delft setting is useful because it provides a well-studied sedimentary
framework. The final scenario is still counterfactual: no dunite body or natural
hydrogen accumulation is inferred at Delft. The body is deliberately inserted
below a simplified representation of the Delft Sandstone Member (DSSM), with
alternative connections through the underlying rocks. Its surface ecology is a
separate natural-landscape analogue.

The distinction is recorded in `data/literature/delft_geology_evidence.csv`.
Published wells, cores and models constrain layer order, approximate depth,
temperature and petrophysical ranges. They do not fix this model's exact
horizontal boundaries, formation surfaces, capillary curves, stress or natural
groundwater gradient. A simplified layered scenario is not the published Delft
digital twin.

## Published constraints

The campus reservoir is the fluvial Lower Cretaceous Delft Sandstone. The
production and injection wells are deviated, so measured depth along the well
cannot be used as vertical depth. Published initial production temperature is
about 80°C. [Vardon et al. (2024)](https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2024/Vardon.pdf).

The end-of-well open-hole figure on PDF page211 covers approximately 2110–2335m
true vertical depth while the along-well interval is around2570–2930m. The logs
show marked vertical heterogeneity. The illustration alone does not justify
assigning a precise member boundary to every log break.
[End-of-well science report](https://doi.org/10.4233/uuid:6ce07471-6986-434e-aa24-ad6e1f6714d9).

The Delft modelling study's Table1 spans porosity0.01–0.264 and permeability
0.004–1113mD. These are model ranges constrained by geological/log information,
not measured hydrogen transport curves. Their use here supplies plausible
contrasts between sandy intervals and baffles.
[Voskov et al. (2024)](https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2024/Voskov.pdf).

The relevant order is Rodenrijs Claystone over Delft Sandstone over Alblasserdam
intervals. This supports testing caprock, reservoir and underlying connectivity
separately. It does not establish a connected pathway to the surface.
[UrbEnLab geological description](https://sd.copernicus.org/articles/35/83/2026/).

## Numerical interpretation

The exact layer table and geometric coordinates are stored in
`hydrogen_research/darts_source/config.json` and repeated in every completed run.
Depth is positive down. The regular horizontal grid and layer-aligned vertical
grid coarsen the geology; they do not reproduce individual fluvial channels or
fault offsets. Layer tops near1900–2800m and internal DSSM divisions are declared
scenario choices. Sensitivity results should therefore be used to identify which
measurements matter, not to claim where hydrogen would be trapped at the campus.

The lower alternative connections are deliberately hypothetical. A sealing
interval tests isolation. A connected corridor tests natural water access and
hydrogen migration from the source to sandstone. The same finite reactive iron
inventory is used when comparing the connections. An imposed hydraulic gradient
represents an externally maintained natural head difference; it is neither a
measured Delft flow nor a engineered pumping schedule. Report pressure, flow and
generation together rather than assuming the open corridor is always beneficial.

The isothermal80°C choice is compatible with the local temperature scale but
does not solve a geothermal profile or reaction heat. At these depths the
original bachelor brine fit's230bar limit matters. The selected fluid formulation
and its independently checked pressure range must be read in the native model
documentation before extending the case deeper.

## Gravity and ecology use the same conceptual case

The gravity forward calculation sums the actual exported source cells. Its static
anomaly is relative to the local host-density model; flat background layers are
removed by definition. A dense nonreactive body can create the same anomaly.
The time-lapse calculation uses changes in fluid mass plus the accounted solid
mass gain, so reaction water transferred into solids is not mistaken for mass
leaving the subsurface. Fixed-volume and prescribed-compliance approximations
limit this diagnostic. No claim of measurable time-lapse gravity is made before
comparing its amplitude with the assumed survey uncertainty.

The grassland/shrubland/woodland analogue changes the surface sampling context,
not the Delft subsurface data. It is deliberately not a model of campus plants.
LiDAR measures structure and terrain; independent H2 and microbial observations
are still required. The physical model produces hydrogen exposure and uptake
scenarios, not an invented species-specific plant-growth response.
