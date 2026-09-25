# Native fields in ParaView

Use `bash run_model.sh export` after completing the native pair, then `bash run_model.sh render --pvpython /absolute/path/to/pvpython`. ParaView 6.1.1 was used for the original images. It runs in its own bundled Python, not in the DARTS environment.

`export` produces `source_model.pvd` and `vtk/snapshot_*.vtu`, and checks each native numeric field, cell centre/width and timestep. `render` creates geology, hydrogen and damage/permeability PNGs, saves `hydrogen_research.pvsm` and reloads it in a fresh process. Keep the PVD and VTU directory together. The state contains local absolute paths and must be regenerated after moving results.

To view another case, use `export --case results/darts/closed --output results/paraview_closed`, followed by `render --output results/paraview_closed`. Coarse pilot exports require `--pilot-label`; they are not scientific reference results.

Positive downward native depth becomes negative VTK elevation. Geometry has no vertical exaggeration. All numeric native CSV columns are cell data; formation and elevation are explicit aliases. The displayed H₂ field is free plus dissolved mass per pore-fluid volume, not gas saturation. The default view hides concentrations below 1% of the maximum for readability, while the reader retains the complete fields. Damage is an empirical continuum proxy, not a resolved fracture network.

The scripts verify file conversion and viewing. Model conservation, numerical convergence and physical validation are separate questions. The ParaView MPI runtime may require local sockets and the rendering backend may require a display/offscreen configuration on your platform.
