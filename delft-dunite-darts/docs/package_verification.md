# Standalone package verification

Checked on macOS arm64 with the existing Python 3.11 native DARTS runtime recorded in `environment-reference.json`.

- Actual engine/API imports, configuration and solubility calculation passed.
- Eleven native unit/conservation tests passed; 21 fluid benchmark rows passed their declared acceptance checks.
- Both complete **standard 3,650-day cases**, at their original 3,648-cell resolution, were rerun through the portable launcher in a separate directory.
- Generated, retained dissolved, retained gas and exported hydrogen agreed with the retained original reference at relative tolerance `1e-10` and absolute tolerance `1e-9 kg`. The native model/configuration fingerprint was unchanged.
- All four accepted spatial snapshots exported to VTK. Every numeric cell array matched exactly; cell geometry and time checks passed.
- ParaView 6.1.1 rendered the 3D fields and reloaded its saved state in a fresh process, checking first/final frames.

This packaging check did not repeat all eight controls or the expensive paired refinement suite: those native files were unchanged, and the original refinement report is retained in `reference/numerical_verification.json`. It explicitly records unmet gas-inventory grid sensitivity targets. No clean dependency installation on a second machine or GitHub-hosted native run is claimed.
