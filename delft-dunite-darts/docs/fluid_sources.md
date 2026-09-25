# Solubility evidence and concentration basis

Verified 25 September 2026. This note supports an idealised model around 353 K and 190–300 bar. It does not establish Delft formation-water chemistry.

[Kerkache et al. (2024)](https://doi.org/10.1016/j.molliq.2024.124497), equations 9–15, give a Krichevsky–Kasarnovsky correlation:

\[
x_{H_2}=\frac{y_{H_2}P\phi_{H_2}}{H^0(T)\gamma(T,m)\Pi},\qquad
\Pi=\exp\!\left[\frac{19\times10^{-6}(P-P_{sat})}{RT}\right].
\]

Use Pa in the Poynting exponent, and the same pressure unit in the numerator and Henry constant. Their polynomial is **ln(H / MPa)**, with T in K:

\[
\ln H=B_0T^4+B_1T^3+B_2T^2+B_3T+B_4,
\quad B=(-4.40347\times10^{-10},8.33536\times10^{-7},-6.05506\times10^{-4},0.190837,-12.7991).
\]

The paper uses pure-H2 fugacity, an IAPWS saturation pressure, and \(y_{H_2}=1-P_{sat}/P\). The 19 cm³/mol partial volume is an **effective fitted quantity**, not a measured high-temperature constant. Its stated domain is 298.15–453.15 K, up to 100 MPa, and 0–6 mol NaCl/kg water. The high-pressure brine extension combines experiments with molecular simulations; it is not all experimentally validated.

Equation 11 is \(\gamma=\exp(K_s m)\). Equation 12 is transcribed from the original PDF text and publisher mathematical text as

\[
K_s=A_0T+A_1/T+A_2T/m+A_3+A_4mT,
\]

with \(A=(0.00159054,263.936,0.0000379192,-1.14604,-0.0000118419)\), \(K_s\) in kg water/mol NaCl and \(m\) in mol NaCl/kg water. The fraction layout could not be visually inspected with the available PDF viewer. The positive-brine expression has a nonzero \(A_2T\) term in \(\ln\gamma\) as \(m\to0\), so it must not silently bridge to the pure-water \(\gamma=1\) branch. A positive-salinity operating range and this limitation must be explicit.

The article distinguishes true liquid mole fraction \(x\) from salt-free \(x'\). Its main text retrieved here does not provide an explicit salt-count denominator. With \(q=\nu m M_w\), a declared convention gives

\[
x'=\frac{x(1+q)}{1+qx},\qquad b_{H_2}=\frac{x'}{M_w(1-x')},\quad M_w=0.01801528\ {m kg/mol}.
\]

Here \(\nu=1\) counts NaCl formula units and \(\nu=2\) counts Na+ and Cl− separately. The molecular description and “true” terminology motivate \(\nu=2\), but agreement with a benchmark alone cannot prove the author's convention. Preserve that qualification and quantify convention sensitivity. DARTS' water/H2 mole fraction should use \(x'\) when salt is a fixed background property rather than an explicit transported component.

## Measured checks

`data/literature/h2_solubility_benchmarks.csv` retains the **reported** units and denominator. Convert only in analysis; retain the original values.

- [Torín-Ollarves and Trusler, ELEGANCY D2.1.6 (2019)](https://www.sintef.no/globalassets/project/elegancy/deliverables/elegancy_d2.1.6_h2-solubility-in-brine_v2.pdf), Table 1, printed page 4, reports bubble pressures for 2.5 mol NaCl/kg water and explicitly salt-free H2 mole fractions. Nine rows span 323.15–423.15 K and 11.64–45.81 MPa. The table also supplies pressure standard uncertainties. The later [2021 journal paper](https://doi.org/10.1016/j.fluid.2021.113025) reports 3% expanded relative solubility uncertainty (coverage factor 2); do not automatically assign this to every preliminary report row.
- [Wolff et al. (2026)](https://www.frontiersin.org/journals/energy-research/articles/10.3389/fenrg.2026.1919332/full), Table 5, provides independent water measurements near 200 bar, in mol H2/kg water. Three water rows are retained. These are checks near 200 bar, not validation at 300–400 bar. The paper's NaCl solutions use **3 M**, not 3 mol/kg water; no unverified molarity-to-molality conversion is included.
- [NIST/IUPAC hydrogen compilation](https://srdata.nist.gov/solubility/IUPAC/SDS-5-6/SDS-5-6.pdf), printed pages 303 and 306, reproduces results attributed to [Wiebe and Gaddy (1934)](https://doi.org/10.1021/ja01316a022). The page-306 mole fractions and page-303 smoothed Kuenen coefficients differ by roughly 12% under an ordinary STP conversion. The cause is unresolved. They are **not** used as an exact acceptance target; no silent selection or recalibration is warranted.

[Chabab et al. (2024)](https://doi.org/10.1016/j.ijhydene.2023.10.290) measures only to 200 bar and explicitly discusses discrepancies between preceding brine studies, reaching 38%. It cannot by itself justify extending the older 230-bar property fit to deeper pressure. Report residuals by source and concentration convention; a plausible trend is not independent validation across the entire correlation domain.

## Bounded empirical alternative

The preliminary ELEGANCY table supplies an unambiguous option for a **fixed 2.5 mol/kg-water NaCl analogue**: linearly interpolate pressure along the 323.15 K and 373.15 K measured isotherms separately, then linearly interpolate temperature between them. Restrict this construction to 323.15–373.15 K and the common pressure support 120.1–394.7 bar. The proposed 353 K and 190–300 bar scenarios fall inside that support. This is an explicitly empirical interpolation of measured salt-free mole fractions. It is neither a new thermodynamic EOS nor an independently validated interpolation at every interior point; it supplies no salinity derivative or justification for another brine composition. Agreement at its own tabulated knots is a transcription check, not independent model validation.

## Delft provenance

The actual formation order is Rodenrijs Claystone above Delft Sandstone above Alblasserdam Claystone. The [2023-well end-of-well report](https://doi.org/10.4233/uuid:6ce07471-6986-434e-aa24-ad6e1f6714d9) and [Bruhn et al. (2026)](https://sd.copernicus.org/articles/35/83/2026/) establish the sedimentary setting. Measured depth along deviated wells is not vertical depth. Coarsened layer boundaries in this model are scenario geometry, not newly interpreted measured member tops. No cited Delft source establishes a dunite body.

[Zaal et al. (2021)](https://doi.org/10.1186/s40517-021-00193-0) predates the drilled doublet and uses a representative model: top 2200 m, thickness 105 m, porosity 5–30%, permeability 1–3400 mD, initial temperature 345.75 K. These are model inputs, not measured 2023 well-log ranges. Its Table 1 overburden permeability lists mutually inconsistent units, `9.9e-15 m²` and `0.001 mD`; do not transfer that row as a verified seal permeability.
