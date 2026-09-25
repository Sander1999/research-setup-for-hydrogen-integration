"""Measured fluid benchmarks; these data were not used to fit our code."""
import json
from pathlib import Path
from ..properties import solubility_x,MW


def benchmarks(output=None):
    rows=[]
    url='https://www.frontiersin.org/journals/energy-research/articles/10.3389/fenrg.2026.1919332/full'
    for t,p,b in [(323.15,195.3,.1380),(348.15,194.4,.1408),(373.15,192.8,.1511)]:
        reference=b/(b+1000/MW[1]);predicted=float(solubility_x(p,t,0))
        rows.append(dict(study='Wolff et al.2026 Table5',temperature_k=t,pressure_bar=p,molality=0.,salt_basis_particles=0,
                         observed_salt_free_x=reference,model_salt_free_x=predicted,relative_error=predicted/reference-1,source_url=url))
    url='https://www.sintef.no/globalassets/project/elegancy/deliverables/elegancy_d2.1.6_h2-solubility-in-brine_v2.pdf'
    for t,p,x in [(323.15,120.1,.001099),(323.15,262.3,.002212),(323.15,394.7,.003410),(373.15,116.4,.001148),(373.15,229.4,.002347),(373.15,277.6,.002873),(373.15,458.1,.004242),(423.15,149.4,.001724),(423.15,387.3,.004557)]:
        for particles in [1.,2.]:
            predicted=float(solubility_x(p,t,2.5,particles))
            rows.append(dict(study='ELEGANCY D2.1.6 Table1',temperature_k=t,pressure_bar=p,molality=2.5,salt_basis_particles=particles,
                             observed_salt_free_x=x,model_salt_free_x=predicted,relative_error=predicted/x-1,source_url=url))
    result={'model':'Kerkache2024 gamma-phi correlation with integrated Lemmon2008 H2 fugacity, IAPWS saturation vapour pressure and explicit salt-free conversion',
        'basis_notice':'The paper distinguishes true and salt-free mole fractions but the retrieved text does not explicitly specify the salt particle-count denominator. The default adopts dissociated Na+/Cl- counting (2); formula-unit counting (1) is included as a model-convention sensitivity. Agreement with measurements is not used as proof of the denominator definition.',
        'equation_transcription':'Eq12 A2*T/ms is consistent across publisher text and original PDF text extraction; positive-brine ln(gamma) is evaluated without division. Pure water is a separate gamma=1 branch.',
        'out_of_case_discrepancy_notice':'All available report points are shown, including ~17–20% mismatch at150C and8.4%at458bar100C. The default80C190–300bar case does not use those conditions. This is model/data uncertainty, not numerical solver error.',
        'model_scope':'High-pressure extension is partly supported by molecular simulation; selected measurements validate only their listed P/T/salinity points. No new coefficients were fitted.',
        'rows':rows}
    if max(abs(r['relative_error']) for r in rows if r['molality']==0)>.04:raise AssertionError('Fresh-water benchmark >4%')
    if max(abs(r['relative_error']) for r in rows if r['salt_basis_particles']==2 and r['temperature_k']<=373.15 and r['pressure_bar']<=400)>.05:raise AssertionError('Brine benchmark within the used T/P neighbourhood >5%')
    if output:
        output=Path(output);output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    return result
