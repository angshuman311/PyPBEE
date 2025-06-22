# PyPBEE — Performance-Based Earthquake Engineering in Python
A modular, high-performance framework that lets researchers and engineers run the full Performance-Based Earthquake Engineering (PBEE) workflow — from site‑specific seismic‑hazard curves to demand, damage and (soon) loss metrics — entirely in Python.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- Add CI, PyPI, or Docs badges here if applicable -->

---

## Key Capabilities
* **End‑to‑end PBEE pipeline**  
  Probabilistic Seismic Hazard Analysis (PSHA), Ground‑Motion Selection (GMS), Non‑linear Time‑History Analysis (NLTHA), Probabilistic Seismic Demand Hazard Analysis (PSDemHA) and Probabilistic Seismic Damage Hazard Analysis (PSDamHA).  
* **Object‑oriented core, plug‑and‑play extensions**  
  Clean template‑method lifecycle (`setup → run → wrap_up`) and abstract entity classes (`IM`, `EDP`, `DS`, `Structure`) make it straightforward to add new intensity measures, demand parameters, fragilities, or analysis stages.  
* **Advanced Uncertainty Quantification**  
  Treat aleatory *and* epistemic sources (random FE parameters, parameter‑estimation uncertainty, model‑form alternatives) with `MultivariateNataf`, `Mixture`, Latin‑Hypercube sampling, etc.  
* **Finite‑element back‑ends out of the box**  
  Adapters for **OpenSeesTcl** and **OpenSeesPy**; extendable to Abaqus, Ansys, …  
* **Scales from laptop to cluster**  
  Local multi‑processing via *pathos* or embarrassingly‑parallel Slurm jobs on HPC/HTC systems with a one‑line `comp_env` switch.  
* **Built‑in visualisation**  
  Hazard curves, conditional spectra, demand / damage hazard surfaces, deaggregation plots — all ready for publication.  

## Installation
```bash
# clone the repo
git clone https://github.com/angshuman311/PyPBEE.git
cd PyPBEE

# Windows users can run the helper script
setup_venv.bat                 # Windows CMD / double‑click
```
> ⚠️  PyPBEE is research software; pinned dependency versions are in `requirements.txt`.

## Quick‑start
```python
from pypbee.structures   import OSB               # Ordinary Standard Bridge
from pypbee.intensity    import AvgSa
from pypbee.edp          import MaxColRebarStrain
from pypbee.analysis     import PrelimAnalysis, PSHA, GMS, NLTHA, PSDemHA, PSDamHA

# 1 ── define model‐space & site info ------------------------------------------------
model_params   = {...}          # dict of deterministic design vars
random_params  = {...}          # dict of {param: (dist, mean, std, ...)}
corr_matrix    = [...]          # numpy array
location_info  = {"lat": 37.7531, "lon": -121.1427, "region": "California"}

# 2 ── create entities ---------------------------------------------------------------
bridge   = OSB(model_params, random_params, corr_matrix, location_info)
im       = AvgSa(bridge)
edp_list = [MaxColRebarStrain(bridge, tag='col_rebar_strain')]

# 3 ── assemble analyses -------------------------------------------------------------
pre  = PrelimAnalysis(bridge)
psha = PSHA(im)
gms  = GMS(im)
nltha   = NLTHA(im, edp_list)
psdemha = PSDemHA(edp_list, im)
psdamha = PSDamHA([], im)       # add DS objects when available

# 4 ── run workflow ------------------------------------------------------------------
for a in (pre, psha, gms, nltha, psdemha, psdamha):
    a.setup()
    a.run()
    a.wrap_up()
```
Full, commented examples live in `examples/`.

## Contributing
1. Fork the repo and create a feature branch.   
2. Open a pull request; describe *why* your change matters. 
3. Bug reports and feature requests are equally welcome — open a GitHub issue.

## Acknowledgements
Development supported by **Caltrans (65A0594, Task 2880)**, **PEER Transportation Systems Research Program (Project #1147‑NCTRTE)**, and the **Reissner Chair, UC San Diego**. HPC tests used TACC resources.

## License
This project is released under the [MIT License](LICENSE.md).
