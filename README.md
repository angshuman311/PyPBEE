# PyPBEE — Performance-Based Earthquake Engineering in Python
A modular, high-performance framework that lets researchers and engineers run the full Performance-Based Earthquake Engineering (PBEE) workflow — from site-specific seismic hazard analysis to demand, damage and (soon) loss hazard analysis — entirely in Python.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Key Capabilities
* **End-to-end PBEE pipeline**  
  Probabilistic Seismic Hazard Analysis (PSHA), Ground Motion Selection (GMS), Nonlinear Time-History Analysis (NLTHA), Probabilistic Seismic Demand Hazard Analysis (PSDemHA) and Probabilistic Seismic Damage Hazard Analysis (PSDamHA).  
* **Object-oriented core, plug-and-play extensions**  
  Clean template-method lifecycle (`setup → run → wrap_up`) and abstract entity classes (`IM`, `EDP`, `DS`, `Structure`) make it straightforward to add new intensity measures, demand parameters, fragilities, or analysis stages.  
* **Advanced Uncertainty Quantification**  
  Treat aleatory *and* epistemic sources (random FE parameters, parameter estimation uncertainty, model-form alternatives) with `MultivariateNataf`, `Mixture`, Latin Hypercube sampling, etc.  
* **Finite-element backends out of the box**  
  Adapters for **OpenSeesTcl** and **OpenSeesPy**; extendable to Abaqus, Ansys, …  
* **Scales from laptop to cluster**  
  Local multiprocessing via *pathos* or embarrassingly-parallel Slurm jobs on HPC/HTC systems.  
* **Built-in visualisation**  
  Hazard curves, conditional spectra, demand / damage hazard surfaces, deaggregation plots — all ready for publication.  

## Quick install (users)

```bash
pip install pypbee
```

`pip` will fetch the latest wheel from PyPI along with NumPy, SciPy, Matplotlib, etc.

---

## Developer setup (work with the source)

```bash
git clone https://github.com/angshuman311/PyPBEE.git
cd PyPBEE
```

### Windows

Simply **double-click** `setup_venv.bat` (or run it in CMD/PowerShell).  
The script:

1. creates a virtual environment `venv\pypbee`;
2. installs everything listed in `requirements.txt`;

Point your IDE’s interpreter to:

```
/path/to/PyPBEE/venv/pypbee/Scripts/python.exe
```

### Linux / macOS

Do the same steps manually:

```bash
python3 -m venv venv/pypbee    
source venv/pypbee/bin/activate

# install requirements
pip install -r requirements.txt
```

Point your IDE’s interpreter to:

```
/path/to/PyPBEE/venv/pypbee/Scripts/python.exe
```

Once the env is active you can edit code under `pypbee/`, run examples in `scripts/`, and submit pull requests.

---

## Conceptual PSEUDO code
```python
from pypbee.structure                               import OSB               # Ordinary Standard Bridge
from pypbee.avg_sa                                  import AvgSa
from pypbee.edp                                     import MaxColRebarStrain
from pypbee.ds                                      import DS
from pypbee.analysis                                import PrelimAnalysis, PSHA, GMS, NLTHA, PSDemHA, PSDamHA
from pypbee.utility                                 import Utility
from pypbee.pygmm_extension.boore_atkinson_2008     import BooreAtkinson2008
from pygmm.baker_jayaram_2008                       import calc_correls

# 1 ── define model‐params & site info ------------------------------------------------
# See files `osb_info_*.py` in `examples/Bridge_*/osb_info_*.py` for model_params and location_info definition
name = 'Bridge_A'
model_files_path = f'path/to/{name}'
model_params = Utility.import_attr_from_module(module_path=model_files_path, module_name=f"osb_info_{name}", attr='model_params')
location_info = Utility.import_attr_from_module(module_path=model_files_path, module_name=f"osb_info_{name}", attr='location_info')
local_python_path = '/path/to/PyPBEE/venv/pypbee/Scripts/python.exe'
model_work_dir_path = f'path/to/a/working/directory/for/results/data/storage/for/{name}'

# 2 ── create entities ---------------------------------------------------------------
structural_analysis_platform = OpenSeesPy(model_files_path, local_python_path)
osb         = OSB(name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform)
im          = AvgSa(osb, gmm=BooreAtkinson2008, correl_func=calc_correls, define_range=['T_1_trans', 'T_1_trans'], range_multiplier=[1, 2.5])
edp_list    = [
                    MaxColRebarStrain(max_what='compression', frame_structure=osb, tag='1', recorder_file_storage='shared'),
                    MaxColRebarStrain(max_what='tension', frame_structure=osb, tag='2', recorder_file_storage='shared'),
                    MaxSpringDeformation(spring_type='shear_key', max_what='compression', frame_structure=osb, tag='4', recorder_file_storage='separate', normalize_with='D3')
              ]
ds_list     = [
                    DS(
                        edp=edp_list[0], predictor=lambda x: 0.004,
                        haz_req={'normalized_fragility_dist': lognorm(0.326, 0, 1.02),
                                'estimation_sample_size': 5},
                        ds_type='col_rebar_strain_damage'
                    )
                    DS(
                        edp=edp_list[1], predictor=lambda x: 0.03 + 700 * x[1] * x[2] / x[3] - 0.1 * x[8] / (x[4] * x[5]),
                        haz_req={'normalized_fragility_dist': lognorm(0.201, 0, 1.05),
                                'estimation_sample_size': 5},
                        ds_type='col_rebar_strain_damage'
                    )
                    DS(
                        edp=edp_list[2], predictor=lambda x: 1.0,
                        haz_req={'normalized_fragility_dist': lognorm(0.11, 0, 1.14),
                                'estimation_sample_size': 5},
                        ds_type='spring_deformation_damage',
                    )
              ]


# 3 ── assemble analyses -------------------------------------------------------------
pre         = PrelimAnalysis(osb, num_modes=8)
psha        = PSHA(im)
gms         = GMS(im)
nltha       = NLTHA(im, edp_list)
psdemha     = PSDemHA(edp_list, im)
psdamha     = PSDamHA(ds_list, im, sol_type='numerical')

# 4 ── run workflow ------------------------------------------------------------------
for a in (pre, psha, gms, nltha, psdemha, psdamha):
    a.setup(...)
    a.run(...)
    a.wrap_up(...)
```
Full, commented examples live in `examples/` and `scripts/`.

## Contributing
1. Fork this repository and clone your fork locally.   
2. Create a feature branch, make your changes (e.g., add a new `IM`, `EDP`, etc.), and push the branch to your fork.  
3. Open a pull request describing *what* you changed and *why* it’s useful.
4. Bug reports and feature requests are equally welcome — open a GitHub issue.

## Acknowledgements
Development supported by **Caltrans (65A0594, Task 2880)**, **PEER Transportation Systems Research Program (Project #1147-NCTRTE)**, and the **Reissner Chair, UC San Diego**. HPC tests used TACC resources.

## License
This project is released under the [MIT License](LICENSE.md).
