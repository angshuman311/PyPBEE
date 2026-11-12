import openseespy.opensees as ops
import numpy as np
import os

# Basic units (output units)
IN   = 1.0      # inch
KIP  = 1.0      # kip
SEC  = 1.0      # second

FT   = 12.0 * IN
G    = 32.2 * FT / SEC**2       # gravitational acceleration

ADD_ZEROS_FOR = 1.0          # seconds of zero-padding after GM
DT_ANALYSIS = 0.001          # analysis time step (sec)


def main(model_input_dict):
    mvals = model_input_dict["model_param_vals"]
    attrs = model_input_dict["model_attributes"]

    k = mvals["k_in_kip_per_in"]
    Fy = mvals["Fy_in_kip"]
    b = mvals["b"]
    m = mvals["W"] / G
    xi = mvals["damp_ratio"]

    ops.wipe()
    ops.model("BasicBuilder", "-ndm", 1, "-ndf", 1)

    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.mass(2, m)

    ops.uniaxialMaterial(attrs["spring_material"], 1, Fy, k, b)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)

    # c = 2 * xi * m * omega
    ops.rayleigh(2.0 * xi * np.sqrt(k / m), 0.0, 0.0, 0.0)


def analysis_gravity(model_input_dict):
    return True


def analysis_prelim(model_input_dict):
    model_info_directory = model_input_dict.get("model_info_directory", ".")

    # default eigen solver needs more than 1 dof
    lam = ops.eigen("-fullGenLapack", 1)
    period = 2.0 * np.pi / np.sqrt(lam[0])

    with open(os.path.join(model_info_directory, "periods.txt"), "w") as f:
        f.write(f"{period:.16g}")

    with open(os.path.join(model_info_directory, "mode_shape_1.txt"), "w") as f:
        f.write(f"2 {ops.nodeEigenvector(2, 1)[0]:.16g}\n")

    print("Modal Analysis COMPLETE!")


def analysis_time_hist(model_input_dict):
    gm_file_loc = model_input_dict.get("gm_file_loc", ".")

    # gm_1.txt : npts, dt, accelerations in g
    gm = np.loadtxt(os.path.join(gm_file_loc, "gm_1.txt")).flatten()
    npts = int(gm[0])
    dt_rec = gm[1]
    acc = np.append(gm[2:2 + npts], np.zeros(int(ADD_ZEROS_FOR / dt_rec)))
    total_time = len(acc) * dt_rec

    ops.timeSeries("Path", 2, "-dt", dt_rec, "-values", *acc.tolist(), "-factor", G)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2)

    ops.test("NormDispIncr", 1e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.50, 0.25)
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.system("FullGeneral")
    ops.analysis("Transient")

    ok = ops.analyze(int(total_time / DT_ANALYSIS), DT_ANALYSIS)

    ops.wipeAnalysis()
    return ok == 0
