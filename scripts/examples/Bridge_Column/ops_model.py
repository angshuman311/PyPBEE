import openseespy.opensees as ops
import numpy as np
import os


def main_column_only(model_input_dict):
    fpc_col_1_secdef_1 = model_input_dict["model_param_vals"]["fpc_col_1_secdef_1"]
    steel_material = model_input_dict["model_attributes"]["steel_material"]

    # Unit weight of concrete [kip/in^3]
    wconc = 143.96 * 0.001 / (12 ** 3)

    # Tributary deck parameters to node 200
    Adeck = 97.546 * (12 ** 2)      # in^2
    L_left  = 1302.96               # in  span 100 → 200
    L_right = 1341.84               # in  span 200 → 300
    W200 = wconc * Adeck * 0.5 * (L_left + L_right)   # kip
    g = 386.4                       # in/s^2
    m200 = W200 / g                 # kip·s^2/in

    # Build reduced model
    ops.wipe()
    ops.model('BasicBuilder', '-ndm', 3, '-ndf', 6)

    # Transformations
    # Column local axes consistent with your original
    ops.geomTransf('Corotational', 2000, -0.5446390350150271, 0.838670567945424, 0.0)
    # Short rigid link from column top to node 200
    ops.geomTransf('Corotational', 3000, 0.838670567945424, 0.5446390350150271, 0.0)

    # --- Column materials (kept from your script; only what the column uses) ---
    # Bottom section set
    ops.uniaxialMaterial('Concrete01', 1001, -8.438902632269224, -0.0028296772213088676,
                         -1.687780526453845, -0.0629141365433386)
    ops.uniaxialMaterial('Concrete01', 2001, -fpc_col_1_secdef_1, -0.0021795371673297504,
                         -1.3, -0.005884823378103793)
    ops.uniaxialMaterial('Elastic', 4001, 6019019.069751658)    # shear Vy,Vz
    ops.uniaxialMaterial('Elastic', 5001, 3654760417.191347)    # torsion

    # Middle section set
    ops.uniaxialMaterial('Concrete01', 1002, -8.438902632269224, -0.0028296772213088676,
                         -1.687780526453845, -0.0629141365433386)
    ops.uniaxialMaterial('Concrete01', 2002, -6.5, -0.0021795371673297504,
                         -1.3, -0.005884823378103793)
    ops.uniaxialMaterial('Elastic', 4002, 6019019.069751658)
    ops.uniaxialMaterial('Elastic', 5002, 3654760417.191347)

    # Top section set
    ops.uniaxialMaterial('Concrete01', 1003, -8.438902632269224, -0.0028296772213088676,
                         -1.687780526453845, -0.0629141365433386)
    ops.uniaxialMaterial('Concrete01', 2003, -6.5, -0.0021795371673297504,
                         -1.3, -0.005884823378103793)
    ops.uniaxialMaterial('Elastic', 4003, 6019019.069751658)
    ops.uniaxialMaterial('Elastic', 5003, 3654760417.191347)

    # Longitudinal reinforcing (kept exactly so section matches your original)
    for tag in range(3001, 3013):
        ops.uniaxialMaterial(steel_material, tag, 69.144, 69.144, 29200.0,
                             0.01, 0.01, 20., 0.925, 0.15)

    # --- Column fiber sections and aggregators (as in your script) ---
    # Bottom fiber section
    ops.section('Fiber', 20001, '-torsion', 5001)
    ops.patch('circ', 2001, 132, 2, 0., 0., 30.24654000098859, 33.06, 0.0, 360.0)
    ops.patch('circ', 1001, 121, 5, 0., 0., 22.28692421125475, 30.24654000098859, 0.0, 360.0)
    ops.patch('circ', 1001, 89, 4, 0., 0., 15.919231579467679, 22.28692421125475, 0.0, 360.0)
    ops.patch('circ', 1001, 64, 3, 0., 0., 11.143462105627375, 15.919231579467679, 0.0, 360.0)
    ops.patch('circ', 1001, 44, 2, 0., 0., 7.9596157897338395, 11.143462105627375, 0.0, 360.0)
    ops.patch('circ', 1001, 32, 2, 0., 0., 4.775769473840303, 7.9596157897338395, 0.0, 360.0)
    ops.patch('circ', 1001, 19, 1, 0., 0., 3.1838463158935357, 4.775769473840303, 0.0, 360.0)
    ops.patch('circ', 1001, 13, 1, 0., 0., 1.5919231579467679, 3.1838463158935357, 0.0, 360.0)
    ops.patch('circ', 1001, 6, 1, 0., 0., 0.0, 1.5919231579467679, 0.0, 360.0)
    ops.layer('circ', 3001, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 0.0, 75.0)
    ops.layer('circ', 3002, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 90.0, 165.0)
    ops.layer('circ', 3003, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 180.0, 255.0)
    ops.layer('circ', 3004, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 270.0, 345.0)
    ops.section('Aggregator', 2001, 4001, 'Vy', 4001, 'Vz', '-section', 20001)

    # Middle fiber section
    ops.section('Fiber', 20002, '-torsion', 5002)
    ops.patch('circ', 2002, 132, 2, 0., 0., 30.24654000098859, 33.06, 0.0, 360.0)
    ops.patch('circ', 1002, 121, 5, 0., 0., 22.28692421125475, 30.24654000098859, 0.0, 360.0)
    ops.patch('circ', 1002, 89, 4, 0., 0., 15.919231579467679, 22.28692421125475, 0.0, 360.0)
    ops.patch('circ', 1002, 64, 3, 0., 0., 11.143462105627375, 15.919231579467679, 0.0, 360.0)
    ops.patch('circ', 1002, 44, 2, 0., 0., 7.9596157897338395, 11.143462105627375, 0.0, 360.0)
    ops.patch('circ', 1002, 32, 2, 0., 0., 4.775769473840303, 7.9596157897338395, 0.0, 360.0)
    ops.patch('circ', 1002, 19, 1, 0., 0., 3.1838463158935357, 4.775769473840303, 0.0, 360.0)
    ops.patch('circ', 1002, 13, 1, 0., 0., 1.5919231579467679, 3.1838463158935357, 0.0, 360.0)
    ops.patch('circ', 1002, 6, 1, 0., 0., 0.0, 1.5919231579467679, 0.0, 360.0)
    ops.layer('circ', 3005, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 0.0, 75.0)
    ops.layer('circ', 3006, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 90.0, 165.0)
    ops.layer('circ', 3007, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 180.0, 255.0)
    ops.layer('circ', 3008, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 270.0, 345.0)
    ops.section('Aggregator', 2002, 4002, 'Vy', 4002, 'Vz', '-section', 20002)

    # Top fiber section
    ops.section('Fiber', 20003, '-torsion', 5003)
    ops.patch('circ', 2003, 132, 2, 0., 0., 30.24654000098859, 33.06, 0.0, 360.0)
    ops.patch('circ', 1003, 121, 5, 0., 0., 22.28692421125475, 30.24654000098859, 0.0, 360.0)
    ops.patch('circ', 1003, 89, 4, 0., 0., 15.919231579467679, 22.28692421125475, 0.0, 360.0)
    ops.patch('circ', 1003, 64, 3, 0., 0., 11.143462105627375, 15.919231579467679, 0.0, 360.0)
    ops.patch('circ', 1003, 44, 2, 0., 0., 7.9596157897338395, 11.143462105627375, 0.0, 360.0)
    ops.patch('circ', 1003, 32, 2, 0., 0., 4.775769473840303, 7.9596157897338395, 0.0, 360.0)
    ops.patch('circ', 1003, 19, 1, 0., 0., 3.1838463158935357, 4.775769473840303, 0.0, 360.0)
    ops.patch('circ', 1003, 13, 1, 0., 0., 1.5919231579467679, 3.1838463158935357, 0.0, 360.0)
    ops.patch('circ', 1003, 6, 1, 0., 0., 0.0, 1.5919231579467679, 0.0, 360.0)
    ops.layer('circ', 3009, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 0.0, 75.0)
    ops.layer('circ', 3010, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 90.0, 165.0)
    ops.layer('circ', 3011, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 180.0, 255.0)
    ops.layer('circ', 3012, 6, 2.8613720136675442, 0., 0., 29.292180006018135, 270.0, 345.0)
    ops.section('Aggregator', 2003, 4003, 'Vy', 4003, 'Vz', '-section', 20003)

    # --- Column element and nodes only ---
    # Column base and top
    ops.node(2001, 1302.96, 0.0,   0.0)
    ops.node(2002, 1302.96, 0.0, 335.0)

    # Superstructure replacement node sitting above the column
    ops.node(200,  1302.96, 0.0, 364.76)

    # Integration rule identical to your original
    ops.beamIntegration('UserDefined', 2001, 8,
                        2001, 2002, 2002, 2002, 2002, 2002, 2002, 2003,
                        0, 0.06413, 0.20415, 0.39535, 0.60465, 0.79585, 0.93587, 1,
                        0.017857, 0.10535, 0.17056, 0.20623, 0.20623, 0.17056, 0.10535, 0.017857)

    # Column element
    ops.element('forceBeamColumn', 2001, 2001, 2002, 2000, 2001, '-iter', 50, 1e-12)

    # Short rigid link to the mass node
    ops.element('elasticBeamColumn', 2002, 2002, 200,
                1.0, 10.e10, 10.e10, 1.0, 1.0, 1.0, 3000)

    # Boundary conditions
    ops.fix(2001, 1, 1, 1, 1, 1, 1)

    # Lumped tributary mass and gravity load at node 200
    ops.mass(200, m200, m200, m200, 0.0, 0.0, 0.0)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(200, 0.0, 0.0, -W200, 0.0, 0.0, 0.0)

    # Optional: keep gravity constant if you plan to run dynamics after this call
    # ops.loadConst('-time', 0.0)

    # For quick verification you can print these two
    print(f"Adeck = {Adeck:.6f} in^2")
    print(f"W200 = wconc * {Adeck * 0.5 * (L_left + L_right):.6f}  -> numeric W200 ≈ {float(W200):.6f} kip")
    print(f"m200 = W200 / g -> numeric m200 ≈ {float(m200):.9f} kip·s^2/in")
