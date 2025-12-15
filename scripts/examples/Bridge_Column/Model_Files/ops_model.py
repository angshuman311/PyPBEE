import openseespy.opensees as ops
import numpy as np
import os
import math

# Basic units (output units)
IN   = 1.0      # inch
KIP  = 1.0      # kip
SEC  = 1.0      # second

# Engineering units and derived units
FT   = 12.0 * IN
KSI  = KIP / IN**2
PSI  = KSI / 1000.0

LBF  = PSI * IN * IN            # pounds force
PCF  = LBF / FT**3              # pounds per cubic foot
PSF  = LBF / FT**2              # pounds per square foot

IN2  = IN**2
IN4  = IN**4

CM   = IN / 2.54
CMSEC2 = CM / SEC**2            # cm/sec^2

M    = CM * 100.0
MM   = CM / 10.0
MM2  = MM**2

KN   = 0.2247 * KIP
N    = 1.0e-3 * KN
MN   = 1.0e6 * N

MPa  = 0.1450 * KSI
GPa  = 1000.0 * MPa

PI   = math.pi
G    = 32.2 * FT / SEC**2       # gravitational acceleration
UBIG = 10.0e10                  # a really large number
USMALL = 1.0 / UBIG             # a really small number


def create_circle_patches(matTag, yCenter, zCenter,
                          radiusInner, radiusOuter,
                          minFiber, maxFiber,
                          startAng, endAng,
                          export_file=None):
    """
    Creates concentric circular patches with varying fiber sizes between
    radiusInner and radiusOuter, from startAng to endAng (degrees).

    If export_file is not None, writes the patch commands to that file.
    """

    PI_local = PI  # use global PI

    # number of concentric rings and their width
    numCirc   = int((radiusOuter - radiusInner) / maxFiber + 0.5)
    if numCirc < 1:
        numCirc = 1
    ringWidth = (radiusOuter - radiusInner) / numCirc

    # ratio of min fiber to max fiber
    fibRatio = minFiber / maxFiber

    jRadiusOuter  = radiusOuter
    jRadiusCutoff = radiusOuter * fibRatio
    jRadiusInner  = radiusInner

    for iRing in range(numCirc, 0, -1):
        iRadius = radiusInner + iRing * ringWidth

        # if inside cutoff, create a patch from jRadiusInner to jRadiusOuter
        if iRadius < jRadiusCutoff:
            jRadiusInner = iRadius

            jNumC = max(int(jRadiusOuter * (endAng - startAng) * (PI_local / 180.0) / maxFiber + 0.5), 1)
            jNumR = max(int((jRadiusOuter - jRadiusInner) / ringWidth + 0.5), 1)

            # openseespy patch circ
            ops.patch("circ", matTag, jNumC, jNumR,
                      yCenter, zCenter,
                      jRadiusInner, jRadiusOuter,
                      startAng, endAng)

            if export_file is not None:
                export_file.write(
                    f"\tpatch circ {matTag} {jNumC} {jNumR} {yCenter} "
                    f"{zCenter} {jRadiusInner} {jRadiusOuter} {startAng} {endAng}\n"
                )

            jRadiusOuter  = jRadiusInner
            jRadiusCutoff = jRadiusOuter * fibRatio

        if iRing == 1:
            jNumC = max(int(jRadiusOuter * (endAng - startAng) * (PI_local / 180.0) / maxFiber + 0.5), 1)
            jNumR = max(int((jRadiusOuter - radiusInner) / ringWidth + 0.5), 1)

            ops.patch("circ", matTag, jNumC, jNumR,
                      yCenter, zCenter,
                      radiusInner, jRadiusOuter,
                      startAng, endAng)

            if export_file is not None:
                export_file.write(
                    f"\tpatch circ {matTag} {jNumC} {jNumR} {yCenter} "
                    f"{zCenter} {radiusInner} {jRadiusOuter} {startAng} {endAng}\n"
                )


def main(model_input_dict):
    """
    Build a reduced OpenSeesPy model of the cantilever bridge column,
    using the parameters defined in model_info.py.

    PyPBEE passes:
      model_input_dict = {
          'model_files_path': ...,
          'model_info_directory': ...,
          'write_model_files': 1 or 0,
          'model_param_vals': model_param_vals,
          'model_attributes': model_attributes,
      }

    Here:
      - model_param_vals contains:
          * random model parameters (from rv_list)
          * primary design parameters
          * other model parameters
      - model_attributes contains selected model-form options

    It is the user's responsibility to keep parameter names consistent
    between model_info.py and this model file.
    """

    ops.wipe()
    ops.model("BasicBuilder", "-ndm", 3, "-ndf", 6)

    # ------------------------------------------------------------------
    # 1. Unpack dictionaries from model_input_dict
    # ------------------------------------------------------------------
    mvals = model_input_dict["model_param_vals"]
    attrs = model_input_dict["model_attributes"]

    write_model_files = model_input_dict.get("write_model_files", 0)
    model_info_directory = model_input_dict.get("model_info_directory", ".")

    # ----------------------------------------------------------------------
    # Inputs from model_param_vals (PyPBEE)
    # ----------------------------------------------------------------------
    all_col_dia_in_ft = mvals["all_col_dia_in_ft"]
    all_rho_long = float(mvals["all_rho_long"])
    num_bar_clusters = int(mvals["num_bar_clusters"])
    num_secdef_per_col = int(mvals["num_secdef_per_col"])

    # ----------------------------------------------------------------------
    # Tributary spans and deck section properties
    # ----------------------------------------------------------------------
    L1 = 108.58 * FT
    L2 = 111.82 * FT
    L_trib = L1 + L2

    Adeck = 97.546 * FT ** 2  # in^2
    Jdeck = 341.442 * FT ** 4  # in^4
    Iydeck = 180.328 * FT ** 4  # in^4
    Izdeck = 3797.9 * FT ** 4  # in^4

    # ----------------------------------------------------------------------
    # Column geometry and reinforcement parameters
    # ----------------------------------------------------------------------

    dcg = 2.48 * FT  # deck centroid to column top distance

    DcolAll = all_col_dia_in_ft * FT  # column diameter in inches

    rho_long_all = all_rho_long
    num_bar_clusters = int(num_bar_clusters)
    num_secdef_per_col = int(num_secdef_per_col)

    # Transverse reinforcement ratio
    rhoTransAll = 0.5 * rho_long_all

    # Column section properties
    Dcol = DcolAll
    Acol = (PI * Dcol ** 2) / 4.0
    Jcol = (PI * (Dcol / 2.0) ** 4) / 2.0
    Izcol = (PI * (Dcol / 2.0) ** 4) / 4.0
    Iycol = Izcol

    rhoLong = rho_long_all
    rhoTrans = rhoTransAll

    numBarCol = 24  # number of longitudinal bars

    # Longitudinal bar diameter and area
    Dbar = math.sqrt(4.0 * rhoLong * Acol / (numBarCol * PI))
    Abar = (PI * Dbar ** 2) / 4.0

    # Cover and hoop spacing
    cover = 2.0 * IN
    sH = 3.3465 * IN

    h2hOuter = Dcol - 2.0 * cover

    # Hoop bar area and diameter
    AbarH = rhoTrans * h2hOuter * sH / 4.0
    DbarH = math.sqrt(4.0 * AbarH / PI)

    sHprime = sH - DbarH
    ds = h2hOuter - DbarH
    Ac = (PI * ds ** 2) / 4.0

    rho_cc = numBarCol * Abar / Ac
    A_cc = Ac * (1.0 - rho_cc)

    if all_col_dia_in_ft < 6.0:
        DbarSingle = 1.410 * IN
    else:
        DbarSingle = 1.693 * IN

    HColumn	= 335. * IN

    # steel material model form
    steel_material = attrs["steel_material"]

    # Open files if requested
    if write_model_files == 1:
        mat_data_path = f"{model_info_directory}/material_data.txt"
        col_rebar_info_path = f"{model_info_directory}/col_rebar_mat_info.txt"
        matDataFileID = open(mat_data_path, "w")
        colRebarMatTagInfoFileID = open(col_rebar_info_path, "w")
    else:
        matDataFileID = None
        colRebarMatTagInfoFileID = None

    # Material tags
    IDconcCore = 1000
    IDconcCover = 2000
    IDSteel = 3000
    IDShear = 4000
    IDTorsion = 5000
    ID_Rigid = 99999

    # Rigid material
    ops.uniaxialMaterial("Elastic", ID_Rigid, UBIG)

    # Concrete unit weight and mass density
    wconc = mvals["wconc_all"]  # kip/in^3
    mconc = wconc / G  # mass per volume

    # Single column and single bent
    nCols = 1
    nBent = 1

    # ----------------------------------------------------------------------
    # Confined and unconfined concrete materials (loop over sections)
    # ----------------------------------------------------------------------
    ctr = 1

    for i in range(1, nCols * nBent + 1):
        for j in range(1, num_secdef_per_col + 1):

            # Predictor files for rebar strain damage (optional)
            if write_model_files == 1:
                if j == 1:
                    pred1_path = f"{model_info_directory}/predictor_info_col_rebar_strain_damage_col_{i}_edge_1.txt"
                    predFileID1 = open(pred1_path, "w")

            # fpc and Ec from model_param_vals (already in ksi)
            fpc = mvals[f"fpc_col_{i}_secdef_{j}"]
            Ec = mvals[f"Ec_col_{i}_secdef_{j}"]
            Uc = 0.2

            # Hoop steel: use bar cluster 1 for confinement
            randCtr = 1
            fyH = mvals[f"fy_col_{i}_secdef_{j}_barcluster_{randCtr}"]
            EsH = mvals[f"Es_col_{i}_secdef_{j}_barcluster_{randCtr}"]

            # Confinement computations
            fl = 2.0 * AbarH * fyH / (ds * sH)
            Ke = (1.0 - sHprime / (2.0 * ds)) / (1.0 - rho_cc)
            fl_prime = fl * Ke
            Kc = 9.0 / math.sqrt(1.0 + 8.0 * fl_prime / fpc) - 2.0
            n_val = 0.8 + fpc / (2500.0 * PSI)
            epc = fpc / Ec * (n_val / (n_val - 1.0))
            r_val = 1.0 / (1.0 - (fpc / epc) / Ec)
            ecu = 0.004
            fcu = r_val * ((ecu / epc) / (r_val - 1.0 + (ecu / epc) ** r_val)) * fpc
            slope_UC = (fpc - fcu) / (epc - ecu)
            ecu_UC = epc - (fpc - 0.2 * fpc) / slope_UC

            # Confined core
            fpcc = fpc + Kc * fl_prime
            # Ec for confined concrete from SDC style formula (psi inside, ksi outside)
            Ecc = (40000.0 * math.sqrt(fpcc / PSI) + 1.0e6) * PSI
            Gcc = Ecc / (2.0 * (1.0 + Uc))

            epcc = epc * (1.0 + 5.0 * ((fpcc / fpc) - 1.0))
            epcc1 = fpcc * epc / fpc
            rc = 1.0 / (1.0 - (fpcc / epcc) / Ecc)
            rho_s = (AbarH * PI * ds) / (PI * (ds ** 2) * sH / 4.0)
            eccu = 0.004 + 3.0 * Ke * rho_s
            fccu = rc * ((eccu / epcc) / (rc - 1.0 + (eccu / epcc) ** rc)) * fpcc
            slope_C = (fpcc - fccu) / (epcc1 - eccu)
            eccu_C = epcc1 - (fpcc - 0.2 * fpcc) / slope_C

            # Core and cover concrete materials
            ops.uniaxialMaterial(
                "Concrete01",
                IDconcCore + ctr,
                -fpcc,
                -epcc1,
                -0.2 * fpcc,
                -eccu_C
            )

            ops.uniaxialMaterial(
                "Concrete01",
                IDconcCover + ctr,
                -fpc,
                -epc,
                -0.2 * fpc,
                -ecu_UC
            )

            if write_model_files == 1:
                matDataFileID.write(
                    f"uniaxialMaterial Concrete01 {IDconcCore + ctr} {-fpcc} {-epcc1} {-0.2 * fpcc} {-eccu_C};\n"
                )
                matDataFileID.write(
                    f"uniaxialMaterial Concrete01 {IDconcCover + ctr} {-fpc} {-epc} {-0.2 * fpc} {-ecu_UC};\n"
                )

            # Shear and torsion materials for section aggregation
            ops.uniaxialMaterial("Elastic", IDShear + ctr, 0.9 * Gcc * Acol)
            ops.uniaxialMaterial("Elastic", IDTorsion + ctr, 1.0 * Gcc * Jcol)

            if write_model_files == 1:
                matDataFileID.write(
                    f"uniaxialMaterial Elastic {IDShear + ctr} {0.9 * Gcc * Acol};\n"
                )
                matDataFileID.write(
                    f"uniaxialMaterial Elastic {IDTorsion + ctr} {1.0 * Gcc * Jcol};\n"
                )

            # Predictor info for plastic hinge edges
            if write_model_files == 1:
                if j == 1:
                    predFileID1.write(f"{rhoTrans}\n")
                    predFileID1.write(f"{fyH}\n")
                    predFileID1.write(f"{EsH}\n")
                    predFileID1.write(f"{fpc}\n")
                    predFileID1.write(f"{(PI * Dcol ** 2) / 4.0}\n")
                    predFileID1.write(f"{numBarCol}\n")
                    predFileID1.close()
            ctr += 1

    # ----------------------------------------------------------------------
    # Reinforcing steel materials
    # ----------------------------------------------------------------------
    ctr = 1
    fy_e = 0.0

    for i in range(1, nCols * nBent + 1):
        for j in range(1, num_secdef_per_col + 1):
            ToverY_max = 0.0
            for k in range(1, num_bar_clusters + 1):

                fy = mvals[f"fy_col_{i}_secdef_{j}_barcluster_{k}"]
                fy_e += fy

                fu = mvals[f"fu_col_{i}_secdef_{j}_barcluster_{k}"]
                ToverY = fu / fy
                if ToverY > ToverY_max:
                    ToverY_max = ToverY

                b_val = mvals[f"b_col_{i}_secdef_{j}_barcluster_{k}"]
                Es = mvals[f"Es_col_{i}_secdef_{j}_barcluster_{k}"]

                Us = 0.2
                Gs = Es / (2.0 * (1.0 + Us))

                Fy = fy
                bs = b_val

                R0 = 20.0
                cR1 = 0.925
                cR2 = 0.15

                if steel_material == "SteelMPF":
                    ops.uniaxialMaterial(
                        "SteelMPF",
                        IDSteel + ctr,
                        Fy,
                        Fy,
                        Es,
                        bs,
                        bs,
                        R0,
                        cR1,
                        cR2
                    )
                elif steel_material == "Steel01":
                    ops.uniaxialMaterial(
                        "Steel01",
                        IDSteel + ctr,
                        Fy,
                        Es,
                        bs
                    )

                if write_model_files == 1:
                    if steel_material == "SteelMPF":
                        matDataFileID.write(
                            f"uniaxialMaterial SteelMPF {IDSteel + ctr} {Fy} {Fy} {Es} {bs} {bs} {R0} {cR1} {cR2};\n"
                        )
                    elif steel_material == "Steel01":
                        matDataFileID.write(
                            f"uniaxialMaterial Steel01 {IDSteel + ctr} {Fy} {Es} {bs};\n"
                        )
                    colRebarMatTagInfoFileID.write(f"{IDSteel + ctr}\n")

                ctr += 1

            if write_model_files == 1:
                if j == 1:
                    pred1_path = f"{model_info_directory}/predictor_info_col_rebar_strain_damage_col_{i}_edge_1.txt"
                    with open(pred1_path, "a") as predFileID1:
                        predFileID1.write(f"{ToverY_max}\n")

    # Mean fy over all bars
    if ctr > 1:
        fy_e /= (ctr - 1)

    if write_model_files == 1:
        matDataFileID.close()
        colRebarMatTagInfoFileID.close()

    # ----------------------------------------------------------------------
    # COLUMN SECTION
    # ----------------------------------------------------------------------

    # Fiber size targets (from Tcl)
    maxFiber = 0.04 * M
    minFiber = 0.03 * M

    ColSecTag = 2000

    yCenter = 0.0
    zCenter = 0.0

    # Radii for core and cover
    intRadCore = 0.0
    extRadCore = Dcol / 2.0 - cover - DbarH
    intRadCover = extRadCore
    extRadCover = Dcol / 2.0

    # Steel layer parameters
    numBar = numBarCol
    areaBar = Abar
    radius = Dcol / 2.0 - cover - DbarH - Dbar / 2.0
    theta = 360.0 / numBar

    ctr = 1
    ctrSec = 1

    for i in range(1, 1 * 1 + 1):  # nCols * nBent but both are 1
        for j in range(1, num_secdef_per_col + 1):

            # Each section definition gets its own file if requested
            if write_model_files == 1:
                sectionID = ColSecTag + ctrSec
                sec_path = f"{model_info_directory}/fib_secdef_sec_{sectionID}.txt"
                ColSecDefFileID = open(sec_path, "w+")
                ColSecDefFileID.write(
                    f"section Fiber {ColSecTag * 10 + ctrSec} "
                    f"-torsion {IDTorsion + ctrSec} {{\n"
                )
            else:
                ColSecDefFileID = None

            # Fiber section
            secTag_fiber = ColSecTag * 10 + ctrSec
            ops.section("Fiber", secTag_fiber, "-torsion", IDTorsion + ctrSec)

            # Cover concrete patches
            create_circle_patches(
                matTag=IDconcCover + ctrSec,
                yCenter=yCenter,
                zCenter=zCenter,
                radiusInner=intRadCover,
                radiusOuter=extRadCover,
                minFiber=minFiber,
                maxFiber=maxFiber,
                startAng=0.0,
                endAng=360.0,
                export_file=ColSecDefFileID,
            )

            # Core concrete patches
            create_circle_patches(
                matTag=IDconcCore + ctrSec,
                yCenter=yCenter,
                zCenter=zCenter,
                radiusInner=intRadCore,
                radiusOuter=extRadCore,
                minFiber=minFiber,
                maxFiber=maxFiber,
                startAng=0.0,
                endAng=360.0,
                export_file=ColSecDefFileID,
            )

            # Steel layers in clusters
            for k in range(1, num_bar_clusters + 1):
                startAng = (k - 1) * ((int(numBar / num_bar_clusters) - 1) * theta + theta)
                endAng = startAng + (int(numBar / num_bar_clusters) - 1) * theta

                ops.layer(
                    "circ",
                    IDSteel + ctr,
                    int(numBar / num_bar_clusters),
                    areaBar,
                    yCenter,
                    zCenter,
                    radius,
                    startAng,
                    endAng,
                )

                if write_model_files == 1:
                    ColSecDefFileID.write(
                        f"\tlayer circ {IDSteel + ctr} "
                        f"{int(numBar / num_bar_clusters)} {areaBar} "
                        f"{yCenter} {zCenter} {radius} {startAng} {endAng};\n"
                    )

                ctr += 1

            # Aggregator section
            ops.section(
                "Aggregator",
                ColSecTag + ctrSec,
                IDShear + ctrSec,
                "Vy",
                IDShear + ctrSec,
                "Vz",
                "-section",
                secTag_fiber,
            )

            if write_model_files == 1:
                ColSecDefFileID.write("}\n")
                ColSecDefFileID.write(
                    f"section Aggregator {ColSecTag + ctrSec} "
                    f"{IDShear + ctrSec} Vy {IDShear + ctrSec} Vz "
                    f"-section {secTag_fiber};\n"
                )
                ColSecDefFileID.close()

            ctrSec += 1

    # ----------------------------------------------------------------------
    # Nodes
    # ----------------------------------------------------------------------
    ops.node(2001, 1302.96, 0.0, 0.0)
    ops.node(2002, 1302.96, 0.0, HColumn)
    ops.node(200, 1302.96, 0.0, HColumn + dcg)

    # ----------------------------------------------------------------------
    # Column element generation based on col_elem_type
    # ----------------------------------------------------------------------

    transfTagCol = 2000
    transfTagColDeckConnection = 3000

    # Choose transformation type, for example "PDelta" or "Corotational"
    gT = "Corotational"

    ops.geomTransf(gT, transfTagCol, 0, 1, 0)
    ops.geomTransf(gT, transfTagColDeckConnection, 0, 1, 0)

    # Read column element type from model_attributes
    col_elem_type = attrs["col_elem_type"]

    # Section tags from fiber section generation
    ColSecTag = 2000
    sec1_tag = ColSecTag + 1  # bottom section definition
    sec2_tag = ColSecTag + 2  # middle section definition
    sec3_tag = ColSecTag + 3  # top section definition

    # Node tags
    baseNode = 2001
    topNode = 2002
    deckNode = 200

    # 8-point Lobatto locations and weights
    locs8 = [0.0, 0.06413, 0.20415, 0.39535, 0.60465, 0.79585, 0.93587, 1.0]
    wts8 = [0.017857, 0.10535, 0.17056, 0.20623, 0.20623, 0.17056, 0.10535, 0.017857]

    # 5-point Lobatto locations and weights
    locs5 = [0.0, 0.17267, 0.5, 0.82733, 1.0]
    wts5 = [0.05, 0.27222, 0.35556, 0.27222, 0.05]

    # ----------------------------------------------------------------------
    # Case 1: col_elem_type == "1"
    # Single force-based beam-column with 8 Lobatto points
    # IP 1 -> secdef 1, IP 8 -> secdef 3, IP 2-7 -> secdef 2
    # ----------------------------------------------------------------------
    if col_elem_type == "1":

        # Section tag sequence along height
        sec_tags_8 = [
            sec1_tag,  # IP 1 (bottom)
            sec2_tag,  # IP 2
            sec2_tag,  # IP 3
            sec2_tag,  # IP 4
            sec2_tag,  # IP 5
            sec2_tag,  # IP 6
            sec2_tag,  # IP 7
            sec3_tag  # IP 8 (top)
        ]

        # Integration tag for this column
        integTag = 1

        # Define UserDefined beam integration for the column
        ops.beamIntegration(
            "UserDefined",
            integTag,
            8,
            *sec_tags_8,
            *locs8,
            *wts8
        )

        # Single forceBeamColumn element from base to top
        ops.element(
            "forceBeamColumn",
            baseNode,
            baseNode,
            topNode,
            transfTagCol,
            integTag,
            "-iter", 500, 1.0e-8
        )

    # ----------------------------------------------------------------------
    # Case 2: col_elem_type == "2"
    # Five dispBeamColumn elements in plastic region (Lp1 + Lp2)
    # One dispBeamColumn element in remaining height.
    #
    # 5-point Lobatto integration in each element.
    # Section assignment:
    #   - First element in plastic region:
    #       IP1 -> secdef 1, IP2-5 -> secdef 2
    #   - Remaining plastic region elements:
    #       all 5 IPs -> secdef 2
    #   - Remaining height element:
    #       IP1-4 -> secdef 2, IP5 -> secdef 3
    # ----------------------------------------------------------------------
    elif col_elem_type == "2":

        Lp1 = 0.08 * HColumn
        Lp2 = 0.022 * (fy_e / MPa) * (DbarSingle / MM) * MM
        Lp_plastic = Lp1 + Lp2

        if Lp_plastic > 0.9 * HColumn:
            Lp_plastic = 0.9 * HColumn

        L_remaining = HColumn - Lp_plastic

        n_plastic_elems = 5
        dz_plastic = Lp_plastic / n_plastic_elems

        plastic_nodes = [baseNode]
        for i_seg in range(1, n_plastic_elems + 1):
            z_i = dz_plastic * i_seg
            node_tag = 2100 + i_seg
            ops.node(2100 + i_seg, 1302.96, 0.0, z_i)
            plastic_nodes.append(node_tag)

        sec_tags_plast_first = [sec1_tag, sec2_tag, sec2_tag, sec2_tag, sec2_tag]
        sec_tags_plast_mid = [sec2_tag] * 5
        sec_tags_remain = [sec2_tag, sec2_tag, sec2_tag, sec2_tag, sec3_tag]

        # Create UserDefined integrations (5-point Lobatto) and reuse them
        integTag_first = 101
        integTag_mid = 102
        integTag_remain = 103

        ops.beamIntegration(
            "UserDefined",
            integTag_first,
            5,
            *sec_tags_plast_first,
            *locs5,
            *wts5
        )

        ops.beamIntegration(
            "UserDefined",
            integTag_mid,
            5,
            *sec_tags_plast_mid,
            *locs5,
            *wts5
        )

        ops.beamIntegration(
            "UserDefined",
            integTag_remain,
            5,
            *sec_tags_remain,
            *locs5,
            *wts5
        )

        # Plastic region elements: element tag = iNode
        for e_idx in range(n_plastic_elems):
            nI = plastic_nodes[e_idx]
            nJ = plastic_nodes[e_idx + 1]

            if e_idx == 0:
                integTag = integTag_first
            else:
                integTag = integTag_mid

            ops.element(
                "dispBeamColumn",
                nI,
                nI,
                nJ,
                transfTagCol,
                integTag
            )

        # Remaining height element: element tag = iNode (2105)
        nI = plastic_nodes[-1]
        nJ = topNode

        ops.element(
            "dispBeamColumn",
            nI,
            nI,
            nJ,
            transfTagCol,
            integTag_remain
        )
    else:
        raise ValueError(f"Unknown col_elem_type: {col_elem_type}")

    if write_model_files == 1:
        sec_info_path = f"{model_info_directory}/col_elem_sec_info_col_1.txt"
        with open(sec_info_path, "w") as fsec:

            if col_elem_type == "1":
                # Single forceBeamColumn element
                colEleTag = baseNode
                sec_tags_8 = [
                    sec1_tag,
                    sec2_tag,
                    sec2_tag,
                    sec2_tag,
                    sec2_tag,
                    sec2_tag,
                    sec2_tag,
                    sec3_tag
                ]
                fsec.write(f"{colEleTag} " + " ".join(str(s) for s in sec_tags_8) + "\n")

            elif col_elem_type == "2":
                # Six dispBeamColumn elements: five in plastic region + one in remaining height
                sec_tags_plast_first = [sec1_tag, sec2_tag, sec2_tag, sec2_tag, sec2_tag]
                sec_tags_plast_mid = [sec2_tag, sec2_tag, sec2_tag, sec2_tag, sec2_tag]
                sec_tags_remain = [sec2_tag, sec2_tag, sec2_tag, sec2_tag, sec3_tag]

                # First five elements in plastic region: element tag = iNode (2001, 2101, 2102, 2103, 2104)
                elem_tags = [baseNode, 2101, 2102, 2103, 2104]

                fsec.write(f"{elem_tags[0]} " + ", ".join(str(s) for s in sec_tags_plast_first) + "\n")
                for ele_tag in elem_tags[1:]:
                    fsec.write(f"{ele_tag} " + ", ".join(str(s) for s in sec_tags_plast_mid) + "\n")

                # Remaining height element: element tag = iNode (2105)
                fsec.write(f"{2105} " + ", ".join(str(s) for s in sec_tags_remain) + "\n")

            else:
                raise ValueError(f"Unknown col_elem_type: {col_elem_type}")

    # ----------------------------------------------------------------------
    # Deck connection (common to both element types)
    # ----------------------------------------------------------------------
    # element elasticBeamColumn $jNodeTag $jNodeTag $deckNode 1. $Ubig $Ubig 1. 1. 1. $transfTagColDeckConnection;

    ops.element(
        "elasticBeamColumn",
        topNode,
        topNode,
        deckNode,
        1.0,
        UBIG,
        UBIG,
        1.0,
        1.0,
        1.0,
        transfTagColDeckConnection
    )

    # ----------------------------------------------------------------------
    # Boundary Condition
    # ----------------------------------------------------------------------
    ops.fix(baseNode, 1, 1, 1, 1, 1, 1)

    # ----------------------------------------------------------------------
    # Mass Assignment
    # ----------------------------------------------------------------------

    rotMassSwitchTorsion = 0.0
    rotMassSwitch = 1.0

    def _add_mass(node_tag, mx, my, mz, mr_x=0.0, mr_y=0.0, mr_z=0.0):
        ops.mass(node_tag, mx, my, mz, mr_x, mr_y, mr_z)

    # ----------------------------------------------------------------------
    # Deck node mass (single deck node at the top)
    # ----------------------------------------------------------------------
    # Deck node is 200, and for this reduced model it gets full L_trib
    m_deck = mconc * Adeck * L_trib
    _add_mass(
        200,
        m_deck, m_deck, m_deck,
        rotMassSwitchTorsion * mconc * L_trib * (Iydeck + Izdeck),
        rotMassSwitch * mconc * L_trib * Iydeck,
        rotMassSwitch * mconc * L_trib * Izdeck,
    )

    # ----------------------------------------------------------------------
    # Column node masses (self-mass of the column)
    # ----------------------------------------------------------------------
    if col_elem_type == "1":
        # One element over full height: lump half to base and half to top
        m_col_total = mconc * Acol * HColumn
        m_end = 0.5 * m_col_total

        _add_mass(
            2001,
            m_end, m_end, m_end,
            rotMassSwitch * 0.5 * mconc * HColumn * Iycol,
            rotMassSwitch * 0.5 * mconc * HColumn * Izcol,
            rotMassSwitchTorsion * 0.5 * mconc * HColumn * (Iycol + Izcol),
        )
        _add_mass(
            2002,
            m_end, m_end, m_end,
            rotMassSwitch * 0.5 * mconc * HColumn * Iycol,
            rotMassSwitch * 0.5 * mconc * HColumn * Izcol,
            rotMassSwitchTorsion * 0.5 * mconc * HColumn * (Iycol + Izcol),
        )

    elif col_elem_type == "2":
        # Six elements: 2001-2101-2102-2103-2104-2105-2002
        # Distribute based on actual segment lengths and lump half to each node.

        node_tags = [2001, 2101, 2102, 2103, 2104, 2105, 2002]
        z_coords = [0.0,
                    dz_plastic * 1.0,
                    dz_plastic * 2.0,
                    dz_plastic * 3.0,
                    dz_plastic * 4.0,
                    dz_plastic * 5.0,
                    HColumn]

        # Accumulate translational and rotational masses per node
        m_node = {tag: 0.0 for tag in node_tags}
        rx_node = {tag: 0.0 for tag in node_tags}
        ry_node = {tag: 0.0 for tag in node_tags}
        rz_node = {tag: 0.0 for tag in node_tags}

        for a in range(len(node_tags) - 1):
            nI = node_tags[a]
            nJ = node_tags[a + 1]
            Lseg = z_coords[a + 1] - z_coords[a]

            m_seg = mconc * Acol * Lseg
            m_half = 0.5 * m_seg

            m_node[nI] += m_half
            m_node[nJ] += m_half

            # Rotational masses are zeroed by switches, but keep the structure
            rx_half = rotMassSwitch * 0.5 * mconc * Lseg * Iycol
            ry_half = rotMassSwitch * 0.5 * mconc * Lseg * Izcol
            rz_half = rotMassSwitchTorsion * 0.5 * mconc * Lseg * (Iycol + Izcol)

            rx_node[nI] += rx_half
            rx_node[nJ] += rx_half
            ry_node[nI] += ry_half
            ry_node[nJ] += ry_half
            rz_node[nI] += rz_half
            rz_node[nJ] += rz_half

        # Apply masses to each column node
        for tag in node_tags:
            _add_mass(tag, m_node[tag], m_node[tag], m_node[tag], rx_node[tag], ry_node[tag], rz_node[tag])

    else:
        raise ValueError(f"Unknown col_elem_type: {col_elem_type}")

    # ----------------------------------------------------------------------
    # Gravity Load Assignment
    # ----------------------------------------------------------------------
    IDGravLoadTag = 1

    ops.timeSeries("Linear", IDGravLoadTag)
    ops.pattern("Plain", IDGravLoadTag, IDGravLoadTag)

    # ----------------------------------------------------------------------
    # Deck node gravity load (single deck node at the top)
    # ----------------------------------------------------------------------
    # load = -wconc * Adeck * L_trib in global Z
    Pdeck = -wconc * Adeck * L_trib
    ops.load(200, 0.0, 0.0, Pdeck, 0.0, 0.0, 0.0)

    # ----------------------------------------------------------------------
    # Column gravity loads (self-weight)
    # ----------------------------------------------------------------------
    if col_elem_type == "1":
        # One element over full height: half weight to base and half to top
        Pcol_total = -wconc * Acol * HColumn
        Pend = 0.5 * Pcol_total

        ops.load(2001, 0.0, 0.0, Pend, 0.0, 0.0, 0.0)
        ops.load(2002, 0.0, 0.0, Pend, 0.0, 0.0, 0.0)

    elif col_elem_type == "2":
        # Six elements: 2001-2101-2102-2103-2104-2105-2002
        node_tags = [2001, 2101, 2102, 2103, 2104, 2105, 2002]
        z_coords = [0.0,
                    dz_plastic * 1.0,
                    dz_plastic * 2.0,
                    dz_plastic * 3.0,
                    dz_plastic * 4.0,
                    dz_plastic * 5.0,
                    HColumn]

        # Accumulate nodal vertical loads
        P_node = {tag: 0.0 for tag in node_tags}

        for a in range(len(node_tags) - 1):
            nI = node_tags[a]
            nJ = node_tags[a + 1]
            Lseg = z_coords[a + 1] - z_coords[a]

            Pseg = -wconc * Acol * Lseg
            Phalf = 0.5 * Pseg

            P_node[nI] += Phalf
            P_node[nJ] += Phalf

        # Apply loads to nodes
        for tag in node_tags:
            ops.load(tag, 0.0, 0.0, P_node[tag], 0.0, 0.0, 0.0)

    else:
        raise ValueError(f"Unknown col_elem_type: {col_elem_type}")


def analysis_gravity(model_input_dict):
    write_model_files = model_input_dict.get("write_model_files", 0)
    model_info_directory = model_input_dict.get("model_info_directory", ".")
    attrs = model_input_dict["model_attributes"]

    col_elem_type = attrs["col_elem_type"]

    if col_elem_type == "1":
        ele_first = 2001
        ele_last = 2001
    elif col_elem_type == "2":
        ele_first = 2001
        ele_last = 2105
    else:
        raise ValueError(f"Unknown col_elem_type: {col_elem_type}")

    nStepsGrav = 10
    ops.test("EnergyIncr", 1e-20, 2000, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nStepsGrav)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.system("ProfileSPD")
    ops.analysis("Static")

    ok = ops.analyze(nStepsGrav)

    ops.loadConst("-time", 0.0)

    if ok == 0:
        print("Gravity Analysis COMPLETE!")
    else:
        print("Gravity Analysis FAILED!")

    # Pull axial force directly from element localForce
    # Tcl used index 1 of the localForce vector
    lf_first = ops.eleResponse(ele_first, "localForce")
    axial_first = abs(lf_first[0])

    # Append to edge 1 predictor file only (cantilever)
    if write_model_files == 1:
        pred_path = os.path.join(
            model_info_directory,
            "predictor_info_col_rebar_strain_damage_col_1_edge_1.txt"
        )
        with open(pred_path, "a") as pf:
            pf.write(f"{axial_first}\n")

    ops.wipeAnalysis()
    return ok == 0

def analysis_prelim(model_input_dict):
    modes = int(model_input_dict["modes"])
    model_info_directory = model_input_dict.get("model_info_directory", ".")

    # Eigenvalues (lambda = omega^2)
    lam = ops.eigen(modes)

    # Convert to numpy array for vectorized ops
    lam = np.asarray(lam, dtype=float)
    omega = np.sqrt(lam)
    periods = (2.0 * np.pi) / omega  # T = 2*pi/sqrt(lambda)

    # Save periods in one line (Tcl wrote "$T" which prints the list on one line)
    os.makedirs(model_info_directory, exist_ok=True)
    with open(os.path.join(model_info_directory, "periods.txt"), "w") as f:
        f.write(" ".join(f"{T:.16g}" for T in periods))

    # Save mode shapes: one file per mode, one line per node: "<node> <evec>"
    all_nodes = ops.getNodeTags()
    for i_mode in range(1, modes + 1):
        out_path = os.path.join(model_info_directory, f"mode_shape_{i_mode}.txt")
        with open(out_path, "w") as f:
            for node in all_nodes:
                # nodeEigenvector returns a list of length ndf for the node
                vec = ops.nodeEigenvector(node, i_mode)
                f.write(f"{node} " + " ".join(f"{v:.16g}" for v in vec) + "\n")

    print("Modal Analysis COMPLETE!")


# ======================================================================
# Global NLTHA control parameters
# ======================================================================

# Ground-motion preprocessing
ADD_ZEROS_FOR = 1.0          # seconds of zero-padding after GM
SCALE_FAC_L = 1.0            # longitudinal scale factor
SCALE_FAC_T = 1.0            # transverse scale factor
SCALE_FAC_V = 1.0            # vertical scale factor
DT_ANALYSIS = 0.001          # analysis time step (sec)
GM_SKEW_DEG = 0.0            # skew angle (degrees)

# Solver and test configuration
ALGORITHM_BASIC = ["NewtonLineSearch"]
# ALGORITHM_BASIC = ["ModifiedNewton", "-initial"]
TEST_BASIC = "RelativeNormDispIncr"
SHOW_TEST_BASIC = 0
SHOW_TEST = 0

# Convergence tolerances
TOL_DYN_BASIC = 1e-3
TOL_DYN_DISP = 1e-3
TOL_DYN_UNB = 1e-2

# Iteration limits
MAX_NUM_ITER_DYN_BASIC = 500
MAX_NUM_ITER_DYN = 2500
MAX_DIM_KRYLOV = 50


def analysis_time_hist(model_input_dict) -> bool:
    """
    Run nonlinear time history analysis (NLTHA) using OpenSeesPy.

    Expected in model_input_dict:
        {
            'write_model_files': 0,
            'model_param_vals': model_param_vals,
            'model_attributes': model_attributes,
            'gm_file_loc': '.',
        }

    Ground motion files (in gm_file_loc):
        gm_1.txt  (longitudinal)
        gm_2.txt  (transverse)
        gm_3.txt  (vertical)

    File format (whitespace-separated):
        npts dt a1 a2 ... anpts
    where accelerations are in units of g (dimensionless).
    """
    attrs = model_input_dict.get("model_attributes", {}) or {}
    gm_file_loc = model_input_dict.get("gm_file_loc", ".")

    # Direction IDs for UniformExcitation patterns
    gm_dirctn_L, gm_dirctn_T, gm_dirctn_V = 1, 2, 3

    # Vertical ground motion switch
    vert = int(attrs.get("vertical_gm", 0))

    # Gravitational acceleration (in/s^2); uses global units if available
    g = globals().get("G", 386.4)

    def _read_gm(path: str) -> tuple[int, float, np.ndarray]:
        with open(path, "r") as f:
            toks = f.read().strip().split()
        if len(toks) < 3:
            raise ValueError(f"GM file {path} does not have enough data.")
        npts = int(float(toks[0]))
        dt = float(toks[1])
        vals = np.array([float(v) for v in toks[2:2 + npts]], dtype=float)
        if vals.size != npts:
            raise ValueError(f"GM file {path} expected {npts} points, got {vals.size}.")
        return npts, dt, vals

    def _append_zeros(arr: np.ndarray, dt_rec: float) -> np.ndarray:
        n0 = int(ADD_ZEROS_FOR / dt_rec)
        if n0 <= 0:
            return arr
        return np.concatenate([arr, np.zeros(n0, dtype=float)])

    # ------------------------------
    # Read and assemble GM vectors
    # ------------------------------
    _, dtRec_L, acc_L = _read_gm(os.path.join(gm_file_loc, "gm_1.txt"))
    _, dtRec_T, acc_T = _read_gm(os.path.join(gm_file_loc, "gm_2.txt"))

    # Vertical file is read only if needed
    if vert == 1:
        _, dtRec_V, acc_V = _read_gm(os.path.join(gm_file_loc, "gm_3.txt"))
    else:
        dtRec_V, acc_V = dtRec_L, np.zeros_like(acc_L)

    # Apply skew rotation between L and T
    c = math.cos(GM_SKEW_DEG * math.pi / 180.0)
    s = math.sin(GM_SKEW_DEG * math.pi / 180.0)
    gmInput_L = acc_L * c + (-acc_T) * s
    gmInput_T = acc_L * s + acc_T * c
    gmInput_V = acc_V.copy()

    # Add trailing zeros
    gmInput_L = _append_zeros(gmInput_L, dtRec_L)
    gmInput_T = _append_zeros(gmInput_T, dtRec_T)
    if vert == 1:
        gmInput_V = _append_zeros(gmInput_V, dtRec_V)

    # Determine analysis duration (use the shortest component duration)
    total_times = [
        len(gmInput_L) * dtRec_L,
        len(gmInput_T) * dtRec_T,
    ]
    if vert == 1:
        total_times.append(len(gmInput_V) * dtRec_V)
    totalAnalysisTime = min(total_times)

    # ------------------------------
    # Define uniform excitation
    # ------------------------------
    gmLoadTag_L, gmLoadTag_T, gmLoadTag_V = 2, 3, 4
    tsTag_L, tsTag_T, tsTag_V = gmLoadTag_L, gmLoadTag_T, gmLoadTag_V

    gmFact_L = g * SCALE_FAC_L
    gmFact_T = g * SCALE_FAC_T
    gmFact_V = g * SCALE_FAC_V

    ops.timeSeries("Path", tsTag_L, "-dt", dtRec_L, "-values", *gmInput_L.tolist(), "-factor", gmFact_L)
    ops.timeSeries("Path", tsTag_T, "-dt", dtRec_T, "-values", *gmInput_T.tolist(), "-factor", gmFact_T)
    if vert == 1:
        ops.timeSeries("Path", tsTag_V, "-dt", dtRec_V, "-values", *gmInput_V.tolist(), "-factor", gmFact_V)

    ops.pattern("UniformExcitation", gmLoadTag_L, gm_dirctn_L, "-accel", tsTag_L)
    ops.pattern("UniformExcitation", gmLoadTag_T, gm_dirctn_T, "-accel", tsTag_T)
    if vert == 1:
        ops.pattern("UniformExcitation", gmLoadTag_V, gm_dirctn_V, "-accel", tsTag_V)

    # ------------------------------
    # Analysis configuration
    # ------------------------------
    ops.test(TEST_BASIC, TOL_DYN_BASIC, int(MAX_NUM_ITER_DYN_BASIC), int(SHOW_TEST_BASIC))
    ops.algorithm(*ALGORITHM_BASIC)
    ops.integrator("Newmark", 0.50, 0.25)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.system("UmfPack")

    try:
        ops.analysis("VariableTransient", "-numSublevels", 10, "-numSubSteps", 5)
    except Exception:
        ops.analysis("Transient")

    # ------------------------------
    # Run analysis (with fallback strategies)
    # ------------------------------
    nSteps = int(totalAnalysisTime / DT_ANALYSIS)
    ok = ops.analyze(nSteps, DT_ANALYSIS)

    if ok != 0:
        tCurrent = ops.getTime()
        ok_step = 0  # 0 means success in OpenSees; nonzero means failure

        while (ok_step == 0) and (tCurrent <= totalAnalysisTime):
            ok_step = ops.analyze(1, DT_ANALYSIS)

            # Fallback 1: KrylovNewton + displacement increment test
            if ok_step != 0:
                ops.test("NormDispIncr", TOL_DYN_DISP, int(MAX_NUM_ITER_DYN), int(SHOW_TEST))
                ops.algorithm("KrylovNewton", "-maxDim", int(MAX_DIM_KRYLOV))
                ok_step = ops.analyze(1, DT_ANALYSIS)
                if ok_step == 0:
                    ops.test(TEST_BASIC, TOL_DYN_BASIC, int(MAX_NUM_ITER_DYN_BASIC), int(SHOW_TEST_BASIC))
                    ops.algorithm(*ALGORITHM_BASIC)

            # Fallback 2: Newton -initial + displacement increment test
            if ok_step != 0:
                ops.test("NormDispIncr", TOL_DYN_DISP, int(MAX_NUM_ITER_DYN), int(SHOW_TEST))
                ops.algorithm("Newton", "-initial")
                ok_step = ops.analyze(1, DT_ANALYSIS)
                if ok_step == 0:
                    ops.test(TEST_BASIC, TOL_DYN_BASIC, int(MAX_NUM_ITER_DYN_BASIC), int(SHOW_TEST_BASIC))
                    ops.algorithm(*ALGORITHM_BASIC)

            # Fallback 3: Newton -Hall + unbalance test
            if ok_step != 0:
                ops.test("NormUnbalance", TOL_DYN_UNB, int(MAX_NUM_ITER_DYN), int(SHOW_TEST))
                ops.algorithm("Newton", "-Hall", 0.5, 0.5)
                ok_step = ops.analyze(1, DT_ANALYSIS)
                if ok_step == 0:
                    ops.test(TEST_BASIC, TOL_DYN_BASIC, int(MAX_NUM_ITER_DYN_BASIC), int(SHOW_TEST_BASIC))
                    ops.algorithm(*ALGORITHM_BASIC)

            # Fallback 4: NewtonLineSearch + unbalance test
            if ok_step != 0:
                ops.test("NormUnbalance", TOL_DYN_UNB, int(MAX_NUM_ITER_DYN), int(SHOW_TEST))
                ops.algorithm("NewtonLineSearch")
                ok_step = ops.analyze(1, DT_ANALYSIS)
                if ok_step == 0:
                    ops.test(TEST_BASIC, TOL_DYN_BASIC, int(MAX_NUM_ITER_DYN_BASIC), int(SHOW_TEST_BASIC))
                    ops.algorithm(*ALGORITHM_BASIC)

            tCurrent = ops.getTime()

        ok = ok_step

    # ------------------------------
    # Final status
    # ------------------------------
    if ok == 0:
        tCurrent = ops.getTime()
        nltha_ok = (tCurrent / totalAnalysisTime) >= 0.85
    else:
        nltha_ok = False

    ops.wipeAnalysis()
    return bool(nltha_ok)