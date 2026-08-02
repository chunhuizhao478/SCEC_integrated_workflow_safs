# SAFS reproduction exercise — results

Target deck: `safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_intermediate_plast_phi30_40_gradedfw_k1p70_nwredM7p8_sefw0_attenuation_deep40km`

## Ladder

| Rung | Verdict | Evidence |
|:--|:--|:--|
| E0 | **PASS** | safalt_0d5Hz_p3_deep40km.puml.h5: 3,210,006 nodes, 15,979,903 tets, 320,560 BC-3 faces, 320,560 of them interior (33.2 s) |
| E1 | **PASS** | 37 CVM slices, 105 CTM slices, CSM present; shipped material nc: 452x361x194, z[-45000,3250] dz=250; CVM slice count 37 vs nc z-levels 194 -- the nc is RESAMPLED onto a uniform axis, so these need not be equal |
| E2 | **PASS** | 37 fields recovered, 0 unknown: [] |
| E2c | **PASS** | shallow freeze measured ON THE FIELD: levels above -500 m are DISTINCT (freeze OFF); the filename said 0.0 |
| E3 | **BLOCKED** | the ['cvm_slices'] reader is not implemented, so the CVM cannot be rebuilt from raw.  Every downstream rung inherits this: E4 needs E3's Sv, E5 needs E3's thermal nc.  This is the honest state -- the reader raises rather than silently producing a wrong field. |
| E4 | **BLOCKED** | inherits E3: E4 consumes a material product |
| E5 | **BLOCKED** | inherits E3: E5 consumes a material product |
| E6 | **BLOCKED** | inherits E3: E6 consumes a material product |

## Recovered parameters

```
    EndTime              : 150.0
    FL                   : 103
    FreqCentral          : 0.5
    FreqRatio            : 100.0
    Plasticity           : 1
    Tv                   : 0.05
    deck                 : safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_intermediate_plast_phi30_40_gradedfw_k1p70_nwredM7p8_sefw0_attenuation_deep40km
    freeze_from_filename : 0.0
    freeze_note          : no _freeze tag in the name -> the shallow freeze is OFF for this deck; VERIFY against the field, do not trust the name
    friction_grid        : x[303000,696000] y[3612000,3901500] z[-21000,200] dx=1500 dz=200
    friction_nc          : safs_friction_thermal_case1.nc
    fw_base              : 0.0
    fw_boundaries_s_km   : [(20.0, 28.0), (150.0, 160.0), (176.0, 182.0), (188.0, 194.0), (214.0, 226.0), (244.0, 252.0)]
    fw_n_regions         : 7
    fw_values            : [0.0, 0.0, 0.045, 0.03, 0.06, 0.0175, 0.05]
    hypocenter_xyz       : (604446.944, 3704576.3853, -10067.9819)
    k_from_filename      : 1.7
    material_grid        : x[28500,705000] y[3543000,4083000] z[-45000,3250] dx=1500 dz=250
    material_nc          : safs_material_cvm.nc
    mesh                 : safalt_0d5Hz_p3_deep40km.puml.h5
    nucleation_amplitude : 75000000.0
    nucleation_radius_m  : 2000.0
    phi_deg_from_field   : [30.0, 40.0]
    phi_note             : read from bulkFriction = tan(phi), NOT from the yaml comment; a deck's prose has been wrong before
    plasticity_grid      : x[28500,705000] y[3543000,4083000] z[-45000,3250] dx=1500 dz=250
    plasticity_nc        : safs_plasticity_phi30_40.nc
    rs_a_constant        : 0.0334
    rs_b_constant        : 0.019
    rs_sl0_constant      : 0.1
    rs_srW_constant      : 1000.0
    stress_box           : {'xmin': 350000.0, 'xmax': 630000.0, 'ymin': 3680000.0, 'ymax': 3850000.0, 'zmin': -20000.0, 'zmax': 0.0, 'dx': 1000.0, 'dz': 250.0}
    stress_grid          : x[350000,630000] y[3680000,3850000] z[-20000,0] dx=1000 dz=250
    stress_nc            : safs_stress_andersonian_k1.7.nc
    strike_azimuth_deg   : 314.0
    strike_origin_xy     : (606971.0, 3707270.0)
    thermal_grid         : x[303000,696000] y[3612000,3901500] z[-21000,200] dx=1500 dz=200
    thermal_nc           : safs_friction_thermal_case1.nc
```

## Blocking issue

`cvm_slices` and `ctm_slices` are not implemented in `deckbuild/material.py`. They raise `MaterialError('not implemented yet')` rather than silently producing a wrong field, which is the right interim state, but it means E3–E6 cannot run and the reproduction claim is **unproven**.

The mesh (E0), the raw-data provenance (E1) and the descriptor recovery (E2) all pass, so the exercise establishes that the deck is readable and its design is fully recoverable — which is the prerequisite for E3–E6, not a substitute for them.
