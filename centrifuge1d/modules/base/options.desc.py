# This file contains the description of base options

CONFIG_OPTIONS = {
    'general': {'g': "Gravitational constant. Default is 981 cm/s. Type:float."},
    'experiment': \
        {'exp_type': ("Identifier of type of experiment do determine which "
                      "module should be used. Type: string"),
         'tube_no': ("Tube number. The given tube number should correspond to "
                     "an entry in the 'constants.ini' file, where starting "
                     "and ending filter proberties should be specified"
                     "Type: integer"),
         'tube_diam': "Diameter of the tube.",
         're': ("Distance (radius) from the centrifuge axis to the end of the "
                "ending filter; on the other hand, value r0 is the distance "
                "the centrifuge axis to the beginning of soil sample, i.e. "
                "re=r0 + L + fl2, where L is soil length (see 'l0') and "
                "'fl2' is the thickness of the ending filter. Either 're' "
                "or 'r0' has to be specified (but not both). Type: float"),
         'r0': "See parameter re",
         'r0_fall': ("When simulating falling head test in gravitational field "
                     "we use centrifuge with a very long arm that rotates at "
                     "slow speed such, that the effect is similar to the g "
                     "level. Corresponding rotational speed is calculated "
                     "automatically. Type: float"),
         'ks': ("Saturated hydraulic conductivity of the soil sample. "
                "Type: float"),
         'l0': "Soil sample (initial) length. Type: float or array of floats.",
         'wt_out': ("Distance from 're' to the basis of the outflow cup. "
                    "Used for computing the force acting on the outflow cup. "
                    "If force is not computed, it can be any value (as it is) "
                    "ignored. Type: float"),
         'wl0': ("Length of water on the inflow (above the soil) in cm. "
                 "Type: float"),
         'ww0': ("Weight of water on the inflow (above the soil) in gramms. "
                 "Only one of 'wl0' and 'ww0' can be specified."
                 "Type: float"),
         'descr': "(optional) Description of the experiment. Type: string"
        },
    'filters': \
        # These values are set based on the value of 'tube_no'. See 'tube_no'
        # for more.
        {'fl1': "Length (thickness) of the starting filter. Type: float",
         'fl2': "Length (thickness) of the ending filter. Type: float",
         'ks1': ("Saturated hydraulic conductivity of starting filter [cm/s]."
                 "Type: float"),
         'ks2': ("Saturated hydraulic conductivity of ending filter [cm/s]."
                 "Type: float"),

        },
    'fluid': \
    # we use water, so these are probably not needed to be changed
        {'density': "Fluid density. For water it is ~1.0 g/cm3. Type: float."},
    'measurements': \
        {'include_acceleration': ("Flag whether to include acceleration and/or "
                                  "deceleration into simulation. If True, "
                                  "acceleration is included and deceleration "
                                  "only if also 'deceleration_duration' > 0. "
                                  "See also 'deceleration_duration'. "
                                  "Type: boolean"),
         'duration': ("Duration of a phase for omega > 0. If the flag "
                      "'include_acceleration' is True, also centrifuge "
                      "acceleration is simulated during this phase."
                      "One iteration consists of an optional centrifugation "
                      "phase followed by an optional deceleration phase and an "
                      "also optional gravitational phase. In each iteration at "
                      "least one phase has to have duration > 0 (otherwise "
                      "there is not much to simulate:D), whereas deceleration "
                      "phase duration > 0 has obviously sense only if a "
                      "centrifugal phase duration is also > 0."
                      "Type: float or array of floats."),
         'deceleration_duration': ("Duration of a deceleration phase.If the "
                                   "flag 'include_acceleration' is True, then "
                                   "also deceleration is simulated. Otherwise "
                                   "this value is ignored. "
                                   "Type: same as 'duration'."),
         'fh_duration': ("Duration of a phase for omega = 0, i.e. only "
                         "gravitational force is applied. Used when sample was "
                         "for example left some time outside the centrifuge "
                         "or for falling-head test simulation. See also option "
                         "'r0_fall'. Type: same as 'duration'."),
         'measurements_keep',: "See 'measurements_remove'.",
         'measurements_remove': \
             ("For filtering measurements 'measurements_keep' (further 'keep') "
              "and 'measurements_remove' (further 'remove') are used. They are "
              "of type dict with keys being the measurements names and values "
              "is a list of indices. If 'keep' is supplied, indices specified "
              "will be preserved whereas the rest is dropped. On the other "
              "if 'remove' is supplied, indices specified in there will be "
              "removed an the rest is preserved. These two parameters are not "
              "complementary, in fact what can be achieved using one, can be "
              "achived also using only the other. But sometimes it is more "
              "convenient to specify it one way than the other."),
         'smoothing': ("Measured data can be 'smoothed', which may improve "
                       "found results. Type: dict with key being measurement "
                       "name and value is one of 'smlin' - linear averaging, "
                       "'smtri' - triangular averaging, 'smgau' - gaussian "
                       "averaging."),
         'l1': ("Soil sample length at measured time. "
                "Type: array of floats or None"),
         'wl1': ("Measured length of water above the soil at given (measured) "
                 "time in cm. Type: array of floats or None"),
         'ww1': ("Measured weight of water above the soil at given (measured) "
                 "time, in gramms. Type: array of floats or None"),
         'wl_out': ("Measured length of water in the outflow chamber at given "
                    "(measured) time. Units: array of floats or None"),
         'gc1': "Measured gravitational center. Units: array of floats or None",
         'rm1': "Measured rotational momentum. Units: array of floats or None",
         'measurements_scale_coefs': \
             ("When running inverse problem, multiple meaurements can be "
              "given, potentially with different scales (e.g. MO ~ 0.5 whereas "
              "GC ~ 2.5, which can cause that solver \"prefers\" the "
              "optimization of GC in favour of MO, because it introduces adds "
              "more error in least-squares. Therefore a scaling is important "
              "to make the data of roughly the same order.) By default the "
              "data is scaled so that the biggest number in measurements of "
              "given type is in <1, 10) interval. See also *_weights, options, "
              "which specify a weight - i.e. importance of given measurement."
              "Value: dict of type: {measurement_name1: scale_coef1, ...}."),
         'gf_mo': ("Measured centrifugal force of the expelled water. More "
                   "precisely it is the force divided by g (= gravitational "
                   "constant)."),
         'gf_mt': ("Measured centrifugal force of the water inside the tube. "
                   "More precisely it is the force divided by g "
                   "(= gravitational constant)."),
         'gf_mo_tara': ("Force implied on the sensor of expelled water by the "
                        "holding aparatus. Two possible values are supported: "
                        "a) list of two items: (omega, gF_tara) where omega is "
                        "the speed (in rpm) at which the force gF_tara was "
                        "measured - force [gF_tara] = 1e-3 kgf (i.e as if "
                        "measured under 1g: gF_tara = F_tara/g with F the "
                        "centrifugal force. When this values is supplied, "
                        "the gF_tara for other speeds is calculated dynamically "
                        "at runtime."
                        "\nb) a dictiory of form: {omega1: gF_tara1, omega2: "
                        "gF_tara2, ...}. In this case forces are taken as "
                        "static, so no dynamic adaptation is supported (which "
                        "implicitly means that gF_tara in the acceleration/"
                        "deceleration phase is wrong)"),
         'gf_mt_tara': ("Force implied on the sensor measuring the water "
                        "inside of the tube. See also 'f_mo_tara'."),
         'gf_mo_calibration_curve':\
             ("Calibration curve for measured f_mo. Useful when there is some "
              "shift in measured data (e.g. in case of calibration with base "
              "not equal to zero). So f_mo:= f_mo - f_mo_calibration_curve. "
              "Calibration curve can be either constant, list (of the same "
              "length as 'f_mo') or a dict of format: "
              "{omega1: base1, omega2: base2, ...}. In the last case the "
              "measurements will be shifted accordin to the rotational speed.\n"
              "Set to None to compute with differences between measurements, "
              "instead of absolute values."),
         'gf_mt_calibration_curve': "See 'f_mo_calibration_curve'.",
         'xvalues': ("Specify the 'xvalues' of given measurement, mostly it's "
                     "it's times at which measurements were taken. 'xvalues' "
                     "are indicated as option '[option_name]_xvalues', e.g. "
                     "'gF_MO_xvalues'. If xvalues are not specified, it is "
                     "implicitly assumed that measurement was taken at the end "
                     "of one run/cycle (consisting of centrifugation phase, "
                     "deceleration phase and g-phase)."
                     "This is kept mainly for backward compatibility."),
         'mo_gc_calibration_curve': \
             ("Calibration curve for the expelled water. Allows to determine "
              "the gravitational centre of the expelled water. Value is a list "
              "of pairs of the format: [(mo1, R1), (mo2, R2), ...], where "
              "mo1 < mo2 < ... and moX (X=1,2,...) is the amount of expelled "
              "water (in cm3) and RX is the radius of GC corresponding to the "
              "moX units of expelled water. The value of R for mo: "
              "moX < mo < moY is linearly interpolated.")
         },
    'solver' : \
      {'atol': ("Set absolute tolerances for the solver. "
                "Type: float or array of floats"),
       'rtol': "Set relative tolerances for the solver. Type: float",
       'first_step_size': ("The size of the first step of the solver (needed "
                           "for the computation of initial condition, as the "
                           "solver can crash if it's too big, but can prolong "
                           "computation time significantly when set too small. "
                           "Usually a value of 1e-8 is sufficient. Type: float"),
        'max_steps': ("Maximal allowed steps the solver will perform. "
                      "Type: integer"),
        'max_step_size': ("Maximal size of one step the solver performs. "
                          "Type: float"),
        'compute_initcond': ("Value of compute_initcond for the solver. This"
                              " can be one of None: don't compute init condition,"
                              " 'yp0': calculates/improves the yp0 initial "
                              "condition (considering y0 to be accurate), or "
                              "'y0' : calculates/improves the y0 initial "
                              "condition (considering yp0 to be accurate)."
                            ),
        'always_restart_solver': ("If the initial condition of the next step "
                                  "is the same as the final state of previous "
                                  "step, the solver just continues - if set to "
                                  "'True',  the solver always calculates the "
                                  "initial condition for the next step. "
                                  "Default is False. Type: boolean")
        },
    'display': \
        {'show_figures': "Draw graphs after computation. Type: boolean",
         'save_figures': "Save graphs to file. Type: boolean",
         'separate_figures': ("Each graph will be drawn in it's own window, "
                              "instead of drawing all of then in one. "
                              "Type: boolean"),
         'save_data': "Saves data to file. Type: boolean",
         'verbosity': ("Specifies how much data will be displayed during "
                       "computation. Default is 1. Value of 0 does not display "
                       "anything. Greater values show more. Type: integer"),
         'params_ref': ("Referencing parameters. The results of the same used "
                        "model but with these changed parameters will be "
                        "displayed in graphs together with computed results."
                        "Type: dictionary")
        },
     'INFORMATIVE': \
     {'measurements_length': ("This option is an integer determining the "
                              "length of measurements and is generated by "
                              "mkini. Is only to inform the user and therefore "
                              "it's ignored.")}
    }

INTERNAL_OPTIONS = {
    'omega2g_fns': ("Specifies the functions used during each phases. See "
                    "'find_omega2g'. This option is used only internally "
                    "and is not directly exposed to the user code. The "
                    "solver internally takes corresponding 'find_omega2g' "
                    "function depending on the current phase. Internally "
                    "it's a dictionary of items phase_name:find_omega2g. "
                    "Type: dict"),
    'find_omega2g': ("Function that returns the value of omega^2/g. When in "
                     "acceleration phase, it should take this into account "
                     "based on the time 't0' when acceleration started and "
                     "'omega_start' the starting speed. This function is "
                     "meant to be used by used inside the code. It takes "
                     "'t' as actual time and 'model' parameters (from which "
                     "it get's all other needed parameters). See also: t0 "
                     "and omega _start. "
                     "Type: function(t, model) -> returns omega^2/g"),
    't0': ("Initial time when given phase started (elapsed time from "
           "experiment start, in seconds). This value is used internaly "
           "by the find_omega2g() functions and is probably not needed for "
           "the user otherwise. It value is set internally by the solver "
           "on each phase start. See also: 'find_omega2g'. Type: float"),
    'omega_start': ("Internal value meant to be used during acceleration "
                    "and deceleration phase. For user this value is probably "
                    "not needed. In the acceleration phase it carries the "
                    "information of 'omega' (rotational) speed of previous "
                    "phase from which we start to accelerate until we reach "
                    "the rotational speed 'omega' desired. By analogy, during "
                    "deceleration phase it holds the information of the "
                    "rotational speed of previous phase, at which we started "
                    "to decelerate. This value is set internally by the "
                    "solver on phase change. Type: float")
  }
