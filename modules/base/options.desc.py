# This file contains the description of base options

CONFIG_OPTIONS = {
    'general': {'g': "Gravitational constant. Default is 981 cm/s. Type:float."},
    'experiment': \
        {'exp_type': ("Identifier of type of experiment do determine which "
                      "module should be used. Type: string")
         'tube_no': ("Tube number. To given tube number should correspong "
                     "an entry in the 'constants.ini' file, where starting "
                     "and ending filter proberties should be specified".
                     "Type: integer"),
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
                 "Type: float")
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
        {'density': "Fluid density. For water it is ~1.0 g/cm3. Type: float."}
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
                      "centrifugal phase duration is also > 0".
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
         'measurements_times': ("Times (in seconds) at which we want to make"
                                "measurement (or at which measurement was "
                                "taken). If not specified it is by default "
                                "taken at the end of one iteration (i.e. after "
                                "one run - consisting of centrifugation phase, "
                                "decelerating phase and a g-phase)."
                                "Type: array of floats or None."),
         'l1': ("Soil sample length at measured time. "
                "Type: array of floats or None"),
         'wl1': ("Measured length of water above the soil at given (measured) "
                 "time. Units: array of floats or None"),
         'wl_out': ("Measured length of water in the outflow chamber at given "
                    "(measured) time. Units: array of floats or None"),
         'gc1': "Measured gravitational center. Units: array of floats or None",
         'rm1': "Measured rotational momentum. Units: array of floats or None",
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
        }
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
                     "'t' as actual time and 'model' parameters (from which
                     it get's all other needed parameters). See also: t0 "
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
