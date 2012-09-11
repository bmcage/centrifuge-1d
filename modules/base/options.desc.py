CONFIG_OPTIONS = {
    'starting-filter': \
        {'n2' : ("(optional) Van Genuchten parameter 'n' of the filter. "
                 "Type: float or list"),
         'gamma2': ("(optional) Van Genuchten parameter 'gamma' of the filter. "
                 "Type: float or list")},
    'soil': \
        {'n' : ("(optional) Van Genuchten parameter 'n' of the soil. "
                 "Type: float or list"),
         'gamma': ("(optional) Van Genuchten parameter 'gamma' of the soil. "
                 "Type: float or list")},
    'ending-filter': \
        {'n2' : ("(optional) Van Genuchten parameter 'n' of the filter. "
                 "Type: float or list"),
         'gamma2': ("(optional) Van Genuchten parameter 'gamma' of the filter. "
                 "Type: float or list")},
    'set-up': \
        {'h_init'    : ("Initial fluid head value (close to 0) for the "
                        "(almost fully saturated) sample.")
                    " soil sample (e.g. end of the starting-filter)."
                    "Type: float or list"),
         'dtype'     : ("Discretization type. Determines how the soil sample "
                        "will be discretized in space. Type: integer"),
         'inner_point': ("How many (inner) point will be used for the sample "
                         "discretization. Type: integer")},
    'experiment': \
        {'show_figures': "Draw graphs after computation. Type: boolean",
         'descr': "(optional) Description of the experiment. Type: string"}
    'solver' : \
      {'atol': "Set absolute tolerances for the solver",
       'rtol': "Set relative tolerances for the solver",
       'first_step_size': ("The size of the first step of the solver (needed "
                           "for the computation of initial condition, as the "
                           "solver can crash if it's too big, but can prolong"
                           "computation time significantly when set too small)",
        'max_steps': "Maximal allowed steps the solver will perform",
        'max_step_size': "Maximal size of one step the solver performs."}
        'always_restart_solver': ("If the initial condition of the next step"
                                  "is the same as the final state of previous"
                                  "step, the solver just continues - if set to"
                                  "'True',  the solver always calculates the"
                                  "initial condition for the next step."
    }
