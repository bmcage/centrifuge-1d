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
                        "(almost fully saturate)d sample.")
                    " soil sample (e.g. end of the starting-filter)."
                    "Type: float or list"),
         'dtype'     : ("Discretization type. Determines how the soil sample "
                        "will be discretized in space. Type: integer"),
         'inner_point': ("How many (inner) point will be used for the sample "
                         "discretization. Type: integer")},
    'experiment': \
        {'draw_graphs': "Draw graphs after computation. Type: boolean",
         'descr': "(optional) Description of the experiment. Type: string"}
    }
