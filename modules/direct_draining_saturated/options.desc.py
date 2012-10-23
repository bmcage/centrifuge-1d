CONFIG_OPTIONS = {
    'soil': \
        {'porosity': "Porosity of the soil sample. Type: float",
         'n': ("Van Genuchten parameter 'n' of the soil. "
               "Type: float or list"),
               'gamma': ("Van Genuchten parameter 'gamma' of the soil. "
                         "Type: float or list")},
    'initial-and-boundary-conditions': \
        {'h_init': ("Initial fluid head value (close to 0) for the "
                    "(almost fully saturated) sample.")
                    " soil sample (e.g. end of the starting-filter)."
                    "Type: float or list"),
         'rb_type': ("Right boundary condition of type:"
                     "\n\t 0 - no outflow [q_last = 0]"
                     "\n\t 1 - free outflow"
                     "\n\t 2 - prescribed pressure (see h_last)"
                     "\nType: integer"),
         'h_last': ("The prescribed pressure on the right boundary - only for "
                    "'rb_type' = 2")
     },
    'solver': \
        {'inner_points': ("How many (inner) point will be used for the sample "
                          "discretization. Type: integer"),
         'dtype': ("Discretization type. Determines how the soil sample "
                   "will be discretized in space. Currently implemented are: "
                   "\n\t 1 - linear discretization"
                   "\n\t 2 - the distance of two points is 'k_dx'*(dist of "
                   "previous two points. See also: k_dx)"
                   "\nType: integer"),


        }
    }
