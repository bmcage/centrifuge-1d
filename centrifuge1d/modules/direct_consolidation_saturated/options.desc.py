CONFIG_OPTIONS = {
    'soil': \
        {'porosity': "Start porosity of the soil sample. Type: float"},
         'l0': "Soil sample (initial) length. Type: float or array of floats.",
         'wl0': ("Length of water on the inflow (above the soil) in cm. "
                 "Type: float"),
         'ww0': ("Weight of water on the inflow (above the soil) in gramms. "
                 "Only one of 'wl0' and 'ww0' can be specified."
                 "Type: float"),
	 'density_s': ("Density of the soil (without void space!!) in g/(cm^3)"),
    'consolidation': \
        {'con_type': ("Type of constitutive law to use for consolidation "
                        "sigma(e) and K(e). slurry=1, preconsolidated=2, freeform=3"),
         'con_max_refine': ("For con_type freeform with refinement, how many "
                         "times are we allowed to refine the grid?"),
         'a' : ("A parameter for the constitutive consolidation laws"),
         'b' : ("A parameter for the constitutive consolidation laws"),
         'c' : ("A parameter for the constitutive consolidation laws"),
         'd' : ("A parameter for the constitutive consolidation laws"),
         'cc' : ("A parameter for the constitutive consolidation laws"),
         'ei' : ("For freeform consolidation grid points in void radio e"),
         'si' : ("For freeform consolidation effective stress sigma values "
                 "with the grid points in void radio e"),
         'ki' : ("For freeform consolidation saturated conductivity values "
                 "with the grid points in void radio e"),
         'eiadd' : ("For freeform consolidation where to add the next "
                     "refinement point in grid points in void radio e"),
        },
    'initial-and-boundary-conditions': \
        {
         'rb_type': ("Right boundary condition of type:"
                     "\n\t 0 - no outflow [q_last = 0]"
                     "\n\t 1 - free outflow"
                     "\n\t 2 - prescribed pressure (see h_last)"
                     "\n\t 3 - sample is dipped in a basin with water level "
                     "         height dip_height"
                     "\nType: integer"),
         'h_last': ("The prescribed pressure on the right boundary - only for "
                    "'rb_type' = 2."),
         'dip_height': ("Maximum height of the water table in cm above "
                         "the base of the sample."),
         'fp1': ("porosity of filter 1, required if fl1>0"),
         'fp2': ("porosity of filter 1, required if fl2>0"),
     },
    'solver': \
        {'inner_points': ("How many (inner) point will be used for the sample "
                          "discretization. Type: integer"),
         'dtype': ("Discretization type. Determines how the soil sample "
                   "will be discretized in space. Currently implemented are: "
                   "\n\t 1 - linear discretization"
                   "\n\t 2 - the distance of two points is 'k_dx'*(dist of "
                   "previous two points.)"
                   "\nType: integer"),
         'k_dx': ("If dtype==2, k_dx value to use in the algorithm, see dtype"),
         'estimate_zp0': ("True or False. If True, the code attempts to help "
                         "the solver by determining good starting value of "
                         "the derivatives of the unknowns ei, dei/dt."),
         'L_atol': ("Default: 1e-8. Length sample is via an algebraic equation."
                     " With L_atol you can set abs tolerance for this. "),
         'numerfact_e0': ("Default: 0.999. Numerical starting value at 0 rpm"
                     " of the void ratio is numerfact_e0 * e0. This because"
                     " staring at e0 leads to problems on first step when value"
                     " of e becomes > e0, and no good penalization can be done." ),
         'e0_overshoot_factor': ("Default: 0.0. Error contribution of solver"
                     " searching solution e>e0. Contributes an error of "
                     "e0_overshoot_factor * (e-e0), forcing solution towards"
                     " e values below e0."),
        }
    }

INTERNAL_OPTIONS = {
  }
