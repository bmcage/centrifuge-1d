# This file contains the description of base options

CONFIG_OPTIONS = {
    'inverse-solver' : \
      {'optimfn': ("Optimization function for inverse problem used."
                   "Values: 'leastsq', 'fmin', 'fmin_powell', 'fmin_cg', "
                   "'fmin_bfgs'")},
    'inverse-solver-options': \
      {'epsn': ("Type: float."
                "Used by: 'leastsq'"),
       'factor': ("Type: float."
                  "Used by: 'leastsq'"),
       'xtol': ("Type: float."
                "Used by: 'leastsq', 'fmin', 'fmin_powell'"),
       'ftol': ("Type: float."
                "\nUsed by: 'leastsq', 'fmin', 'fmin_powell'"),
       'gtol': ("Type: float."
                "\nUsed by: 'fmin_cg', 'fmin_bfgs'"),
       'max_fev': ("Maximal number of direct function evaluations."
                   "\nUsed by: 'fmin', 'fmin_powell"),
       'max_inv_iter': ("Maximal number of iterations of the inverse solver."
                        "\nUsed by: 'fmin', 'fmin_powell, 'fmin_cg', "
                        "'fmin_bfgs'"),
       'disp_inv_conv': ("Type: boolean."
                         "\nUsed by: 'fmin', 'fmin_powell, 'fmin_cg', "
                         "'fmin_bfgs'")
      }
}

INTERNAL_OPTIONS = {}

PROVIDE_OPTIONS {
    'inv_init_params': ("Optimized parameters are set by the inverse solver, "
                        "so none of them should be specified in the inifile.")}
