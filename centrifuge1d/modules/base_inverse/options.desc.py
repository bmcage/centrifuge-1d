# This file contains the description of base options

CONFIG_OPTIONS = {
    'inverse-solver' : \
      {'optimfn': ("Optimization function for inverse problem used."
                   "Values: 'leastsq', 'fmin', 'fmin_powell', 'fmin_cg', "
                   "'fmin_bfgs', 'fmin_slsqp'")},
    'inverse-solver-options': \
      {'epsfcn': ("Type: float."
                "Used by: 'leastsq', 'fmin_slsqp'"),
       'factor': ("Type: float."
                  "Used by: 'leastsq'"),
       'xtol': ("Type: float."
                "Used by: 'leastsq', 'fmin', 'fmin_powell', 'fmin_slsqp'"),
       'ftol': ("Type: float."
                "\nUsed by: 'leastsq', 'fmin', 'fmin_powell'"),
       'gtol': ("Type: float."
                "\nUsed by: 'fmin_cg', 'fmin_bfgs'"),
       'max_fev': ("Maximal number of direct function evaluations."
                   "\nUsed by: 'fmin', 'fmin_powell"),
       'max_inv_iter': ("Maximal number of iterations of the inverse solver."
                        "\nUsed by: 'fmin', 'fmin_powell, 'fmin_cg', "
                        "'fmin_bfgs', 'fmin_slsqp'"),
       'disp_inv_conv': ("Type: boolean."
                         "\nUsed by: 'fmin', 'fmin_powell, 'fmin_cg', "
                         "'fmin_bfgs', 'fmin_slsqp'")
      },
    'additional': \
      {'transform_params':("Specifies whether to transform optimized "
                           "parameters as during inverse problem computation "
                           "unsupported values may appear."
                           "Type: boolean. Default: True"),
       'untransformed_cov': ("If optimized parameters are transformed, "
                             "computed 'cov' is corresponds to these "
                             "transfomed parameters. If 'untransformd_cov' "
                             "is set to True, the covariance corresponding "
                             "to untransformed parameters is computed.")
    }
}

INTERNAL_OPTIONS = {}

PROVIDE_OPTIONS = {
    'inv_init_params': ("Optimized parameters are set by the inverse solver, "
                        "so none of them should be specified in the inifile.")}
