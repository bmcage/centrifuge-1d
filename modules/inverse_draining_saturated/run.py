from modules.direct_draining_saturated.run import solve as solve_direct, \
     display_graphs, multiple_solves, get_refencing_models
from modules.shared.functions import measurements_time
from modules.shared.solver import simulate_inverse
from modules.shared.vangenuchten import h2u
from modules.shared.show import disp_status as display_status
from numpy import alen

def update_runtime_variables(model):
    if model.dynamic_h_init:
            model.h_init = min(model.c_gammah / model.gamma, model.h_init_max)
            if model.verbosity > 1:
                print('\nh_init: ', model.h_init,
                      'u_init', h2u(model.h_init, model.n,
                                    model.m, model.gamma), '\n')

def solve(model):

    no_measurements = []

    def ip_direct_drainage(model):
        update_runtime_variables(model)
        (flag, t, z, gc1, rm1, u, wm, wm_in_tube) = solve_direct(model)

        contains_data = (alen(t) > 1)

        # we discard values at t=0 (for given measurement)
        if model.calc_wl_out and contains_data:
            wl_out = z[1:, model.mass_out_idx].transpose()
        else:
            wl_out = no_measurements

        if model.calc_wl_in and contains_data:
            wl_in = z[1:, model.mass_in_idx].transpose()
        else:
            wl_in = no_measurements

        if model.calc_gc and contains_data:
            gc1 = gc1[1:]
        else:
            gc1 = no_measurements

        if model.calc_rm and contains_data:
            rm1 = rm1[1:]
        else:
            rm1 = no_measurements

        return (flag, t, wl_in, wl_out, gc1, rm1)

    t_meas = measurements_time(model)

    (inv_params, cov) = \
      simulate_inverse(t_meas, ip_direct_drainage, model, model.inv_init_params,
                       wl_in_meas  = model.wl1,
                       wl_out_meas = model.wl_out,
                       gc_meas     = model.gc1,
                       rm_meas     = model.rm1,
                       wl_in_weights  = model.wl1_weights,
                       wl_out_weights = model.wl_out_weights,
                       gc_weights     = model.gc1_weights,
                       rm_weights     = model.rm1_weights,
                       optimfn=model.optimfn)

    return (inv_params, cov)

def make_status_plots(results, annotation, model):
    status_items = []
    for (data_id, data_c) in zip(annotation, results):

        if not data_id in ['MI', 'MO', 'GC', 'RM']: continue

        if data_id == 'MI': model_id = 'wl1'
        elif data_id == 'MO': model_id = 'wl_out'
        else: model_id = data_id

        if not model_id in model: continue

        data_m = getattr(model, model_id)
        if data_m: status_items.append(mk_status_item(data_id, data_c, data_m))

    return status_items

def run(model):
    (inv_params, cov) = solve(model)

    # DISPLAY RESULTS:
    if inv_params:
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0
        referencing_models = get_refencing_models(model)
        (results, annotation) = multiple_solves(model, referencing_models)

        dg_options = {'save_figures': model.save_figures,
                      'separate_figures': model.separate_figures,
                      'save_as_text': model.save_as_text,
                      'show_figures': model.show_figures,
                      'experiment_info': model.experiment_information}
        display_graphs(model, results, annotation, dg_options)
        display_status(data_plots=make_status_plots(results, annotation, model),
                       params=inv_params, cov=cov)
        model.verbosity = model_verbosity # restore verbosity

    return inv_params
