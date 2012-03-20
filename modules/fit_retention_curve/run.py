import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import modules.base.run as base
from config import merge_cfgs
from shared_functions import h2u

PARAMETERS = {'inverse': ['inv_init_params']}

CFG_ADDITIONAL_PARAMETERS = {}

def base_cfg():
    return CFG_ADDITIONAL_PARAMETERS

def adjust_cfg(flattened_cfg):

    required_parameters = ['inv_init_params', 'p', 'theta']

    for param in required_parameters:
        if not param  in flattened_cfg:
            print('CFG:check: Missing paramter in configuration file(s): %s' 
                  % param)
            exit(1)

    theta_s_p = 'theta_s' in flattened_cfg
	theta_r_p = 'theta_r' in flattened_cfg
	inv_init_params_len = len(flattened_cfg['inv_init_params'])

    if ((theta_s_p and theta_r_p and inv_init_params == 2)
		or (theta_s_p and (not theta_r_p) and inv_init_params == 3)
		or ((not theta_s_p) and theta_r_p and inv_init_params == 3)
		or ((not theta_s) and (not theta_r) and inv_init_params == 4)):

        # check correctness of input:
        # if 2 initial guesses given, optimize only n, gamma
		# if 3 given, we optimize also either theta_s or theta_r
		# if 4 we optimize also theta_s and theta_r
		# theta_s and theta_r;
		pass
	else:
		th_s_str = ''
		th_r_str = ''
		if theta_s_p: th_s_str = ' = ' + str(flattened_cfg['theta_s'])
		if theta_r_p: th_s_str = ' = ' + str(flattened_cfg['theta_r'])

        print('Inconsistent initial guesses inv_init_params = %s'
			  % flattened_cfg['inv_init_params'])
		print("with 'theta_s'%s and 'theta_r'%s" % th_s_str, th_r_str)
		exit(1)

	flattened_cfg['theta_s_p'] = theta_s_p

def experiments_files(first_experiment, last_experiment, tubes):
    files = []
    identifiers = []

    for exp_no in range(first_experiment, last_experiment+1):
            inifilename = 'experiment_' + str(exp_no) +'.ini'
            identifier  = 'experiment ' + str(exp_no)

            files.append(inifilename)
            identifiers.append(identifier)

    return (identifiers, files)


def solve(model):
    def lsq_fn(xdata, *optim_args):

        optim_args = len(optim_args)

        if optim_args == 4
			h = xdata
			(n, gamma, theta_s, theta_r) = optim_args
		elif optim_args == 3:
			(h, theta_s_p, theta_sr) = xdata
			if theta_s_p:
				theta_s = theta_sr
				(n, gamma, theta_r) = optim_args
			else:
				theta_r = theta_sr
				(n, gamma, theta_s) = optim_args
		else:
			(h, theta_s, theta_r) = xdata
			(n, gamma) = optim_args

        m = 1. - 1./n

        u = h2u(h, n, m, gamma)

        theta = theta_r + (theta_s - theta_r) * u

        return theta

    inv_params_init = model.inv_init_params
	inv_init_params_len = len(model.inv_init_params)

	if inv_init_params_len == 4:
		xdata = model.h
	elif inv_init_params_len == 3:
		theta_s_p = has_attr('theta_s', model)
		if theta_s_p: theta_sr = model.theta_s
		else: theta_sr = model.theta_r

		xdata = (model.h, theta_s_p, theta_sr)
	 else:
		 theta_s = model.theta_s
		 theta_r = model.theta_r

		 xdata = (model.h, theta_s, theta_r)

    data_measured   = model.theta

    inv_params, cov_inv = curve_fit(lsq_fn, xdata,
                                   data_measured, p0 = inv_params_init)

    if inv_init_params_len == 4:
		(n, gamma, theta_s, theta_r) = inv_params
		print(' n     found: %s\n gamma found: %s' % tuple(inv_params))
	elif inv_init_params_len == 3:
		if theta_s_p:
			(n, gamma, theta_r) = inv_params
			print(' n:\t %s\n gamma:\t %s\n theta_r:\t' % tuple(inv_params))
		else:
			(n, gamma, theta_s) = inv_params
			print(' n:\t %s\n gamma:\t %s\n theta_s:\t' % tuple(inv_params))
	else:
		print(' n:\t %s\n gamma:\t %s\n theta_s:\t \n theta_r:\t'
			  % tuple(inv_params))
    print(' Cov: %s' % cov_inv)

    u_inv = lsq_fn(xdata, *inv_params)

    draw_graphs(1, model.h, model.theta, n, gamma, theta_s, theta_r)

    return inv_params

def draw_graphs(fignum, h_measured, u_measured, inv_found):
    import matplotlib.pyplot as plt

    plt.figure(fignum, figsize=(8, 4.5))

    [n_found, gamma_found] = inv_found
    h_calc = -np.arange(0, 1600, 1)
    u_calc = h2u(h_calc, n_found, 1.-1./n_found, gamma_found)

    plt.plot(u_calc, -h_calc, '-', u_measured, -h_measured, 'x',)
    plt.yscale('log')
    plt.ylabel('Hydraulic pressure ''h'' [$Pa$]')
    plt.xlabel('Relative saturation ''u'' ')

    plt.show()
