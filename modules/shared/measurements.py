from __future__ import division, print_function

import numpy as np
from collections import OrderedDict
from shared import get_directories, flatten
from os import makedirs, path
from modules.shared.functions import rpm2radps, compare_data, \
    smoothing_gaussian, smoothing_triangle, smoothing_linear


# MEASUREMENTS_NAMES are the mapping between the internal
# denotation of measured data (used during computation and for
# plotting) and corresponding options that are to be found in
# configuration inifile(s).
MEASUREMENTS_NAMES = {'MI': 'wl1', 'MO': 'wl_out', 'GC': 'gc1', 'RM': 'rm1',
                      'gF_MO': 'gf_mo', 'gF_MT': 'gf_mt',
                      'dgF_MO': None, 'dgF_MT': None,
                      'theta': 'theta'}

def determine_scaling_factor(v):
    return np.power(10, -np.floor(np.log10(np.max(np.abs(v)))))

def same_value(item):
    if np.isscalar(item): return item

    ref = item[0]
    for value in item:
        if not value == ref: return False

    return ref

def determine_measurements_scans(cfg):
    """
      Determine scans of times at which measurements were taken.
      This also determines when measurements should be computed.
    """
    from modules.shared.functions import phases_end_times

    phases_scans = \
        phases_end_times(cfg.get_value('duration', not_found=None),
                         cfg.get_value('deceleration_duration', not_found=None),
                         cfg.get_value('fh_duration', not_found=None),
                         cfg.get_value('include_acceleration', not_found=True))

    if phases_scans is None:
        print("\nPhases scand could not be determined. Were "
              "'duration'/'fh_duration'/'deceleration_duration' "
              "set properly?\n Cannot continue, exiting...")
        exit(1)

    if cfg.get_value('measurements_times'):
        scans = np.asarray(flatten(cfg.get_value('measurements_times')),
                           dtype=float)
        if not scans[0] == 0.0:
            scans = np.concatenate(([0.0], scans))
    else:
        scans = phases_scans

    cfg.del_value('measurements_times') # could have been set to "None"

    scan_span = float(cfg.get_value('scan_span', not_found=1.0))

    # adjust durations accordingly to scan_span
    if scan_span != 1.0:
        for name in ('duration', 'fh_duration', 'deceleration_duration'):
            duration = scan_span * np.asarray(cfg.get_value(name),
                                              dtype=float)
            if not np.alen(duration) > 0:
                duration = list(duration)
            cfg.set_value(name, duration)

    return (scans, scan_span)

def determine_weights(measurements, measurements_diff, cfg):
    measurements_weights = {}

    same_weights = True
    common_weight  = None

    for name in measurements.keys():
        iname_w = MEASUREMENTS_NAMES[name] + '_weights'
        if name in measurements_diff:
            iname_w = 'd' + iname_w

        value_w = cfg.get_value(iname_w, not_found=1.0)
        if not np.isscalar(value_w):
            value_w = flatten(value_w)
        cfg.del_value(iname_w)

        if same_weights: # are all weights the same float/integer?
            # (if yes, we don't have to preallocate array for each)
            w_ref = same_value(value_w)
            if w_ref is False:
                same_weights = False
            elif common_weight is None:
                common_weight = w_ref
            elif common_weight != w_ref:
                same_weights = False

        measurements_weights[name] = value_w

    if same_weights:
        measurements_weights = None
    else:
        for name in measurements_weights.keys():
            value_w = measurements_weights[name]
            value_m = measurements[name]

            m_length = np.alen(value_m)
            if name in measurements_diff:
                m_length -= 1

            if np.isscalar(value_w):
                value_w = value_w * np.ones((1, m_length), dtype=float)
            else:
                value_w = np.asarray(value_w, dtype=float)

                if not np.alen(value_w) == m_length:
                    if name in measurements_diff:
                        name = 'd' + name
                        value_m = value_m[1:] - value_m[:-1]
                    print('Length of measurements array and it''s weight '
                          'has to be the same:\nMeasurement name: ', name,
                          '\nMeasurements:\n  ', value_m,
                          '\nWeights:\n  ', value_w,
                          '\nLenghts: Measurements: ', m_length,
                          '  Weights: ', np.alen(value_w))
                    exit(1)

            measurements_weights[name] = value_w

    return measurements_weights

class Measurements():
    """
    This class is for handling with measurements of all types - measured data
    just as computed data. It reads needed informations from a configuration
    object and/or from saved data to files.
    For handling the (external) measurements it transforms the data to needed
    format and stores also additional information (weights of individual
    measurements i.e. the importance just as scaling them to avoid too wide
    range of the measured data that may lead to unintentional "preference" of
    larger-value measuremts).
    For the computed data it contains methods for computing needing measured
    value and storing the value, together with easy re-setting (meant for
    running an inverse problem to avoid unnecessary re-allocation).
    """
    def __init__(self):
        # Stored computed measurements
        self._computed = {} # values
        self._indexes  = {} # index of next measurement to be stored
        self._times = None  # times at which all the measurements were computed

        # User supplied (physically done) measurements
        self._measurements = OrderedDict() # values
        self._original_measurements = OrderedDict() # untransformed values
        self._measurements_weights = OrderedDict() # weights of each of measurement's' values
        self._measurements_xvalues = OrderedDict() # x-axes values
        self._original_measurements_xvalues = OrderedDict() # x-axes value of untransformed measurements
        self._measurements_times   = None # times at which measurements were taken
        self._measurements_diff    = []

        self._original_already_shifted = False # 'original' weren't adjusted WR_tara
        # Maximal number of values stored per measurement
        self._measurements_nr = -1 # i.e. length(self._times)

        # Radius * weight of tara (R is at the gravitational centre)
        self.WR_tara = {} # {'gF_MO': None, 'gF_MT': None}

        # Measurements scales
        self._scales_coefs = None
        self._scales = None

        # Measurements weights
        self._weights = None

        # Measurements arrays (stores all measurements in 1D array)
        #     The array self._{computations,error}_array is allocated
        #      when the size of self._measurements_array is determined
        self._measurements_array = None
        self._computations_array = None
        self._error_array        = None

        # Do we run a direct or inverse problem?
        self._run_inverse_p = None

    def read(self, cfg):
        """
        Read and transform data from configuration object.
        Data are the (external) measurements. Read and save also additional
        information - measurements times, weights and scaling.
        """


        measurements = self._measurements
        measurements_xvalues = self._measurements_xvalues
        measurements_diff    = self._measurements_diff

        # 0. Determine whether we run direct or inverse problem
        self.run_inverse_problem_p(cfg.get_value('inv_init_params'))

        # 1. determine measurements times
        [scans, scan_span] = determine_measurements_scans(cfg)

        m_first_idx = scans[1]  # scan[0] == 0.0
        m_last_idx  = scans[-1]
        t = scans * scan_span   # actual times of measurements
        t_meas = t[1:]
        scans_meas = scans[1:]

        # 2. a) determine measured data
        smoothing = cfg.get_value('smoothing', not_found={})
        if smoothing and (not type(smoothing) == dict):
            print('Smoothing has do be a dict where key is measurement'
                  'name and value is smoothing algorithm to be used.',
                  '\nCurrent value: ', smoothing,
                  '\n\nCannot continue, exiting...')
            exit(0)
        elif not smoothing:
            smoothing = {} # make sure it's dict

        for (name, iname) in MEASUREMENTS_NAMES.items():
            value = cfg.get_value(iname, not_found=None)
            if value == None: continue

            if type(value) in [int, float]: value = (value, )

            value = np.asarray(value, dtype=float)

            if name in smoothing:
                sm_opts = smoothing[name]

                if type(sm_opts) == str:
                    sm_opts = {'name': sm_opts}
                elif not type(sm_opts) is dict:
                    print('The value in smoothing has to be either string ',
                          'of used method, or dict. Current value: ', sm_opts,
                          '\nExiting...')
                    exit(0)

                sm_alg    = None
                sm_degree = 5

                for (sm_key, sm_value) in sm_opts.items():
                    if sm_key == 'name':
                        sm_alg = sm_value
                    elif sm_key == 'degree':
                        sm_degree = sm_value
                    else:
                        print('Unknown option for smoothing:', sm_key,
                              'Exiting...')
                        exit(0)

                if not sm_alg: continue

                # store untransformed value
                self._original_measurements[name] = value
                self._original_measurements_xvalues[name] = \
                  scan_span * np.arange(0, np.alen(value))

                if sm_alg == 'smlin': # linear smoothing
                    value = smoothing_linear(value, sm_degree)
                elif sm_alg == 'smtri': # triangular smoothing
                    value = smoothing_triangle(value, sm_degree)
                elif sm_alg == 'smgau': # gaussian smoothing
                    value = smoothing_gaussian(value, sm_degree)
                else:
                    print('Unknown smoothing value:', sm_alg,  'for key:', name)
                    exit(0)

            measurements[name] = value
            measurements_xvalues[name] = t_meas

            cfg.del_value(iname) # remove from cfg

        #    b) postprocessing MO: MO is a cumulative value
        if 'MO' in measurements:
            measurements['MO'] = np.cumsum(measurements['MO'])

        #    c) postprocessing gF_MO, gF_MT:
        #           - filter out values outside the measured times
        #           - if calibration curve present, use directly force as
        #              measurement, otherwise use difference between two
        #              subsequent force measurements
        g = cfg.get_value('g')
        if not scans_meas is None:
            filter_idxs = np.asarray(scans_meas, dtype=int)
        for F_name in ('gF_MT', 'gF_MO'):
            if F_name in measurements:

                # Determine the influence of tara
                gF_tara_calibration = cfg.get_value(F_name.lower() + '_tara')
                if not gF_tara_calibration is None:
                    (omega_rpm, gF_tara) = gF_tara_calibration
                    omega_radps = rpm2radps(omega_rpm)
                    # Centrifugal force = F = M omega^2 r,
                    # sensor gives us m kg, so F = m g, with m measurement
                    # M r = m g/ omega^2
                    WR_tara = gF_tara * g / omega_radps / omega_radps

                    if F_name == 'gF_MT':
                        # as extra we need to subtract the water inside the tube
                        extra_values = {}
                        for name in ('porosity', 'fl1', 'fl2', 'fp1', 'fp2'):
                            extra_values[name] = cfg.get_value(name)
                            if not np.isscalar(extra_values[name]):
                                print(name + ' has to be scalar to compute '
                                      'gF_MT. Aborting.')
                                exit(1)
                        for name in ('wl0', 'l0', 're'):
                            value = cfg.get_value(name)
                            if value is None:
                                extra_values[name] = 0.0
                            elif np.isscalar(value):
                                extra_values[name] = value
                            else:
                                extra_values[name] = value[0]

                        (rE, l0)   = (extra_values['re'], extra_values['l0'])
                        (fl1, fl2) = (extra_values['fl1'], extra_values['fl2'])
                        r0 = rE - fl2 - l0 - fl1
                        wl0 = extra_values['wl0']

                        liq_dens = cfg.get_value('density')
                        WR_fluid = (extra_values['porosity']
                                    * calc_sat_force(r0+fl1, rE-fl2, liq_dens))
                        if wl0 > 0.0:
                            WR_fluid += wl0 * (r0 - wl0/2.0)
                        if fl1 > 0.0:
                            WR_fluid += (extra_values['fp1']
                                         * calc_sat_force(r0, r0+fl1, liq_dens))
                        if fl2 > 0.0:
                            WR_fluid += (extra_values['fp2']
                                         * calc_sat_force(rE-fl2, rE, liq_dens))
                        tube_area = np.power(cfg.get_value('tube_diam')/2, 2) * np.pi
                        WR_fluid *= tube_area
                        # WR_fluid =  gF_fluid * g/(omega_radps^2)
                        # Hence:
                        # WR_tara = (gF_tara - gF_fluid) * g / (omega_radps^2)
                        WR_tara = WR_tara - WR_fluid

                    self.WR_tara[F_name] = WR_tara

                # Process the force measurements
                F  = measurements[F_name]

                # Leave only values at desired point (t_meas)
                F_filter = F[filter_idxs]

                if F_name in self._original_measurements:
                    orig_value  = self._original_measurements[F_name]
                    orig_xvalue = self._original_measurements_xvalues[F_name]

                    self._original_measurements[F_name] = \
                      orig_value[m_first_idx:m_last_idx+1]
                    self._original_measurements_xvalues[F_name] = \
                      orig_xvalue[m_first_idx:m_last_idx+1]

                calibration_curve = cfg.get_value(MEASUREMENTS_NAMES[F_name]
                                                  + '_calibration_curve')
                if calibration_curve is None:
                    measurements_diff.append(F_name)

                elif np.isscalar(calibration_curve):
                    F_filter -= calibration_curve
                    F = F_filter

                elif type(calibration_curve) is dict:
                    first_idx = 0
                    omegas = cfg.get_value('omega')

                    # phases_scans = [0.0, phase1_end, phase2_end, ...]
                    for (phase_end, omega) in zip(phases_scans[1:], omegas):
                        if not omega in calibration_curve:
                            print('Value of omega=' + str(omega) + ' not '
                                  'found in ' + F_name + '_calibration'
                                  '_curve. Cannot proceed, exiting.')
                            exit(1)

                        base_weight = calibration_curve[omega]
                        F[first_idx: phase_end+1] -= base_weight

                        first_idx = phase_end+1

                    F_filter = F[filter_idxs]
                else:
                    # type(calibration_curve) == list:
                    calibration_curve = np.asarray(calibration_curve,
                                                   dtype=float)
                    F_filter -= calibration_curve

                measurements[F_name] = F_filter

        #    d) theta uses pressure as x-value (and not time)
        if 'theta' in measurements:
            measurements_xvalues['theta'] = np.asarray(cfg.get_value('p'),
                                                       dtype=float)

        # 3. determine weights of measurements
        self._measurements_weights = determine_weights(measurements,
                                                       measurements_diff, cfg)

        # 4. set remaining data to internal variables
        self._measurements_times   = t_meas
        self._times                = t
        self._measurements_nr      = np.alen(t)
        self._weights              = None # weights as numpy array

        # 5. scaling measurements
        self._scales_coefs = cfg.get_value('measurements_scale_coefs')
        if not self._scales_coefs: self._scales_coefs = {}

        cfg.del_value('measurements_scale_coefs')

        # 6. convert omega from rpm to radians/s
        #    NOTE: we do it here because *_calibration_curve is expressed
        #          in terms of omega(rpm)
        for key in ['omega']:
            value = cfg.get_value(key, not_found=None)
            if value is None:
                continue
            elif type(value) == list:
                cfg.set_value(key, [rpm2radps(omega) for omega in value])
            else:
                cfg.set_value(key, rpm2radps(value))

    def dump(self):
        """ Collect stored data to simple format for saving. """
        return (self._measurements, self._measurements_weights,
                self._measurements_xvalues, self._measurements_times,
                self._original_measurements, self._original_measurements_xvalues,
                self._times, self._measurements_nr, self._computed,
                self._indexes, self._scales_coefs)

    def load(self, raw_data):
        """ Restore data in simple format from saved file. """
        (self._measurements, self._measurements_weights,
         self._measurements_xvalues, self._measurements_times,
         self._original_measurements, self._original_measurements_xvalues,
         self._times, self._measurements_nr, self._computed,
         self._indexes, self._scales_coefs) = raw_data

    def get_times(self):
        """ Return time of measurements. """
        return self._times

    def get_names(self):
        """ Return names of (external) measurements that are stored. """
        return list(self._measurements.keys())

    def get_values(self):
        """ Return values of (external) measurements that are stored. """
        return list(self._measurements.values())

    def _get_scales(self):
        """
        Return scaling coeficients of (external) measurements that are stored.
        """

        if self._scales is None:
            scales_coefs = self._scales_coefs

            scales = []
            for (name, measurement) in self._measurements.items():
                if not name in scales_coefs:
                    #scales_coefs[name] = determine_scaling_factor(measurement)
                    scales_coefs[name] = 1.0

                length = np.alen(measurement)

                if name in self._measurements_diff:
                    length -= 1

                scales.append(scales_coefs[name]
                              * np.ones((1, length), dtype=float))
            self._scales = np.concatenate(scales, axis=1)[0]

        return self._scales

    def _get_weights(self):
        """ Return weights of (external) measurements that are stored. """
        if self._weights is None:
            weights = []

            if not self._measurements_weights is None:
                for weight in self._measurements_weights.values():
                    weights.append(weight)

            if not weights:
                self._weights = 1.0
            else:
                self._weights = np.concatenate(weights)

        return self._weights

    def _get_measurements(self):
        if self._measurements_array is None:
            values = []
            for (name, value) in self._measurements.items():
                if name in ('gF_MO', 'gF_MT'):
                    # For the forces we need to determine the influence of tara
                    # and subtract it from the measured values of force; but as
                    # we don't have the value from measurement, we subtract the
                    # force of tara found during computation - but this
                    # information is only available at runtime
                    if (name + '_tara' in self._computed):
                        F_tara = self._computed[name + '_tara']

                        self._measurements[name][:] -= F_tara[1:]

                if name in self._measurements_diff:
                    values.append(value[1:] - value[:-1])
                else:
                    values.append(value)

            self._measurements_array = np.concatenate(values)
            # preallocate arrays
            self._computations_array = \
              np.empty(self._measurements_array.shape, dtype=float)
            self._error_array = \
              np.empty(self._measurements_array.shape, dtype=float)

        return self._measurements_array

    def _get_computations(self):
        computations = self._computations_array
        data = self._computed

        iS = 0
        for name in self._measurements.keys():
            if name == 'theta':
                # theta is in x and not t
                value = data[name]
            else:
                value = data[name][1:]

            if name in self._measurements_diff:
                iE = iS + np.alen(value) - 1
                computations[iS:iE] = value[1:] - value[:-1]
            else:
                iE = iS + np.alen(value)
                computations[iS:iE] = value
            iS = iE

        return computations

    def get_xvalues(self):
        """ Return x-axis values of (external) measurements that are stored. """
        return list(self._measurements_xvalues.values())

    def store_calc_measurement(self, meas_id, value):
        """ Store (calculated) value of measurement. """
        if not meas_id in self._computed:
            self._computed[meas_id] = np.empty([self._measurements_nr, ],
                                               dtype=float)
            self._indexes[meas_id]  = 0

        idx = self._indexes[meas_id]
        self._computed[meas_id][idx] = value
        self._indexes[meas_id] += 1

        return value

    def reset_calc_measurements(self):
        """ Reset all calculated values of measurements. """
        indexes = self._indexes
        for calc_id in indexes.keys():
            indexes[calc_id] = 0

    def iterate_calc_measurements(self):
        """ Iterate over stored values of calculated measurements. """
        t = self.get_times()
        measurements_diff = self._measurements_diff

        for (name, yvalue) in self._computed.items():
            if name in ('x'): continue

            if name in ('h', 'u'):
                xvalue = self._computed['x']
            elif name == 'theta':
                xvalue = self._measurements_xvalues[name]
            else:
                xvalue = t

            yield (name, xvalue, yvalue)

            if name in measurements_diff:
                yield ('d' + name, xvalue[1:], yvalue[1:] - yvalue[:-1])

    def iterate_meas_measurements(self, untransformed=False, model=None):
        """
          Iterate over measured values. If untransformed==True, original
          (not in any way transformed) measurements are returned.
        """
        measurements_diff = self._measurements_diff

        if untransformed:
            measurements = self._original_measurements
            xvalues      = self._original_measurements_xvalues

            if self._original_already_shifted:
                pass # nothing to do, measurements were already adjusted
            elif not model:
                print("Warning: 'original' line will not be shifted as 'model' "
                      "variable was not set.")
            else:
                for name in ('gF_MO', 'gF_MT'):
                    if name in measurements:
                        assert name in self.WR_tara
                        omega2g = determine_omega2g(model, xvalues[name])
                        measurements[name] -= omega2g * self.WR_tara[name]
                self._original_already_shifted = True
        else:
            measurements = self._measurements
            xvalues      = self._measurements_xvalues

        for (name, yvalue) in measurements.items():
            xvalue = xvalues[name]

            if name in measurements_diff:
                if not untransformed:
                    yield ('d' + name, xvalue[1:], yvalue[1:] - yvalue[:-1])
            else:
                yield (name, xvalue, yvalue)

    def calc_measurement_p(self, measurement_name):
        if measurement_name == 'RM':
            return False

        return ((not self._run_inverse_p)
                or (measurement_name in self._measurements))

    def run_inverse_problem_p(self, flag):
        self._run_inverse_p = bool(flag)

    def store_calc_theta(self, h, SC, theta_s, theta_r, rho, g):
        if not 'theta' in self._computed:
            self._computed['theta'] = np.empty(h.shape, dtype=float)
            self._indexes['theta'] = 0

        idx = self._indexes['theta']
        size = np.alen(h)
        value = SC.retention_curve(theta_s, rho, g, theta_r, p=None, h=h,
                                   find_p=False)[1]

        self._computed['theta'][idx:idx+size] = value
        self._indexes['theta'] += size

        return value

    def store_calc_u(self, x, h, SC, method_name='h2u'):
        """
        Calculate and store relative saturation u.

        Arguments:
        x  - x-coordinates
        h  - pressure head
        SC - saturation curve object
        method_name - accessor function for the computation of the saturation
        """
        if not 'u' in self._computed:
            # preallocate and initialize index
            for meas_id in ('x', 'h', 'u'):
                self._computed[meas_id] = \
                  np.empty([self._measurements_nr, np.alen(h)], dtype=float)
                self._indexes[meas_id]  = 0

        # store x
        self._computed['x'][self._indexes['x'], :] = x
        self._indexes['x'] += 1

        # store h
        self._computed['h'][self._indexes['h'], :] = h
        self._indexes['h'] += 1

        # store u
        value = self._computed['u'][self._indexes['u'], :] # reference

        h2u = getattr(SC, method_name)
        h2u(h, value)
        self._indexes['u'] += 1

        return value

    def store_calc_wm(self, u, dy, s1, s2, mass_in, saturated_length,
                      free_fluid_length, porosity, fl2, fp2, tube_area=1.0,
                      store=True):
        """
          Calculate and store the value of water contained in the experiment.
          The amount is the sum of water in inflow chamber (mass_in), outflow
          chamber + any other source like water in basin (free_fluid_length)
          and saturated zone (saturated_length) just as unsaturated zone
          (which is between s1 and s2). The soil is characterized by it's
          porosity (porosity).
          Also water in filter of length fl2 with porosity fp2 is considered,
          but only for fully saturated filter.

          Arguments:
          u - relative saturation
          dy - difference of y: dy = diff(y) = y(i+1) - y(i)
          s1, s2 - interfaces (s1 < s2); between them is the unsaturated part
          mass_in - water in the inflow chamber
          saturated_length - length of saturated part inside the sample
          free_fluid_length - length of mass outside the sample (measured
                     to the 1D length relatively to the tube diameter)
          soil_porosity - porosity of the soil sample
          fl2, fp2 - ending filter length, porosity and distance from
                     sample beginning (to filter's beginning)
          tube_area - the cross-section area of the tube containing the sample
          from_end - (optional) if specified, computed GC will be returned as
                     distance from the "from_end" point
          tube_area - the cross-section area of the tube containing the sample
          store - if set to False, value will not be saved
        """
        # Water mass
        ds = s2 - s1
        unsat = ds/2  * (dy[0]* u[0] + dy[-1]*u[-1]
                         + np.sum((dy[:-1] + dy[1:])*u[1:-1]))

        WM_in_tube = (mass_in + porosity*(saturated_length + unsat))
        WM_in_tube += fp2 * fl2
        WM_in_tube *= tube_area

        WM_total   = WM_in_tube + free_fluid_length * tube_area

        if store:
            self.store_calc_measurement('WM', WM_total)
            self.store_calc_measurement('WM_in_tube', WM_in_tube)

        return (WM_total, WM_in_tube)

    def store_calc_gc(self, u, y, dy, s1, s2, mass_in, d_sat_s, d_sat_e,
                      soil_porosity, fl2, fp2, fr2, WM_in_tube, fluid_density,
                      tube_area=1.0, from_end=None):
        """
          Calculate and store the value of the gravitational center of water
          contained in the sample - the water on the inflow is taken into
          account (but not the water in the outflow chamber).
          Computed GC is measured from the BEGINNING of the soil sample
          (in the direction from centrifuge axis)

          Arguments:
          u - relative saturation
          y - transformed interval where holds: u(x) = u(fl1 + s1 + (s2-s1)y)
          dy - difference of y: dy = diff(y) = y(i+1) - y(i)
          s1, s2 - interfaces (s1 < s2); between them is the unsaturated part
          mass_in - water in the inflow chamber
          d_sat_s - distance of the beginning of saturated part from r0
          d_sat_e - distance of the end of saturated part from r0
          soil_porosity - porosity of the soil sample
          fl2, fp2, fr2 - ending filter length, porosity and distance from
                          sample beginning (to filter's beginning)
          WM_in_tube - amount of water contained inside the tube
          fluid_density - density of used fluid (in g/cm3)
          tube_area - the cross-section area of the tube containing the sample
          from_end - (optional) if specified, computed GC will be returned as
                     distance from the "from_end" point
        """
        r0 = 0.0

        gc =  (calc_force(u, y, dy, r0, s1, s2, mass_in, d_sat_s, d_sat_e,
                          soil_porosity, fl2, fp2, fr2, fluid_density)
                / WM_in_tube * tube_area)

        if not from_end is None: gc = from_end - gc

        return self.store_calc_measurement('GC', gc)

    def store_calc_gf_mo(self, omega2g, mo_1d, MO_calibration_curve,
                         fluid_density, tube_area):
        """
          Calculate and store the value of the calculated centrifugal force
          caused by the expelled water in water chamber. More preciselly it's
          the weight that would be measured in 1g environment, i.e. W=F/g with
          F being the centrifugal force and g the normal Gravitational
          acceleration.

          Arguments:
          omega2g - value of omega^2/g, where omega is the rotational speed and
                    g is the gravitation acceleration ~981 cm/s^2
          mo      - amount of expelled water in gramms
          MO_calibration_curve - the curve converting the the amount of expelled
                    water to the value of the gravitational center of the
                    expelled water.
          tube_area - the cross-sectional area of the tube containing the sample
                    (should not be contained in the mo?)
          fluid_density - density of used fluid (in g/cm3)

          Returns the weight W.
        """
        calibration = MO_calibration_curve

        mo = mo_1d * tube_area * fluid_density
        (mo0, gc0) = calibration[0]

        if mo < mo0:
            print('Amount of expelled fluid (in gramms): ', mo,
                  ' is less than the amount specified as the first point '
                  'of the MO_calibration curve: ', calibration[0][0],
                  '. Cannot proceed, exiting...')
            exit(1)

        GC = -1.0

        for (mo1, gc1) in calibration[1:]:
            if mo < mo1:
                GC = gc0 + (mo - mo0)/(mo1 - mo0) * (gc1 - gc0)
                break
            else:
                mo0 = mo1
                gc0 = gc1

        if GC < 0.0:
            print('Amount of expelled water: ', mo,  ' is more than the amount '
                  'at the last point on the calibration curve: ',
                  calibration[-1][0], '. Cannot proceed, exiting...')
            exit(1)

        F = omega2g * GC * mo

        if 'gF_MO' in self.WR_tara:
            self.store_calc_measurement('gF_MO_tara',
                                        omega2g * self.WR_tara['gF_MO'])

        return self.store_calc_measurement('gF_MO', F)

    def store_calc_gcf_mo(self, omega2g, rS, rE, fluid_density, chamber_area):
        """
          Computes the centrifugal force acting on the expelled mass based on
          the assumption of mass being uniformly distributed

          Arguments:
          fluid_density - density of used fluid (in g/cm3)
          chamber_area - the cross-section area of the outflow chamber
        """
        raise NotImplementedError('Broken. For calculation the amount of'
                                  'expelled water needs to be recalculated'
                                  '- chamber area vs. tube area.')
        F = omega2g * calc_sat_force(rS, rE, fluid_density) * chamber_area
        return self.store_calc_measurement('gCF_MO', F)

    def store_calc_gf_mt(self, omega2g, u, y, dy, r0, s1, s2, mass_in,
                         d_sat_s, d_sat_e, soil_porosity, fl2, fp2, fr2, l0,
                         fluid_density, tube_area):
        """
          Calculate and store the value of the centrifugal force of water in
          the sample. The water on the inflow is taken into account, but not
          the water on the outflow. Computed centrifugal force is measured
          from the BEGINNING of the soil sample (in the direction from
          centrifuge axis)

          Arguments:
          omega2g - value of (omega^2/g) at given time
          u - relative saturation
          y - transformed interval where holds: u(x) = u(fl1 + s1 + (s2-s1)y)
          dy - difference of y: dy = diff(y) = y(i+1) - y(i)
          r0 - distance from the centrifuge axis to soil sample beginning
          s1, s2 - interfaces (s1 < s2); between them is the unsaturated part
          mass_in - water in the inflow chamber
          d_sat_s - distance of the beginning of saturated part from r0
          d_sat_e - distance of the end of saturated part from r0
          soil_porosity - porosity of the soil sample
          fl2, fp2, fr2 - ending filter length, porosity and distance from
                          sample beginning (to filter's beginning)
          l0 - length of soil sample
          fluid_density - density of used fluid (in g/cm3)
          tube_area - the cross-sectional area of the tube
        """
        F = (calc_force(u, y, dy, r0, s1, s2, mass_in, d_sat_s, d_sat_e,
                        soil_porosity, fl2, fp2, fr2, fluid_density)
             * omega2g * tube_area)

        if 'gF_MT' in self.WR_tara:
            self.store_calc_measurement('gF_MT_tara',
                                        omega2g * self.WR_tara['gF_MT'])

        return self.store_calc_measurement('gF_MT', F)

    def store_calc_rm(t, u, mass_in, mass_out, s1, s2, model):

        raise NotImplementedError('Calculation of rotational momentum is not '
                                  'verified and therefore not provided.')

        porosity = model.porosity
        y  = model.y
        dy = model.dy
        L  = model.l0
        l0_out = L + model.wt_out
        l_out  = L + model.wt_out - mass_out

        ds = s2 - s1

        P = np.pi * model.d / 4
        omega2g = find_omega2g(t, model.omega, model)

        # Rotational momentum
        r0 = model.r0
        rm_unsat = (porosity * 1/4 * model.density * ds
                    * (np.power(r0 + s1, 2)*dy[0]*u[0]
                       + np.power(r0 + s2, 2)*dy[-1]*u[-1]
                       + np.sum((dy[:-1]+dy[1:]) * u[1:-1]
                                * np.power(r0 + s1 + ds*y[1:-1], 2))))
        rm_sat   = (1/6 * model.density
                    * (porosity * (np.power(r0 + s1, 3) - np.power(r0, 3))
                       + (np.power(r0, 3) - np.power(r0 - mass_in, 3))
                       + (np.power(r0 + l0_out, 3) - np.power(r0 + l_out, 3))))

        RM = omega2g * P * (rm_unsat + rm_sat)

        return RM

    def get_error(self):
        """ Returns the difference between computed and measured values """

        scale   = self._get_scales()
        weights = self._get_weights()

        measurements = self._get_measurements()
        computations = self._get_computations()
        #error        = self._error_array
        # inverse solver needs to have freshly allocated
        error         = np.empty(measurements.shape, dtype=float)

        error[:] = scale * (computations - measurements)
        if not np.isscalar(weights):
            error[:] *= weights

        return error

    def display_error(self, stream=None):
        self._get_measurements() # ensure all runtime data is set

        status_items = []
        measurements_diff = self._measurements_diff

        computed = self._computed
        measured = self._measurements

        for (name, measured_value) in self._measurements.items():
            if name == 'theta':
                # theta is based on pressures in x, not t
                computed_value = computed[name]
            else:
                computed_value = computed[name][1:]

            if name in measurements_diff:
                compare_data('d'+name, computed_value[1:] - computed_value[:-1],
                             measured_value[1:] - measured_value[:-1],
                             stream=stream)
            else:
                compare_data(name, computed_value, measured_value, stream=stream)

    def get_penalized(self, penalization, scalar=False):
        """
          Return penalization for given measurement.
          If 'scalar'=True, return it as a single scalar value.
        """
        scale   = self._get_scales()
        weights = self._get_weights()
        if self._measurements_array is None:
            print('ERROR: Cannot determine penalization if measurements'
                  'are not specified. Exiting.')
            exit(1)
        elif scalar:
            return penalization * np.alen(self._measurements_array)
        else:
            error =  (penalization + self._measurements_array)
            error[:] *= scale
            if not np.isscalar(weights):
                error[:] *= weights
            return error

##################################################################
#                     Auxiliary functions                        #
##################################################################

def calc_unsat_force(r0, u, s1, s2, y, dy, soil_porosity, fluid_density):
    """
    Calculate gram-force of the water in the unsaturated part of the sample
    Multply with g for actual Newton!
    """
    ds = s2 - s1
    if ds == 0.:
        #support calling by saturated experiment
        return 0.

    r = r0 + s1 + ds*y
    F_unsat = (soil_porosity * fluid_density * ds/2
               * (dy[0]*u[0]*r[0] + dy[-1]*u[-1]*r[-1]
                  + np.sum((dy[:-1] + dy[1:])*u[1:-1]*r[1:-1])))

    return F_unsat

def calc_sat_force(rS, rE, fluid_density):
    """
    Calculate gram-force of the water in the saturated part of the sample
    Multply with g for actual Newton!
    """
    return 1/2 * fluid_density * (np.power(rE, 2) - np.power(rS, 2))

def calc_force(u, y, dy, r0, s1, s2, mass_in, d_sat_s, d_sat_e,
               soil_porosity, fl2, fp2, fr2, fluid_density):
    """
    Calculate gram-force of the water in the sample
    Multply with g for actual Newton!
    """

    F_unsat = calc_unsat_force(r0, u, s1, s2, y, dy, soil_porosity,
                               fluid_density)

    F_sat = (soil_porosity * calc_sat_force(r0+d_sat_s, r0+d_sat_e,
                                            fluid_density)
                + calc_sat_force(r0-mass_in, r0, fluid_density))
    if fl2 > 0.0:
        F_sat += fp2 * calc_sat_force(r0+fr2, r0+fr2+fl2, fluid_density)

    return F_unsat + F_sat

def determine_omega2g(model, measurements_times):
    from scikits.odes.sundials.ida import IDA_RhsFunction
    from modules.shared.solver import simulate_direct
    from modules.shared.measurements import Measurements

    class _omega2g_residual(IDA_RhsFunction):
        def evaluate(self, t, z, zdot, result, model):
            result[:] = zdot[:]
            return 0

    residual_fn = _omega2g_residual()

    def _omega2g_on_measurement(t, z, model, measurements):
        measurements.store_calc_measurement('omega2g', model.find_omega2g(t))
    def _omega2g_initialize_z0(z0, model):
        z0[:] = -1.0

    measurements = Measurements()
    measurements._times = measurements_times
    measurements._measurements_nr = np.alen(measurements_times)

    z_size_backup = model.z_size
    model.z_size = 1

    (flag, t, z, i) = simulate_direct(_omega2g_initialize_z0, model,
                                      measurements, residual_fn,
                                      root_fn = None, nr_rootfns=None,
                                      on_measurement=_omega2g_on_measurement)

    model.z_size = z_size_backup

    if not flag:
        print('Computation of omega2g was not successfull.\nDetails:',
              't=', t,'\nz=', z, '\nindex=', i, '\Exiting.')
        exit(1)

    return measurements._computed['omega2g']
