from __future__ import division

import numpy as np
from shared import get_directories, flatten
from os import makedirs, path
from modules.shared.vangenuchten import h2u
from modules.shared.functions import rpm2radps

# MEASUREMENTS_NAMES are the mapping between the internal
# denotation of measured data (used during computation and for
# plotting) and corresponding options that are to be found in
# configuration inifile(s).
MEASUREMENTS_NAMES = {'MI': 'wl1', 'MO': 'wl_out', 'GC': 'gc1', 'RM': 'rm1',
                      'F_MO': 'f_mo', 'F_MT': 'f_mt',
                      'dF_MO': None, 'dF_MT': None,
                      'theta': 'theta'}

def calc_f_tara(omega2g, WR_tara):
    if WR_tara is None:
        raise NotImplementedError

    return omega2g * WR_tara

def determine_scaling_factor(v):
    return np.power(10, -np.floor(np.log10(np.max(np.abs(v)))))

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
        self._measurements_weights = OrderedDict() # weights of each of measurement's' values
        self._measurements_xvalues = OrderedDict() # x-axes values
        self._measurements_times   = None # times at which measurements were taken

        # Maximal number of values stored per measurement
        self._measurements_nr = -1 # i.e. length(self._times)

        # Radius * weight of tara (R is at the gravitational centre)
        self.WR_mt_tara = None
        self.WR_mo_tara = None

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

    def read(self, cfg):
        """
        Read and transform data from configuration object.
        Data are the (external) measurements. Read and save also additional
        information - measurements times, weights and scaling.
        """

        from modules.shared.functions import phases_end_times
        from collections import OrderedDict

        measurements = self._measurements
        measurements_xvalues = self._measurements_xvalues
        measurements_weights = self._measurements_weights

        # 1. determine measurements times
        phases_scans = \
          phases_end_times(cfg.get_value('duration', not_found=None),
                           cfg.get_value('deceleration_duration', not_found=None),
                           cfg.get_value('fh_duration', not_found=None),
                           cfg.get_value('include_acceleration', not_found=True))

        if cfg.get_value('measurements_times'):
            scans = np.asarray(flatten(cfg.get_value('measurements_times')),
                           dtype=float)
            if not scans[0] == 0.0:
                scans = np.concatenate(([0.0], scans))
        else:
            scans = phases_end_times

        scan_span = float(cfg.get_value('scan_span', not_found=None))
        if scan_span != 1.0:
            for name in ('duration', 'fh_duration', 'deceleration_duration'):
                duration = scan_span * np.asarray(cfg.get_value(name),
                                                  dtype=float)
                if not np.alen(duration) > 0:
                    duration = list(duration)
                cfg.set_value(name, duration)

        if scans is None:
            t = None
            t_meas = None
        else:
            t = scans * scan_span
            t_meas = t[1:]
            scans_meas = scans[1:]

        cfg.del_value('measurements_times')

        # 2. a) determine measured data
        same_weights = True
        common_weight  = None

        for (name, iname) in MEASUREMENTS_NAMES.items():
            value = cfg.get_value(iname, not_found=None)
            if value == None: continue

            if type(value) in [int, float]: value = (value, )

            value = np.asarray(value, dtype=float)

            measurements[name] = value
            measurements_xvalues[name] = t_meas

            cfg.del_value(iname) # remove from cfg

        #    b) postprocessing MO: MO is a cumulative value
        if 'MO' in measurements:
            measurements['MO'] = np.cumsum(measurements['MO'])

        #    c) postprocessing F_MO, F_MT:
        #           - filter out values outside the measured times
        #           - if calibration curve present, use directly force as
        #              measurement, otherwise use difference between two
        #              subsequent force measurements
        filter = np.asarray(scans_meas, dtype=int)
        for F_name in ('F_MT', 'F_MO'):
            if F_name in measurements:

                # Determine the influence of tara
                F_tara = cfg.get_value(F_name.lower() + '_tara')
                if not F_tara is None:
                    (omega_rpm, weight_tara) = F_tara
                    g = cfg.get_value('g')
                    omega_radps = rpm2radps(omega_rpm)
                    WR_tara = weight_tara * g / omega_radps / omega_radps
                    setattr(self, 'WR' + F_name[1:].lower() + '_tara', WR_tara)

                # Process the force measurements
                F  = measurements[F_name]

                # Leave only values at desired point (t_meas)
                F_filter = F[filter]

                calibration_curve = cfg.get_value(MEASUREMENTS_NAMES[F_name]
                                                  + '_calibration_curve')
                if calibration_curve is None:
                    del measurements[F_name]
                    del measurements_xvalues[F_name]

                    dF = F_filter[1:] - F_filter[:-1] # forces difference

                    measurements['d'+F_name] = dF
                    measurements_xvalues[F_name] = t_meas[1:]

                    continue

                if np.isscalar(calibration_curve):
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

                    F_filter = F[filter]
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

        # 3. a) determine weights of measurements
        for (name, value) in measurements.items():
            if name in ['dF_MT', 'dF_MO']:
                # we use the same weights for forces differences (except first)
                iname_w = name[1:].lower() + '_weights'
            else:
                iname_w = MEASUREMENTS_NAMES[name] + '_weights'
            value_w = cfg.get_value(iname_w, not_found=None)
            if value_w == None:
                value_w = 1.0
            else:
                cfg.del_value(iname_w)

            if same_weights: # are all weights the same float/integer?
                # (if yes, we don't have to preallocate array for each)
                if not type(value_w) in [float, int]:
                    same_weights = False
                elif common_weight is None:
                    common_weight = value_w
                elif common_weight != value_w:
                    same_weights = False

            if type(value_w) in [float, int]:
                value_w = value_w * np.ones(value.shape, dtype=float)
            else:
                value_w = np.asarray(value_w, dtype=float)

            if not np.alen(value) == np.alen(value_w):
                print('Length of measurements array and it''s weight has to be '
                      'the same:\nMeasurement name: ', name,
                      '\nMeasurements:\n  ', value,
                      '\nWeights:\n  ', value_w,
                      '\nLenghts: Measurements: ', np.alen(value),
                      '  Weights: ', np.alen(value_w))
                exit(1)

            if name in ['dF_MT', 'dF_MO']:
                value_w = value_w[1:] # cut off the first value
            measurements_weights[name] = value_w

        #    b) weights postpocessing - if all weights are the same, drop them
        if same_weights: measurements_weights = OrderedDict()

        # 4. set remaining data to internal variables
        self._measurements_times   = t_meas
        self._times                = t
        self._measurements_nr      = np.alen(t)
        self._weights              = None # weights as numpy array

        # 5. scaling measurements
        self._scales_coefs = cfg.get_value('measurements_scale_coefs',
                                           not_found={})
        cfg.del_value('measurements_scale_coefs')

        # 6. convert omega from rpm to radians/s
        #    NOTE: we do it here because *_calibration_curve is expressed
        #          in terms of omega(rpm)
        for key in ['omega']:
            value = cfg.get_value(key)
            if type(value) == list:
                cfg.set_value(key, [rpm2radps(omega) for omega in value])
            else:
                cfg.set_value(key, rpm2radps(value))

    def dump(self):
        """ Collect stored data to simple format for saving. """
        return (self._measurements, self._measurements_weights,
                self._measurements_xvalues, self._measurements_times,
                self._times, self._measurements_nr, self._computed,
                self._indexes, self._scales_coefs)

    def load(self, raw_data):
        """ Restore data in simple format from saved file. """
        (self._measurements, self._measurements_weights,
         self._measurements_xvalues, self._measurements_times,
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
            scale_coef = self._scales_coefs

            scales = []

            for (name, measurement) in self._measurements.items():
                if not name in measurements_scales:
                    scale_coef[name] = determine_scaling_factor(measurement)

                scales.append(scale_coef[name]
                              * np.ones(measurement.shape, dtype=float))

            self._scales = np.concatenate(scales)

        return self._scales

    def _get_weights(self):
        """ Return weights of (external) measurements that are stored. """
        if self._weights is None:
            weights = self._measurements_weights.values()
            if weights:
                self._weights = np.concatenate(weights)
            else:
                self._weights = np.asarray([], dtype=float)

        return self._weights

    def _get_measurements(self):
        if self._measurements_array is None:
            for name in ('F_MO', 'F_MT'):
                # For the forces we need to determine the influence of tara and
                # subtract it from the measured values of force; but as we don't
                # have the value from measurement, we subtract the force of tara
                # found during computation - but this information is only
                # available at runtime
                if (name + '_tara' in self._computed):
                    F_tara = self._computed[name + '_tara']
                    self._measurements[name][:] -= F_tara
            self._measurements_array = \
              np.concatenate(self._measurements.values())
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
            value = data[name]
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

    def get_calc_measurement(self, meas_id, not_found=None):
        """
            Return stored times and value of calculated measurements with
            measurement ID 'meas_id'.
        """
        if meas_id in self._computed:
            times  = self._times
            values = self._computed[meas_id]

        elif meas_id in ('dF_MT', 'dF_MO'):
            mid = meas_id[1:]
            if mid in self._computed:
                F = self._computed[mid]

                times  = self._times[1:]
                values = F[1:] - F[:-1]
            else:
                return not_found
        else:
            return not_found

        return (times, values)

    def iterate_calc_measurements(self):
        """ Iterate over stored values of calculated measurements. """
        for meas_id in self._computed.keys():
            (times, values) = self.get_calc_measurement(meas_id)
            yield (meas_id, times, values)

        for meas_id in ('dF_MT', 'dF_MO'):
            calc_data = self.get_calc_measurement(meas_id)
            if not calc_data is None:
                yield (meas_id, calc_data[0], calc_data[1])

    def store_calc_u(self, h, n, m, gamma):
        """
        Calculate and store relative saturation u.

        Arguments:
        h - hydraulic head
        n, m, gamma - corresponding van Genuchten parameters
        """
        meas_id = 'u'
        if not meas_id in self._computed:
            self._computed[meas_id] = \
              np.empty([self._measurements_nr, np.alen(h)], dtype=float)
            self._indexes[meas_id]  = 0

        value = self._computed[meas_id][self._indexes[meas_id], :] # reference

        h2u(h, n, m, gamma, value)
        self._indexes[meas_id] += 1

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
          fluid_density - density of used fluid (in g/cm3)
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

    def store_calc_f_mo(self, omega2g, mo_1d, MO_calibration_curve, tube_area):
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

          Returns the weight W.
        """
        calibration = MO_calibration_curve

        mo = mo_1d * tube_area
        (mo0, gc0) = calibration[0]

        if mo < mo0:
            print('Amount of water: ', mo, ' is less than the amount at the '
                  'first point on the calibration curve: ', calibration[0][0],
                  '. Cannot proceed, exiting...')
            exit(1)

        GC = -1.0

        for (mo1, gc1) in calibration[1:]:
            if mo < mo1:
                GC = gc0 + (mo - mo0)/(mo1 - mo0) * (gc1 - gc0)

        if GC < 0.0:
            print('Amount of expelled water: ', mo,  ' is more than the amount '
                  'at the last point on the calibration curve: ',
                  calibration[-1][0], '. Cannot proceed, exiting...')
            exit(1)

        F = omega2g * GC * mo

        if not self.WR_mo_tara is None:
            self.store_calc_measurement('F_MO_tara',
                                        calc_f_tara(omega2g, self.WR_mo_tara))

        return self.store_calc_measurement('F_MO', F)

    def store_calc_cf_mo(self, omega2g, rS, rE, fluid_density, chamber_area):
        """
          Computes the centrifugal force acting on the expelled mass based on
          the assumption of mass being uniformly distributed

          Arguments:
          chamber_area - the cross-section area of the outflow chamber
        """
        raise NotImplementedError('Broken. For calculation the amount of'
                                  'expelled water needs to be recalculated'
                                  '- chamber area vs. tube area.')
        F = omega2g * calc_sat_force(rS, rE, fluid_density) * chamber_area
        return self.store_calc_measurement('CF_MO', F)

    def store_calc_f_mt(self, omega2g, u, y, dy, r0, s1, s2, mass_in,
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
          tube_area - the cross-sectional area of the tube
        """
        F = (calc_force(u, y, dy, r0, s1, s2, mass_in, dsat_s, d_sat_e,
                        soil_porosity, fl2, fp2, fr2, fluid_density)
             * omega2g * tube_area)

        if not self.WR_mt_tara is None:
            self.store_calc_measurement('F_MT_tara',
                                        calc_f_tara(omega2g, self.WR_mt_tara))

        return self.store_calc_measurement('F_MT', F)

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
        error        = self._error_array

        error[:] = scale * (computation - measurements)

        if weights:
            error[:] *= weights

        return error

##################################################################
#                     Auxiliary functions                        #
##################################################################

def calc_unsat_force(r0, u, s1, s2, y, dy, soil_porosity, fluid_density):
    ds = s2 - s1
    r = r0 + s1 + ds*y
    F_unsat = (soil_porosity * fluid_density * ds/2
               * (dy[0]*u[0]*r[0] + dy[-1]*u[-1]*r[-1]
                  + np.sum((dy[:-1] + dy[1:])*u[1:-1]*r[1:-1])))

    return F_unsat

def calc_sat_force(rS, rE, fluid_density):
    return 1/2 * fluid_density * (np.power(rE, 2) - np.power(rS, 2))

def calc_force(u, y, dy, r0, s1, s2, mass_in, d_sat_s, d_sat_e,
               soil_porosity, fl2, fp2, fr2, fluid_density):

    F_unsat = calc_unsat_force(r0, u, s1, s2, y, dy, soil_porosity,
                               fluid_density)

    F_sat = (soil_porosity * calc_sat_force(r0+d_sat_s, r0+d_sat_e,
                                            fluid_density)
                + calc_sat_force(r0-mass_in, r0, fluid_density))
    if fl2 > 0.0:
        F_sat += fp2 * calc_sat_force(r0+fr2, r0+fr2+fl2, fluid_density)

    return F_unsat + F_sat
