from __future__ import division

import numpy as np
from modules.shared.functions import determine_scaling_factor
from shared import get_directories, flatten
from os import makedirs, path
from modules.shared.vangenuchten import h2u

# MEASUREMENTS_NAMES are the mapping between the internal
# denotation of measured data (used during computation and for
# plotting) and corresponding options that are to be found in
# configuration inifile(s).
MEASUREMENTS_NAMES = {'MI': 'wl1', 'MO': 'wl_out', 'GC': 'gc1', 'RM': 'rm1',
                      'F_MO': 'f_mo', 'F_MT': 'f_mt', 'theta': 'theta'}

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
        self._measurements = {}         # values
        self._measurements_weights = {} # weights of each of measurement's' values
        self._measurements_xvalues = {} # x-axes values
        self._measurements_times   = None # times at which measurements were taken

        # Maximal number of values stored per measurement
        self._measurements_nr = -1 # i.e. length(self._times)

    def read(self, cfg):
        """
        Read and transform data from configuration object.
        Data are the (external) measurements. Read and save also additional
        information - measurements times, weights and scaling.
        """

        from modules.shared.functions import phases_end_times
        from collections import OrderedDict

        measurements = OrderedDict()
        measurements_xvalues = OrderedDict()
        measurements_weights = OrderedDict()

        # 1. determine measurements times
        if cfg.get_value('measurements_times'):
            t = np.asarray(flatten(cfg.get_value('measurements_times')),
                           dtype=float)
            if not t[0] == 0.0:
                t= np.concatenate(([0.0], t))
        else:
            t = phases_end_times(cfg.get_value('duration', not_found=None),
                                 cfg.get_value('deceleration_duration',
                                               not_found=None),
                                 cfg.get_value('fh_duration', not_found=None),
                                 cfg.get_value('include_acceleration',
                                               not_found=True))

        if t is None:
            t_meas = None
        else:
            t_meas = t[1:]
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
        #           - correct data usind calibration curve
        #           - filter out values outside the measured times
        F_filter = None
        for F_name in ('F_MT', 'F_MO'):
            if F_name in measurements:
                calibration_curve = \
                  np.asarray(cfg.get_value(MEASUREMENTS_NAMES[F_name]
                                           + '_calibration_curve'),
                             dtype=float)

                if F_filter is None: F_filter = np.asarray(t_meas, dtype=int)
                measurements[F_name] -= calibration_curve

                measurements[F_name]  = measurements[F_name][F_filter]

        #    d) theta uses pressure as x-value (and not time)
        if 'theta' in measurements:
            measurements_xvalues['theta'] = np.asarray(cfg.get_value('p'),
                                                       dtype=float)

        # 3. a) determine weights of measurements
        for (name, value) in measurements.items():
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

            measurements_weights[name] = value_w

        #    b) weights postpocessing - if all weights are the same, drop them
        if same_weights: measurements_weights = OrderedDict()

        # 4. set all data to internal variables
        self._measurements = measurements
        self._measurements_weights = measurements_weights
        self._measurements_xvalues = measurements_xvalues
        self._measurements_times   = t_meas
        self._times                = t
        self._measurements_nr      = np.alen(t)

    def dump(self):
        """ Collect stored data to simple format for saving. """
        return (self._measurements, self._measurements_weights,
                self._measurements_xvalues, self._measurements_times,
                self._times, self._measurements_nr, self._computed,
                self._indexes)

    def load(self, raw_data):
        """ Restore data in simple format from saved file. """
        (self._measurements, self._measurements_weights,
         self._measurements_xvalues, self._measurements_times,
         self._times, self._measurements_nr, self._computed,
         self._indexes) = raw_data

    def get_times(self):
        """ Return time of measurements. """
        return self._times

    def get_names(self):
        """ Return names of (external) measurements that are stored. """
        return list(self._measurements.keys())

    def get_values(self, scaling=True):
        """ Return values of (external) measurements that are stored. """
        return list(self._measurements.values())

    def get_scales(self, scaling_coefs={}):
        """
        Return scaling coeficients of (external) measurements that are stored.
        """
        scales = []

        for (name, measurement) in self._measurements.items():

            if name in scaling_coefs:
                scaling_coef = scaling_coefs[name]
            else:
                scaling_coef = determine_scaling_factor(measurement)

            data_scale = scaling_coef * np.ones(measurement.shape, dtype=float)

            scales.append(data_scale)

        return scales

    def get_weights(self):
        """ Return weights of (external) measurements that are stored. """
        return list(self._measurements_weights.values())

    def get_xvalues(self):
        """ Return x-axis values of (external) measurements that are stored. """
        return list(self._measurements_xvalues.values())

    def store_calc_measurement(self, meas_id, value):
        """ Store (calculated) value of measurement. """
        if not meas_id in self._computed:
            self._computed[meas_id] = np.empty([self._measurements_nr, ],
                                          dtype=float)
            self._indexes[meas_id]  = 0

        self._computed[meas_id][self._indexes[meas_id]] = value
        self._indexes[meas_id] += 1

        return value

    def reset_calc_measurements(self):
        """ Reset all calculated values of measurements. """
        indexes = self._indexes
        for calc_id in indexes.keys():
            indexes[calc_id] = 0

    def get_calc_measurement(self, meas_id, not_found=None):
        """ Return stored value of calculated measurements with ID meas_id. """
        if meas_id in self._computed: return self._computed[meas_id]
        else: return not_found

    def iterate_calc_measurements(self):
        """ Iterate over stored values of calculated measurements. """
        for (meas_id, value) in self._computed.items():
            yield (meas_id, value)

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
            print('Amount of water: ', mo, ' is less than first point on the '
                  'calibration curve: ', calibration[0],'. Cannot proceed, '
                  'exiting...')
            exit(1)

        GC = -1.0

        for (mo1, gc1) in calibration[1:]:
            if mo < mo1:
                GC = gc0 + (mo - mo0)/(mo1 - mo0) * (gc1 - gc0)

        if GC < 0.0:
            print('Amount of expelled water: ', mo,  ' is more than the last '
                  'point on the calibration curve: ', calibration[0],'. Cannot '
                  'proceed, exiting...')
            exit(1)

        F = omega2g * GC * mo

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
