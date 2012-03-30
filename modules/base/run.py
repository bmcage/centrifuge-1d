
# unsaturated_cfg = {
#            'fluid': {'s1_0': 0.1, 's2_0': 0.2 },
#   'discretization': {'inner_points': 80, 'first_point_offset': 80.0,
#                      'dtype': 1, 'percent_in_saturation': 40.0,
#                      'approximation_type': 5, 'mb_epsilon': 1e-5}

def adjust_cfg(cfg):
    """
      This method is called after the configuration is read from a file
      (and was validated). Allows to process configuration data supplied
       by configuration file(s), e.g. allocate the discretized interval
       based on the discretization type and number of inner points.
    """
    pass

def generate_tubes_suffixes(tubes_numbers):
    if not tubes_numbers:
        print('CFG: Error: No tubes numbers specified: "%s"' % tubes_numbers)
        exit(1)

    suffixes = ['-tube' + str(tube_no) for tube_no in tubes_numbers]
    identifiers = [', tube ' + str(tube_no) for tube_no in tubes_numbers]

    return suffixes, identifiers
