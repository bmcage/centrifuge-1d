
def test_related_options(cfg, options_names, mode='exact1'):
    """
      Test the presence of nonempty options valeus in configuration 'cfg'
      specified by 'options_names' with respect to 'mode'.
      Modes:
        exact1  - exactly one option must contain data
        atmost1 - 0 or 1 options must contain data
    """
    found = False
    for name in options_names:
        value = cfg.get_value(name)

        if value:
            if not found:
                found = True
            else:
                print('Only one of given options can be specified: ',
                      ', '.join(options_names))
                return False

    if (not found) and (mode == 'exact1'):
        print('Exactly one of following options has to be specified :',
              ', '.join(options_names))
        return False

    return True
