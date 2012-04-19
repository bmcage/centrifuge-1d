CONFIG_OPTIONS = {
    'set-up': \
        {'rb_type': ("(required) Right boundary type.\n"
                     "\t 0 - no outflow [q_last = 0]"
                     "\t 1 - free outflow"
                     "\t 2 - prescribed pressure (see h_last)"),
         'h_last': ("(dependent) The prescribed pressure on the right"
                    " boundary\n\t if 'rb_type' = 2")
        }
    }
