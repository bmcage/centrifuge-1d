# direct_saturated:
def characteristics(t, mass_in, mass_out, model):
    porosity = model.porosity
    r0 = model.r0
    L  = model.l
    l0_out = L + model.l0_out
    l_out  = L + model.l0_out - mass_out

    GC = np.empty(t.shape, float)
    RM = np.empty(t.shape, float)

    P = np.pi * model.d / 4

    WM = P * (mass_in + porosity * L + mass_out)

    for i in range(len(t)):
        omega2g = find_omega2g(t, model)
        
        # Gravitational center
        gc_sat   = (1/2 * model.density
                    * (porosity * (np.power(r0 + L, 2) - np.power(r0, 2))
                       + (np.power(r0, 2) - np.power(r0 - mass_in[i], 2))
                       + (np.power(r0 + l0_out, 2) - np.power(r0 + l_out[i], 2))))
        GC[i] =  P * gc_sat / WM[i]

        # Rotational momentum
        rm_sat   = (1/6 * model.density
                    * (porosity * (np.power(r0 + L, 3) - np.power(r0, 3))
                        + (np.power(r0, 3) - np.power(r0 - mass_in[i], 3))
                        + (np.power(r0 + l0_out, 3) - np.power(r0 + l_out[i], 3))))

        RM[i] = omega2g * P * rm_sat

    return GC, RM, WM


def total_water_volume(model):
    #TODO: we assume that we start with no outspelled water
    return (model.l0_in + model.l * model.porosity)

def extract_saturated_characteristics(t, z, model):
    GC, RM = characteristics(model.tspan, z[:, mass_in_idx],
                             z[:, mass_out_idx], model)[:2]

    return GC, RM

#=====================================================================
# inverse_saturated:

def ip_direct_saturated_characteristics(model, Ks):
    model.ks = Ks
    if model.debugging:
        print(Ks)
    _flag, t, z = solve_direct_saturated_problem(model)
    GC, RM = extract_saturated_characteristics(t, z, model)
    return np.concatenate((GC, RM))
