from auxiliaryfunctions import save_data
from direct_saturated import direct_saturated_problem, extract_saturated_characteristics

def generate_data(filename_out, model):
    _flag, t, z = direct_saturated_problem(model)
    GC, RM = extract_saturated_characteristics(t, z, model)
    save_data(filename_out, {'t': t, 'GC': GC, 'RM': RM})

if __name__ == "__main__":
    from direct_saturated import run_direct

    model, GC, RM = run_direct(False)

    if not model.inverse_data_filename:
        raise ValueError('Output data file ''inverse_data_filename'' for generated values not specified !')
    save_data(model.inverse_data_filename, {'t': model.tspan, 'GC': GC, 'RM': RM})
