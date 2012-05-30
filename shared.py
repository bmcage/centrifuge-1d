import numpy as np

def generate_tubes_suffixes(tubes_numbers):
    if not tubes_numbers:
        print('generate_tubes_suffixes error: No tubes (numbers) were '
              'specified - input values was: "%s"' % tubes_numbers)
        exit(1)

    suffixes = ['-tube' + str(tube_no) for tube_no in tubes_numbers]
    identifiers = [', tube ' + str(tube_no) for tube_no in tubes_numbers]

    return suffixes, identifiers


def print_by_tube(tube_number, tube_data):
    print('Tube number: ', tube_number)
    for tdata in tube_data:
        print(tdata[0] / 100.)

def make_collector(tubes_numbers):
    by_tube_collection = {tube_no:[] for tube_no in tubes_numbers}
    fifo_collection = []

    def collection(command, data = None, tube_no = None, print_format_fn=print):
        if command == 'collect':
            if tube_no is None:
                fifo_collection.append(data)
            else:
                by_tube_collection[tube_no].append(data)

        elif command in ['print', 'print-by-tube', 'print-fifo']:
            if command == 'print-by-tube':
                for (tube_no, tube_data) in by_tube_collection.items():
                    print('Tube:', tube_no)
                    for value in tube_data:
                        print_format_fn(value)
            else:
                for value in fifo_collection:
                    print_format_fn(value)

        elif command == 'get':
            return (fifo_collection, by_tube_collection)
        else:
            raise ValueError('Collector: Unknown command: %s.\n Valid commans'
                             ' are: "collect", "print", "print-by-tube",'
                             ' "print-fifo" or "get".' % command)
    return collection
