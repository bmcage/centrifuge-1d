import numpy as np

def generate_tubes_suffixes(tubes_numbers):
    if not tubes_numbers:
        print('generate_tubes_suffixes error: No tubes (numbers) were '
              'specified - input values was: "%s"' % tubes_numbers)
        exit(1)

    numbers = tubes_numbers.split(',')

    suffixes = ['-tube' + str(tube_no) for tube_no in numbers]
    identifiers = [', tube ' + str(tube_no) for tube_no in numbers]

    return suffixes, identifiers


def print_by_tube(tube_number, tube_data):
    print('Tube number: ', tube_number)
    for tdata in tube_data:
        print(tdata[0] / 100.)

def make_collector(tubes_numbers):
    by_tube_collection = {tube_no:[] for tube_no in tubes_numbers}
    fifo_collection = []

    def collection(command, data = None, tube_no = None):
        if command == 'collect':
            if tube_no:
                by_tube_collection[tube_no].append(data)
            else:
                raise ValueError('Collector:collect - tube number or data'
                                 ' not specified')
        elif command in ['print', 'print-by-tube', 'print-fifo']:
            if data:
                if command == 'print-by-tube':
                    for (tube_no, tube_data) in by_tube_collection.items():
                        data(tube_no, tube_data)
                else:
                    data(fifo_collection)
            else:
                raise ValueError('Collector:print data not specified. Data'
                                 ' has to contain a print-data function.')
        elif command == 'get':
            return (fifo_collection, by_tube_collection)
        else:
            raise ValueError('Collector: Unknown command: %s.\n Valid commans'
                             ' are: "collect", "print", "print-by-tube",'
                             ' "print-fifo" or "get".' % command)
    return collection
