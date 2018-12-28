from argparse import ArgumentParser
import numpy as np
import itertools


def main():
    '''given a list of arguments and values, generates all
    permutations of the parameter list to standard output.
    Example: Parameter_Sweeper.py -p "a;0:1:2" -p "b;0,2" ->
    -a 0 -b 0
    -a 1 -b 0
    -a 2 -b 0
    -a 0 -b 2
    -a 1 -b 2
    -a 2 -b 2
    '''

    args, formats = get_arguments()
    args = get_permutations(args, formats)
    for arg in args:
        print('{}'.format(arg))


def get_arguments(args=None):
    parser = ArgumentParser(description="Get list of arguments")
    parser.add_argument("-p", "--parameter",
                        dest="params",
                        action='append',
                        nargs='+',
                        help="provide parameter name and list of values "
                        "as PARAM;range,range;FORMAT. FORMAT is any valid "
                        "format notation, e.g ':.2f' or ':.2e'. Default is "
                        "no added formatting, taken for the first argument")
    options = parser.parse_args(args)

    if options.params is None:
        parser.print_help()
        raise ValueError("No arguments provided!")

    result = {}
    formats = {}
    for param in options.params:
        param = param[0]  # weird argument parser thing
        # if fmt is not provided it defaults to empty string
        key, args, *fmt = param.split(';')
        fmt = fmt[0] if fmt else ''

        values = []
        for arg in args.split(','):
            rangeVals = arg.split(':')
            if len(rangeVals) == 1:
                try:
                    values.append(float(arg))
                    continue
                except ValueError:
                    raise ValueError('Ill formed number in {}'.format(param))

            if len(rangeVals) > 3:
                raise SyntaxError("Too many range values for {}".format(param))

            try:
                rangeVals = [float(rv) for rv in rangeVals]
            except ValueError:
                raise ValueError('Ill formed number in {}'.format(param))

            if len(rangeVals) == 2:
                # use default 1 spacing
                rangeVals.insert(1, 1)

            start = rangeVals[0]
            stop = rangeVals[2]
            step = rangeVals[1]

            numsteps = int(np.ceil((stop - start) / step)) + 1

            values += np.linspace(start, stop, numsteps,
                                  endpoint=True).tolist()

        if key not in result:
            formats[key] = fmt
            result[key] = []

        result[key] += values

    return result, formats


def get_permutations(values, formats):
    result = []

    for key, args in values.items():
        if len(key) > 1:
            flag = '--' + key
        else:
            flag = '-' + key

        line_args = [('{} {' + formats[key] + '}')
                     .format(flag, arg) for arg in args]

        if len(result) == 0:
            result = line_args
        else:
            result = ['{} {}'.format(prod[0], prod[1]) for
                      prod in itertools.product(result, line_args)]

    return result


if __name__ == "__main__":
    main()
