import pandas as pd
import click
from collections import OrderedDict


@click.command()
@click.option('-s', '--simulated', type=click.File('r'))
@click.option('-p', '--pi', type=click.File('r'))
@click.option('-a', '--admixed', type=click.File('r'))
@click.option('-o', '--output', type=click.File('w'))
def main(simulated, pi, admixed, output):
    outputs = OrderedDict()

    if simulated:
        # handle desert distribution
        simulated = pd.read_csv(simulated, sep='\t',
                                usecols=[1, 3],
                                names=['winsize', 'sim_prop'])
        for row in simulated.itertuples(index=False):
            outputs[f'desert-{row.winsize}'] = row.sim_prop

    if pi:
        # add in pi results
        keys = pi.readline().split()
        values = pi.readline().split()

        for k, v in zip(keys, values):
            outputs[f'pi-{k}'] = v

    if admixed:
        for line in admixed:
            tokens = line.split()
            # remove trailing ':' on key
            outputs[f'admix-{tokens[0][:-1]}'] = tokens[1]

    items = outputs.items()
    output.write('\t'.join([item[0] for item in items]) + '\n')
    output.write('\t'.join([str(item[1]) for item in items]) + '\n')


if __name__ == '__main__':
    main()
