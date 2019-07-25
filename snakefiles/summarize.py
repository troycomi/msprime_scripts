import pandas as pd
import click
from collections import OrderedDict


@click.command()
@click.option('-t', '--observed', type=click.File('r'))
@click.option('-s', '--simulated', type=click.File('r'))
@click.option('-p', '--pi', type=click.File('r'))
@click.option('-a', '--admixed', type=click.File('r'))
@click.option('-o', '--output', type=click.File('w'))
def main(observed, simulated, pi, admixed, output):
    outputs = OrderedDict()

    if observed and simulated:
        # handle desert distribution
        observed = pd.read_csv(observed, sep='\t',
                               usecols=[0, 3],
                               index_col=0,
                               names=['winsize', 'ob_prop'])
        simulated = pd.read_csv(simulated, sep='\t',
                                usecols=[1, 3],
                                index_col=0,  # this is relative to usecols
                                names=['winsize', 'sim_prop'])
        joined = simulated.join(observed, how='left')
        outputs['desert-MSE'] = ((joined['sim_prop'] -
                                  joined['ob_prop'])**2).mean()

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
