import elfi
from elfi.clients import multiprocessing
import click
import numpy as np
import os
import shutil
import subprocess


base_output = '/tigress/tcomi/abwolf_abc/'
summary_output = '/tigress/tcomi/abwolf_abc/SplitPop/summary.txt'
@click.command()
@click.option('-b', '--base', default="",
              help='Root directory to store temporary output')
@click.option('-o', '--output', default="",
              help='Output file containing parameters and summary statistics')
@click.option('-n', '--number-processes', default=5,
              help='Number of simulations to run in this instance')
def main(base, output, number_processes):
    # generate null (should only execute once)
    subprocess.call(["snakemake", "null", "--profile", "elfi_profile"])
    if base:
        global base_output
        base_output = base
    if output:
        global summary_output
        summary_output = output
    # setup priors
    elfi.set_client(multiprocessing.Client(num_processes=number_processes))
    model = elfi.ElfiModel(name='msprime')
    elfi.Prior('uniform', 0, 0.5, model=model, name='n1')
    elfi.Prior('uniform', 0, 0.5, model=model, name='n2')
    elfi.Prior('uniform', 0, 1, model=model, name='split_prop')
    elfi.Prior('uniform', 1, 1861, model=model, name='split_size')
    snake_sim = elfi.tools.external_operation(command,
                                              prepare_inputs=prepare_inputs,
                                              process_result=process_result,
                                              stdout=False)
    vec_snake = elfi.tools.vectorize(snake_sim)
    empirical = np.array([0,  # desert MSE
                          0.00047, 0.00035, 0.00035,  # pi
                          0.1, 0.1, 0.1,  # fst
                          0.03, 0.037])  # admix
    snake_node = elfi.Simulator(vec_snake, model['n1'],
                                model['n2'],
                                model['split_prop'],
                                model['split_size'],
                                observed=empirical,
                                name='snake_sim')
    snake_node.uses_meta = True
    distance = elfi.Distance('euclidean', snake_node, name='dist')
    pool = elfi.ArrayPool(['n1', 'n2', 'split_prop',
                           'split_size', 'snake_sim'],
                          name='bolfi_pool')

    bolfi = elfi.BOLFI(distance, batch_size=1,
                       initial_evidence=6, update_interval=3,
                       bounds={'n1': (0, 0.3),
                               'n2': (0, 0.5),
                               'split_prop': (0, 1),
                               'split_size': (1, 1861)},
                       pool=pool)
    post = bolfi.fit(n_evidence=10, bar=False)
    pool.save()
    post.save()

    return

    print(snake_sim(0.1, 0.1, 0.1, 100,
                    seed=2, meta={
                        'model_name': 'test',
                        'batch_index': 1,
                        'submission_index': 2}))
    return

    rej = elfi.Rejection(distance, model['n1'], batch_size=1)
    rej.sample(5)
    return


command = 'snakemake --profile elfi_profile --configfile {config_file}'


def prepare_inputs(*inputs, **kwinputs):
    n1, n2, split_prop, split_size = inputs
    split_size = int(split_size)
    meta = kwinputs['meta']
    seed = kwinputs['seed']

    admix_base = '{model_name}_{batch_index}_{submission_index}'
    if 'index_in_batch' in meta:
        admix_base += '_{index_in_batch}'
    admix_base = admix_base.format(**meta)

    yml_file = f'{admix_base}.yaml'

    with open(yml_file, 'w') as writer:
        writer.write(f'''---
paths:
    base_output: "{base_output}"
    admixed_dir: "__BASE_OUTPUT__/{admix_base}"
    summary: "__ADMIXED_DIR__/summary.txt"

msprime:
    base_seed: {seed}
    n1: {n1}
    n2: {n2}
    split_population_proportion: {split_prop}
    split_population_size: {split_size}
''')

    kwinputs['config_file'] = yml_file
    kwinputs['output_root'] = f'{base_output}SplitPop/{admix_base}'
    kwinputs['output_file'] = f'{base_output}SplitPop/{admix_base}/summary.txt'

    return inputs, kwinputs


def process_result(completed_process, *inputs, **kwinputs):
    output_file = kwinputs['output_file']

    results = np.loadtxt(output_file, skiprows=1)

    if not os.path.exists(summary_output):
        # write header
        with open(summary_output, 'w') as writer:
            writer.write(
                '\t'.join(('seed n1 n2 split_prop split_size '
                           'desert-mse pi-AF pi-EU pi-AS '
                           'pi-AF-EU pi-AF-AS pi-EU-AS '
                           'admix-ASN admix-EUR').split()))
            writer.write('\n')

    with open(summary_output, 'a') as writer:
        writer.write(f'{kwinputs["seed"]}\t')
        writer.write(f'{inputs[0]}\t')
        writer.write(f'{inputs[1]}\t')
        writer.write(f'{inputs[2]}\t')
        writer.write(f'{inputs[3]}\t')
        writer.write('\t'.join([str(r) for r in results]) + '\n')

    # clean up files
    shutil.rmtree(kwinputs['output_root'])
    os.remove(kwinputs['config_file'])

    return results


if __name__ == "__main__":
    main()
