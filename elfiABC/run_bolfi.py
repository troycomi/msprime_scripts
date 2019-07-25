import elfi
from elfi.clients import multiprocessing
import click
import numpy as np
import os


base_output = '/tigress/tcomi/abwolf_abc/'
summary_output = '/tigress/tcomi/abwolf_abc/Sriram/summary.txt'
@click.command()
@click.option('-b', '--base', default="",
              help='Root directory to store temporary output')
@click.option('-o', '--output', default="",
              help='Output file containing parameters and summary statistics')
@click.option('-n', '--number-processes', default=5,
              help='Number of simulations to run in this instance')
def main(base, output, number_processes):
    # TODO either generate null as a separate option or add it in
    if base:
        global base_output
        base_output = base
    if output:
        global summary_output
        summary_output = output
    # setup priors
    elfi.set_client(multiprocessing.Client(num_processes=number_processes))
    model = elfi.ElfiModel(name='msprime')
    elfi.Prior('uniform', 0, 0.3, model=model, name='n1')
    snake_sim = elfi.tools.external_operation(command,
                                              prepare_inputs=prepare_inputs,
                                              process_result=process_result,
                                              stdout=False)
    print(snake_sim(0.1, seed=2, meta={
        'model_name': 'test',
        'batch_index': 1,
        'submission_index': 2}))
    return
    vec_snake = elfi.tools.vectorize(snake_sim)
    empirical = np.array([0, 0.00047, 0.00035, 0.00035, 0.1, 0.1, 0.03, 0.037])
    snake_node = elfi.Simulator(vec_snake, model['n1'], observed=empirical)
    snake_node.uses_meta = True
    distance = elfi.Distance('euclidean', snake_node)

    rej = elfi.Rejection(distance, model['n1'], batch_size=1)
    rej.sample(5)
    return

    bolfi = elfi.BOLFI(snake_node, batch_size=1,
                       initial_evidence=6, update_interval=3,
                       bounds={'n1': (0, 0.3)})
    bolfi.fit(n_evidence=10)
    return


command = 'snakemake --profile elfi_profile --configfile {config_file}'


def prepare_inputs(*inputs, **kwinputs):
    n1, = inputs
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
''')

    kwinputs['config_file'] = yml_file
    kwinputs['output_root'] = f'{base_output}Sriram/{admix_base}'
    kwinputs['output_file'] = f'{base_output}Sriram/{admix_base}/summary.txt'

    return inputs, kwinputs


def process_result(completed_process, *inputs, **kwinputs):
    output_file = kwinputs['output_file']

    results = np.loadtxt(output_file, skiprows=1)
    print(results.shape)
    return results

    if not os.path.exists(summary_output):
        # write header
        with open(summary_output, 'w') as writer:
            writer.write(
                '\t'.join(('seed n1 desert-mse pi-AF pi-EU pi-AS '
                           'pi-AF-EU pi-AF-AS pi-EU-AS admix-ASN admix-EUR').split()))
            writer.write('\n')

    with open(summary_output, 'a') as writer:
        writer.write(f'{kwinputs["seed"]}\t')
        writer.write(f'{inputs[0]}\t')
        writer.write('\t'.join([str(r) for r in results]) + '\n')

    # TODO clean up
    # os.remove...

    return results


if __name__ == "__main__":
    main()
