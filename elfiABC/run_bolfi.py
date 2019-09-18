import elfi
from elfi.clients import multiprocessing
from scipy.spatial.distance import euclidean
import click
import numpy as np
import os
import shutil
import subprocess


base_output = '/Genomics/akeylab/abwolf/SimulatedDemographic/ABC/{model}/bolfi'
admixed_size = 500
summary_output = ('/Genomics/akeylab/abwolf/SimulatedDemographic/ABC/'
                  f'SplitPop/bolfi/bolfi_summary.txt')
model = 'SplitPop'
empirical_order = (
    'desert-5 desert-6 desert-7 desert-8 desert-9 desert-10 '
    'pi-AF pi-EU pi-AS pi-AF-EU pi-AF-AS pi-EU-AS admix-ASN admix-EUR'
).split()

empirical = np.array([0.0375339, 0.274105, 0.0214289, # desert 5-7
                      0.0176096, 0.0164607, 0.0136497,  # desert 8-10
                      0.291761, 0.328705, 0.335145,  # pi
                      0.12985, 0.156539, 0.0921547,  # fst
                      0.016, 0.019])  # admix
command = ('snakemake '
           '--configfile {config_file} '
           '--quiet ')

@click.command()
@click.option('-b', '--base', default="",
              help='Root directory to store temporary output')
@click.option('-o', '--output', default="",
              help='Output file containing parameters and summary statistics')
@click.option('-n', '--number-processes', default=10,
              help='Number of simulations to run in this instance')
def main(base, output, number_processes):
    main_helper(base, output, number_processes)


def main_helper(base, output, number_processes):
    global base_output
    use = 'BOLFI'

    # write null yaml to override base output
    with open('null.yml', 'w') as writer:
        writer.write(f'''---
paths:
    base_output: "{base_output}"
''')

    # generate null (should only execute once)
    # click.secho('Generating null dataset', fg='yellow')
    # subprocess.check_call(["snakemake", "null",
    #                        "--profile", "elfi_profile",
    #                        "--configfile", "null.yml"])
    # click.secho('Creating environments', fg='yellow')
    # subprocess.check_call(["snakemake",
    #                  "--profile", "elfi_profile",
    #                  "--create-envs-only"])
    # click.secho('Starting BOLFI', fg='yellow')

    if base:
        base_output = base
    if output:
        global summary_output
        summary_output = output
    global command
    if use == 'LOCAL':
        command += '--profile elfi_single '
    else:
        command += '--profile elfi_profile '
    # setup priors
    elfi.set_client(multiprocessing.Client(num_processes=number_processes))
    model = elfi.ElfiModel(name='msprime')
    elfi.Prior('uniform', 0, 0.5, model=model, name='n1')
    # elfi.Prior('uniform', 0, 0.5, model=model, name='n2')
    elfi.Constant(0, model=model, name='n2')
    elfi.Prior('uniform', 0, 1, model=model, name='split_prop')
    elfi.Prior('uniform', 1, 2758, model=model, name='split_size')
    snake_sim = elfi.tools.external_operation(command,
                                              prepare_inputs=prepare_inputs,
                                              process_result=process_result,
                                              stdout=False)
    vec_snake = elfi.tools.vectorize(snake_sim)
    snake_node = elfi.Simulator(vec_snake, model['n1'],
                                model['n2'],
                                model['split_prop'],
                                model['split_size'],
                                observed=empirical,
                                name='snake_sim')
    snake_node.uses_meta = True
    distance = elfi.Distance('euclidean', snake_node, name='dist')
    if use == 'SINGLE' or use == 'LOCAL':
        if use == 'LOCAL':
            number_processes = 1
        from multiprocessing import Pool
        with Pool(number_processes) as pool:
            print(pool.starmap(
                vec_snake,
                [
                    (0.2, 0.0, 0.1, 1000),
                    (0.1, 0.0, 0.1, 1000),
                    (0.0, 0.0, 0.1, 1000),
                    (0.2, 0.0, 0.0, 1000),
                    # (0.1, 0.0, 0.0, 1000),
                    (0.0, 0.0, 0.0, 1000),
                    (0.2, 0.1, 0.1, 1000),
                    (0.1, 0.1, 0.1, 1000),
                    (0.0, 0.1, 0.1, 1000),
                    (0.2, 0.1, 0.0, 1000),
                    (0.1, 0.1, 0.0, 1000),
                    (0.0, 0.1, 0.0, 1000),
                ]))
        return

    elif use == 'BOLFI':
        try:
            pool = elfi.ArrayPool.open(name='bolfi_pool')
            click.secho('Opened existing pool', fg='green')
        except:
            pool = elfi.ArrayPool(['n1', 'split_prop',
                                   'split_size', 'snake_sim'],
                                  name='bolfi_pool')
            click.secho('Creating new pool', fg='yellow')

        bolfi = elfi.BOLFI(distance,
                           batch_size=1,
                           #initial_evidence=500,
                           initial_evidence=read_previous_results(summary_output),
                           update_interval=10,
                           bounds={'n1': (0, 0.5),
                                   # 'n2': (0, 0.5),
                                   'split_prop': (0, 1),
                                   'split_size': (1, 2758)},
                           pool=pool)
        elfi.set_client('native')
        #post = bolfi.fit(n_evidence=525)
        # click.secho('Saving results', fg='yellow')
        # pool.save()
        # pool.close()

        post = bolfi.fit(n_evidence=100)
        result = bolfi.sample(500, n_chains=4)
        result.save('sample500.pkl')
        print(result)
        click.secho('Done', fg='green')


def prepare_inputs(*inputs, **kwinputs):
    n1, n2, split_prop, split_size = inputs
    split_size = int(split_size)

    if 'seed' in kwinputs:
        seed = kwinputs['seed']
    else:
        seed = 2
        kwinputs['seed'] = seed

    if 'meta' in kwinputs:
        meta = kwinputs['meta']
        admix_base = '{model_name}_{batch_index}_{submission_index}'
        if 'index_in_batch' in meta:
            admix_base += '_{index_in_batch}'
        admix_base = admix_base.format(**meta)

    else:
        admix_base = f'single_{n1}_{n2}_{split_prop}'

    yml_file = f'{admix_base}.yaml'

    global base_output
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
    admixed_simulations: {admixed_size}
''')

    kwinputs['config_file'] = yml_file
    temp_output = base_output.format(model=model)
    kwinputs['output_root'] = f'{temp_output}/{admix_base}'
    kwinputs['output_file'] = f'{temp_output}/{admix_base}/summary.txt'

    return inputs, kwinputs


def process_result(completed_process, *inputs, **kwinputs):
    output_file = kwinputs['output_file']

    results = np.loadtxt(output_file, skiprows=1)

    if not os.path.exists(summary_output):
        # write header
        with open(summary_output, 'w') as writer:
            writer.write(
                '\t'.join(('seed n1 n2 split_prop split_size '
                           'desert-5 desert-6 desert-7 '
                           'desert-8 desert-9 desert-10 '
                           'pi-AF pi-EU pi-AS '
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


def read_previous_results(fname):
    result = {}
    with open(fname, 'r') as infile:
        headers = infile.readline().split()
        for h in headers:
            result[h] = []
        for line in infile:
            for h, v in zip(headers, line.split()):
                result[h].append(float(v))

    # create distance measure from empirical values
    result['dist'] = []
    for i in range(len(result[headers[0]])):
        obs = np.array([result[eo][i] for eo in empirical_order])
        result['dist'].append(euclidean(empirical, obs))

    for k in result:
        result[k] = np.array(result[k][:50])
    # TODO return only first 50 results!
    return result


if __name__ == "__main__":
    main()
