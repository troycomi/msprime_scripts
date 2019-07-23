from elfi import ElfiModel, Prior, tools
import click
import numpy as np


base_output = '/tigress/tcomi/abwolf_abc/'
@click.command()
@click.option('-b', '--base', default="",
              help='Root directory to store temporary output')
def main(base):
    # TODO eithe generate null as a separate option or add it in
    if base:
        global base_output
        base_output = base
    # setup priors
    model = ElfiModel(name='msprime')
    n1 = Prior('uniform', 0, 0.3, model=model, name='n1')
    snake_sim = tools.external_operation(command,
                                         prepare_inputs=prepare_inputs,
                                         process_result=process_result,
                                         stdout=False)
    print(snake_sim(0.1, seed=2, meta={
        'model_name': 'test',
        'batch_index': 1,
        'submission_index': 2}))
    # bolfi = elfi.BOLFI()


command = 'snakemake --profile elfi_profile --configfile {config_file}'


def prepare_inputs(*inputs, **kwinputs):
    n1, = inputs
    meta = kwinputs['meta']
    seed = kwinputs['seed']

    output_base = '{model_name}_{batch_index}_{submission_index}'.format(**meta)
    yml_file = f'{output_base}.yaml'

    with open(yml_file, 'w') as writer:
        writer.write(f'''---
paths:
    base_output: "{base_output}"
    admixed_dir: "__BASE_OUTPUT__/{output_base}"
    summary: "__BASE_OUTPUT__/summary.txt"

msprime:
    base_seed: {seed}
    n1: {n1}
''')

    kwinputs['config_file'] = yml_file
    kwinputs['output_file'] = f'{base_output}/summary.txt'

    return inputs, kwinputs


def process_result(completed_process, *inputs, **kwinputs):
    output_file = kwinputs['output_file']
    results = np.loadtxt(output_file, skiprows=1)

    # TODO determine how to save results, clean up
    # os.remove...

    return results


if __name__ == "__main__":
    main()
