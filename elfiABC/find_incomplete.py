import click


@click.command()
@click.option('-j', '--job-output', type=click.File('r'))
def main(job_output):
    jobs = {}
    last_output = None
    last_jobid = None
    for line in job_output:
        line = line.strip()
        if line.startswith('output:'):
            last_output = line.split(': ')[1].split()
        if line.startswith('jobid:'):
            jobs[line.split(': ')[1]] = last_output

        if line.startswith('Finished job '):
            jobs.pop(line.split()[-1][:-1])

    for k, v in jobs.items():
        #print(f'{k:<10}')
        for val in v:
            print(f'{val.strip(",")}')
            #print(f'{"":<10}{val}')


if __name__ == '__main__':
    main()
