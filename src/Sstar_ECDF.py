import click
from typing import List
import os
import pandas as pd
import numpy as np


@click.group()
def main():
    pass


@main.command()
@click.option('--region-dir', default='',
              help='directory containing null windowcalc files')
@click.option('--chr-list',
              help='File location of chromosome list')
@click.option('--outfile',
              help='File location for pandas dataframe pickle database')
def build_null_db(region_dir, chr_list, outfile):
    '''
    Generate the null database for producing ECDF estimates
    '''
    db = Null_DB()
    for filename in region_files(region_dir, get_chromosomes(chr_list)):
        click.echo(filename)
        db.read_windowcalc(filename)

    click.echo(f'Finished reading, writing to {outfile}')
    db.save(outfile)
    click.echo('Done')


@main.command()
@click.option('--outfile',
              help='File location for combined pandas dataframe')
@click.argument('input-dbs', nargs=-1)
def combine_null_dbs(outfile, input_dbs):
    '''
    Combine several null databases into a single pickle
    '''
    db = Null_DB()
    db.load(input_dbs[0])
    click.echo(input_dbs[0])
    for filename in input_dbs[1:]:
        db2 = Null_DB()
        db2.load(filename)
        db.combine(db2)
        click.echo(filename)

    db.save(outfile)
    click.echo('Done')


def region_files(region_dir: str,
                 chromosomes: List[str],
                 basename: str = '{chrom}.windowcalc.gz'):
    '''
    yield filenames of the region directory and chromosomes
    will join region dir with basename, formatting basename with chrom=chrom
    '''
    for chrom in chromosomes:
        yield os.path.join(region_dir, basename.format(chrom=chrom))


def get_chromosomes(chromosome_list: str):
    '''
    read in chromosomes from a chromosome list, yielding results
    '''
    with open(chromosome_list, 'r') as infile:
        for line in infile:
            yield line.strip()


@main.command()
@click.option('--null-db', required=True,
              help='DB pickle created by build-null-db')
@click.option('--sstar-pval', default=0.01,
              help='pvalue for filtering on sstar ECDF')
@click.option('--match-pval', default=0.05,
              help='pvalue for filtering on match statistic')
@click.option('--region-dir', default='',
              help='directory containing admixed windowcalc files')
@click.option('--tsv-dir', default='',
              help='directory containing admixed tsv files from match pct')
@click.option('--chr-list',
              help='File location of chromosome list')
@click.option('--chrom', default=-1,
              help='Perform analysis on one chromosome')
@click.option('--outfile',
              help='Output file with \'{chrom}\' wildcard')
def generate_bed(null_db,
                 sstar_pval,
                 match_pval,
                 region_dir,
                 tsv_dir,
                 chr_list,
                 chrom,
                 outfile):
    '''
    Generate bed file
    '''
    # load null db to use throughout
    null = Null_DB()
    null.load(null_db)

    # for each chromosome in regionfiles
    if chrom == -1:
        chroms = list(get_chromosomes(chr_list))
    else:
        chroms = [chrom]
    for window_file, match_file, chrom in zip(
            region_files(region_dir, chroms),
            region_files(tsv_dir, chroms, '{chrom}.tsv.gz'),
            chroms
    ):
        # open windowcalc
        window = process_windowcalc(window_file)

        window = filter_by_match(window, match_file, match_pval)

        window = filter_by_sstar(window, null, sstar_pval)

        if '{chrom}' in outfile:
            output = outfile.format(chrom=chrom)
        else:
            output = outfile
        window.to_csv(output, sep='\t', index=False)


def process_windowcalc(filename, required_hap_frac=0.8):
    '''
    Read in the windowcalc file and process until other files are needed
    Returns a pandas dataframe
    '''
    result = pd.read_csv(filename, sep='\t', usecols=[
        'chrom', 'winstart', 'winend', 'n_region_ind_snps', 'ind_id',
        'pop', 's_star', 'num_s_star_snps',
        'n_s_star_snps_hap1', 'n_s_star_snps_hap2'
    ])

    # filter out ambiguous haplotypes
    result.n_s_star_snps_hap1 /= result.num_s_star_snps
    result.n_s_star_snps_hap2 /= result.num_s_star_snps
    result = result.loc[
        (result.n_s_star_snps_hap1 >= required_hap_frac) !=  # XOR
        (result.n_s_star_snps_hap2 >= required_hap_frac)]

    # determine haplotype
    result['haplotype'] = 1
    # at this point, this should be the only check to perform
    result.loc[result.n_s_star_snps_hap1 <
               result.n_s_star_snps_hap2, 'haplotype'] = 2

    # build msp_ID
    result.haplotype = result.ind_id + ':' + result['haplotype'].map(str)
    result['msp_ID'] = result.haplotype + '_' + result['chrom'].map(str)

    # rename columns
    result.rename(columns={'winstart': 'start',
                           'winend': 'end'}, inplace=True)

    # remove uneeded columns
    return result.drop(columns=['n_s_star_snps_hap1', 'n_s_star_snps_hap2',
                                'ind_id', 'chrom', 'num_s_star_snps'])


def filter_by_match(window: pd.DataFrame,
                    match_file, match_pval) -> pd.DataFrame:
    '''
    Remove rows of the window dataframe with haplotype entries in match_file
    with pvalue <= match_pval
    Returns the filtered dataframe
    '''
    match = pd.read_csv(match_file, sep='\t', usecols=[
        'start', 'end', 'haplotype', 'population', 'pvalue'])
    match = match.loc[match.pvalue <= match_pval]
    match = match.rename(columns={'population': 'pop'})

    result = pd.merge(window, match,
                      on=['start', 'end', 'haplotype', 'pop'])
    return result.drop(columns=['haplotype', 'pvalue'])


def filter_by_sstar(window, null_db, sstar_pval) -> pd.DataFrame:
    '''
    Remove rows failing the ecdf cutoff of sstar values in null_db
    '''
    ecdf_val = round(1 - sstar_pval, 4)
    if window.empty:
        return pd.DataFrame({
            'msp_ID': [],
            'start': [],
            'end': []
        })
    return pd.concat([
        table.loc[
            np.around(null_db.ecdf(k[0], k[1], table.s_star), 4) >= ecdf_val,
            ('msp_ID', 'start', 'end')
        ]
        for k, table in window.groupby(['pop', 'n_region_ind_snps'])
    ], ignore_index=True)


class Null_DB():
    def __init__(self):
        self.DB = pd.Series(dtype='int64',
                            index=pd.MultiIndex.from_tuples(
                                [],
                                names=['pop', 'n_region_ind_snps', 's_star']))

    def read_windowcalc(self, filename):
        '''
        Read the windowcalc file and add its counts to the database
        '''
        df = pd.read_csv(filename, sep='\t', usecols=[
            'winstart', 'winend', 'n_region_ind_snps', 'pop', 's_star'
        ])
        hist = df.loc[df.s_star != 0].groupby(['pop',
                                               'n_region_ind_snps',
                                               's_star'])['winstart'].count()

        self.DB = self.DB.combine(hist, lambda s1, s2: s1+s2, fill_value=0)

    def combine(self, other):
        '''
        Combine this database with another null DB object
        '''
        self.DB = self.DB.combine(other.DB, lambda s1, s2: s1+s2, fill_value=0)

    def ecdf(self, pop: str, n_region_ind_snps: int, values) -> np.array:
        '''
        Produce an ECDF object for the provided population and n_region
        DB must be sorted prior to usage (save/load to sort normally)
        '''
        sstar = self.get_sstar(pop, n_region_ind_snps)
        # calculate ecdf as cumsum/sum with 0 at index 0
        ecdf = pd.concat([pd.Series([0]), sstar.cumsum() / sstar.sum()])
        # get only closest 'right' value of ecdf on sstar index
        return ecdf.iloc[sstar.index.searchsorted(values, side='right')].values

    def get_sstar(self, pop: str, n_snps: int) -> pd.Series:
        '''
        Return the counts for the provided pop and n_region, performing
        interpolation and correction for missing values of n_region
        '''
        if pop not in self.DB:
            raise ValueError(f'Population "{pop}" not found in null database')

        pop_vals = self.DB[pop]
        n_vals = pop_vals.index.get_level_values('n_region_ind_snps')
        if n_snps in n_vals:
            return pop_vals[n_snps]

        if n_snps < n_vals.min():
            return pop_vals[n_vals.min()]

        if n_snps > n_vals.max():
            return pop_vals[n_vals.max()]

        # interpolate
        else:
            pop_vals = self.DB[pop].astype(float)
            ind = n_vals.searchsorted(n_snps)
            low_val = n_vals[ind-1]
            high_val = n_vals[ind]
            low_w = (high_val - n_snps) / (high_val - low_val)
            high_w = (n_snps - low_val) / (high_val - low_val)
            result = pop_vals[low_val].combine(
                pop_vals[high_val],
                lambda s1, s2: s1 * low_w + s2 * high_w,
                fill_value=0)
            return result

    def save(self, filename):
        self.DB.sort_index().to_pickle(filename)

    def load(self, filename):
        self.DB = pd.read_pickle(filename)


if __name__ == "__main__":
    main()
