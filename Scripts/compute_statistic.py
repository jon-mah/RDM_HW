"""
Computes expected, approximate, and actual number of doubletons by n-ton count.

JCM 20201214
"""

import math
import sys
import os
import argparse
from scipy import stats

import pandas as pd
import numpy as np


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""

    def error(self, message):
        """Print error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


class ComputeHardyWeinbergDepartureFromExpectation():
    """Wrapper class to allow functions to reference each other."""

    def ExistingFile(self, fname):
        """Return *fname* if existing file, otherwise raise ValueError."""
        if os.path.isfile(fname):
            return fname
        else:
            raise ValueError("%s must specify a valid file name" % fname)

    def IntGreaterThanZero(self, n):
        """If *n* is an integer > 1, returns it, otherwise an error."""
        try:
            n = int(n)
        except ValueError:
            sys.exit("%s is not an integer" % n)
        if n <= 0:
            raise ValueError("%d is not > 1" % n)
        else:
            return n

    def compute_hw_parser(self):
        """Return *ArgumentParser* for ``compute_HW_expectation.py``."""
        parser = ArgumentParserNoArgHelp(
            description=(
                'Given a number of individuals, this script computes '
                'the expected number of allele frequencies from that sample'),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument(
            'input_vcf', type=self.ExistingFile,
            help=('Variant Call Format file containing Single-nucleotide-'
                  'polymorphism data.'))
        return parser

    def compute_fisher_p(self, num_ind, allele_count_a, allele_count_b,
                         obs_homo_a, obs_hetero, obs_homo_b):
        """Return p-value for Fisher's exact test."""
        # Fisher's exact test for HWE
        fisher_numerator = math.factorial(num_ind) * \
            math.factorial(allele_count_a) * \
            math.factorial(allele_count_b) * \
            2 ** obs_hetero
        fisher_denominator = math.factorial(obs_homo_a) * \
            math.factorial(obs_hetero) * \
            math.factorial(obs_homo_b) * \
            math.factorial(2 * num_ind)
        return fisher_numerator / fisher_denominator

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.compute_hw_parser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        input_vcf = args['input_vcf']

        # Initialize empty lists for dataframe
        chrom = []
        pos = []
        allele_count = []
        num_heterozygotes = []
        num_homozygotes = []
        fisher_pr_values = []
        snp_count = 0
        allele_count_dict = {}
        n_homo_dict = {}
        exp_homo_dict = {}

        f = open(input_vcf, 'r')

        for line in f:
            if ";MT=1;" in line:
                obs_homo_a = 0
                obs_hetero = 0
                obs_homo_b = 0
                obs_homo_a += line.count('0|0')
                obs_hetero += line.count('1|0')
                obs_hetero += line.count('0|1')
                obs_homo_b += line.count('1|1')
                n_ton = obs_hetero + 2 * obs_homo_b
                allele_count_a = 2 * obs_homo_a + obs_hetero
                allele_count_b = 2 * obs_homo_b + obs_hetero
                total_allele_count = allele_count_a + allele_count_b
                num_ind = int(total_allele_count / 2)
                if n_ton >= 1000:
                    approx_homo_b = (n_ton / (2 * num_ind)) ** 2 * num_ind
                    if n_ton in n_homo_dict.keys():
                        n_homo_dict[n_ton].append(obs_homo_b)
                    else:
                        n_homo_dict[n_ton] = [obs_homo_b]
                        exp_homo_dict[n_ton] = approx_homo_b
                if n_ton >= 2 and n_ton <= 999:
                    # Compute summary statistics for detecing HWE departure
                    snp_count += 1

                    # Compute summary information for dataframe
                    chrom.append(line.split()[0])
                    pos.append(line.split()[1])
                    this_allele_count = obs_hetero + 2 * obs_homo_b
                    allele_count.append(this_allele_count)
                    num_heterozygotes.append(obs_hetero)
                    num_homozygotes.append(obs_homo_b)
                    if this_allele_count not in allele_count_dict:
                        allele_count_dict[this_allele_count] = \
                            np.array([obs_hetero, obs_homo_b])
                    else:
                        allele_count_dict[this_allele_count] = \
                            allele_count_dict[this_allele_count] + \
                            np.array([obs_hetero, obs_homo_b])

                    # Compute Hardy-Weinberg expectation
                    this_combination = ([obs_homo_a,
                                         obs_hetero, obs_homo_b])
                    n_combinations = int(n_ton / 2) + 1
                    this_combination_vectors = []
                    this_fisher_pr_values = []
                    for i in range(n_combinations):
                        if num_ind - n_ton + i >= 0:
                            combination_vector = [num_ind - n_ton + i,
                                                  n_ton - 2 * i, i]
                            this_combination_vectors.append(
                                combination_vector)
                            fisher_pr = self.compute_fisher_p(
                                num_ind, 2 * num_ind - n_ton, n_ton,
                                num_ind - n_ton + i, n_ton - 2 * i, i)
                            this_fisher_pr_values.append(fisher_pr)

                    this_index = this_combination_vectors.index(
                        this_combination)

                    fisher_pr_values.append(this_fisher_pr_values[this_index])

                    exp_homo_b = 0
                    for i in range(len(this_combination_vectors)):
                        exp_homo_b += (
                            this_combination_vectors[i])[2] * \
                            this_fisher_pr_values[i]

                    if n_ton in n_homo_dict.keys():
                        n_homo_dict[n_ton].append(obs_homo_b)
                    else:
                        n_homo_dict[n_ton] = [obs_homo_b]
                        exp_homo_dict[n_ton] = exp_homo_b
        f.close()

        n_ton_list = []
        avg_homo_b_list = []
        exp_homo_b_list = []
        for key in sorted(n_homo_dict.keys()):
            n_ton_list.append(key)
            avg_homo_b_list.append(
                sum(n_homo_dict[key]) / len(n_homo_dict[key]))
            exp_homo_b_list.append(exp_homo_dict[key])

        exp_minus_avg = np.subtract(exp_homo_b_list, avg_homo_b_list)
        n_homo_df = pd.DataFrame({
            'n_ton': n_ton_list,
            'avg_homo_b': avg_homo_b_list,
            'exp_homo_b': exp_homo_b_list,
            'exp_minus_avg': exp_minus_avg})

        n_homo_csv = input_vcf.replace('.vcf', '_avg_vs_exp_homo.csv')
        n_homo_df.to_csv(path_or_buf=n_homo_csv, index=True)

        summary_statistics = input_vcf.replace(
            '.vcf', '_summary_statistics.txt')
        two_sample_ks = stats.ks_2samp(
            avg_homo_b_list, exp_homo_b_list)
        with open(summary_statistics, 'w') as f:
            f.write('The KS test statistic is: ' +
                    str(two_sample_ks[0]) + '.\n')
            f.write('The KS test p-value is: ' +
                    str(two_sample_ks[1]) + '.\n')
            if (two_sample_ks[1] >= two_sample_ks[0]):
                f.write('We cannot reject the null hypothesis since '
                        'the p-value is greater than or equal to the '
                        'KS statistic.')
            else:
                f.write('We can reject the null hypothesis since '
                        'the p-value is less than the KS statistic.')


if __name__ == '__main__':
    ComputeHardyWeinbergDepartureFromExpectation().main()
