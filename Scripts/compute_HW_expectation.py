"""
Computes allele and genotype expectation under HW equilibrium.

JCM 20201109
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


class ComputeHardyWeinbergExpectation():
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
            'num_ind', type=self.IntGreaterThanZero,
            help=('Number of individuals in calculation.'))
        parser.add_argument(
            'n_ton', type=self.IntGreaterThanZero,
            help=('Allele frequency class whose expectation is computed.'))
        parser.add_argument(
            '--compute_by_hand', type=bool, default=False,
            help=('Set this option to true if you want several calculations '
                  'by hand.'))
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

    def compute_chi_p(self, num_ind, allele_count_a, allele_count_b,
                      obs_homo_a, obs_hetero, obs_homo_b):
        """Return probability for Chi-squared test."""
        total_allele_count = allele_count_a + allele_count_b
        allele_freq_a = allele_count_a / total_allele_count
        allele_freq_b = allele_count_b / total_allele_count
        # chi-squared
        exp_homo_a = allele_freq_a ** 2 * total_allele_count / 2
        exp_hetero = allele_freq_a * allele_freq_b * \
            total_allele_count
        exp_homo_b = allele_freq_b ** 2 * total_allele_count / 2
        f_obs = [obs_homo_a, obs_hetero, obs_homo_b]
        f_exp = [exp_homo_a, exp_hetero, exp_homo_b]
        return stats.chisquare(f_obs, f_exp)[1]

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.compute_hw_parser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        num_ind = args['num_ind']
        n_ton = args['n_ton']
        compute_by_hand = args['compute_by_hand']

        if n_ton > num_ind:
            raise Exception('Sorry, the allele class must be less than the '
                            'number of individuals at this time.')

        n_combinations = int(n_ton / 2) + 1

        combination_vectors = []
        fisher_pr_values = []
        chi_pr_values = []
        for i in range(n_combinations):
            combination_vector = [num_ind - n_ton + i, n_ton - 2 * i, i]
            combination_vectors.append(combination_vector)
            chi_pr = self.compute_chi_p(num_ind, 2 * num_ind - n_ton, n_ton,
                                        num_ind - n_ton + i, n_ton - 2 * i, i)
            chi_pr_values.append(chi_pr)
            fisher_pr = self.compute_fisher_p(num_ind,
                                              2 * num_ind - n_ton,
                                              n_ton, num_ind - n_ton + i,
                                              n_ton - 2 * i, i)
            fisher_pr_values.append(fisher_pr)

        # One sided p-value with less than or equal homozygotes.
        chi_p_n_homo_leq = []
        fisher_p_n_homo_leq = []
        for i in range(len(fisher_pr_values)):
            timer = i
            temp_chi = 0
            temp_fisher = 0
            while timer >= 0:
                temp_chi += chi_pr_values[timer]
                temp_fisher += fisher_pr_values[timer]
                timer -= 1
            chi_p_n_homo_leq.append(temp_chi)
            fisher_p_n_homo_leq.append(temp_fisher)

        # One sided p-value with greater than or equal homozygotes.
        chi_p_n_homo_geq = []
        fisher_p_n_homo_geq = []
        for i in range(len(fisher_pr_values)):
            timer = i
            temp_chi = 0
            temp_fisher = 0
            while timer <= ((len(fisher_pr_values)) - 1):
                temp_chi += chi_pr_values[timer]
                temp_fisher += fisher_pr_values[timer]
                timer += 1
            chi_p_n_homo_geq.append(temp_chi)
            fisher_p_n_homo_geq.append(temp_fisher)

        # Two sided p-value with every more extreme case.
        chi_p_two_sided = []
        fisher_p_two_sided = []
        for i in range(len(fisher_pr_values)):
            temp_chi = 0
            temp_fisher = 0
            for j in range(len(fisher_pr_values)):
                if chi_pr_values[j] <= chi_pr_values[i]:
                    temp_chi += chi_pr_values[j]
                if fisher_pr_values[j] <= fisher_pr_values[i]:
                    temp_fisher += fisher_pr_values[j]
            chi_p_two_sided.append(temp_chi)
            fisher_p_two_sided.append(temp_fisher)

        summary_df = pd.DataFrame({'combination': combination_vectors,
                                   'chi_pr_values': chi_pr_values,
                                   # 'chi_p_n_homo_leq': chi_p_n_homo_leq,
                                   # 'chi_p_n_homo_geq': chi_p_n_homo_geq,
                                   # 'chi_p_two_sided': chi_p_two_sided,
                                   'fisher_pr_values': fisher_pr_values,
                                   'fisher_p_n_homo_leq': fisher_p_n_homo_leq,
                                   'fisher_p_n_homo_geq': fisher_p_n_homo_geq,
                                   'fisher_p_two_sided': fisher_p_two_sided})
        print(summary_df)

        if (compute_by_hand):
            # Doubletons
            print('Computing expectation for doubletons.')
            print('The probability of two heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 2, 2,
                                            num_ind - 2, 2, 0)) + '.')
            print('The probability of a single homozygote is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 2, 2,
                                            num_ind - 1, 0, 1)) + '.')

            # Tripletons
            print('Computing expectation for tripletons.')
            print('The probability of three heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 3, 3,
                                            num_ind - 3, 3, 0)) + '.')
            print('The probability of one homozygote and one heterozygote '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 3, 3,
                              num_ind - 2, 1, 1)) + '.')

            # 4-tons
            print('Computing expectation for 4-tons.')
            print('The probability of four heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 4, 4,
                                            num_ind - 4, 4, 0)) + '.')
            print('The probability of one homozygote and two heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 4, 4,
                              num_ind - 3, 2, 1)) + '.')
            print('The probability of two homozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 4, 4,
                                            num_ind - 2, 0, 2)) + '.')

            # 5-tons
            print('Computing expectation for 5-tons.')
            print('The probability of five heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 5, 5,
                                            num_ind - 5, 5, 0)) + '.')
            print('The probability of one homozygote and three heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 5, 5,
                              num_ind - 4, 3, 1)) + '.')
            print('The probability of two homozygotes and one heterozygote '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 5, 5,
                              num_ind - 3, 1, 2)) + '.')

            # 6-tons
            print('Computing expectation for 6-tons.')
            print('The probability of six heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 6, 6,
                                            num_ind - 6, 6, 0)) + '.')
            print('The probability of one homozygote and 4 heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 6, 6,
                                            num_ind - 5, 4, 1)) + '.')
            print('The probability of two homozygote and 4 heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 6, 6,
                                            num_ind - 4, 2, 2)) + '.')
            print('The probability of three homozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 6, 6,
                                            num_ind - 3, 0, 3)) + '.')
            # 7-tons
            print('Computing expectation for 7-tons.')
            print('The probability of seven heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 7, 7,
                                            num_ind - 7, 7, 0)) + '.')
            print('The probability of one homozygote and five heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 7, 7,
                              num_ind - 6, 5, 1)) + '.')
            print('The probability of two homozygotes and three '
                  'heterozygotes is ' + str(self.compute_fisher_p(
                                            num_ind, 2 * num_ind - 7, 7,
                                            num_ind - 5, 3, 2)) + '.')
            print('The probability of three homozygotes and one heterozygote '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 7, 7,
                              num_ind - 4, 1, 3)) + '.')

            # 8-tons
            print('Computing expectation for 8-tons.')
            print('The probability of eight heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 8, 8,
                                            num_ind - 8, 8, 0)) + '.')
            print('The probability of one homozygote and six heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 8, 8,
                              num_ind - 7, 6, 1)) + '.')
            print('The probability of two homozygote and four heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 8, 8,
                              num_ind - 6, 4, 2)) + '.')
            print('The probability of three homozygotes and two '
                  'heterozygotes is ' + str(self.compute_fisher_p(
                                            num_ind, 2 * num_ind - 8, 8,
                                            num_ind - 5, 2, 3)) + '.')
            print('The probability of four homozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 8, 8,
                                            num_ind - 4, 0, 4)) + '.')

            # 9-tons
            print('Computing expectation for 9-tons.')
            print('The probability of seven heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 9, 9,
                                            num_ind - 9, 9, 0)) + '.')
            print('The probability of one homozygote and seven heterozygotes '
                  'is ' + str(self.compute_fisher_p(num_ind, 2 * num_ind - 9,
                                                    9, num_ind - 8, 7,
                                                    1)) + '.')
            print('The probability of two homozygotes and five '
                  'heterozygotes ' + str(self.compute_fisher_p(
                                         num_ind, 2 * num_ind - 9, 9,
                                         num_ind - 7, 5, 2)) + '.')
            print('The probability of three homozygotes and three '
                  'heterozygote is ' + str(self.compute_fisher_p(
                                           num_ind, 2 * num_ind - 9, 9,
                                           num_ind - 6, 3, 3)) + '.')
            print('The probability of four homozygotes and one heterozygote '
                  'is ' + str(self.compute_fisher_p(num_ind, 2 * num_ind - 9,
                                                    9, num_ind - 5, 1,
                                                    4)) + '.')

            # 10-tons
            print('Computing expectation for 10-tons.')
            print('The probability of six heterozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 10, 10,
                                            num_ind - 10, 10, 0)) + '.')
            print('The probability of one homozygote and eight heterozygotes '
                  'is ' + str(self.compute_fisher_p(num_ind, 2 * num_ind - 10,
                                                    10, num_ind - 9, 8,
                                                    1)) + '.')
            print('The probability of two homozygote and six heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 10, 10,
                              num_ind - 8, 6, 2)) + '.')
            print('The probability of three homozygotes and four '
                  'heterozygotes is ' + str(self.compute_fisher_p(
                                            num_ind, 2 * num_ind - 10, 10,
                                            num_ind - 7, 4, 3)) + '.')
            print('The probability of four homozygotes and two heterozygotes '
                  'is ' + str(self.compute_fisher_p(
                              num_ind, 2 * num_ind - 10, 10,
                              num_ind - 6, 2, 4)) + '.')
            print('The probability of five homozygotes is ' +
                  str(self.compute_fisher_p(num_ind, 2 * num_ind - 10, 10,
                                            num_ind - 5, 0, 5)) + '.')
        return summary_df


if __name__ == '__main__':
    ComputeHardyWeinbergExpectation().main()
