"""
Computes allele and genotype freuency, as well as deviations from HW.

JCM 20201007
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


class ComputeHardyWeinbergDeparture():
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
        """Return *argparse.ArgumentParser* for ``fitdadi_infer_DFE.py``."""
        parser = ArgumentParserNoArgHelp(
            description=(
                'Given a `.vcf` file, this script computes the genotype and '
                'allele frequences, as well as computes the departure from '
                'Hardy-Weinberg equilibrium, if any.'),
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
        # if allele_count_a < 0:
        #     print(allele_count_a)
        # if allele_count_b < 0:
        #     print(allele_count_a)
        if obs_homo_a < 0:
            print(obs_homo_a)
        # if obs_hetero < 0:
        #     print(obs_hetero)
        # if obs_homo_b < 0:
        #     print(obs_homo_b)
        # if (2 * num_ind) <= 0:
        #     print(2 * num_ind)
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
        input_vcf = args['input_vcf']

        # Initialize empty lists for dataframe
        chrom = []
        pos = []
        allele_count = []
        num_heterozygotes = []
        num_homozygotes = []
        chi_pr_values = []
        chi_p_n_homo_leq = []
        chi_p_n_homo_geq = []
        chi_p_two_sided = []
        fisher_pr_values = []
        fisher_p_n_homo_leq = []
        fisher_p_n_homo_geq = []
        fisher_p_two_sided = []
        snp_count = 0
        f = open(input_vcf, 'r')
        allele_count_dict = {}

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
                if n_ton >= 2:
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
                    this_combination = str([obs_homo_a,
                                            obs_hetero, obs_homo_b])
                    n_combinations = int(n_ton / 2) + 1
                    this_combination_vectors = []
                    this_fisher_pr_values = []
                    this_chi_pr_values = []
                    for i in range(n_combinations):
                        if num_ind - n_ton + i >= 0:
                            combination_vector = str([num_ind - n_ton + i,
                                                      n_ton - 2 * i, i])
                            this_combination_vectors.append(
                                combination_vector)
                            chi_pr = self.compute_chi_p(
                                num_ind, 2 * num_ind - n_ton, n_ton,
                                num_ind - n_ton + i, n_ton - 2 * i, i)
                            fisher_pr = self.compute_fisher_p(
                                num_ind, 2 * num_ind - n_ton, n_ton,
                                num_ind - n_ton + i, n_ton - 2 * i, i)
                            this_chi_pr_values.append(chi_pr)
                            this_fisher_pr_values.append(fisher_pr)

                    # One sided p-value with less than or equal homozygotes.
                    this_chi_p_n_homo_leq = []
                    this_fisher_p_n_homo_leq = []
                    for i in range(len(this_fisher_pr_values)):
                        timer = i
                        temp_chi = 0
                        temp_fisher = 0
                        while timer >= 0:
                            temp_chi += this_chi_pr_values[timer]
                            temp_fisher += this_fisher_pr_values[timer]
                            timer -= 1
                        this_chi_p_n_homo_leq.append(temp_chi)
                        this_fisher_p_n_homo_leq.append(temp_fisher)

                    # One sided p-value with greater than or equal homozygotes.
                    this_chi_p_n_homo_geq = []
                    this_fisher_p_n_homo_geq = []
                    for i in range(len(this_fisher_pr_values)):
                        timer = i
                        temp_chi = 0
                        temp_fisher = 0
                        while timer <= ((len(this_fisher_pr_values)) - 1):
                            temp_chi += this_chi_pr_values[timer]
                            temp_fisher += this_fisher_pr_values[timer]
                            timer += 1
                        this_chi_p_n_homo_geq.append(temp_chi)
                        this_fisher_p_n_homo_geq.append(temp_fisher)

                    # Two sided p-value with every more extreme case.
                    this_chi_p_two_sided = []
                    this_fisher_p_two_sided = []
                    for i in range(len(this_fisher_pr_values)):
                        temp_chi = 0
                        temp_fisher = 0
                        for j in range(len(this_fisher_pr_values)):
                            if this_chi_pr_values[j] <= this_chi_pr_values[i]:
                                temp_chi += this_chi_pr_values[j]
                            if this_fisher_pr_values[j] <= \
                               this_fisher_pr_values[i]:
                                temp_fisher += this_fisher_pr_values[j]
                        this_chi_p_two_sided.append(temp_chi)
                        this_fisher_p_two_sided.append(temp_fisher)

                    this_index = this_combination_vectors.index(
                        this_combination)
                    chi_pr_values.append(this_chi_pr_values[this_index])
                    chi_p_n_homo_leq.append(this_chi_p_n_homo_leq[this_index])
                    chi_p_n_homo_geq.append(this_chi_p_n_homo_geq[this_index])
                    chi_p_two_sided.append(this_chi_p_two_sided[this_index])

                    fisher_pr_values.append(this_fisher_pr_values[this_index])
                    fisher_p_n_homo_leq.append(
                        this_fisher_p_n_homo_leq[this_index])
                    fisher_p_n_homo_geq.append(
                        this_fisher_p_n_homo_geq[this_index])
                    fisher_p_two_sided.append(this_fisher_p_two_sided[this_index])

                    ###
                    # allele_count_a = 2 * obs_homo_a + obs_hetero
                    # allele_count_b = 2 * obs_homo_b + obs_hetero
                    # total_allele_count = allele_count_a + allele_count_b
                    # num_ind = total_allele_count / 2

                    # chi-squared
                    # allele_freq_a = allele_count_a / total_allele_count
                    # allele_freq_b = allele_count_b / total_allele_count
                    # exp_homo_a = allele_freq_a ** 2 * total_allele_count / 2
                    # exp_hetero = allele_freq_a * allele_freq_b * \
                    #     total_allele_count
                    # exp_homo_b = allele_freq_b ** 2 * total_allele_count / 2
                    # f_obs = [obs_homo_a, obs_hetero, obs_homo_b]
                    # f_exp = [exp_homo_a, exp_hetero, exp_homo_b]
                    # chi_pr_values.append(stats.chisquare(f_obs, f_exp)[1])

                    # Fisher's exact test for HWE
                    # fisher_numerator = math.factorial(num_ind) * \
                    #     math.factorial(allele_count_a) * \
                    #     math.factorial(allele_count_b) * \
                    #     2 ** obs_hetero
                    # fisher_denominator = math.factorial(obs_homo_a) * \
                    #    math.factorial(obs_hetero) * \
                    #     math.factorial(obs_homo_b) * \
                    #     math.factorial(2 * num_ind)
                    # fisher_pr_values.append(
                    #     fisher_numerator / fisher_denominator)
        f.close()

        allele_count_df = pd.DataFrame.from_dict(
            allele_count_dict,
            orient='index',
            columns=['heterozygotes', 'homozygotes'])
        allele_count_df.rename_axis('allele_count', axis='columns')

        summary_df = pd.DataFrame({
            'chrom': chrom,
            'pos': pos,
            'allele_count': allele_count,
            'num_heterozygotes': num_heterozygotes,
            'num_homozygotes': num_homozygotes,
            'chi_pr_values': chi_pr_values,
            # 'chi_p_n_homo_leq': chi_p_n_homo_leq,
            # 'chi_p_n_homo_geq': chi_p_n_homo_geq,
            # 'chi_p_two_sided': chi_p_two_sided,
            'fisher_pr_values': fisher_pr_values,
            'fisher_p_n_homo_leq': fisher_p_n_homo_leq,
            'fisher_p_n_homo_geq': fisher_p_n_homo_geq,
            'fisher_p_two_sided': fisher_p_two_sided})

        # Chi-squared test for HWE
        chi_pr_values.sort()
        low_chi_p_count = 0
        sig_chi_p_count = 0
        for val in chi_pr_values:
            if val <= 0.5:
                low_chi_p_count += 1
            if val <= 0.05:
                sig_chi_p_count += 1
        low_chi_p_freq = low_chi_p_count / len(chi_pr_values)
        sig_chi_p_freq = sig_chi_p_count / len(chi_pr_values)

        fisher_p_n_homo_geq.sort()
        # print(fisher_pr_values[0])
        low_fisher_p_count = 0
        sig_fisher_p_count = 0
        for val in fisher_p_n_homo_geq:
            if val <= 0.5:
                low_fisher_p_count += 1
            if val <= 0.05:
                sig_fisher_p_count += 1
        low_fisher_p_freq = low_fisher_p_count / len(fisher_p_n_homo_geq)
        sig_fisher_p_freq = sig_fisher_p_count / len(fisher_p_n_homo_geq)

        with open('../Data/output.txt', 'a') as f:
            f.write('Outputting summary results for ' + str(input_vcf) + '.\n')
            f.write('There are ' + str(snp_count) + ' SNPs in this sample.\n')
            f.write('The minimum chi-squared p-value is ' +
                    str(chi_pr_values[0]) + '.\n')
            f.write('The minimum Fisher p-value is ' +
                    str(fisher_p_n_homo_geq[0]) + '.\n')
            f.write('The proportion of chi-squared p-values below 0.5 is ' +
                    str(low_chi_p_freq) + '.\n')
            f.write('The proportion of chi-squared p-values below 0.05 is ' +
                    str(sig_chi_p_freq) + '\n.')
            f.write('The proportion of Fisher p-values below 0.5 is ' +
                    str(low_fisher_p_freq) + '.\n')
            f.write('The proportion of Fisher p-values below 0.05 is ' +
                    str(sig_fisher_p_freq) + '\n.')
            f.write('\n')

        table_df = pd.DataFrame.from_records([
            {'sample_size': num_ind,
             'snp_count': snp_count,
             'minimum_chi_p': chi_pr_values[0],
             'minimum_fisher_p': fisher_p_n_homo_geq[0],
             'chi < 0.5': low_chi_p_freq,
             'chi < 0.05': sig_chi_p_freq,
             'fisher < 0.5': low_fisher_p_freq,
             'fisher < 0.05': sig_fisher_p_freq}])
        table_csv = input_vcf.replace('.vcf', '_table.csv')
        table_df.to_csv(path_or_buf=table_csv, index=False)

        output_csv = input_vcf.replace('.vcf', '_summary.csv')
        output_allele_count = input_vcf.replace('.vcf', '_allele_count.csv')
        summary_df.to_csv(path_or_buf=output_csv, index=False)
        allele_count_df.to_csv(path_or_buf=output_allele_count, index=True)


if __name__ == '__main__':
    ComputeHardyWeinbergDeparture().main()
