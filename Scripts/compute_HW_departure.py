"""
Computes allele and genotype freuency, as well as deviations from HW.

JCM 20201007
"""

import sys
import os
import argparse
from scipy import stats


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

    def compute_chi_sq_iter(self, obs, exp):
        """Return chi-square summation iteration."""
        return (obs - exp) ** 2 / exp

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

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.compute_hw_parser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        input_vcf = args['input_vcf']

        count_homo_a = 0
        count_hetero = 0
        count_homo_b = 0
        f = open(input_vcf, 'r')
        p_values = []
        chi_sq_values = []
        for line in f:
            if ";MT=1;" in line:
                count_homo_a += line.count('0|0')
                count_hetero += line.count('1|0')
                count_hetero += line.count('0|1')
                count_homo_b += line.count('1|1')
                allele_count_a = 2 * count_homo_a + count_hetero
                allele_count_b = 2 * count_homo_b + count_hetero
                total_allele_count = allele_count_a + allele_count_b

                obs_homo_a = count_homo_a / (total_allele_count / 2)
                obs_hetero = count_hetero / (total_allele_count / 2)
                obs_homo_b = count_homo_b / (total_allele_count / 2)
                allele_freq_a = allele_count_a / total_allele_count
                allele_freq_b = allele_count_b / total_allele_count
                exp_homo_a = allele_freq_a ** 2
                exp_hetero = 2 * allele_freq_a * allele_freq_b
                exp_homo_b = allele_freq_b ** 2
                f_obs = [obs_homo_a, obs_hetero, obs_homo_b]
                f_exp = [exp_homo_a, exp_hetero, exp_homo_b]
                p_values.append(stats.chisquare(f_obs, f_exp)[1])
                chi_sq_values.append(stats.chisquare(f_obs, f_exp)[0])

                chi_sq_value = self.compute_chi_sq_iter(
                        obs_homo_a, exp_homo_a) + \
                    self.compute_chi_sq_iter(
                        obs_hetero, exp_hetero) + \
                    self.compute_chi_sq_iter(
                        obs_homo_b, exp_homo_b)
        f.close()
        p_values.sort()
        chi_sq_values.sort()
        print(p_values[0])
        # print(chi_sq_values)


if __name__ == '__main__':
    ComputeHardyWeinbergDeparture().main()
