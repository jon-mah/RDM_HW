"""
Computes allele and genotype freuency, as well as deviations from HW.

JCM 20201007
"""

import sys
import os
import argparse


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

    def compute_chi_square_iteration(self, observed, expected):
        """Return chi-square summation iteration."""
        return (observed - expected) ** 2 / expected

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
        for line in f:
            count_homo_a += line.count('0|0')
            count_hetero += line.count('1|0')
            count_hetero += line.count('0|1')
            count_homo_b += line.count('1|1')
        f.close()
        allele_count_a = 2 * count_homo_a + count_hetero
        allele_count_b = 2 * count_homo_b + count_hetero
        total_allele_count = allele_count_a + allele_count_b
        print(total_allele_count)

        observed_homo_a = count_homo_a / (total_allele_count / 2)
        observed_hetero = count_hetero / (total_allele_count / 2)
        observed_homo_b = count_homo_b / (total_allele_count / 2)
        allele_frequency_a = allele_count_a / total_allele_count
        allele_frequency_b = allele_count_b / total_allele_count
        expected_homo_a = allele_frequency_a ** 2
        expected_hetero = 2 * allele_frequency_a * allele_frequency_b
        expected_homo_b = allele_frequency_b ** 2
        chi_square_statistic = self.compute_chi_square_iteration(
                observed_homo_a, expected_homo_a) + \
            self.compute_chi_square_iteration(
                observed_hetero, expected_hetero) + \
            self.compute_chi_square_iteration(
                observed_homo_b, expected_homo_b)

        print('The observed frequency of A homozygotes is ' +
              str(observed_homo_a))
        print('The observed frequency of heterozygotes is ' +
              str(observed_hetero))
        print('The observed frequency of B homozygotes is ' +
              str(observed_homo_b))
        print('The expected frequency of A homozygotes is ' +
              str(expected_homo_a))
        print('The expected frequency of heterozygotes is ' +
              str(expected_hetero))
        print('The expected frequency of B homozygotes is ' +
              str(expected_homo_b))
        print('The computed chi-square statstic is ' +
              str(chi_square_statistic))


if __name__ == '__main__':
    ComputeHardyWeinbergDeparture().main()
