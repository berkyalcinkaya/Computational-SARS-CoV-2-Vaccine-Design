#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import tempfile
import logging
logging.basicConfig(level=logging.WARNING, format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', datefmt='%Y-%m-%d:%H:%M:%S',)

# adding all methods to the python path
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir )
methods_dir = script_dir + '/../method'
methods_base_dirs = (
    'allele-info',
    'iedbtools-utilities',
    'netmhc-3.4-executable',
    'netmhc-4.0-executable',
    'netmhccons-1.1-executable',
    'netmhcpan-2.8-executable',
    'netmhcstabpan-1.0-executable',
    'pickpocket-1.1-executable',
    'netmhcpan-4.0-executable',
    'netmhcpan-4.1-executable',
    'mhci-netmhcpan-4.0-percentile-data',
    'mhci-netmhcpan-4.0-el-percentile-data-',
    'mhci-netmhcpan-4.1-ba-percentile-data',
    'mhci-netmhcpan-4.1-el-percentile-data',
    'mhci-ann-predictor-percentile-data',
    'mhci-netmhcstabpan-percentile-data',
    'mhci-netmhccons-percentile-data'
)
for method_base_dir in methods_base_dirs:
    sys.path.append(methods_dir + '/' + method_base_dir)

from allele_info import is_user_defined_allele

from seqpredictor import MHCBindingPredictions, MHCIPercentilesCalculator, SMMMatrix, PredictorSet, SetupInfo
from util import InputError, UnexpectedInputError, PredictorError, PeptideSequenceInput, get_species, InputData, get_peptides, MHCIAlleleData, MethodSet, get_mhc_list
from setupinfo import SetupInfo

from optparse import OptionParser
from functools import reduce

from netmhc_4_0_executable import pep_score_predict_from_peptide_file as predict_from_peptide_file_netmhc_4_0
from netmhcstabpan_1_0_executable import pep_score_predict_from_peptide_file as predict_from_peptide_file_netmhcstabpan

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_percentile_for_score(score, allele, peptide, method) :
    ''' Returns the percentile scores for the raw scores passed.
    '''
    from mhci_ann_predictor_percentile_data import score_distributions as netmhc_4_score_distributions
    from mhci_netmhcstabpan_predictor_percentile_data import score_distributions as netmhcstabpan_score_distributions
    if method == 'ann':
        percentiles_calculator = MHCIPercentilesCalculator(netmhc_4_score_distributions)
        allele = allele.replace("*","")
    elif method == 'netmhcstabpan':
        percentiles_calculator = MHCIPercentilesCalculator(netmhcstabpan_score_distributions)


    binding_length = len(peptide)
    try:
        percentile = percentiles_calculator.get_percentile_scores(
                    [score,], method, allele, binding_length)[0]
    except ValueError:
        raise

    return percentile

def get_percentiles_for_scores_netmhcpan(raw_scores, alleles, peptides, el=False):
    ''' Returns the percentile scores for the raw scores passed.
    '''
    from mhci_netmhcpan_4_1_el_percentile_data import score_distributions as netmhcpan_4_el_score_distributions
    from mhci_netmhcpan_4_1_ba_percentile_data import score_distributions as netmhcpan_4_score_distributions
    if el:
        percentiles_calculator = MHCIPercentilesCalculator(netmhcpan_4_el_score_distributions)
        method_name = 'netmhcpan_el'
    else:
        percentiles_calculator = MHCIPercentilesCalculator(netmhcpan_4_score_distributions)
        method_name = 'netmhcpan_ba'
    percentiles = []
    for score, allele, peptide in zip(raw_scores, alleles, peptides):
        binding_length = len(peptide)
        try:
            percentile = percentiles_calculator.get_percentile_scores(
                        [score,], method_name, allele, binding_length)[0]
        except ValueError:
            raise
        percentiles.append(percentile)
    return percentiles

def read_peptides(fname):
	with open(fname, 'r') as r_file:
		peptides = r_file.readlines()
		peptides = [row.strip() for row in peptides if row.strip()]
		return peptides

def group_peptides_by_length(peptide_list):
	peptide_groups_by_length = []
	lengths = set(map(len, peptide_list))
	for length in lengths:
		peptide_groups_by_length.append([pep for pep in peptide_list if len(pep) == length])
	return peptide_groups_by_length

class Prediction:

    def __init__(self):
        self.version = '20130222'
        self.row_data = []

    @staticmethod
    def is_valid_file(arg):
        """ Check if arg is a valid file that already exists on the file system. """
        arg = os.path.abspath(arg)
        if not os.path.exists(arg):
            sys.stderr.write("The file {} does not exist!\n".format(arg))
            exit(1)
        else:
            return arg

    @staticmethod
    def read_protein(fname):
        file_contents = open(fname, 'r').read()
        protein = PeptideSequenceInput(file_contents)
        return protein

    @staticmethod
    def insert_dash(method, actual_methods_used, score_list):
        scores = []
        consensus_methods = ['ann', 'smm', 'comblib_sidney2008']
        m = map(lambda v: v in actual_methods_used, consensus_methods)
        m = list(m)
        dashes = ['-', '-']
        for score in map(list, score_list):
            if not m[0]:
                score.extend(dashes)
            if not m[1]:
                score.extend(dashes)
            if not m[2]:
                score.extend(dashes)
            if method == 'IEDB_recommended':
                score.extend(dashes)
            scores.append(tuple(score))
        return scores

    def commandline_input_peptide_file(self, args):
        """ Make predictions given user provided list of sequences. The input sequence is in peptides format. """
        method, input_allele, input_length, fname = args
        # for method netmhcpan_ba, recommended
        method = method.replace('netmhcpan_ba', 'netmhcpan').replace('IEDB_recommended', 'netmhcpan_el')

        allele_list = input_allele.split(",")

        if method == 'netmhcpan_el':
            el = True
            score_unit = 'Score'
        else:
            el = False
            score_unit = 'IC50'

        from netmhcpan_4_1_executable import predict_many_peptides_file
        scores_by_peptides = predict_many_peptides_file(fname, allele_list, el=el)

        alleles, peptides = zip(*scores_by_peptides.keys())
        scores, cores, icores = zip(*scores_by_peptides.values())
        percentiles = get_percentiles_for_scores_netmhcpan(scores, alleles, peptides, el=el)

        DNA_sequence_input = False # Changed to that any of the input sequence is DNA means True

        self.is_valid_file(fname)

        combined_table_rows = list(zip(alleles,peptides,cores,icores,scores,percentiles))

        header = ('allele','peptide','core','icore',score_unit,'rank')
        combined_table_rows.sort(key=lambda tup: tup[4])

        print('\t'.join(header))

        if method in ['netmhcstabpan', 'netmhcpan_el',]:
            combined_table_rows.reverse()
        for row in combined_table_rows:
            print( '\t'.join(map(str, row)))
        if DNA_sequence_input:
            eprint( '# Warning: Potential DNA sequence(s) found! This tool is intended to predict for amino acid sequences. Please double check your input fasta file.')


    def commandline_input_smm(self, args):
        """ Make predictions given user provided list of sequences. The input sequence is in peptides format. """
        (method, input_allele, input_length, fname) = args


        allele_list = input_allele.split(",")

        if method == 'smm':
            score_unit = 'IC50'

        combined_table_rows = []
        for allele in allele_list:
            setupinfo = SetupInfo(version=self.version)
            predictor_set = PredictorSet(setupinfo)
            smm_predictor = predictor_set.get_predictor(method)
            #smm_predictor = SMMMatrix()
            smm_predictor.initialize(input_allele, int(input_length))

            for peptide, score, percentile in smm_predictor.predict_peptides_file(fname):
                combined_table_rows.append((allele, peptide, score, percentile))

        DNA_sequence_input = False # Changed to that any of the input sequence is DNA means True

        self.is_valid_file(fname)

        header = ('allele','peptide', score_unit,'rank')
        combined_table_rows.sort(key=lambda tup: tup[2])

        print('\t'.join(header))

        if method in ['netmhcstabpan', 'netmhcpan_el',]:
            combined_table_rows.reverse()
        for row in combined_table_rows:
            print( '\t'.join(map(str, row)))
        if DNA_sequence_input:
            eprint( '# Warning: Potential DNA sequence(s) found! This tool is intended to predict for amino acid sequences. Please double check your input fasta file.')


    def commandline_input_netmhcstabpan(self, args):
        """ Make predictions given user provided list of sequences. The input sequence is in peptides format. """
        (method, input_allele, input_length, fname) = args


        allele_list = input_allele.split(",")

        if method == 'netmhcstabpan':
            predict_many_peptides_file = predict_from_peptide_file_netmhcstabpan
            score_unit = 'Score'

        to_delete = []

        combined_table_rows = []
        peptide_list = read_peptides(fname)
        for peptides in group_peptides_by_length(peptide_list):
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_peptides_file:
                to_delete.append(tmp_peptides_file.name)
                tmp_peptides_file.write('\n'.join(peptides))
            for allele in allele_list:
                peptide_score_tuple = predict_many_peptides_file(allele, 9, tmp_peptides_file.name)
                for peptide, score in peptide_score_tuple:
                    score = float(score)
                    percentile = get_percentile_for_score(score, allele, peptide, method)
                    combined_table_rows.append((allele, peptide, score, percentile))
        # to delete tmp file
        for tmp_file in to_delete:
            os.remove(tmp_file)

        DNA_sequence_input = False # Changed to that any of the input sequence is DNA means True

        self.is_valid_file(fname)

        header = ('allele','peptide', score_unit,'rank')
        combined_table_rows.sort(key=lambda tup: tup[2])

        print('\t'.join(header))

        if method in ['netmhcstabpan', 'netmhcpan_el',]:
            combined_table_rows.reverse()
        for row in combined_table_rows:
            print( '\t'.join(map(str, row)))
        if DNA_sequence_input:
            eprint( '# Warning: Potential DNA sequence(s) found! This tool is intended to predict for amino acid sequences. Please double check your input fasta file.')


    def commandline_input_peptide(self, args):
        """ Make predictions given user provided list of sequences. The input sequence is in peptides format. """
        (method, input_allele, input_length, fname, seq_file_type) = args
        # for method netmhcpan_ba, recommended
        method = method.replace('netmhcpan_ba', 'netmhcpan').replace('IEDB_recommended', 'netmhcpan_el')

        if method in ['netmhcpan_ba', 'netmhcpan', 'netmhcpan_el']:
            return self.commandline_input_peptide_file((method, input_allele, input_length, fname))

        elif method not in ['ann', 'netmhcstabpan', 'smm']:
            return self.commandline_input((method, input_allele, input_length, fname))

        elif method in ['netmhcstabpan']:
            peptide_list = read_peptides(fname)
            # if not all peptides have same length
            if len(set(map(len, peptide_list))) != 1:
                return self.commandline_input_netmhcstabpan((method, input_allele, input_length, fname))

        elif method in ['smm']:
            return self.commandline_input_smm((method, input_allele, input_length, fname))

        allele_list = input_allele.split(",")

        if method == 'netmhcstabpan':
            predict_many_peptides_file = predict_from_peptide_file_netmhcstabpan
            score_unit = 'Score'
        elif method == 'ann':
            predict_many_peptides_file = predict_from_peptide_file_netmhc_4_0
            score_unit = 'IC50'

        combined_table_rows = []
        for allele in allele_list:
            peptide_score_tuple = predict_many_peptides_file(allele, 9, fname)
            for peptide, score in peptide_score_tuple:
                score = float(score)
                percentile = get_percentile_for_score(score, allele, peptide, method)
                combined_table_rows.append((allele, peptide, score, percentile))

        DNA_sequence_input = False # Changed to that any of the input sequence is DNA means True

        self.is_valid_file(fname)

        header = ('allele','peptide', score_unit,'rank')
        combined_table_rows.sort(key=lambda tup: tup[2])

        print('\t'.join(header))

        if method in ['netmhcstabpan', 'netmhcpan_el',]:
            combined_table_rows.reverse()
        for row in combined_table_rows:
            print( '\t'.join(map(str, row)))
        if DNA_sequence_input:
            eprint( '# Warning: Potential DNA sequence(s) found! This tool is intended to predict for amino acid sequences. Please double check your input fasta file.')

    def commandline_input(self, args):
        """ Make predictions given user provided list of sequences. The input sequence is in fasta format. """
        (method, input_allele, input_length, fname) = args
        # for method netmhcpan_ba, recommended
        method = method.replace('netmhcpan_ba', 'netmhcpan').replace('IEDB_recommended', 'netmhcpan_el')

        alleles = input_allele.split(",")
        lengths = input_length.split(",")
        DNA_sequence_input = False # Changed to that any of the input sequence is DNA means True

        # check if number of alleles and lengths are same
        if len(alleles) != len(lengths):
            sys.stderr.write("ERROR: Number of alleles and corresponding lengths are not equal!\n")
            exit(1)

        self.is_valid_file(fname)

        species = [get_species(allele) for allele in alleles]
        negative_inputs = self.check_for_negative_inputs(method, alleles, lengths, species) if not method == "IEDB_recommended" else []
        if negative_inputs:
            for negative_input in negative_inputs:
                allele, length, species = negative_input
                sys.stderr.write("ERROR: length '{}' for allele '{}' doesn't exist!\n".format(length, allele))
                exit(1)

        combined_table_rows = []
        for allele, length in zip(alleles, lengths):
            proteins = self.read_protein(fname)
            hla_seq = ''
            for sequences in proteins.as_amino_acid_text():
                if len(sequences) < int(''.join(length)):
                    sys.stderr.write("ERROR: Input sequence is too short!\n")
                    exit(1)
                for amino_acid in sequences.strip():
                    if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                        sys.stderr.write("Sequence: '%s' contains an invalid character: '%c' at position %d.\n" % (sequences, amino_acid, sequences.find(amino_acid)))
                        exit(1)

                # Check if string is DNA sequence
                if DNA_sequence_input or re.match('^[ACGT]+$', sequences.upper()):
                    DNA_sequence_input = True
                else:
                    DNA_sequence_input = False
            use_cutoff = cutoff_value = None
            input = InputData(self.version, method, allele, hla_seq, length, proteins, species)
            mhc_predictor = MHCBindingPredictions(input)
            mhc_scores = mhc_predictor.predict(input.input_protein.as_amino_acid_text())
            logging.debug('mhc_scores:%s' % str(mhc_scores))
            table_rows = self.format_binding(input, mhc_scores, method, use_cutoff, cutoff_value)
            logging.debug('table_rows:%s' % str(table_rows))
            method_used = ','.join(mhc_predictor.get_method_set_selected(method))
            table_rows.sort(key=lambda tup: tup[6])
            table_rows = self.add_method_used(table_rows, method)
            logging.debug('table_rows:%s' % str(table_rows))
            combined_table_rows.extend(table_rows)
        combined_table_rows.sort(key=lambda tup: tup[2])
        logging.debug('combined_table_rows:%s' % str(combined_table_rows))
        # headers for different methods
        if method == 'IEDB_recommended':
            header = ('allele','seq_num','start','end','length','peptide','method',mhc_predictor.get_score_unit(),'ann_ic50','ann_rank','smm_ic50','smm_rank','comblib_sidney2008_score','comblib_sidney2008_rank','netmhcpan_el_score','netmhcpan_rank')
            combined_table_rows.sort(key=lambda tup: tup[7])
        elif method == 'consensus':
            header = ('allele','seq_num','start','end','length','peptide',mhc_predictor.get_score_unit(),'ann_ic50','ann_rank','smm_ic50','smm_rank','comblib_sidney2008_score','comblib_sidney2008_rank')
            combined_table_rows.sort(key=lambda tup: tup[6])
        elif method in ['netmhcpan','netmhcpan_el', 'netmhcpan_ba']:
            header = ('allele','seq_num','start','end','length','peptide','core','icore',mhc_predictor.get_score_unit(),'rank')
            combined_table_rows.sort(key=lambda tup: tup[8])
        else:
            header = ('allele','seq_num','start','end','length','peptide',mhc_predictor.get_score_unit(),'rank')
            combined_table_rows.sort(key=lambda tup: tup[6])
        print('\t'.join(header))

        if method in ['netmhcstabpan', 'netmhcpan_el',]:
            combined_table_rows.reverse()
        for row in combined_table_rows:
            print( '\t'.join(map(str, row)))
        if DNA_sequence_input:
            eprint( '# Warning: Potential DNA sequence(s) found! This tool is intended to predict for amino acid sequences. Please double check your input fasta file.')


    def modify(self, lst):
        return[tuple(self.flatten(x)) for x in lst]

    @staticmethod
    def flatten(tup):
        from itertools import chain
        return list(chain(*(i if isinstance(i, tuple) else (i,) for i in tup)))

    def format_binding(self, proteins, results, method, cutoff, value):
        for length, allele, score, methods in results:
            actual_methods_used = methods.split(",")
            if method == 'consensus' or method == 'IEDB_recommended':
                if any(m in actual_methods_used for m in ['ann', 'smm', 'comblib_sidney2008']):
                    logging.debug("any(m in actual_methods_used for m in ['ann', 'smm', 'comblib_sidney2008'])")
                    score_list = []
                    for s in score:

                        ranks_scores = reduce(lambda x, y: x + y, s[1])
                        logging.debug('ranks_scores: %s' % str(ranks_scores))
                        scores = list(zip(s[0], list(zip(*ranks_scores))))
                        logging.debug('s[0]:%s' % str(s[0]))
                        logging.debug('scores:%s' % scores)
                        scores = self.insert_dash(method, actual_methods_used, self.modify(scores))
                        score_list.append(scores)
                elif all(m in actual_methods_used for m in ['netmhcpan']):
                    score_list = []
                    for results in score:
                        # results is a tuple of the form:
                        #  (<score>, <percentile>) where <percentile tuple> and
                        # IC50 scores and their percentiles for a single sequence/allele-length
                        # prediction.
                        print(repr(results))
                        score_row = [
                            (p, '-', '-', '-', '-', '-', '-', s, p) for core, icore, s, p in results
                        ]
                        score_list.append(score_row)
                else:
                    logging.warning('actual_methods_used:%s' % actual_methods_used)
                self.add_rows_binding(allele, length, proteins, score_list, actual_methods_used, cutoff, value)
            else:
                self.add_rows_binding(allele, length, proteins, score, actual_methods_used, cutoff, value)
        return self.modify(self.row_data)

    def add_rows_binding(self, allele, length, proteins, score_list, actual_methods_used, cutoff, value):
        for (i,(sequence, predictions)) in enumerate(zip(proteins.input_protein.as_amino_acid_text(), score_list)):
            for (k, prediction) in enumerate(predictions):
                sequence_source = "%s" %(i+1)
                sequence_start = "%s" %(k + 1)
                sequence_end = "%s" %(k + int(length))
                scanned_sequence = sequence[k : k + length]
                self.row_data.append((allele, sequence_source, sequence_start, sequence_end, length, scanned_sequence, prediction, '-'.join(actual_methods_used)))

    @staticmethod
    def cons_netmhcpan(scores):
        score_list = []
        for score in scores:
            lis = list(score)
            del lis[-1]
            item2 = list(score[1])
            item2.append(score[0])
            lis.append(tuple(item2))
            score_list.append(tuple(lis))
        return score_list

    @staticmethod
    def add_method_used(table_rows, method):
        formated_data = []
        for row in table_rows:
            lis = list(row)
            if method == 'IEDB_recommended':
                if '-' not in lis[-1]:
                    lis.insert(6, lis[-1])
                else:
                    lis.insert(6, "Consensus ("+lis[-1].replace("-","/")+")")
                del lis[-1]
                formated_data.append(tuple(lis))
            else:
                del lis[-1]
                formated_data.append(tuple(lis))
        return formated_data

    def commandline_input_mhc(self, parser, fname):
        """ This version takes a file containing an mhc sequence as input.
        Make predictions given user provided list of sequences. The input sequence is in fasta format. """
        (options, args) = parser.parse_args()
        DNA_sequence_input = False
        if len(args) == 3:
            (method, length, fname) = args
        else: (method, length) = args
        self.is_valid_file(fname)

        # for method netmhcpan_ba
        method = method.replace('netmhcpan_ba', 'netmhcpan').replace('IEDB_recommended', 'netmhcpan_el')

        if method not in [ 'netmhcpan', 'netmhcpan_el']:
            parser.error('Only netmhcpan has the option to take user-provided mhc sequence.')
        mhc = 'User defined'

        mhc_fh = open(options.filename_mhc, 'r')
        hla_seq = ''.join(mhc_fh.readlines())
        allele = [''.join(mhc_fh.readlines())]

        self.check_fasta(hla_seq)

        length = length.split()
        proteins = self.read_protein(fname)
        species = ['human']

        for sequences in proteins.as_amino_acid_text():
            if len(sequences) < int(''.join(length)):
                sys.stderr.write("ERROR: Input sequence is too short!\n")
                exit(1)
            for amino_acid in sequences.strip():
                if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                    sys.stderr.write("Sequence: '%s' contains an invalid character: '%c' at position %d.\n" % (sequences, amino_acid, sequences.find(amino_acid)))
                    exit(1)

            # Check if string is DNA sequence
            if DNA_sequence_input or re.match('^[ACGT]+$', sequences.upper()):
                DNA_sequence_input = True
            else:
                DNA_sequence_input = False

        input_data = InputData(self.version, method, allele, hla_seq, length, proteins, species)

        length = length[len(length) - 1]

        mhc_predictor = MHCBindingPredictions(input_data)
        mhc_scores = mhc_predictor.predict(input_data.input_protein.as_amino_acid_text())

        (results_peptide_length, results_allele, scores, method_used) = mhc_scores[0]

        header = ['allele', 'length', 'peptide', 'core', 'icore', mhc_predictor.get_score_unit()]
        print( '\t'.join(header))

        score_list = []
        core_list = []
        icore_list = []
        for seq_index in range(len(scores)):
            seq = proteins.as_amino_acid_text()[seq_index]
            peptide_list = get_peptides(seq, int(length))
            seq_scores = scores[seq_index]

            # Flatten list of tuples to list
            seq_scores = [item for x in seq_scores for item in x]

            for ind in range(0, len(seq_scores)-3, 3) :
                core_list.append( seq_scores[ind] )
                icore_list.append( seq_scores[ind + 1] )
                score_list.append( seq_scores[ind + 2] )

            # Have to add last set
            core_list.append(seq_scores[-3])
            icore_list.append(seq_scores[-2])
            score_list.append(seq_scores[-1])

            for (peptide, score, core, icore) in zip(peptide_list, score_list, core_list, icore_list):
                row = [mhc, length, peptide, core, icore, score]
                print( '\t'.join(map(str,row)))
        if DNA_sequence_input:
            eprint('# Warning: Potential DNA sequences found! This tool is intended to predict for amino acid sequences. Please double check your input fasta file.')


    @staticmethod
    def check_fasta(fasta):
        seq_list = filter(None, [x.strip() for x in fasta.split('>')])
        seq_list = list(seq_list)
        if len(seq_list) > 1:
            print( "File must contain a single MHC sequence in fasta format.")
            sys.exit(1)

        for seq in seq_list:
            sequence = seq.split("\n")[1]
            for amino_acid in sequence.strip():
                if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                    print( "Sequence: '%s' contains an invalid character: '%c' at position %d." %(sequence, amino_acid, sequence.find(amino_acid)))
                    sys.exit(1)

    def check_for_negative_inputs(self, _method_name, allele_list, length_list, species_list):
        negatives = []

        length_list = map(int, length_list)

        # user-defined alleles don't have length-species info like allele names
        if not any([is_user_defined_allele(a) for a in allele_list]):
            miad = MHCIAlleleData()
            for allele, length, species in zip(allele_list, length_list, species_list):
                length_list = miad.get_allowed_peptide_lengths(method_name=_method_name.replace('netmhcpan_el','netmhcpan'), allele_name=allele)
                if length not in length_list:
                    negatives.append((allele, length, species))
        return negatives

    @staticmethod
    def query_yes_no(question, default="yes"):
        """Ask a yes/no question via raw_input() and return their answer.

        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

        The "answer" return value is one of "yes" or "no".
        """
        valid = {"yes": True, "y": True,  "ye": True,
                 "no": False, "n": False}
        if default == None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)

    def commandline_mhc(self, argv):
        """Return all available MHC molecules against which predictions can be made."""
        (method, mhc) = argv
        # for method netmhcpan_ba
        method = method.replace('netmhcpan_ba', 'netmhcpan').replace('IEDB_recommended', 'netmhcpan_el')

        ms = MethodSet()
        method_index = ms.get_method_index(method)

        if method_index == 3:
            print( "'arb' has been removed from the list of available methods.")
            sys.exit(1)

        mhc_list = get_mhc_list(method_index)
        print("List of available (MHC,PeptideLength) for "+method)

        if method == 'netmhcpan' or method == 'IEDB_recommended':
            answer = self.query_yes_no("The list is very long. Do you still want to print it?")
            if answer == False:
                exit
            else:
                print("Species", "\t", "MHC", "\t", "PeptideLength")
                for (species,mhc, peptide_length) in mhc_list:
                    if mhc is not None:
                        print( species, "\t", ""+mhc+"", "\t", peptide_length)
        else:
            print("Species", "\t", "MHC", "\t", "PeptideLength")
            for (species,mhc, peptide_length) in mhc_list:
                if mhc is not None:
                    print(species, "\t", ""+mhc+"", "\t", peptide_length)

    def commandline_method(self):
        """Return all available prediction methods."""
        method_list = ['ann', 'comblib_sidney2008', 'consensus', 'IEDB_recommended', 'netmhcpan_ba', 'netmhcpan_el', 'smm', 'smmpmbec', 'pickpocket', 'netmhccons', 'netmhcstabpan', ]
        print("MHC-I prediction methods:")
        print("-------------------------")
        for method in method_list:
            print(method)
        print()

    @staticmethod
    def print_method_versions():
        version_info = (
            "method\tversion\n"
            "------\t-------\n"
            "IEDB_recommended\t2020.09\n"
            "consensus\t2.18\n"
            "ann\t4.0\n"
            "netmhcpan_ba\t4.1\n"
            "netmhcpan_el\t4.1\n"
            "pickpocket\t1.1\n"
            "netmhccons\t1.1\n"
            "netmhcstabpan\t1.0\n"
            "comblib_sidney2008\t1.0\n"
            "smm\t1.0\n"
            "smmpmbec\t1.0\n"
        )
        print(version_info)

    @staticmethod
    def commandline_help():
        print (" _______________________________________________________________________________________________________________________")
        print ("|***********************************************************************************************************************|")
        print ("| * List all available commands.                                                                                        |")
        print ("| ./src/predict_binding.py                                                                                              |")
        print ("|_______________________________________________________________________________________________________________________|")
        print ("| * List all available MHC-I prediction methods.                                                                        |")
        print ("| ./src/predict_binding.py method                                                                                       |")
        print ("|_______________________________________________________________________________________________________________________|")
        print ("| * List all available (MHC,peptide_length) for a given method.                                                         |")
        print ("| ./src/predict_binding [method] mhc                                                                                    |")
        print ("| Example: ./src/predict_binding.py ann mhc                                                                             |")
        print ("|_______________________________________________________________________________________________________________________|")
        print ("| * Make predictions given a file containing a list of sequences.                                                       |")
        print ("| ./src/predict_binding [method] [mhc] [peptide_length] [input_file]                                                    |")
        print ("| Example: ./src/predict_binding.py ann HLA-A*02:01 9 ./examples/input_sequence.fasta                                   |")
        print ("|_______________________________________________________________________________________________________________________|")
        print ("| * Make predictions given a file containing a list of sequences AND user-provided MHC sequence.                        |")
        print ("| ** Only netmhcpan has this option.                                                                                    |")
        print ("| ./src/predict_binding [method] -m [input_file_mhc] [peptide_length] [input_file]                                      |")
        print ("| Example: ./src/predict_binding.py netmhcpan -m ./examples/protein_mhc_B0702.fasta 9 ./examples/input_sequence.fasta   |")
        print ("|_______________________________________________________________________________________________________________________|")
        print ("| * You may also redirect (pipe) the input file into the script.                                                        |")
        print ("| Examples:                                                                                                             |")
        print ("| echo -e ./examples/input_sequence.fasta | ./src/predict_binding.py ann HLA-A*02:01 9                                  |")
        print ("| echo -e ./examples/input_sequence.fasta | ./src/predict_binding.py netmhcpan -m ./examples/protein_mhc_B0702.fasta 9  |")
        print ("|_______________________________________________________________________________________________________________________|")


    def main(self):
        import select

        try:
            usage = "usage: %prog method allele or [options] arg length\n---\n\
Following are the available choices - \n   \
method: ann, comblib_sidney2008, consensus, IEDB_recommended, netmhcpan_ba, netmhcpan_el, smm, smmpmbec, pickpocket, netmhccons, netmhcstabpan\n   \
allele: an allele name\n   \
length: a length"

            parser = OptionParser(usage=usage, version="%prog {}".format(self.version))
            parser.add_option("-v", "--versions",
                              action="store_true",
                              dest="version_flag",
                              default=False,
                              help="print specific methods and their versions.")
            parser.add_option("-m", dest="filename_mhc",
                              help="FILE containing a single MHC sequence in fasta format.", metavar="FILE")

            (options, args) = parser.parse_args()

            if options.version_flag:
                self.print_method_versions()
                exit(0)

            if not sys.stdin.isatty():
                stdin = sys.stdin.readline().strip()
                args.append(stdin)

            args = list(filter(None, args))
            if len(args) == 0:
                self.commandline_help()
            elif (len(args) == 1) and (args[0] == 'method'):
                self.commandline_method()
            elif (len(args) == 2) and (args[1] == 'mhc'):
                self.commandline_mhc(args)
            elif len(args) == 3:
                self.commandline_input_mhc(parser, args[2])
            elif len(args) == 4:
                self.commandline_input(args)
            elif (len(args) == 5) and (args[4] == 'peptides'):
                self.commandline_input_peptide(args)
            else:
                parser.error("incorrect number of arguments")
                self.commandline_help()

        except UnexpectedInputError as e:
            print( str(e))

if __name__ == '__main__':
    Prediction().main()
