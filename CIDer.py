#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# 2016-11-02
# 2017-03-29
#
#
# Chromosome-InDel-merger (CIDer)
#
# by Sven T. Bitters
# sven.bitters@gmail.com
#
# MIT License


import os
import sys
import collections
import pandas as pd
import numpy as np
import regex as re
from multiprocessing import Pool, cpu_count
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from datetime import datetime
from traceback import print_exception
from natsort import natsorted

version_info = "0.2.5"

help_text = '\n' + "CIDer v" + version_info + """

OVERVIEW:
    Chromosome-InDel-merger (CIDer) is a Python program which introduces InDels
    and SNVs/MNVs to a DNA sequence based on information read from a VCF file.
    Additionally, CIDer updates the GFF file associated with the DNA sequence
    by reflecting the InDels introduced to the sequence on the gene annotation.
    For this, CIDer is supplied with a FASTA, a GFF, and a VCF file containing
    DNA sequences, annotations, and InDels (respectively) and returns a single
    FASTA and one GFF file which both have had the VCF file's information
    introduced to them.
    Thus, CIDer can perform a complete "lift over" procedure by itself.

ADDITIONAL INFORMATION:
    Depending on the length of the individual sequences and the number of
    InDels associated with them, one CIDer run can take several hours to
    complete. Especially the validation process takes a huge amount of time
    - which is why you can disable validation of CIDer's results by setting
    the "--noValidation" flag when you start CIDer.

    In order to speed the InDel merging process up, CIDer can run in parallel
    on multiple CPUs with each CPU handling one sequence in the FASTA file.
    However, CIDer requires quite a lot of memory for every single CPU process.
    Thus, please keep in mind the limitations of the system you are running
    CIDer on when deciding to utilize CIDer's multiprocessing capabilities.

    CIDer requires VCF files to follow VCFv4.1 specifications
    (https://samtools.github.io/hts-specs/VCFv4.1.pdf),
    GFF files to follow GFF3 specifications
    (http://www.sequenceontology.org/gff3.shtml), and FASTA files to be
    compatible with the FASTA format as defined for BioPerl
    (http://bioperl.org/formats/sequence_formats/FASTA_sequence_format).

    CIDer has been tested on Ubuntu 14.04 (trusty) with Python 3.5.2,
    biopython 1.68, pandas 0.19.0, natsort 5.0.1, and VCFv4.1 files generated
    by GATK HaplotypeCaller (https://software.broadinstitute.org/gatk).

USAGE:
    $ CIDer [--out /path/to/output] [--ex list,of,IDs] [--cpus ##]
      [--type GFFtype] </path/to/FASTA/file.fa> </path/to/GFF/file.gff>
      </path/to/VCF/file.vcf>

    "merge mode"
    $ CIDer [--h] [--l] [--out /path/to/output]
      --merge </path/to/cache/directory>

    [] = optional input
    <> = mandatory input

    You do not have to follow a defined order of inputs. However, the option
    identifiers must be followed directly by their respective values/arguments.

OPTIONS:
    --h      Show this help message and exit.

    --l      Show license information and exit.

    --out    Specify an output directory (full path). By default, the output
             directory is placed in the directory where the input FASTA file
             is located.
             Usage:  --out /full/path/to/output/directory

    --ex     Specify whether certain sequences from the FASTA file you have
             supplied as input should be excluded from the process.  A sequence
             to be excluded should be identified by its ID (see e.g. the
             GenBank format specifications for more information about sequence
             IDs). If you want to exclude multiple sequences from your task,
             you can supply a list of IDs separated by commas without any blank
             spaces between list items and commas.
             Usage:  --ex SeqID1,SeqID2,SeqID3,SeqID4

    --type   Specify an annotation type for which to filter the GFF file upon
             import. Only the specified annotation type (refers to the "type"
             column in the GFF file) will be used for further processing, all
             other types will be disregarded. The output GFF file's name will
             be appended by the type you had specified at the beginning.
             Usage:  --type CDS

    --cpus   CIDer can use multiple CPUs for your InDel merge task. If you do
             not specifiy a number of CPUs to use for parallel processesing,
             the number will default to 1. If your computer has multiple CPUs
             available, you can speed up the InDel merge process by supplying
             a value for --cpus. Keep in mind that each process also requires
             memory and that by starting too many parallel processes you can
             easily occupy substantial amounts of the memory available to your
             system.
             Usage:  --cpus 42

    --merge  CIDer will save the results of each single process in the "cache"
             directory. These cache files will be automatically merged once all
             sequences have been processed. Using the --merge option, you can
             enter "merge mode" in which these part files can be merged
             manually in order to generate a FASTA and a GFF file which then
             contain the information of all cached files.
             When running CIDer in "merge mode" you can still specify an output
             directory using --out.
             Usage:  --merge /full/path/to/directory/containing/cache/files

    --FASTAonly
             If you set the --FASTAonly flag, CIDer will only merge the FASTA
             file you put in with the VCF file's contents and will not attempt
             to update the annotations of your input genome. You do not have to
             specify a path to a GFF file when you set the --FASTAonly flag.
             Usage:  --FASTAonly

    --noValidation
             If you set the --noValidation flag, CIDer will not validate the
             annotations after they have been updated.
"""


license_text = """
MIT License

Copyright (c) 2017 Sven T. Bitters

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


def process_input(inputs, help_text, license_text):
    ## Read in input from the command line

    # Direct exit options:
    if len(inputs) == 1 or '--h' in inputs:
        print(help_text)
        sys.exit()

    if '--l' in inputs:
        print(license_text)
        sys.exit()

    # Figure out whether existing part files shall be merged
    try:
        merge_dir = inputs[inputs.index('--merge') + 1]

        if merge_dir[-1] != '/':
            merge_dir = ''.join([merge_dir, '/'])

        path_output_root = merge_dir
    except:
        merge_dir = ''

    # Find out whether CIDer shall focus on only one "type" of annotations in the GFF file
    try:
        annot_type = inputs[inputs.index('--type') + 1]
    except:
        annot_type = ''

    # Set number of CPUs to use for multiprocessing
    try:
        num_cpus = int(inputs[inputs.index('--cpus') + 1])
        if num_cpus > cpu_count():
            num_cpus = cpu_count()
    except:
        num_cpus = 1

    # Should certain chromosomes be excluded?
    try:
        exclude_id = inputs[inputs.index('--ex') + 1]
    except:
        exclude_id = list()

    # Should CIDer only introduce InDels into the DNA sequence and not update the GFF file?
    try:
        fa_only = inputs[inputs.index('--FASTAonly')]
    except:
        fa_only = False

    # Should CIDer not validate the annotations' positions?
    try:
        no_validation = inputs[inputs.index('--noValidation')]
    except:
        no_validation = False

    # Find the VCF, GFF, and FASTA file for processing
    try:
        # In order to make this process more convenient for the user, I just assume that the files will have proper
        # filename extensions by which they can be detected in the input arguments
        path_vcf = [item for item in inputs if '.vcf' in item.lower()]
        path_fasta = [item for item in inputs if '.fa' in item.lower()]

        if not fa_only:
            path_gff = [item for item in inputs if '.gff' in item.lower()]
            paths = [path_vcf, path_fasta, path_gff]
        else:
            paths = [path_vcf, path_fasta]

        full_paths = list()
        cwd = os.getcwd()
        for file_path in paths:
            if file_path:
                file_path = file_path[0]
                if file_path[0] != '/':
                    full_paths.append(''.join([cwd, '/', file_path]))
                else:
                    full_paths.append(file_path)

        for filepath in full_paths:
            if not os.path.isfile(filepath):
                print('\n > ERROR')
                print(' > ' + str(filepath))
                print(' > does not exist or is not a file!')
                print(' > Please provide a valid file as input.')
                print(' > Aborting CIDer...\n')
                sys.exit()

        if not fa_only:
            path_vcf = full_paths[0]
            path_fasta = full_paths[1]
            path_gff = full_paths[2]
        else:
            path_vcf = full_paths[0]
            path_fasta = full_paths[1]

        path_output_root = path_fasta[:(-(path_fasta[::-1].find('/')) + len(path_fasta))]

    except:
        # If the VCF, GFF, and FA file are not handed over to CIDer and "merge" mode was not selected, an error
        # is thrown and the program is ended
        if not merge_dir:
            print('\n > ERROR')
            print(' > Check your input!')
            print(' > You must at least provide a .fa, a .gff, and a .vcf file!')
            print(' > Aborting CIDer...\n')
            sys.exit()

    # Find file format version of VCF file
    ii = 0
    vcf_version = 0
    with open(path_vcf, 'r') as fio_indel_firstline:
        for line in fio_indel_firstline:
            ii += 1
            match = re.search(r'(?<=##fileformat=VCFv)[0-9].+', line)
            if match:
                vcf_version = float(match.group())
                break

            if ii == 100:
                break

    if not fa_only:
        # Find file format version of GFF file
        ii = 0
        gff_version = 0
        with open(path_gff, 'r') as fio_gff_firstline:
            for line in fio_gff_firstline:
                ii += 1
                match = re.search(r'(?<=##gff-version\ )[0-9]', line)
                if match:
                    gff_version = int(match.group())
                    break

                if ii == 100:
                    break

    # Find out whether an output directory is specified by the user
    try:
        path_output = inputs[inputs.index('--out') + 1]
        if path_output[-1] != '/':
            path_output = path_output + '/'
    except:
        # If no output directory is specified, use the directory where the FASTA file is located in for outputting
        # results
        if not merge_dir:
            path_output = path_output_root + 'CIDer_output/'
        else:
            path_output = path_output_root + 'CIDer_merge_output/'

    # Based on the output directory, specify the cache output directory
    cache_path_output = path_output + 'cache/'

    # If "merge" mode was selected initially, run "merge" mode - otherwise continue with the rest of the program
    if merge_dir:
        mergefiles(merge_dir, path_output)

    # Check wheter data from a previous CIDer run is located in the output directory
    if os.path.isfile(path_output + 'CIDer.finished') or os.path.isfile(cache_path_output + 'comments.txt'):
        print('\n > WARNING')
        print(' > You are about to overwrite data from a previous CIDer run!')
        print(' > In order to keep your data, please specify a different output location')
        print(' > using the --out option when starting CIDer.')
        print(' >')
        print(' > IF YOU CONTINUE YOUR DATA WILL BE IRRECOVERABLY LOST!')

        # This is basically my method for making sure nobody deletes "/" and all directories within by accident - i.e. I do
        # not actually delete directories but rather just the files in them
        overwrite = ''
        while overwrite not in ['Y', 'N', 'y', 'n']:
            overwrite = input(' > Really continue? [y/n] \n >> ')
        if overwrite == 'n' or overwrite == 'N':
            print('\n Aborting CIDer...\n')
            sys.exit()
        elif overwrite == 'y' or overwrite == 'Y':
            try:
                for file in os.listdir(path_output):
                    if os.path.isfile(path_output + file):
                        os.remove(path_output + file)

                for file in os.listdir(cache_path_output):
                    if os.path.isfile(cache_path_output + file):
                        os.remove(cache_path_output + file)
            except FileNotFoundError:
                pass

    # Create output directories
    try:
        if not os.path.isdir(path_output):
            os.makedirs(path_output)
        if not os.path.isdir(cache_path_output):
            os.makedirs(cache_path_output)
    except PermissionError:
        print('\n > ERROR')
        print(' > Cannot create output directory!')
        print(" > CIDer is lacking write permissions in the directory where you want to create")
        print(" > CIDer's output directory:")
        print(" > " + path_output)
        print(" > Please change r/w permissions or choose a different directory using --out.")
        print(' > Aborting CIDer...\n')
        sys.exit()

    if fa_only:
        path_gff = ''
        gff_version = ''

    input_set = (path_vcf, path_gff, path_fasta, path_output, cache_path_output, num_cpus, exclude_id, vcf_version, gff_version, annot_type, fa_only, no_validation)

    return input_set


def read_gff(path_gff, annot_type):
    ## Read in the GFF file and store it in memory in the form of a pandas data frame
    with open(path_gff, 'r') as genes_file:
        genes_comments = ''
        genes_input = ''

        # Some GFFs contain comments - these must be removed before the file's content is read into
        # a pandas data frame because otherwise the read_csv parser will not be able to handle the file
        for line in genes_file:
            if not line.startswith("#"):
                genes_input += line
            else:
                genes_comments += line

        try:
            genes_input = StringIO(genes_input)
            genes_df = pd.read_csv(genes_input, sep='\t', header=None)
            genes_df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        except:
            print('\n > CRITICAL ERROR')
            print(' > Your GFF file does not follow GFF3 file format definitions!')
            print(' > You can look up GFF3 file format definitions at')
            print(' > http://www.sequenceontology.org/gff3.shtml')
            print(' > Aborting CIDer...\n')
            sys.exit()

    # If there had been comments in the GFF file, these will be stored in a seperate file which will later be included
    # in the final InDel_GFF file
    if genes_comments:
        with open((cache_path_output + 'comments.txt'), 'w') as handle_output_comment:
            handle_output_comment.write(genes_comments)

    # Filter the GFF file for a specific annotation type (e.g. CDS or gene)
    if annot_type:
        try:
            genes_df = genes_df[genes_df["type"] == str(annot_type)]
            genes_df = genes_df.reset_index(drop=True)

            if not len(genes_df.index):
                raise ValueError

        except ValueError:
            print('\n > ERROR')
            print(' > The annotation type you specified in --type does not exist in your GFF file!')
            print(' > Aborting CIDer...\n')
            sys.exit()

    return genes_df


def make_indel_df(path_vcf, out_path):
    ## Create the InDel data frame from the VCF file, format, and clean it.

    rightnow = str(datetime.now().strftime("%Y%m%d-%H%M%S"))
    tempcopy_path = out_path + '.' + rightnow + '_vcfcopy.tmp'

    # Regexes used for finding the individual elements of the VCF file
    format_fields_regex = re.compile(r'(?<=\#\#FORMAT\=\<ID\=)[a-zA-Z0-9]+?(?=\,)')
    info_fields_regex = re.compile(r'(?<=\#\#INFO\=\<ID\=)[a-zA-Z0-9]+?(?=\,)')
    names_regex = re.compile(r'(?<=[\t;])\w+?(?==)')
    value_regex = re.compile(r'(?<==)[0-9\.-]+?(?=[,;\t])')
    format_regex = re.compile(r'(?<=\t)[a-zA-Z]{2}(:[a-zA-Z]{2,3})+(?=\t)')
    rgsm_regex = re.compile(r'(?<=\t)[a-zA-Z0-9\/_-]+(:[0-9,\.\/-]+)+')

    print(' reading VCF file...')
    # Find out how many lines the VCF file has. Also, append a newline to the VCF file and write a temporary copy of the
    # modified file to disk
    with open(path_vcf, 'r') as indels_file:
        with open(tempcopy_path, 'w') as tempcopy_file:
            for num_rows, l in enumerate(indels_file):
                tempcopy_file.write(l)
            tempcopy_file.write('\n')
            num_rows += 1

    # Go through the VCF file line by line and find all elements in the INFO column; additionally find the fields in the
    # FORMAT column. Since the INFO elements should be declared in the beginning of the file, this information should be
    # found rather quickly. However the FORMAT elements need to be retrieved from the actual VCF file contents.
    print(' retrieving format information from VCF file...')
    format_layout = []
    fullinfo_label = []
    colformat_label = ''
    with open(tempcopy_path, 'r') as indels_file:
        for line in indels_file:

            try:
                info_match = info_fields_regex.search(line)
                fullinfo_label.append(info_match.group())
            except:
                pass

            # While the fields in the FORMAT column could also be retrieved from the beginning of the VCF file, this
            # declaration does not necessarily follow the order found later in the actual VCF file rows - thus, it is
            # simpler to just find this information in an actual data row and then breaking this loop.
            try:
                match = re.search(format_regex, line)
                colformat_label = re.sub(r':', '\t', match.group())
            except:
                pass

            try:
                format_match = format_fields_regex.search(line)
                format_layout.append(format_match.group())
            except:
                pass

            if not line.startswith('##') and colformat_label and len(colformat_label.split("\t")) == len(format_layout):
                break

    colformat_dict = {colname: "." for colname in colformat_label.split("\t")}
    full_colformat_label = colformat_label.split("\t")

    colformat_label = re.sub('AD', 'nREF\tnALT', colformat_label)

    # Regexes used for identifying certain special ALT variants and for finding parts of the VCF file which will be
    # replaced with restructured information (contains \t)
    info_regex_output = re.compile(fullinfo_label[0] + "=[0-9]\S+")
    format_regex_output = re.compile(r'(?<=\t)[a-zA-Z]{2}(:[a-zA-Z]{2,3})+\t')
    altnucs_regex = re.compile(r'(?<=\t)[GATC]+(\,[GATC]+)+(?=\t)')
    altcount_regex = re.compile(r'(?<=\:)[0-9]+(\,[0-9]+){2,10}(?=\:)')

    # This dict contains all identifiers used in the INFO column
    colname_dict = {key: '.' for key in fullinfo_label}

    write_buffer = list()
    buffer_size = 750000
    with open(tempcopy_path, 'r') as indels_file:
        slash_pos = re.search(r'(?<=\/)[a-zA-Z\_0-9.]+(?=\.)', path_vcf)
        out_path = out_path + path_vcf[slash_pos.start():] + '_table'
        with open(out_path, 'w') as output_file:
            ii = 0
            for line in indels_file:

                progress_indel(ii, num_rows, 'formatting VCF file', nums=3, decs=0, inline=True)
                ii += 1

                # Replace the INFO label by the names of the elements in the INFO column in order to create an
                # individual column for each element
                if not line.startswith('##') and line.rstrip():
                    if line.startswith("#"):
                        line = line.replace('#', '', 1)
                        line = line.replace('INFO', '\t'.join(fullinfo_label))
                        line = line.replace('FORMAT', str(colformat_label))
                        line = line.replace('RGSM', '')
                        line = line.rstrip()
                        line = line + '\n'

                    else:
                        # Find ALT entries, where there is more than one alternative sequence. Delete all but the most
                        # abundant.
                        line = find_doublealt(line, altnucs_regex, altcount_regex)

                        # Find parameter names and values corresponding to FORMAT
                        colformat, colrgsm_list = find_format(line, format_regex, rgsm_regex, colformat_dict, full_colformat_label)

                        # Find parameter names and values corresponding to INFO
                        colname_list, colvalue_list = find_info(line, names_regex, value_regex)

                        # Replace the placeholders in the colname_dict with the INFO column values of the current line
                        # and extract them again from the dictionary in order to put them in the correct order in the
                        # rowvalues list.
                        rowdict = colname_dict.copy()
                        for jj in range(0, len(colvalue_list)):
                            rowdict[colname_list[jj]] = colvalue_list[jj]
                            rowvalues = [rowdict[name] for name in fullinfo_label]

                        # Replace the INFO column with the contents of the rowvalues list, delete the FORMAT column,
                        # and adjust the RGSM column in the current line
                        parameter_substitution = '\t'.join(rowvalues)
                        line = re.sub(info_regex_output, parameter_substitution, line)
                        line = re.sub(format_regex_output, '', line)
                        line = re.sub(rgsm_regex, "\t".join(colrgsm_list), line)

                    write_buffer.append(line)

                if len(write_buffer) == buffer_size:
                    write_buffer = ''.join(write_buffer)
                    output_file.write(write_buffer)
                    write_buffer = list()

            # write remaining write_buffer lines to file
            write_buffer = ''.join(write_buffer)
            output_file.write(write_buffer)

    # After the indel_df has been written to file, import it again as a data frame for further processing
    print('')
    print(' creating InDel data frame...')
    try:
        indel_df = pd.read_table(out_path, sep='\t', na_values= '.', dtype={'CHROM': 'str',
                                                                            'POS': np.int64,
                                                                            'REF': 'str',
                                                                            'ALT': 'str',
                                                                            'QUAL': np.float32})
    except:
        try:
            os.remove(out_path)
        except FileNotFoundError:
            pass

        print('\n> CRITICAL ERROR')
        print('> Data frame could not be created!')
        print('> See Traceback for more information:\n')

        try:
            exc_info = sys.exc_info()
            try:
                print_exception(*exc_info)
            except:
                pass
            del exc_info
        except:
            pass

        print('\n> Aborting CIDer...\n')
        sys.exit()

    try:
        os.remove(tempcopy_path)
    except FileNotFoundError:
        pass

    return indel_df


def progress_indel(dividend, divisor, text, nums=3, decs=0, inline=False):
    ## Display the progress of a process

    percent = abs(round((dividend / divisor) * 100, decs))
    if percent != round(((dividend-1) / divisor) * 100, decs):
        if inline:
            format_str = '\r {} - progress: {:>' + str(nums) + '.' + str(decs) + 'f} % '
            print(format_str.format(text, percent), end='')


def find_doublealt(line, altnucs_regex, altcount_regex):
    # It might happen that there are two or more alternative ALT sequences in the VCF file for one entry - find the one
    # that has more reads associated to it and delete the other one in order to get one REF and one ALT sequence.
    altnucs_match = re.search(altnucs_regex, line)

    try:
        altnucs = altnucs_match.group().split(',')
        count_match = re.search(altcount_regex, line)
        count = count_match.group().split(',')

        preferred_count = [count[0]]
        count.pop(0)
        preferred_count.append(max(count))
        preferred_alt = altnucs[count.index(max(count))]

        line = re.sub(altnucs_regex, preferred_alt, line)
        line = re.sub(altcount_regex, ','.join(preferred_count), line)
    except AttributeError:
        # An AttributeError will be raised when there are not multiple sequences in ALT. Just pass it.
        pass

    return line


def find_format(line, format_regex, rgsm_regex, colformat_dict, full_colformat_label):
    try:
        rgsm_dict = colformat_dict.copy()

        # Find parameter names corresponding to FORMAT
        match = re.search(format_regex, line)
        cols = match.group()
        cols = cols.split(":")
        colformat = re.sub(r':', '\t', match.group())
        colformat = re.sub('AD', 'nREF\tnALT', colformat)

        # Find parameter values corresponding to FORMAT
        colrgsm_list = list()
        for match in rgsm_regex.finditer(line):
            colvals = match.group().split(":")

            vallistcount = 0
            for colname in cols:
                rgsm_dict[colname] = colvals[vallistcount]
                vallistcount += 1

            if rgsm_dict["AD"] != ".":
                rgsm_dict["AD"] = re.sub(r',', '\t', rgsm_dict["AD"])
            else:
                rgsm_dict["AD"] = ".\t."

            for key in full_colformat_label:
                colrgsm_list.append(rgsm_dict[key])

        return colformat, colrgsm_list

    except:
        pass


def find_info(line, names_regex, value_regex):
    try:
        # Find parameter names in INFO
        colname = [match.group() for match in names_regex.finditer(line)]

        # Find parameter values in INFO
        colvalue = [match.group() for match in value_regex.finditer(line)]

        return colname, colvalue

    except:
        pass


def make_new_genome(chrom_num_id, description, orig_DNA_sequence, chromosome_indel, chromosome_gff,
                      process_id, cache_output, fa_only, no_validation):
    ## THIS IS WHERE THE MAGIC HAPPENS
    ## This procedure is the starting platform from which the introduction of InDels to the DNA sequence and the
    ## update of the GFF file are launched

    print(' ' + chrom_num_id + ' - process initialized')

    no_gff = False
    validation_result = False
    deletedsize = 0
    val_prob_no = 0
    serious_errors = 0

    try:
        general_path = cache_output + str(process_id) + '_' + chrom_num_id

        fio_error_path = general_path + '.error'
        fio_error = open(fio_error_path, 'w')

        fio_run_path = general_path + '.running'
        fio_run = open(fio_run_path, 'w')

        write_text_to_file(fio_run, 'process initialized')

        if len(chromosome_indel) > 0:
            # If there are InDels associated with the current DNA sequence proceed here
            write_text_to_file(fio_run, 'removal of bad VCF file entries in progress...')

            clean_chromosome_indel, dropped_indel_df, cache_indel_path, deletedsize = clean_indels(chromosome_indel, fio_run, general_path)

            if len(clean_chromosome_indel) > 0:
                # If there are InDels left in the VCF file after the cleanup, proceed here
                write_text_to_file(fio_run, 'DNA sequence merging in progress...')
                new_DNA_sequence, gff_change_list, serious_errors = update_fa(orig_DNA_sequence, clean_chromosome_indel, fio_run, fio_error, serious_errors)
                write_text_to_file(fio_run, 'DNA sequence merging completed')
                print(' ' + chrom_num_id + ' - DNA completed')

                if len(chromosome_gff) > 0 and not fa_only:
                    write_text_to_file(fio_run, 'GFF update in progress...')
                    dropped_indel_df = get_overlaps(dropped_indel_df)

                    if len(dropped_indel_df) > 0:
                        chromosome_gff = pd.concat([chromosome_gff, dropped_indel_df], ignore_index=True)

                    chromosome_gff = chromosome_gff.reset_index(drop=True)
                    orig_chromosome_gff = chromosome_gff.copy()
                    new_chromosome_gff, serious_errors = update_gff(chromosome_gff, gff_change_list, new_DNA_sequence, fio_run, fio_error, serious_errors)
                    write_text_to_file(fio_run, 'GFF update completed')
                    print(' ' + chrom_num_id + ' - GFF completed')

                    if not no_validation:
                        write_text_to_file(fio_run, 'validation of updated GFF in progress...')

                        validation_result, val_prob_no = validate(orig_chromosome_gff, orig_DNA_sequence, new_chromosome_gff, new_DNA_sequence, clean_chromosome_indel, general_path, fio_run)

                        if validation_result is 'fail':
                            write_text_to_file(fio_run, 'validation completed but there have been ' + str(val_prob_no) + ' disagreements')
                            print(' ' + chrom_num_id + ' - validation disagreement')

                        else:
                            write_text_to_file(fio_run, 'validation completed without disagreements')
                            print(' ' + chrom_num_id + ' - validation completed')
                    else:
                        validation_result = 'noValidation'
                        val_prob_no = -1

                else:
                    write_text_to_file(fio_run, 'there is no GFF file associated with this sequence')
                    no_gff = True
                    validation_result = 'skipped'
                    val_prob_no = -1

                if serious_errors > 0:
                    write_text_to_file(fio_run, 'InDel-merging completed but there have been ' + str(serious_errors) + ' serious errors')
                    print(' ' + chrom_num_id + ' - serious errors')
                    serious_errors = 1

            else:
                # Continue here if all InDels were removed during cleanup
                print(' ' + chrom_num_id + ' - all InDels were bad')
                write_text_to_file(fio_run, 'all InDels were removed during cleanup of the VCF file entries - no further action needs to be taken')
                os.rename(fio_run_path, general_path + '.noInDels')
                fio_run_path = general_path + '.noInDels'

                new_DNA_sequence = orig_DNA_sequence

                if not fa_only:
                    new_chromosome_gff = chromosome_gff
                else:
                    no_gff = True

        else:
            # Continue here if no InDels are associated with the current DNA sequence
            print(' ' + chrom_num_id + ' - no InDels listed')
            write_text_to_file(fio_run, 'there are no InDels associated with this DNA sequence')
            os.rename(fio_run_path, general_path + '.noInDels')
            fio_run_path = general_path + '.noInDels'

            new_DNA_sequence = orig_DNA_sequence

            if not fa_only:
                new_chromosome_gff = chromosome_gff
            else:
                no_gff = True

        # Write FASTA of this chromosome to file
        write_text_to_file(fio_run, 'writing FASTA')
        with open((general_path + '.fa'), 'w') as handle_output_fa:
            SeqIO.write(SeqRecord(new_DNA_sequence, id=description, name='', description=''), handle_output_fa, 'fasta')

        if not no_gff:
            # Write GFF of this chromosome to file
            write_text_to_file(fio_run, 'writing GFF')
            new_chromosome_gff = new_chromosome_gff[new_chromosome_gff["type"] != "!row_del!"]

            with open((general_path + '.gff'), 'w') as handle_output_gff:
                new_chromosome_gff.to_csv(handle_output_gff, sep='\t', header=False, index=False)

        write_text_to_file(fio_run, 'finished process')
        print(' --> ', chrom_num_id, ' - completed')

        try:
            fio_run.close()
            os.rename(fio_run_path, general_path + '.finished')
        except FileNotFoundError:
            pass

        try:
            fio_error.close()
            if os.path.getsize(fio_error_path) == 0.00:
                os.remove(fio_error_path)
        except FileNotFoundError:
            pass

        return chrom_num_id, 'success', validation_result, deletedsize, val_prob_no, serious_errors

    except:
        error_msg = 'CRITICAL ERROR\n' \
                    'unexpected error in {} - shutting down this process\n' \
                    'other processes continue, however.\n'.format(
                        str(chrom_num_id))
        print(error_msg)

        try:
            try:
                write_text_to_file(fio_error, error_msg)
            except:
                pass
            try:
                fio_run.close()
            except:
                pass

            try:
                exc_info = sys.exc_info()
                try:
                    fio_error.write('See traceback for more information:\n')
                    print_exception(*exc_info, file=fio_error)
                except:
                    pass
                del exc_info
            except:
                pass

            try:
                fio_error.close()
            except:
                pass

            return chrom_num_id, 'error', validation_result, deletedsize, val_prob_no, serious_errors

        except:
            return chrom_num_id, 'error', validation_result, deletedsize, val_prob_no, serious_errors


def clean_indels(vcf_df, fio_run, cache_path):
    ## Elaborate filtering and cleaning of the InDel data frame.
    df_orig_size = vcf_df.shape[0]

    # Certain entries in the InDel table can be disregarded
    # If Reference and Allele sequence are the same, then this InDel table entry contains no useful information and can
    # be disregarded without any worries.
    vcf_df = vcf_df[vcf_df["REF"] != vcf_df["ALT"]]

    # If there is no ALT sequence, there is nothing to do for CIDer
    vcf_df = vcf_df[(vcf_df["ALT"] != '-') & (vcf_df["ALT"] != '.')]

    # If the Phred quality score normalized against the read depth - QD - of an InDel is too low, it is unreliable and
    # should be disregarded. If read depth DP is too low, the InDel could be unreliable - too high DP are probably
    # mapping artefacts, though. High absolute scores for ReadPosRankSum indicate that there is a bias towards the end
    # of the read to either call REF or ALT. If the reference sequence does not have a much higher frequency in the
    # genome assembly than the alternative sequence, just use the alternative sequence - otherwise keep the reference
    # sequence.
    # This choice was made since our plants are pure-bred and all alleles that are relevant for our
    # phenotype should be homozygous - thus, the relevance of heterozygous alleles for our research is
    # uncertain.
    vcf_df = vcf_df[(((vcf_df["QD"] >= 20) & ((vcf_df["DP"] >= 5) & (vcf_df["DP"] <= 200)) & (vcf_df["MQ"] >= 40) & (~(vcf_df["nREF"] > vcf_df["nALT"] * 3))) |
                        ((vcf_df["AF"] == 1) & (vcf_df["QD"] >= 20) & ((vcf_df["DP"] >= 5) & (vcf_df["DP"] <= 200)) & (vcf_df["MQ"] >= 40)))]
    vcf_df = vcf_df.reset_index(drop=True)

    rows_delete = list()
    for index, row in vcf_df.iterrows():
        # Note: In contrast to the filtering above where values that will be IN the data frame were selected, here
        # values that will NOT be in the data frame are selected.

        if row["ReadPosRankSum"] != np.nan:
            if abs(float(row["ReadPosRankSum"])) > 2:
                rows_delete.append(index)

        if row["FS"] != ".":
            if row["FS"] > 60:
                rows_delete.append(index)

    rows_delete = list(set(rows_delete))

    vcf_df = vcf_df.drop(vcf_df.index[rows_delete])
    vcf_df = vcf_df.reset_index(drop=True)

    df_bad_params = df_orig_size - vcf_df.shape[0]
    percent_bad = round((df_bad_params / df_orig_size)*100, 1)
    write_text_to_file(fio_run, 'removed ' + str(df_bad_params) + ' of ' + str(df_orig_size) + ' entries (' + str(percent_bad) + ' %) due to low quality')

    try:
        del vcf_df["index"]
    except:
        pass

    # This section makes sure that overlapping InDels are removed from the data frame because it is difficult to decide
    # which overlapping InDel is more relevant for the genome sequence than the other(s). Since we are working with
    # pure-bred plants which should be homozygous for all alleles relevant for the phenotype, the relevance of
    # overlapping heterozygous alleles cannot be easily determined. Thus, by default, the reference sequence will be
    # kept in the new genome sequence.
    delete_rows = list()
    index_dict = collections.OrderedDict()
    dict_maxsize = 2500
    for index, row in vcf_df.iterrows():
        for ii in range(int(row.POS), (int(row.POS) + len(row.REF))):

            index_dict[ii] = index_dict.get(ii, []) + [index]

            if len(index_dict) == dict_maxsize:
                sorted_dict = collections.OrderedDict(sorted(index_dict.items(), key=lambda x: len(x[1]), reverse=True))
                for value in sorted_dict.values():
                    if len(value) > 1:
                        delete_rows += value
                    else:
                        break
                while len(index_dict) != dict_maxsize/2:
                    index_dict.popitem(last=False)

    for value in index_dict.values():
        if len(value) > 1:
            delete_rows += value

    # Get unique data frame indexes and drop these rows
    delete_rows = list(set(delete_rows))
    dropped_vcf_df = vcf_df.copy()
    dropped_vcf_df = dropped_vcf_df.loc[delete_rows]
    vcf_df = vcf_df.drop(vcf_df.index[delete_rows])
    vcf_df = vcf_df.reset_index(drop=True)

    df_overlapping = (df_orig_size - df_bad_params) - vcf_df.shape[0]
    percent_overlapping = round((df_overlapping / df_orig_size)*100, 1)
    write_text_to_file(fio_run, 'removed ' + str(df_overlapping) + ' of ' + str(df_orig_size) + ' entries (' + str(percent_overlapping) + ' %) due to overlapping')

    # Reverse the order of the vcf_df. This is relevant for the updating of the GFF file because in a reversed
    # vcf_df there is no index shifting which could impair the quality of the updated GFF.
    vcf_df = vcf_df.reindex(index=vcf_df.index[::-1])

    percent_remaining = round((vcf_df.shape[0] / df_orig_size)*100, 1)
    write_text_to_file(fio_run, 'final size of the cleaned VCF file is ' + str(vcf_df.shape[0]) + ' entries (' + str(percent_remaining) + ' %)')

    cache_indel_path = cache_path + ".indel"
    with open(cache_indel_path, "w") as fio_indel:
        vcf_df.to_csv(fio_indel, sep="\t")

    return vcf_df, dropped_vcf_df, cache_indel_path, (df_overlapping + df_bad_params)


def update_fa(sequence, chromosome_indel, fio_run, fio_error, serious_errors):
    ## This procedure introduces the InDels from the VCF file to the DNA sequence and at the same time prepares the
    ## GFF file update.

    # This dictionary will be filled with key: value pairs where keys are index positions in the FASTA file and
    # values are the InDels corresponding to these positions.
    chrom_seq_dict = dict()

    # For speed reasons, every sequence is handled as a list of nucleotides until it has been fully processed.
    chrom_seq_list = list(sequence)

    len_chromosome_indel = len(chromosome_indel)

    aligndiff_regex = re.compile(r'.-+')

    gff_change_list = list()
    previous_percent_progress = 0
    aa = 1
    # Go through the InDel table line by line from start to end
    for indel_index, indel_line in chromosome_indel.iterrows():
        percent_progress = round((aa / len_chromosome_indel) * 100, 1)
        if percent_progress != previous_percent_progress and percent_progress % 2.5 == 0:
            write_text_to_file(fio_run, 'progress - {:>5.1f} %'.format(percent_progress))
            previous_percent_progress = percent_progress
        aa += 1

        ref_seq = indel_line["REF"]
        alt_seq = indel_line["ALT"]
        pos_indel = indel_line["POS"]

        # General notice:
        # Since Python lists/strings start at index position 0 and FASTA sequences start with index position 1,
        # all indexes in the InDel table have to be shifted by 1 to the left (-1) in order to match the DNA
        # sequence which uses Python notation.
        index_1 = int(pos_indel) - 1
        index_2 = index_1 + len(ref_seq)

        sequence_in_chrom = ''.join(chrom_seq_list[index_1: index_2])

        try:
            if ref_seq.isalpha() and alt_seq.isalpha():

                if sequence_in_chrom.upper() == ref_seq.upper():

                    # SNPs (or MNPs) have same length REFs and ALTs
                    if len(ref_seq) == len(alt_seq):
                        if len(ref_seq) == 1:
                            # Handle SNPs
                            chrom_seq_dict[index_1] = alt_seq
                            gff_change_list.append([index_1, 0, 1])
                        else:
                            # Handle MNPs by basically splitting them in SNPs
                            jj = 0
                            for pos in range(index_1, index_2):
                                chrom_seq_dict[pos] = alt_seq[jj]
                                jj += 1
                                gff_change_list.append([pos, 0, 1])

                    # Insertions have smaller REFs than ALTs
                    elif len(ref_seq) < len(alt_seq):
                        if len(ref_seq) == 1:
                            chrom_seq_dict[index_1] = alt_seq
                            lendiff = len(alt_seq) - len(ref_seq)
                            gff_change_list.append([index_1, lendiff, 2])

                        else:
                            # Handle more complex insertions by reducing longer insertion stretches to simple insertions
                            # and converting SNPs contained in the insertion to simple SNPs.
                            # This process is based on creating gap-less alignments of REF and ALT and then finding
                            # simple insertions in the complex ALT "statement"
                            alt_ref_alignment = pairwise2.align.localms(alt_seq, ref_seq, 5, -4, -3, -.1)
                            alt_align_res = alt_ref_alignment[0][0]
                            ref_align_res = alt_ref_alignment[0][1]

                            while ref_align_res.startswith("-"):
                                add_seq = sequence[index_1 - 1]
                                alt_seq = add_seq + alt_seq
                                ref_seq = add_seq + ref_seq
                                index_1 -= 1

                                alt_ref_alignment = pairwise2.align.localms(alt_seq, ref_seq, 5, -4, -3, -.1)
                                alt_align_res = alt_ref_alignment[0][0]
                                ref_align_res = alt_ref_alignment[0][1]

                            ref_diff = re.finditer(aligndiff_regex, ref_align_res)

                            diff_pos = list()
                            for diff in ref_diff:
                                subinsert = alt_align_res[diff.start():diff.end()]

                                until_subinsert = ref_align_res[0:diff.start()]
                                shift = len(until_subinsert) - until_subinsert.count("-")
                                chrom_seq_dict[(index_1 + shift)] = subinsert

                                orig_element = re.sub("-", "", diff.group())
                                diff_length = len(subinsert) - len(orig_element)
                                gff_change_list.append([index_1 + shift, diff_length, 2])

                                diff_pos = diff_pos + list(range(diff.start(), diff.end()))

                            for nucleotide_pos in range(0, len(ref_align_res)):
                                ref_nuc = ref_align_res[nucleotide_pos]
                                alt_nuc = alt_align_res[nucleotide_pos]

                                if ref_nuc != alt_nuc:
                                    # Handle SNPs contained in complex insertions
                                    if not nucleotide_pos in diff_pos:
                                        until_diff = ref_align_res[0:nucleotide_pos]
                                        shift = len(until_diff) - until_diff.count("-")
                                        chrom_seq_dict[(index_1 + shift)] = alt_nuc

                                        gff_change_list.append([(index_1 + shift), 0, 1])

                    # Deletions have longer REFs than ALTs
                    elif len(ref_seq) > len(alt_seq):
                        if len(alt_seq) == 1:
                            chrom_seq_dict[index_1] = alt_seq
                            # For deletions, chrom_seq_dict will be filled with zero-length strings on every deleted
                            # index position
                            for ii in range(index_1 + 1, index_2):
                                chrom_seq_dict[ii] = ''
                            lendiff = len(alt_seq) - len(ref_seq)
                            gff_change_list.append([index_1, lendiff, 3])
                        else:
                            # Handle more complex deletions by splitting them in single nucleotide deletions
                            alt_ref_alignment = pairwise2.align.localms(alt_seq, ref_seq, 5, -4, -3, -.1)
                            alt_align_res = alt_ref_alignment[0][0]
                            ref_align_res = alt_ref_alignment[0][1]

                            while alt_align_res.startswith("-") or ref_align_res.startswith("-"):
                                add_seq = sequence[index_1 - 1]
                                alt_seq = add_seq + alt_seq
                                ref_seq = add_seq + ref_seq
                                index_1 -= 1

                                alt_ref_alignment = pairwise2.align.localms(alt_seq, ref_seq, 2, -1, -3, -.1)
                                alt_align_res = alt_ref_alignment[0][0]
                                ref_align_res = alt_ref_alignment[0][1]

                            ii = 0
                            for nt in alt_align_res:
                                print(nt)
                                if nt == '-':
                                    chrom_seq_dict[(index_1 + ii)] = ''
                                else:
                                    chrom_seq_dict[(index_1 + ii)] = nt
                                ii += 1

                            alt_diff = re.finditer(aligndiff_regex, alt_align_res)

                            diff_pos = list()
                            for diff in alt_diff:
                                subdeletion = alt_align_res[diff.start():diff.end()]
                                until_subdeletion = alt_align_res[0:diff.start()]

                                shift = len(until_subdeletion)

                                """
                                ii = diff.start()
                                for nt in subdeletion:
                                    print(nt)
                                    if nt == '-':
                                        chrom_seq_dict[(index_1 + ii)] = ''
                                    else:
                                        chrom_seq_dict[(index_1 + ii)] = nt
                                    ii += 1

                                print(chrom_seq_dict)
                                """
                                orig_element = re.sub("-", "", diff.group())
                                diff_length = len(orig_element) - len(subdeletion)
                                gff_change_list.append([index_1 + shift, diff_length, 3])

                            """
                                diff_pos = diff_pos + list(range(diff.start(), diff.end()))

                            for nucleotide_pos in range(0, len(ref_align_res)):
                                ref_nuc = ref_align_res[nucleotide_pos]
                                alt_nuc = alt_align_res[nucleotide_pos]

                                if ref_nuc != alt_nuc:
                                    if not nucleotide_pos in diff_pos:
                                        if alt_nuc == "-":
                                            chrom_seq_dict[(index_1 + nucleotide_pos)] = ''
                                            gff_change_list.append([(index_1 + nucleotide_pos), -1, 3])
                                        else:
                                            chrom_seq_dict[(index_1 + nucleotide_pos)] = alt_nuc
                            """

                    else:
                        error_msg = 'SERIOUS ERROR\n' \
                                    'Could not handle\n{}\n' \
                                    'For support, please contact me via sven.bitters@gmail.com.\n\n'.format(indel_line)
                        write_text_to_file(fio_error, error_msg)
                        serious_errors += 1

                elif sequence_in_chrom.upper() != ref_seq.upper():
                    error_msg = 'WARNING\n' \
                                'Sequence and VCF file do not match starting from position {}!\n' \
                                'REF:    {}\nFound:  {}\nALT:    {}\n' \
                                'This VCF file entry will be disregarded!\n\n'.format(
                        (index_1 + 1), ref_seq, sequence_in_chrom, alt_seq)
                    write_text_to_file(fio_error, error_msg)

            else:
                error_msg = 'SERIOUS ERROR\n' \
                            'Could not handle\n{}\n' \
                            'For support, please contact me via sven.bitters@gmail.com.\n\n'.format(indel_line)
                write_text_to_file(fio_error, error_msg)
                serious_errors += 1

        except TypeError:
            error_msg = 'ERROR\n' \
                        'Something seems to be wrong with your VCF file in line {}!\n' \
                        '{}\n' \
                        'This VCF file entry will be disregarded!' \
                        'For support, please contact me via sven.bitters@gmail.com.\n\n'.format(
                indel_index + 1, indel_line)
            write_text_to_file(fio_error, error_msg)

    # Take the sequence list and replace the list items at the positions where an InDel or S/MNP occured with the
    # corresponding new sequence. Thus, for every key in the chrom_seq_dict put the corresponding value in
    # chrom_seq_list.
    chrom_seq_list = list(sequence)
    for key in chrom_seq_dict.keys():
        chrom_seq_list[key] = chrom_seq_dict[key]
    chrom_seq_str = ''.join(chrom_seq_list)
    chrom_seq = Seq(chrom_seq_str)

    return chrom_seq, gff_change_list, serious_errors


def update_gff(gff_df, input_listoflists ,new_DNA_sequence, fio_run, fio_error, serious_errors):
    ## This procedure updates the GFF file according to the changes made to the sequence in make_new_genome / update_fa

    gff_len = len(gff_df)
    previous_percent_progress = 0
    counter = 1

    # Convert to Python index because the elements of input_listoflists are in Python indexes
    gff_df["start"] = gff_df["start"] - 1
    gff_df["end"] = gff_df["end"] - 1

    # gff_dict basically contains all GFF file lines in the form {index: line}. Storing the GFF file in a dictionary
    # is much quicker than looping through the pandas data frame
    gff_dict = {index: [str(content.seqid), str(content.source), str(content.type), int(content.start), int(content.end), str(content.score), str(content.strand), str(content.phase), str(content.attributes)] for index, content in gff_df.iterrows()}

    # Loop through every entry in the dictionary
    for gff_index in range(0, len(gff_df.index)):
        percent_progress = round((counter / gff_len) * 100, 1)
        if percent_progress != previous_percent_progress and percent_progress % 2.5 == 0:
            write_text_to_file(fio_run, 'progress - {:>5.1f} %'.format(percent_progress))
            previous_percent_progress = percent_progress
        counter += 1

        # There exists some problem that the region annotation is longer than the DNA sequence in the end
        # While this is not a proper solution, I use this bypass in order to get something I can work with for now!
        #if gff_dict[gff_index][2] == "region":
        #    gff_dict[gff_index][4] = len(new_DNA_sequence) - 1
        #    continue

        annotation_start = int(gff_dict[gff_index][3])
        annotation_end = int(gff_dict[gff_index][4])

        # First, handle all indels before annotation_start
        indels_before = [indel for indel in input_listoflists if (int(indel[0]) < annotation_start and int(indel[2]) != 3) or (int(indel[0]) + abs(int(indel[1])) < annotation_start and indel[2] == 3)]

        before_lendiff = sum([lendiff[1] for lendiff in indels_before])

        indels_inside = [indel for indel in input_listoflists if (annotation_start <= int(indel[0]) <= annotation_end and int(indel[2]) != 3) or (annotation_start <= int(indel[0]) + abs(int(indel[1])) and int(indel[0]) <= annotation_end and int(indel[2]) == 3) or ((int(indel[0]) <= annotation_start and annotation_end <= int(indel[0]) + abs(int(indel[1]))) and int(indel[2]) == 3)]
        indels_inside = sorted(indels_inside, reverse=True)

        if indels_inside:
            for sublist in indels_inside:
                # indexpos is a Python index!
                indexpos = int(sublist[0])
                # lendiff = len(ALT) - len(REF) --> lendiff < 0 for deletions, lendiff > 0 for insertions
                # thus, indexpos + abs(lendiff) --> endpos of the InDel
                lendiff = int(sublist[1])
                endpos = indexpos + abs(lendiff)
                # SNP        -->  indel_flag = 1
                # insertion  -->  indel_flag = 2
                # deletion   -->  indel_flag = 3
                indel_flag = int(sublist[2])

                if gff_dict[gff_index][2] == "!row_del":
                    continue

                change_att = False
                row_del = False

                # S/MNVs
                if indel_flag == 1:
                    if indexpos in range(int(gff_dict[gff_index][3]), int(gff_dict[gff_index][4]) + 1):
                        change_att = True

                    #if indexpos in range(int(gff_dict[gff_index][3]), int(gff_dict[gff_index][3]) + 3):
                    #    change_att = True

                    #elif indexpos in range(int(gff_dict[gff_index][4]) - 2, int(gff_dict[gff_index][4]) + 1):
                    #    change_att = True

                    else:
                        error_msg = 'SERIOUS ERROR\n' \
                                    'InDel data in {} could not be handled properly.\n' \
                                    'This is likely very bad!\n' \
                                    'Annotation_start: {}\n' \
                                    'Annotation_end:   {}\n' \
                                    'InDels_inside:    {}\n' \
                                    'Annotation:       {}\n' \
                                    'For support, please contact me via sven.bitters@gmail.com.\n\n'.format(str(sublist), str(gff_dict[gff_index][3]), str(gff_dict[gff_index][4]), str(indels_inside), gff_dict[gff_index])
                        write_text_to_file(fio_error, error_msg)
                        serious_errors += 1

                # Insertions
                elif indel_flag == 2:
                    if (int(gff_dict[gff_index][3])) <= indexpos < int(gff_dict[gff_index][4]):
                        # Insertion directly at the start or somewhere in the middle of an annotated element:
                        gff_dict[gff_index][4] = gff_dict[gff_index][4] + lendiff
                        change_att = True

                    elif indexpos == (int(gff_dict[gff_index][4])):
                        # Insertion on the last nucleotide of the annotated element:
                        if gff_dict[gff_index][2] == "region":
                            gff_dict[gff_index][4] = gff_dict[gff_index][4] + lendiff
                        else:
                            pass

                    else:
                        error_msg = 'SERIOUS ERROR\n' \
                                    'InDel data in {} could not be handled properly.\n' \
                                    'This is likely very bad!\n' \
                                    'Annotation_start: {}\n' \
                                    'Annotation_end:   {}\n' \
                                    'InDels_inside:    {}\n' \
                                    'Annotation:       {}\n' \
                                    'For support, please contact me via sven.bitters@gmail.com.\n\n'.format(str(sublist), str(gff_dict[gff_index][3]), str(gff_dict[gff_index][4]), str(indels_inside), gff_dict[gff_index])
                        write_text_to_file(fio_error, error_msg)
                        serious_errors += 1

                # Deletions
                elif indel_flag == 3:

                    if indexpos < int(gff_dict[gff_index][3]) and int(gff_dict[gff_index][3]) <= endpos < int(gff_dict[gff_index][4]):
                        # If a deletion only partially overlaps with the start of an annotation:
                        gff_dict[gff_index][3] = indexpos + 1
                        gff_dict[gff_index][4] = gff_dict[gff_index][4] + lendiff
                        change_att = True

                    elif indexpos == int(gff_dict[gff_index][3]) and int(gff_dict[gff_index][3]) < endpos < int(gff_dict[gff_index][4]):
                        # If a deletion coincides with the start of an annotation:
                        gff_dict[gff_index][4] = gff_dict[gff_index][4] + lendiff
                        change_att = True

                    elif int(gff_dict[gff_index][3]) < indexpos and endpos == int(gff_dict[gff_index][4]):
                        # If a deletion coincides with the end of an annotation:
                        gff_dict[gff_index][4] = indexpos
                        change_att = True

                    elif int(gff_dict[gff_index][3]) < indexpos < int(gff_dict[gff_index][4]) < endpos:
                        # If a deletion only partially overlaps with the end of an annotation:
                        gff_dict[gff_index][4] = indexpos
                        change_att = True

                    elif indexpos <= int(gff_dict[gff_index][3]) and int(gff_dict[gff_index][4]) <= endpos:
                        # If a deletion stretches accross a whole annotated element, remove that element
                        gff_dict[gff_index][2] = "!row_del!"
                        change_att = False
                        row_del = True

                    elif indexpos == int(gff_dict[gff_index][4]):
                        # If the start of a deletion coincides with the end of an annotation
                        gff_dict[gff_index][4] = gff_dict[gff_index][4]
                        change_att = True

                    elif (int(gff_dict[gff_index][3])) < indexpos and endpos < (int(gff_dict[gff_index][4])):
                        # If an InDel is located inside an annotation:
                        gff_dict[gff_index][4] = gff_dict[gff_index][4] + lendiff
                        change_att = True

                    else:
                        error_msg = 'SERIOUS ERROR\n' \
                                    'InDel data in {} could not be handled properly.\n' \
                                    'This is likely very bad!\n' \
                                    'Annotation_start: {}\n' \
                                    'Annotation_end:   {}\n' \
                                    'InDels_inside:    {}\n' \
                                    'Annotation:       {}\n' \
                                    'For support, please contact me via sven.bitters@gmail.com.\n\n'.format(str(sublist), str(gff_dict[gff_index][3]), str(gff_dict[gff_index][4]), str(indels_inside), gff_dict[gff_index])
                        write_text_to_file(fio_error, error_msg)
                        serious_errors += 1

                if not row_del and int(gff_dict[gff_index][4]) - int(gff_dict[gff_index][3]) < 1:
                    # It might happen that some InDels slowly decrease the size of an annotated element until there is
                    # nothing left of it. In this case, this element should be deleted.
                    gff_dict[gff_index][2] = "!row_del!"
                    change_att = False

                if change_att:
                    # If an InDel is inside an annotation, mark the annotation by adding "_![x]" to its 'Name' or 'ID'
                    #orig_attribute = gff_dict[gff_index][8]

                    #if "Name" in orig_attribute:
                    #    att_pattern = "Name"
                    #elif "ID" in orig_attribute:
                    #    att_pattern = "ID"
                    #else:
                    #    continue

                    #attribute_regex = re.compile(att_pattern + r'=.+?((?=$)|(?=\;))')
                    #attribute_match = re.search(attribute_regex, orig_attribute)
                    #attribute_group = attribute_match.group()

                    #if change_att not in attribute_group:
                    #    new_attribute_part = attribute_group + change_att

                        # Replace the gene's attribute with the updated version containing "_![x]"
                     #   mod_attribute = re.sub(attribute_regex, new_attribute_part, orig_attribute)
                     #   gff_dict[gff_index][8] = mod_attribute
                    gff_dict[gff_index][1] = "CIDer"

        gff_dict[gff_index][3] = gff_dict[gff_index][3] + before_lendiff
        gff_dict[gff_index][4] = gff_dict[gff_index][4] + before_lendiff

    # Create a new data frame from the gff_dict
    gff_df = gff_df.from_dict(gff_dict, orient='index')
    gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    gff_df["start"] = gff_df["start"] + 1
    gff_df["end"] = gff_df["end"] + 1

    # Return the GFF file
    return gff_df, serious_errors


def validate(orig_chromosome_gff, orig_chromosome_seq, new_chromosome_gff, new_chromosome_seq, indel_df, fio_path, fio_run):
    ## Validate the GFF file by checking whether the DNA sequence in the beginning of each annotation is the same in the
    ## original and the updated version of the genome. If InDels are located in this sequence, the original and the
    ## updated sequence for the annotations will not match - thus, the validate procedure will use the information in
    ## the InDel table to update the original DNA sequence in order to recapitulate the initial update_fa process which
    ## resulted in the updated sequence. Then, this "new" original sequence will be compared again with the updated
    ## sequence. If the make_new_genome process worked as expected, there should not remain any disagreements between
    ## original and updated sequence per annotation.

    #try:
    fio_report_path = fio_path + '.validation_report'
    fio_report = open(fio_report_path, 'w')

    slash_pos = re.search(r'(?<=\/)[0-9]+\_[a-zA-Z\_0-9.]+(?=$)', fio_path)
    fio_dummy_path = fio_path[:slash_pos.start()] + '.' + fio_path[slash_pos.start():] + '.null'
    fio_dummy = open(fio_dummy_path, 'w')

    write_text_to_file(fio_report, 'starting validation of annotations')

    comparison_length = 99
    validation_problem_no = 0
    serious_errors = 0

    # GFF file entries with the 'variant' type will not be evaluated
    orig_chromosome_gff = orig_chromosome_gff[orig_chromosome_gff["type"] != 'variant']
    new_chromosome_gff = new_chromosome_gff[new_chromosome_gff["type"] != 'variant']

    orig_seq_list = list()
    new_seq_list = list()

    # Store all annotations from the original and the updated GFF set in lists
    for index, orig_gff_line in orig_chromosome_gff.iterrows():
        orig_seq_list.append([orig_gff_line, orig_chromosome_seq[(int(orig_gff_line.start) - 1): (int(orig_gff_line.start) + comparison_length)]])

    for index, new_gff_line in new_chromosome_gff.iterrows():
        new_seq_list.append([new_gff_line, new_chromosome_seq[(int(new_gff_line.start) - 1): (int(new_gff_line.start) + comparison_length)]])

    orig_seq_list_len = len(orig_seq_list)
    previous_percent_progress = 0
    counter = 1

    # Loop through these lists element by element, i.e. one annotation after the other
    for ii in range(0, len(orig_chromosome_gff.index)):

        percent_progress = round((counter / orig_seq_list_len) * 100, 1)
        if percent_progress != previous_percent_progress and percent_progress % 2.5 == 0:
            write_text_to_file(fio_run, 'progress - {:>5.1f} %'.format(percent_progress))
            previous_percent_progress = percent_progress
        counter += 1

        # Compare the original sequence of a given annotation with the new sequence of the updated annotation
        if orig_seq_list[ii][1].lower() != new_seq_list[ii][1].lower():
            # If the sequences do not match, try to recapitulate the update process of the original sequence which
            # resulted in the updated sequence

            if new_seq_list[ii][0].type == '!row_del!':
                rowdel_msg = 'WARNING' \
                             'Found an annotation which was deleted from the new GFF file\n' \
                             'Original sequence:\n' + str(orig_seq_list[ii][1]) + '\n' \
                             'Original annotation:\n' + str(orig_seq_list[ii][0]) + '\n\n'
                write_text_to_file(fio_report, rowdel_msg)

            # GFF file entries with the 'variant' type will not be evaluated
            elif new_seq_list[ii][0].type != 'variant':
                orig_GFF = orig_seq_list[ii][0]
                orig_DNA = orig_seq_list[ii][1]
                new_GFF = new_seq_list[ii][0]
                new_DNA = new_seq_list[ii][1]

                relevant_indels = list()
                indel_startpos = list()
                indel_endpos = list()

                startpos = int(orig_GFF.start)
                indel_startpos.append(startpos)

                indel_endpos.append((startpos + comparison_length))

                # Get all InDels which occur before the currently investigated annotation
                all_indels_lower = indel_df[indel_df["POS"] <= (startpos + comparison_length)]

                # Find all InDels which are located (completely or partially) in the query stretch of the
                # investigated annotation. Only the InDels which span across the first N nucleotides of the
                # annotation will be evaluated later.
                for index, indel_row in all_indels_lower.iterrows():
                    if ((startpos <= int(indel_row.POS) + len(indel_row.REF) - 1) and int(indel_row.POS) <= startpos + comparison_length) or (int(indel_row.POS) <= startpos and (startpos + comparison_length) <= (int(indel_row.POS) + len(indel_row.REF) - 1)):
                        relevant_indels.append(index)
                        indel_startpos.append(indel_row.POS)

                        # Potentially, the InDel ends after the annotation's end - in that case, shift the
                        # end of the comparison sequence to the end of the InDel.
                        if indel_row.POS + len(indel_row.REF) - 1 > startpos + comparison_length:
                            indel_endpos.append(int(indel_row.POS) + len(indel_row.REF) - 1)

                relevant_indels_table = indel_df.ix[relevant_indels]
                indels_table = relevant_indels_table.copy()

                if indel_startpos:
                    shift = 0
                    # Similarly, the InDel might start before the annotation's start - in that case, shift the
                    # start of the comparison sequence to the start of the InDel.
                    if min(indel_startpos) < startpos:
                        comparison_startpos = min(indel_startpos)
                        shift = startpos - comparison_startpos
                    else:
                        comparison_startpos = startpos

                    comparison_endpos = max(indel_endpos)
                    endpos_shift = comparison_endpos - (startpos + comparison_length)

                    # Now that the start and end position of the sequence in the original FASTA file (which will be
                    # compared with the updated sequence) have been determined, extract this test sequence
                    test_sequence = orig_chromosome_seq[comparison_startpos-1: comparison_endpos]

                    # Based on the start position of the test sequence, adjust the POS of all relevant InDels.
                    # Afterwards, perform an update of the test sequence using this InDel table in order to
                    # generate th updated original sequence.
                    indels_table["POS"] = indels_table["POS"] - comparison_startpos + 1
                    updated_orig_seq, val_gff_lol, null = update_fa(test_sequence, indels_table, fio_dummy, fio_dummy, serious_errors)

                    cummulative_lendiff = sum([element[1] for element in val_gff_lol])

                    adjusted_new_DNA = ''

                    indels_table = indels_table[indels_table["POS"] == 1]
                    if len(indels_table) > 0:
                        for index, entry in indels_table.iterrows():
                            if len(entry["REF"]) > len(entry["ALT"]):
                                if shift != 0:
                                    updated_orig_seq = updated_orig_seq[1:]

                                adjusted_new_DNA = new_chromosome_seq[new_seq_list[ii][0].start - 1: new_seq_list[ii][0].start + comparison_length + endpos_shift + cummulative_lendiff]

                                if len(updated_orig_seq) > len(adjusted_new_DNA):
                                    updated_orig_seq = updated_orig_seq[:len(adjusted_new_DNA)]
                    else:
                        adjusted_new_DNA = new_chromosome_seq[new_seq_list[ii][0].start - 1: new_seq_list[ii][0].start + comparison_length + endpos_shift + cummulative_lendiff]

                    if not adjusted_new_DNA:
                        adjusted_new_DNA = new_chromosome_seq[new_seq_list[ii][0].start - 1: new_seq_list[ii][0].start + comparison_length + endpos_shift + cummulative_lendiff]

                    if updated_orig_seq.lower() != adjusted_new_DNA.lower():
                        adjusted_new_DNA = new_chromosome_seq[new_seq_list[ii][0].start - 1 + shift: new_seq_list[ii][0].start + comparison_length + endpos_shift + cummulative_lendiff]

                        if updated_orig_seq.lower() != adjusted_new_DNA.lower():

                            error_msg = 'WARNING' \
                                        "Found some disagreement between old and new annotations and sequences which cannot be explained using the VCF file's information\n" \
                                        'Original sequence:\n' + str(orig_DNA) + '\n' \
                                        'Original annotation:\n' + str(orig_GFF) + '\n\n' \
                                        'New sequence:\n' + str(new_DNA) + '\n' \
                                        'New annotation:\n' + str(new_GFF) + '\n\n'
                            write_text_to_file(fio_report, error_msg)

                            indels_table.to_csv(fio_report, sep="\t")

                            test_msg = '\nCS:       ' + str(comparison_startpos) + '\n' \
                                       'Shift:    ' + str(shift) + '\n' \
                                       'TestSeq:  ' + str(test_sequence) + '\n' \
                                       'Updated:  ' + str(updated_orig_seq) + '\n' \
                                       'newShort: ' + str(adjusted_new_DNA) + '\n\n\n'
                            write_text_to_file(fio_report, test_msg)

                            """
                            print(indels_table)
                            print('CS:       ' + str(comparison_startpos))
                            print('Shift:    ' + str(shift))
                            print('newStart: ' + str(find_new_start))
                            print('TestSeq:  ' + test_sequence)
                            print('Updated1: ' + updated_orig_seq)
                            print('Updated:  ' + updated_seq)
                            print('NewDNA:   ' + new_DNA_short)
                            """

                            validation_problem_no += 1

                else:
                    print('Hit ValueError!')
                    error_msg = "ValueError case\n" \
                                'original sequence:\n' + str(orig_DNA) + '\n' \
                                'original annotation:\n' + str(orig_GFF) + '\n\n' \
                                'new sequence:\n' + str(new_DNA) + '\n' \
                                'new annotation:\n' + str(new_GFF) + '\n\n'
                    write_text_to_file(fio_report, error_msg)
                    print(indel_row.POS in range(startpos, startpos + comparison_length) or indel_row.POS + len(indel_row.REF) in range(startpos, startpos + comparison_length))

    if validation_problem_no == 0:
        write_text_to_file(fio_report, 'validation completed without disagreement')
    else:
        write_text_to_file(fio_report, 'validation failed - encountered ' + str(validation_problem_no) + ' disagreements during validation of shifted sequence annotations')

    try:
        os.remove(fio_dummy_path)
    except FileNotFoundError:
        pass

    try:
        fio_report.close()
    except:
        pass

    try:
        if validation_problem_no == 0:
            os.rename(fio_report_path, fio_path + '.validation_success')
        else:
            os.rename(fio_report_path, fio_path + '.validation_fail')
    except:
        pass

    if validation_problem_no == 0:
        validation_status = 'success'
    else:
        validation_status = 'fail'
    
    return validation_status, validation_problem_no

    """
    except:
        try:
            os.remove(fio_dummy_path)
        except FileNotFoundError:
            pass

        try:
            exc_info = sys.exc_info()
            try:
                print_exception(*exc_info, file=fio_report_path)
            except:
                pass
            del exc_info
        except:
            pass

        try:
            fio_report.close()
        except:
            pass

        try:
            os.rename(fio_report_path, fio_path + '.validation_error')
        except:
            pass

        return 'error'
    """


def get_overlaps(indel_drop_df):
    ## Retrieve the overlapping InDels from the InDel data frame and store them in GFF format as type 'variant' in their
    ## own data frame. This data frame will later be appended to the big GFF file.

    line_list = [['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    line_no = 0
    for index, line in indel_drop_df.iterrows():
        seqid = line["CHROM"]
        source = 'CIDer_Var'
        type = 'variant'
        start = line["POS"]
        end = line["POS"] + len(line["REF"])
        score = '.'
        strand = '+'
        phase = '.'
        attributes = 'Name=Variant' + str(line_no) + ';alternative_sequence=' + line["ALT"]
        line_list.append([seqid, source, type, start, end, score, strand, phase, attributes])

        line_no += 1

    overlap_df = pd.DataFrame(line_list[1:], columns=line_list[0])

    return overlap_df


def mergefiles(merge_dir, path_output, annot_type=None, fio_CIDerrun=None, from_within=False, gff_name='', fasta_name='', fasta_only=False):
    ## If the user has chosen to merge cache files, this procedure will handle the merge process

    # Create output folder
    if not os.path.isdir(path_output):
        os.makedirs(path_output)

    if from_within is False:
        print('')
        print('CIDer +++ merge mode')
        print('')
        print(' inputs')
        print(' input directory:  ' + merge_dir)
        print(' output directory: ' + path_output)
        print('')
        print('')

    ## First, CIDer merge looks for all FA and GFF files in the specified input directory. Then, the appropriate FA and
    ## GFF files are read and their contents are written to CIDer_merge.fa and CIDer_merge.gff, respectively.

    if not fasta_only:
        extension_list = ['.fa', '.gff']
    else:
        extension_list = ['.fa']

    for extension in extension_list:
        if extension == '.gff':
            write_text_to_file(fio_CIDerrun, 'writing the updated annotation list to GFF')
            print(' writing the updated annotation list to GFF...', end=' ')
        else:
            write_text_to_file(fio_CIDerrun, 'writing the updated FASTA file')
            print(' writing the updated FASTA file...', end=' ')
        fio_CIDerrun.flush()

        # Find files
        if extension == '.gff':
            try:
                comment_file = [file for file in os.listdir(merge_dir) if file == 'comments.txt']
                comment_file = comment_file[0]
                with open((merge_dir + comment_file), 'r') as comment:
                    comment_cont = comment.read()
            except:
                comment_cont = ''

        ext_files = [file for file in os.listdir(merge_dir) if file.endswith(extension)]
        ext_files = natsorted(ext_files)

        if gff_name and fasta_name:
            if extension == '.fa':
                output_path = path_output + fasta_name

                new_file_dummy = open(output_path, 'w')
                new_file_dummy.close()

            else:
                if annot_type:
                    pos = re.search('.gff', gff_name.lower())
                    output_path = path_output + gff_name[:pos.start()] + '_' + annot_type + gff_name[pos.start():]
                else:
                    output_path = path_output + gff_name

        if extension == '.gff':
            with open(output_path, 'w') as fio_dummy:
                fio_dummy.write(comment_cont)

        for file in ext_files:
            input_path = merge_dir + file
            with open(input_path, 'r') as input_file:
                with open(output_path, 'a') as dest_file:
                    input_content = input_file.read()
                    dest_file.write(input_content)
        print('done!')


def write_text_to_file(file, text):
    ## Makes writing to a file easier.

    file.write(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '\t' + text + '\n')
    file.flush()


def check_input_version(version, type):
    ## Check the version of VCF and GFF files and display a warning message if necessary

    if type == 'VCF':
        if version == 0:
            print(' WARNING'
                  ' The file format of your VCF file could not be determined.\n'
                  ' The recommended VCF file specification version is VCFv4.1\n'
                  ' CIDer might still work with your VCF file version, though.\n\n')
        elif 0 < version < 4.1:
            print(' WARNING'
                  ' Your VCF file follows VCFv' + str(version) + ' specifications.\n'
                  ' The recommended VCF file specification version is VCFv4.1\n'
                  ' CIDer might still work with this VCF version, though.\n\n')

    elif type == 'GFF':
        if version == 0:
            print(' WARNING'
                  ' The file format of your GFF file could not be determined.\n'
                  ' The recommended GFF file specification is GFF3\n'
                  ' CIDer might still work with your GFF file, though.\n\n')


## MAIN ###############################################################################################################

## Input handling
# Read in the command line inputs and evaluate
input_set = process_input(sys.argv[:], help_text, license_text)

path_vcf, path_gff, path_fasta, path_output, cache_path_output, num_cpus, exclude_id, vcf_version, gff_version, annot_type, fa_only, no_validation = input_set

# Create a file - CIDer.running - which tracks overall CIDer progress
fio_CIDerrun_path = path_output + 'CIDer.running'
fio_CIDerrun = open(fio_CIDerrun_path, 'w')

print('\n\nChromosome-InDel-merger\n')

if not fa_only:
    inputs = ' VCF file:\t' + path_vcf + '\n' \
             ' GFF file:\t' + path_gff + '\n' \
             ' FASTA file:\t' + path_fasta + '\n\n' \
             ' output directory:\t' + path_output[:-1] + '\n' \
             ' cache directory:\t' + cache_path_output[:-1] + '\n'
else:
    inputs = ' VCF file:\t' + path_vcf + '\n' \
             ' FASTA file:\t' + path_fasta + '\n\n' \
             ' output directory:\t' + path_output[:-1] + '\n' \
             ' cache directory:\t' + cache_path_output[:-1] + '\n'

print(inputs)

options = '\noptions\n' \
          ' CPUs:\t' + str(num_cpus) + '\n' \
          ' Excluded:\t' + str(exclude_id) + '\n' \
          ' Type:\t' + str(annot_type) + '\n' \
          ' VCFv:\t' + str(vcf_version) + '\n' \
          ' GFFv:\t' + str(gff_version) + '\n' \
          ' FASTAonly:\t' + str(fa_only) + '\n' \
          ' noValidation:\t' + str(no_validation) + '\n\n'

# Write basic input information to CIDer.running
start_log = ' CIDer v' + version_info + '\ninputs\n' + inputs + options
fio_CIDerrun.write(start_log)
write_text_to_file(fio_CIDerrun, 'starting CIDer')

# Check whether the input VCF and GFF files follow the required file type specifications
check_input_version(float(vcf_version), 'VCF')

if not fa_only:
    check_input_version(int(gff_version), 'GFF')

    # Read in input GFF file and store contents in a pandas data frame
    write_text_to_file(fio_CIDerrun, 'reading GFF file')
    print('\n reading GFF file...')
    genes_df = read_gff(path_gff, annot_type)
    print(' done!\n')

# Read in the InDel table and store it in a pandas data frame
write_text_to_file(fio_CIDerrun, 'processing VCF file')
indels_df = make_indel_df(path_vcf, cache_path_output)
size_vcf = len(indels_df)
print(' done!\n')

# Parse the chromosomes/sequences from the FASTA file one by one; each sequence is stored in form of a Biopython Seq.
# Subsequently, several variables
#
# * chromosome.id - contains the ID of the sequence extracted from the FASTA file
# * chromosome.description - contains the description of the sequence extracted from the FASTA file
# * chromosome.seq - contains the sequence extracted from the FASTA file
# * chromosome_indel - holds the InDel table
# * chromosome_gff - contains the gene table
# * progress_list - is needed to hand over the progress_list which is a list associated with multiprocess.Manager)
# * index - is the number of sequence in the FASTA file
# * cache_path_output - is the path to the directory, where cache files shall be written to
#
# are handed over to the make_new_genome procedure. In parallel, a subprocess is initiated in which make_new_genome
# will run. The make_new_genome procedure uses the information stored in the VCF file in order to introduce InDels
# into chromosome.seq. Based on these sequence changes, the GFF file is updated in order to in the end hold the
# correct positions for each annotation.
# While the make_new_genome subprocesses are running, the main program is halted and waits until all subprocesses are
# completed. Once the make_new_genome process for a given sequence is completed, the results of this run are saved to
# the cache directory, and the subprocess is terminated.
# After all make_new_genome subprocesses are finished, the GFF and FASTA files saved in the cache directory are merged
# in order to create one complete GFF and one complete FASTA file which then hold data from all chromosomes.
chrdesc_list = list()
chrnewseq_list = list()
chrnewgff_list = list()
processes = list()
indel_output = list()
genome_result_list = list()
validation_result_list = list()
serious_error_list = list()
chromosome_id_list = list()
result_msg_list = list()
deletedsize = 0
validation_problem_number = 0

write_text_to_file(fio_CIDerrun, 'accessing FASTA file')
print(' accessing FASTA file...')

index = 0
with open(path_fasta, 'r') as full_genome_file:
    with Pool(processes=int(num_cpus)) as pool:

        if num_cpus > 1:
            write_text_to_file(fio_CIDerrun, 'spawning subprocesses')
            print(' spawning subprocesses...')
        else:
            write_text_to_file(fio_CIDerrun, 'starting construction of new genome')
            print(' construction of new genome in progress...')

        for chromosome in SeqIO.parse(full_genome_file, 'fasta'):
            if chromosome.id not in exclude_id:
                chrdesc_list.append(chromosome.description)

                chromosome_indel = indels_df[indels_df.CHROM == chromosome.id]

                if not fa_only:
                    chromosome_gff = genes_df[genes_df.seqid == chromosome.id]
                else:
                    chromosome_gff = pd.DataFrame()

                ## THIS IS WHERE THE MAGIC HAPPENS
                indel_output.append(pool.apply_async(make_new_genome, (chromosome.id, chromosome.description,
                                                     chromosome.seq, chromosome_indel,
                                                     chromosome_gff, index, cache_path_output, fa_only, no_validation)))

                # delete the part of indels_df that has just been processed
                indels_df = indels_df[indels_df.CHROM != chromosome.id]

                index += 1

        # make_new_genome returns a set containing the chromosome ID, 'success' or 'error' depending on whether the
        # subprocess could be completed without a critical error or not, and a return value of the validation process
        # which is either 'success', 'error', or 'disagreement' depending on whether validation was successful, had a
        # critical error, or some disagreements between sequence and annotation were detected
        for result in indel_output:
            result_set = result.get()

            chromosome_id, genome_result, validation_result, vcf_deleted_size, val_prob_no, serious_errors = result_set
            chromosome_id_list.append(chromosome_id)
            genome_result_list.append(genome_result)
            genome_result_list.append(genome_result)
            validation_result_list.append(validation_result)
            serious_error_list.append(serious_errors)
            deletedsize += vcf_deleted_size
            validation_problem_number += val_prob_no

            if genome_result != 'success':
                genome_msg = chromosome_id + ' - genome: ' + genome_result
                result_msg_list.append(genome_msg)

            if validation_result == 'disagreement' or validation_result == 'error':
                validation_msg = chromosome_id + ' - validation: ' + validation_result
                result_msg_list.append(validation_msg)

if not no_validation or not fa_only:
    validation_msg = 'number of validation problems\t' + str(validation_problem_number)
    result_msg_list.append(validation_msg)

deleted_msg = 'deleted VCF rows\t' + str(deletedsize) + '\noriginal size\t' + str(size_vcf)
result_msg_list.append(deleted_msg)
print(' done!\n')

## OUTPUT
# Determine the names of the new FASTA and GFF files based on the original ones
name_new_genes = 'CIDer_' + path_gff[(-(path_gff[::-1].find('/'))+len(path_gff)):]
name_new_chromosomes = 'CIDer_' + path_fasta[(-(path_fasta[::-1].find('/'))+len(path_fasta)):]

# Start to merge the individual FASTA and GFF files in the cache directory in order to obtain a single genome FASTA and
# a single genome GFF file
write_text_to_file(fio_CIDerrun, 'writing output files')
mergefiles(cache_path_output, path_output, annot_type, fio_CIDerrun, from_within=True, gff_name=name_new_genes, fasta_name=name_new_chromosomes, fasta_only=fa_only)

# Do some evaluations of CIDer's overall performance in order to show the user some informations about how confident
# they can be about CIDer's results and then end the program
num_genome_errors = genome_result_list.count('error')
num_serious_errors = sum(serious_error_list)
num_validation_errors = validation_result_list.count('error')
num_validation_problem_nos = validation_result_list.count('disagreement')
result_detail = '\n'.join(result_msg_list)

write_text_to_file(fio_CIDerrun, 'CIDer completed')

if not no_validation or not fa_only:
    end_report = '\noverview\n' \
                 'encountered ' + str(num_genome_errors) + ' critical error(s) during genome construction\n' \
                 'encountered ' + str(num_serious_errors) + ' chromosome(s) with serious errors during genome construction\n' \
                 'encountered ' + str(num_validation_errors) + ' chromosome(s) with validation errors\n' \
                 'encountered ' + str(num_validation_problem_nos) + ' chromosome(s) with validation disagreements\n'\
                 '\ndetails\n' + result_detail
else:
    end_report = '\noverview\n' \
                 'encountered ' + str(num_genome_errors) + ' critical error(s) during genome construction\n' \
                 'encountered ' + str(num_serious_errors) + ' chromosome(s) with serious errors during genome construction\n' \
                 '\ndetails\n' + result_detail

fio_CIDerrun.write(end_report)
fio_CIDerrun.close()
os.rename(fio_CIDerrun_path, path_output + 'CIDer.finished')

print('')
if num_genome_errors == 0 and num_validation_errors == 0 and num_validation_problem_nos == 0:
    print(' completed without errors!')
else:
    if not no_validation or not fa_only:
        print(' completed but encountered\n'
              ' - ' + str(num_genome_errors) + ' critical error(s)\n'
              ' - ' + str(num_serious_errors) + ' serious error(s)\n'
              ' - ' + str(num_validation_errors) + ' validation error(s)\n'
              ' - ' + str(num_validation_problem_nos) + ' validation disagreement(s)')
    else:
        print(' completed but encountered\n'
              ' - ' + str(num_genome_errors) + ' critical error(s)\n'
              ' - ' + str(num_serious_errors) + ' serious error(s)\n')

print('')
