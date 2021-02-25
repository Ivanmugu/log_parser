#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on Sat Feb 20, 2021
@author: ivanmugu
"""

import sys
import argparse
import textwrap
import os
import re
import csv
from itertools import zip_longest


def user_input():
    """
    Gets input from user to run the program, parse the arguments and check
    correct usage.

    Returns
    -------
    argparse object
        Contains all the needed information of the arguments provided by the
        user, i.e. the path to the input and output directories.
    """
    epilog_msg = (
        """\
        To run this program, you need a folder containing subfolders with
        unicycler results. For example, you can have a folder named results
        that will be the input folder and will contain the SW0001 and SW0002
        unicycler result subfolders as follows:

        results/
            SW0001/
                unicycler.log
            SW0002/
                unicycler.log

        The program will first search for unicycler.log in the SW0001
        subfolder, parse the file and create two tables as csv files with
        relevant information about the assembly. The first column of each row
        will contain the name of the subfolder, in this case SW0001. The
        program will continue with the next subfolders, SW0002, and append the
        extracted information to the csv files. The generated csv files will be
        named molecules_summary and assemblies_summary. The subfolders (SWXXXX)
        can have other files or folders, the program will only analyze
        unicycler.log files.

        Usage examples:
        python3 log_parser.py -i ~/Documents/results -o ~/Documents/results

        python3 log_file_parser -i . -o .

        python3 log_file_parser -o ~/Desktop -i ~/Documents/assembly/results
        """)
    # Parsing aguments and providing help
    parser = argparse.ArgumentParser(
        prog='log_parser.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Parser of unicycler.log files.",
        epilog=textwrap.dedent(epilog_msg))
    parser.add_argument('-i', '--input', help="path to input directory")
    parser.add_argument('-o', '--output', help="path to output directory")
    args = parser.parse_args()
    # Checking if user provided mandatory Argument
    if (args.input is None) and (args.output is None):
        parser.exit(1, message=textwrap.dedent("""\
        error: missing arguments, you did not provide paths to input directory
        nor to output directory\n"""))
    if (vars(args)).get("input") is None:
        parser.exit(1, message=textwrap.dedent("""\
        error: missing argument, you did not provide the path to input
        directory\n"""))
    if (vars(args)).get("output") is None:
        parser.exit(1, message=textwrap.dedent("""\
        error: missing argument, you did not provide the path to output
        directory\n"""))
    # Checking if output and input folder exist
    if not  os.path.exists(args.input):
        sys.exit("error: input directory does not exist")
    if not  os.path.exists(args.output):
        sys.exit("error: output directory does not exist")

    return args

def dir_names_addresses(input_directory, file_name):
    """
    Iterates over subdirectories of a given directory (input_directory) looking
    for files with the same name (file_name). For the purpose of our lab,
    file_name is "unicycler.log". Returns the path to these file_names.
    Additionaly, returns the names af all the subfolders cantaining file_name.

    Parameters
    ----------
    input_directory : str
        Name of the directory to analyze.
    file_name : str
        Name of the file which path is needed

    Returns
    -------
    tuple
        Index 0 has a list of the addresses.
        Index 1 has a list of the directories' names.
    """
    # List of files' addresses
    file_addresses = []
    # List of directory's names
    dir_names = []
    # Get all files' and folders' names in the indicated directory
    file_dictectory_names = os.listdir(input_directory)
    # Iterate over all the files and folders contained in input_directory
    for filename in file_dictectory_names:
        # Check if the currect object is a folder or not
        if os.path.isdir(os.path.join(input_directory, filename)):
            # Checking if folder contains file_name
            if not os.path.exists(
                    os.path.join(input_directory, filename, file_name)):
                print("folder " + filename + " does not have " + file_name)
                continue
            # Getting folder's name
            dir_names.append(filename)
            # Getting path of 'unicycler.log'
            file_addresses.append(
                os.path.join(input_directory, filename, file_name))
    # print(file_addresses, dir_names)
    return (file_addresses, dir_names)

def extractor(table, headers):
    """
    Reads rows of a table from an input file and convert it into a dictionary.
    It uses the values of the headers argument as keys for the dictionary. The
    table argument get input form a file that was opened for reading. This
    function was designed to be used in the molecularSummary function, but can
    be addapted for other uses.

    Parameters
    ----------
    table : iterable object
        File opened for reading that will be processed.
    headers : list object
        Contains a list of the table's headers

    Returns
    -------
    dictionary
        A dictionary of dictionaries. The key of the primary dictionary is the
        Lenght of the molecule that correspond to every row of the analyzed
        table.
        For example, a possible output is this:
        {'5,179,485': {'Component': '1', 'Length': '5,179,485', 'Status': 'complete'},
         '131,127': {'Component': '2', 'Length': '131,127', 'Status': 'complete'},
         '4,074': {'Component': '3', 'Length': '4,074', 'Status': 'complete'}}
    """
    extracted_table = []
    # Iterate over file (table) to extract data
    for row in table:
        # Break when reach the end of the table. log files separeates the
        # end of a table with '\n'
        if row == '\n':
            break
        # If 'none found' in row replace with 'none_found' to facilitate parsing
        if 'none found' in row:
            row = re.sub('none found', 'none_found', row)
        # Format lines by replacing line's spaces with tabs and convert line into a list
        line_list = re.sub('\\s+', '\t', row.strip()).split('\t')
        # If data in first column is not numeric don't get info
        if not line_list[0].isnumeric():
            continue
        # Convert the list created above (line_list) into dictionary. Use the values of
        # 'headers' (provided as argument in the function) as keys. Then, append the
        # dictionary to the variable extracted_table of type list. This will create a
        # list of dictionaries. For example, a possible output is this:
        # [{'Component': 1, 'Length': 5179485}, {'Component': 2, 'Length': 131127}]
        extracted_table.append(dict(zip_longest(headers, line_list)))
    # Convert the list of dictionaries into a dictionary of dictionaries
    final = {}
    for index in extracted_table:
        final[index.get('Length')] = index
    return final

def concatenate_molecules_summary(
        status, depth, input_folder_name, output_file, table_headers):
    """
    Concatenates the extracted results from status and depth tables at the end
    of the output_file.

    Parameters
    ----------
    status : dictionary
        Dictionary of dictionaries containg the information of the status table
        extracted from the log file using the extractor function
    depth : dictionary
        Dictionary of dictionaries containg the information of the depth table
        extracted from the log file using the extractor function
    input_folder_name : string
        Name of directory that contains the unicycler.log file being analyzed
    output_file : string
        Path to file used to concatenate the results of the status and depth
        tables.
    table_headers : list
        Headers of the table molecules_summary

    Returns
    -------
    bool
        True for success, otherwise False.
    """
    # Opening outfile
    with open(output_file, 'a') as outfile:
        # Creting an object to write the csv file
        writer = csv.DictWriter(outfile, fieldnames=table_headers)
        # Both status and depth are dictionaries of dictionaries, and the keys
        # are the lenght of the molecules. Check the documentation of the extractor
        # function for more details.
        for key in status:
            # Some of the information in status and depth tables is redundant, like
            # segment or length. Additionaly, some molecules don't have information
            # in the depth table. Therefore, if the length of the molecule (key) is
            # present in both depth and status tables, the program will extract
            # information from depth and status. Oherwise, the columns
            # corresponding to the part of the depth table will be fill with None.
            #
            # If Length (key) from status table is in depth table, get information
            # from both status and depth tables.
            if key in depth:
                relevant = {
                    'Folder_name': input_folder_name,
                    'Component': status.get(key).get('Component'),
                    'Segments': status.get(key).get('Segments'),
                    'Links': status.get(key).get('Links'),
                    'Length': status.get(key).get('Length'),
                    'N50': status.get(key).get('N50'),
                    'Longest_segment': status.get(key).get('Longest_segment'),
                    'Status': status.get(key).get('Status'),
                    'Depth': depth.get(key).get('Depth'),
                    'Starting_gene': depth.get(key).get('Starting_gene'),
                    'Position': depth.get(key).get('Position'),
                    'Strand': depth.get(key).get('Strand'),
                    'Identity': depth.get(key).get('Identity'),
                    'Coverage': depth.get(key).get('Coverage')}
            # Otherwhise put 'None' in Depth
            else:
                relevant = {
                    'Folder_name': input_folder_name,
                    'Component': status.get(key).get('Component'),
                    'Segments': status.get(key).get('Segments'),
                    'Links': status.get(key).get('Links'),
                    'Length': status.get(key).get('Length'),
                    'N50': status.get(key).get('N50'),
                    'Longest_segment': status.get(key).get('Longest_segment'),
                    'Status': status.get(key).get('Status'),
                    'Depth': None,
                    'Starting_gene': None,
                    'Position': None,
                    'Strand': None,
                    'Identity': None,
                    'Coverage': None}
            # Append relevant information in the outfile
            writer.writerow(relevant)
    return True


def molecules_summary(file_addresses, dir_names, output_folder):
    """
    Creates a csv file with relevant information retrieved from all the
    unicycler.log files contained in the primary subfolders of a given
    directory. The tables to parse are the status and depth tables.

    Parameters
    ----------
    file_addresses : list object
        List of addresses that will be used to open the unicycler.log files

    dir_names : list object
        List of directory names that have the unicycler.log files

    output_folder : string
        Path to the directory that will be used to save the output
        molecules_summary.csv file

    Returns
    -------
    bool
        True for success, otherwise False.
        If the the path to output_folder doesn't exist returns False
    """
    # Checking if output_folder exists, otherwise return False
    if os.path.exists(output_folder):
        output_path = os.path.join(output_folder, "molecules_summary.csv")
    else:
        return False
    # Opening the outfile (molecules_summary) to save the summary and make the
    # headers of the table
    with open(output_path, 'a') as outfile:
        # Headers of the molecules_summary table. These headers correspond to the ones
        # found in unicycler.log
        fieldnames = [
            'Folder_name','Component', 'Segments', 'Links', 'Length', 'N50',
            'Longest_segment', 'Status', 'Depth', 'Starting_gene', 'Position',
            'Strand', 'Identity', 'Coverage']
        # Creting an object to write the csv file
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        # Write header of table in outfile
        writer.writeheader()
    # Iterate over each directory
    for i, address in enumerate(file_addresses):
        # Dictionaries to save the status and depth tables
        status = {}
        depth = {}
        # Opening log file
        with open(address, 'r') as log_file:
            # Iterate over log file
            for line in log_file:
                # If 'Component' and 'Status' are found in line, then extract table status
                if re.search('^Component.*Status', line):
                    # Convert header 'Longest segment' into 'Longest_segment'
                    headers = re.sub('Longest segment', 'Longest_segment', line)
                    # Replace line's spaces with tabs and convert headers into a list
                    headers = re.sub('\\s+', '\t', headers.strip()).split('\t')
                    # Extract table status using the extractor function
                    status = extractor(log_file, headers)
                # If 'Segment' and 'Depth' are found in line extract table depth
                if re.search('^Segment.*Depth', line):
                    # Convert header 'Starting gene' into 'Starting_gene'
                    headers = re.sub('Starting gene', 'Starting_gene', line)
                    # Replace line's spaces with tabs and convert headers into a list
                    headers = re.sub('\\s+', '\t', headers.strip()).split('\t')
                    # Extract table depth using the extractor function
                    depth = extractor(log_file, headers)
            # Saving (concatenate) info from status and depth variables into outfile
            concatenate_molecules_summary(
                status, depth, dir_names[i], output_path, fieldnames)
    return True


def extract_best_k_mer(table):
    """
    Read rows of a table from an input file until it finds best K-mer and
    converts it into a list

    Parameters
    ----------
    table : iteratable object
        File opened for reading that will be processed

    Returns
    -------
    list
        List containing the row with the best K-mer. The order of the headers
        are the following: [0]->K-mer, [1]->Contings, [2]->Dead_ends,
        [3]->Score
    """
    # Looking for the best in table
    for best_line in table:
        if 'best' in best_line:
            # Get the row, replace row's spaces with tabs, and
            # convert row into a list
            best = re.sub('\\s+', '\t',
                          best_line.strip()).split('\t')
            break
    return best

def extract_alignment_summary(table):
    """
    Read rows of the Read alignment summary table from a unicycler.log file and
    convert them into a list

    Parameters
    ----------
    table : iteratable object
        File opened for reading that will be processed

    Returns
    ------
    list
        List containing the information of the table. The order of the list is
        as follows: [0]->Total_read_count, [1]->Fully_aligned_reads,
        [2]->Partially_aligned_reads, [3]->Unaligned_reads,
        [4]->Total_bases_aligned, [5]->Mean_alignment_identity
    """
    # List to save info
    alignment_summary_list = []
    for alignment_summary in table:
        # Break when reach the end of the table. log files separeates the
        # end of a table with '\n'
        if alignment_summary == '\n':
            break
        # If it find a row with '--' ignore and continue
        if '--' in alignment_summary:
            continue
        # Replace single line's spaces with '_'
        alignment_summary = re.sub(r'([^\s])(\s)([^\s])',
                                   r'\1_\3', alignment_summary)
        # Replace multiple line's spaces with '\t' and conver line in list
        alignment_summary = (re.sub(r'\s+', r'\t',
                                    alignment_summary)).split('\t')
        # Extract relevant data. The second column is the one with therefore
        # values; therefore the program uses index 1 of aligment_summary
        alignment_summary_list.append(alignment_summary[1])
    return alignment_summary_list

def concatenate_assemblies_summary(
        best, alignment_summary_list, input_folder_name, output_file,
        table_headers):
    """
    Concatenate the extracted results from K-mer and Read alignment summary
    tables at the end of the output_file.

    Parameters
    ----------
    best : list
        List containing the best K-mer parameter.
    alignment_summary_list : list
        List containing the results of the Read alignment summary table.
    input_folder_name : string
        Name of directory that contains the unicycler.log file being analyzed
    output_file : string
        Path to file used to concatenate the results of the K-mer and Read
        alignment summary tables.
    table_headers : list
        Headers of the table assemblies_summary
    """
    # Opening output_file
    with open(output_file, "a") as outfile:
        # Creating an object to write the csv file
        writer = csv.DictWriter(outfile, fieldnames=table_headers)
        # Compaling information
        relevant = {'Folder_name': input_folder_name,
                    'K-mer_best': best[0],
                    'Contigs_best': best[1],
                    'Dead_ends_best': best[2],
                    'Score_best': best[3],
                    'Total_read_count': alignment_summary_list[0],
                    'Fully_aligned_reads': alignment_summary_list[1],
                    'Partially_aligned_reads': alignment_summary_list[2],
                    'Unaligned_reads': alignment_summary_list[3],
                    'Total_bases_aligned': alignment_summary_list[4],
                    'Mean_alignment_identity': alignment_summary_list[5]}
        # Write relevant info in outfile
        writer.writerow(relevant)

def assemblies_summary(file_addresses, dir_names, output_folder):
    """
    Creates a csv file with relevant information retrieved from all the unicycler.log
    files contained in the primary subfolders of a given directory. The tables to parse
    are the K-mer and Read alignment summary tables.

    Parameters
    ----------
    file_addresses : list object
        List of addresses that will be used to open the unicycler.log files

    dir_names : list object
        List of directory names that have the unicycler.log files

    output_folder : string
        Path to save the output assemblies_summary.csv file

    Returns
    -------
    bool
        True for success, otherwise False.
        If the the path to output_folder doesn't exist returns False.
    """
    # Checking if output_folder exists, otherwise return False
    if os.path.exists(output_folder):
        output_path = os.path.join(output_folder, "assemblies_summary.csv")
    else:
        return False
    # Opening the outfile to save the summary
    with open(output_path, 'a') as outfile:
        # Headers of table assemblies summary
        fieldnames = ['Folder_name', 'K-mer_best','Contigs_best', 'Dead_ends_best',
                      'Score_best', 'Total_read_count', 'Fully_aligned_reads',
                      'Partially_aligned_reads', 'Unaligned_reads',
                      'Total_bases_aligned', 'Mean_alignment_identity']
        # Creting an object to write the csv file
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        # Write header of table in outfile
        writer.writeheader()
    # Iterate over each directory
    for i, address in enumerate(file_addresses):
        # List to save best K-mer
        best = []
        # Opening log file
        with open(address, 'r') as log_file:
            # Iterate over log file
            for line in log_file:
                # If 'K-mer', 'Contigs', 'Dead ends' and 'Score' are found in line extract table
                if re.search('^K-mer.*Contigs.*Dead ends.*Score', line):
                    best = extract_best_k_mer(log_file)
                # If 'Read alignment summary' in line extract table
                if re.search('Read alignment summary', line):
                    # List to save info
                    alignment_summary_list = extract_alignment_summary(log_file)
        # If the leng of best is zero, it means that the file doesn't have the
        # table k-mer. Therefore, we don't need the info of that unicycler.log file
        if len(best) == 0:
            continue
        concatenate_assemblies_summary(
            best, alignment_summary_list, dir_names[i], output_path,
            fieldnames)
    return True

def main():
    """
    Main function
    """
    # Getting input from user
    args = user_input()
    # Getting path to the directory that carries the subdirectories to be analyzed
    input_directory = args.input
    # Getting path to output directory
    output_directory = args.output
    # Getting a list of the file addresses and directory names
    file_addresses = dir_names_addresses(input_directory, 'unicycler.log')[0]
    dir_names = dir_names_addresses(input_directory, 'unicycler.log')[1]
    # Making new files with extracted information from unicycler.log files
    molecules_summary(file_addresses, dir_names, output_directory)
    assemblies_summary(file_addresses, dir_names, output_directory)
    # If everything went OK
    print(
        "The assemblies_summary.csv and molecules_summary.csv files are in:")
    print(f"{os.path.abspath(output_directory)}")
    print("log_parser.py is done!")
    sys.exit(0)

if __name__ == '__main__':
    main()
