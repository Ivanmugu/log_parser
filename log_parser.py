#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File name: log_parser.py
Author: Ivan Munoz-Gutierrez
Date created: 02/20/2021
Date last modified: 03/18/2021
Python version: 3.9
"""

import sys
import os
import argparse
import textwrap
import re
import csv
from itertools import zip_longest


def user_input():
    """
    Parse command line arguments provided by the user and check correct usage.

    Returns
    -------
    argparse object (.input and .output)
        .input : holds the path to the input directory
        .output : holds the path to the output directory
    """
    epilog_msg = ("""
    To run this program, you need a folder containing subfolders with Unicycler
    results. For example, you can have a hypothetical folder named results that
    could be the input folder and could contain the SW0001 and SW0002 Unicycler
    result subfolders as follows:

    ~/Documents/results/
                    SW0001/
                        unicycler.log
                    SW0002/
                        unicycler.log

    The program will first search for unicycler.log in the SW0001 subfolder,
    parse the file and create two tables as csv files with relevant information
    about the assembly. The first column of each row will contain the name of
    the subfolder, in this case SW0001. The program will continue with the next
    subfolders, SW0002, and append the extracted information to the csv files.
    The generated csv files will be named molecules_summary and
    assemblies_summary. The subfolders (SWXXXX) can have other files or
    folders, the program will only analyze unicycler.log files.

    Example of usage
    ----------------
    python3 log_parser.py -i ~/Documents/results -o ~/Documents/results

    python3 log_file_parser --input . --output .

    python3 log_file_parser -o ~/Desktop -i ~/Documents/assembly/results
    """)
    # Parsing aguments and providing help.
    parser = argparse.ArgumentParser(
        prog='log_parser.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Parser of unicycler.log files.",
        epilog=textwrap.dedent(epilog_msg))
    parser.add_argument('-i', '--input', help="path to input directory")
    parser.add_argument('-o', '--output', help="path to output directory")
    args = parser.parse_args()
    # Checking if user provided mandatory argument.
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
    # Checking if output and input folder exist.
    if not  os.path.exists(args.input):
        sys.exit("error: input directory does not exist")
    if not  os.path.exists(args.output):
        sys.exit("error: output directory does not exist")

    return args

def get_file_paths(input_directory, file_name):
    """
    Get the paths of files with the same name that are contained in primary
    subdirectories.

    Iterate over subdirectories of a given directory (input_folder) looking
    for files with the same name (file_name). Return the paths to these files.
    This function was designed to get the path of the unicycler.log files
    generated by Unicycler.

    Parameters
    ----------
    file_name : str
        Name of the file which path is needed, i.e. unicycler.log.
    input_folder : str
        Path to the directory to analyze.

    Returns
    -------
    file_addresses : list
        List of paths to file_name.

    For example, if you have the hypothetical tree:

    ~/assemblies/
            SW0001/
                unicycler.log
            SW0002/
                unicycler.log

    You will get the following:
    ["~/assemblies/SW0001/unicycler.log",
     "~/assemblies/SW0002/unicycler.log"]
    """
    # List of files' addresses.
    file_addresses = []
    # Getting all files' and folders' names in from input_directory.
    file_dictectory_names = os.listdir(input_directory)
    # Iterating over all the files and folders contained in input_directory.
    for folder in file_dictectory_names:
        # Check if the currect object is a folder or not.
        if os.path.isdir(os.path.join(input_directory, folder)):
            # If folder contains file_name, get path and append it to
            # file_addresses. Otherwise, print an error message and continue.
            if not os.path.exists(
                    os.path.join(input_directory, folder, file_name)):
                print("folder " + folder + " does not have " + file_name)
                continue
            file_addresses.append(
                os.path.join(input_directory, folder, file_name))

    return file_addresses

def extractor(table, headers):
    """
    Read rows of a table from an input file and convert them in a dictionary.

    The headers parameter provides the keys for the dictionary. The table
    parameter gets input form a file that was opened for reading.

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

        Example of possible output:
        {'5,179,485': {
            'Component': '1', 'Length': '5,179,485', 'Status': 'complete'},
         '131,127': {
             'Component': '2', 'Length': '131,127', 'Status': 'complete'},
         '4,074': {
             'Component': '3', 'Length': '4,074', 'Status': 'complete'}}
    """
    extracted_table = []
    # Iterating over file (table) to extract data.
    for row in table:
        # If the end of the table is reached, break. log files separeates the
        # end of a table with '\n'.
        if row == '\n':
            break
        # If 'none found' in row, replace it with 'none_found' to facilitate
        # parsing.
        if 'none found' in row:
            row = re.sub('none found', 'none_found', row)
        # Formating lines by replacing line's spaces with tabs and convert line
        # into a list.
        line_list = re.sub('\\s+', '\t', row.strip()).split('\t')
        # If data in first column is not numeric don't get info.
        if not line_list[0].isnumeric():
            continue
        # Converting the list created above (line_list) into dictionary.
        # 
        # Use the values of 'headers' (provided as argument in the function) as
        # keys. Then, append the dictionary to the variable extracted_table of
        # type list. This will create a list of dictionaries. For example, a
        # possible output is this:
        # [{'Component': 1, 'Length': 5179485},
        #  {'Component': 2, 'Length': 131127}]
        extracted_table.append(dict(zip_longest(headers, line_list)))
    # Converting the list of dictionaries into a dictionary of dictionaries.
    final = {}
    for index in extracted_table:
        final[index.get('Length')] = index

    return final

def concatenate_molecules_summary(
        status, depth, input_folder_name, output_file, table_headers):
    """
    Concatenate the extracted results from the status and depth tables at the
    end of the output_file.

    Parameters
    ----------
    status : dictionary
        Dictionary of dictionaries containg the information of the status table
        extracted from the log file using the extractor function.
    depth : dictionary
        Dictionary of dictionaries containg the information of the depth table
        extracted from the log file using the extractor function.
    input_folder_name : string
        Name of directory that contains the unicycler.log file being analyzed.
    output_file : string
        Path to file used to concatenate the results of the status and depth
        tables.
    table_headers : list
        Headers of the table molecules_summary.

    Returns
    -------
    bool
        True for success, otherwise False.
    """
    # Opening outfile.
    with open(output_file, 'a') as outfile:
        # Creting an object to write the csv file.
        writer = csv.DictWriter(outfile, fieldnames=table_headers)
        # Both status and depth are dictionaries of dictionaries, and the keys
        # are the lenght of the molecules. status and depth were created with
        # the extractor function. Check the documentation of the extractor
        # function for more details.
        for key in status:
            # Some of the information in the status and depth tables is
            # redundant, like segment or length. Additionaly, some molecules
            # don't have information in the depth table. Therefore, if the
            # length of the molecule (key) is present in both depth and status
            # tables, the program will extract information from depth and
            # status. Oherwise, the columns corresponding to the part of the
            # depth table will be fill with None.
            #
            # If Length (key) from status table is in depth table, get
            # information from both status and depth tables. Otherwise, put
            # 'None' in Depth.
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
            # Appending relevant information in the outfile.
            writer.writerow(relevant)

    return True


def molecules_summary(file_addresses, output_folder):
    """
    Make a summary using the status and depth tables of unicycler.log files.

    Creates a csv file with relevant information retrieved from all the
    unicycler.log files contained in the primary subfolders of a given
    directory.

    Parameters
    ----------
    file_addresses : list object
        List of addresses that will be used to open the unicycler.log files.

    output_folder : string
        Path to the directory that will be used to save the output
        molecules_summary.csv file.

    Returns
    -------
    bool
        True for success, otherwise False.
        If the the path to output_folder doesn't exist returns False.
    """
    # Checking if output_folder exists, otherwise return False.
    if os.path.exists(output_folder):
        output_path = os.path.join(output_folder, "molecules_summary.csv")
    else:
        return False
    # Opening the outfile (molecules_summary) to save the summary and make the
    # headers of the table.
    with open(output_path, 'a') as outfile:
        # Headers of the molecules_summary table. These headers correspond to
        # the ones found in unicycler.log.
        fieldnames = [
            'Folder_name','Component', 'Segments', 'Links', 'Length', 'N50',
            'Longest_segment', 'Status', 'Depth', 'Starting_gene', 'Position',
            'Strand', 'Identity', 'Coverage']
        # Creting an object to write the csv file.
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        # Writing header's table in outfile.
        writer.writeheader()
    # Iterating over each directory.
    for _, address in enumerate(file_addresses):
        # Dictionaries to save the status and depth tables.
        status = {}
        depth = {}
        # Getting path to folder containg input file
        folder_path = os.path.dirname(address)
        # Getting name of folder containing the input file
        folder_name = os.path.basename(folder_path)
        # Opening log file.
        with open(address, 'r') as log_file:
            # Iterating over log file.
            for line in log_file:
                # If 'Component' and 'Status' are found in line, then extract
                # table status.
                if re.search('^Component.*Status', line):
                    # Convert header 'Longest segment' into 'Longest_segment'.
                    headers = re.sub(
                        'Longest segment', 'Longest_segment', line)
                    # Replace line's spaces with tabs and convert headers into
                    # a list.
                    headers = re.sub(
                        '\\s+', '\t', headers.strip()).split('\t')
                    # Extract table status using the extractor function.
                    status = extractor(log_file, headers)
                # If 'Segment' and 'Depth' are found in line extract table
                # depth.
                if re.search('^Segment.*Depth', line):
                    # Convert header 'Starting gene' into 'Starting_gene'.
                    headers = re.sub(
                        'Starting gene', 'Starting_gene', line)
                    # Replace line's spaces with tabs and convert headers into
                    # a list.
                    headers = re.sub(
                        '\\s+', '\t', headers.strip()).split('\t')
                    # Extract table depth using the extractor function.
                    depth = extractor(log_file, headers)
            # Saving (concatenate) info from status and depth variables into
            # outfile.
            concatenate_molecules_summary(
                status, depth, folder_name, output_path, fieldnames)

    return True


def extract_best_k_mer(table):
    """
    Read rows of a table from an input file until it finds best K-mer and
    converts it into a list.

    Parameters
    ----------
    table : iteratable object
        File opened for reading that will be processed

    Returns
    -------
    list (k-mer, Contigs, Dead_ends, Score)
        List of strings containing the row with the best K-mer. 
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
    convert them into a list.

    Parameters
    ----------
    table : iteratable object
        File opened for reading that will be processed.

    Returns
    ------
    list (
        Total_read_count, Fully_aligned_reads, Partially_aligned_reads,
        Unaligned_reads, Total_bases_aligned, Mean_alignment_identity)
        List containing the information of the table.
    """
    # List to save info.
    alignment_summary_list = []
    for alignment_summary in table:
        # Break when reach the end of the table. log files separeates the
        # end of a table with '\n'.
        if alignment_summary == '\n':
            break
        # If it find a row with '--' ignore and continue.
        if '--' in alignment_summary:
            continue
        # Replacing single line's spaces with '_'.
        alignment_summary = re.sub(r'([^\s])(\s)([^\s])',
                                   r'\1_\3', alignment_summary)
        # Replacing multiple line's spaces with '\t' and conver line in list.
        alignment_summary = (re.sub(r'\s+', r'\t',
                                    alignment_summary)).split('\t')
        # Extracting relevant data. The second column is the one with values;
        # therefore use index 1 of aligment_summary
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
        Name of directory that contains the unicycler.log file being analyzed.
    output_file : string
        Path to file used to concatenate the results of the K-mer and Read
        alignment summary tables.
    table_headers : list
        Headers of the table assemblies_summary.
    """
    # Opening output_file.
    with open(output_file, "a") as outfile:
        # Creating an object to write the csv file.
        writer = csv.DictWriter(outfile, fieldnames=table_headers)
        # Compaling information.
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
        # Write relevant info in outfile.
        writer.writerow(relevant)

def assemblies_summary(file_addresses, output_folder):
    """
    Make a summary using the K-mer and alignment tables of unicycler.log files.

    Creates a csv file with relevant information retrieved from all the
    unicycler.log files contained in the primary subfolders of a given
    directory.

    Parameters
    ----------
    file_addresses : list object
        List of addresses that will be used to open the unicycler.log files.
    output_folder : string
        Path to save the output assemblies_summary.csv file.

    Returns
    -------
    bool
        True for success, otherwise False.
        If the the path to output_folder doesn't exist returns False.
    """
    # Checking if output_folder exists, otherwise return False.
    if os.path.exists(output_folder):
        output_path = os.path.join(output_folder, "assemblies_summary.csv")
    else:
        return False
    # Opening the outfile to save the summary.
    with open(output_path, 'a') as outfile:
        # Headers of table assemblies summary.
        fieldnames = ['Folder_name', 'K-mer_best','Contigs_best', 'Dead_ends_best',
                      'Score_best', 'Total_read_count', 'Fully_aligned_reads',
                      'Partially_aligned_reads', 'Unaligned_reads',
                      'Total_bases_aligned', 'Mean_alignment_identity']
        # Creting an object to write the csv file.
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        # Writing header of table in outfile.
        writer.writeheader()
    # Iterating over each directory.
    for _, address in enumerate(file_addresses):
        # Getting path to folder containg input file.
        folder_path = os.path.dirname(address)
        # Getting name of folder containing the input file.
        folder_name = os.path.basename(folder_path)
        # List to save best K-mer.
        best = []
        # Opening log file.
        with open(address, 'r') as log_file:
            # Iterating over log file.
            for line in log_file:
                # If 'K-mer', 'Contigs', 'Dead ends' and 'Score' are found in
                # line, extract table.
                if re.search('^K-mer.*Contigs.*Dead ends.*Score', line):
                    best = extract_best_k_mer(log_file)
                # If 'Read alignment summary' in line extract table.
                if re.search('Read alignment summary', line):
                    # List to save info.
                    alignment_summary_list = extract_alignment_summary(
                        log_file)
        # If the leng of best is zero, it means that the file doesn't have the
        # table k-mer. Therefore, we don't need the info of that file.
        if len(best) == 0:
            continue
        concatenate_assemblies_summary(
            best, alignment_summary_list, folder_name, output_path,
            fieldnames)

    return True

def main():
    """Parse and extract information from unicycler.log files"""
    # Getting input from user.
    args = user_input()
    # Getting path to the directory that carries the subdirectories to be
    # analyzed.
    input_directory = args.input
    # Getting path to output directory.
    output_directory = args.output
    # Getting a list of the file addresses.
    file_addresses = get_file_paths(input_directory, 'unicycler.log')
    # Making new files with extracted information from unicycler.log files.
    molecules_summary(file_addresses, output_directory)
    assemblies_summary(file_addresses, output_directory)
    # If everything went OK print message.
    print(
        "The assemblies_summary.csv and molecules_summary.csv files are in:")
    print(f"{os.path.abspath(output_directory)}")
    print("log_parser.py is done!")

    sys.exit(0)

if __name__ == '__main__':
    main()
