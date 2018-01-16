#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Run CNMF-E recursively on all of detected data files.

This script can be used to process many image stacks using CNMF-E
in MATLAB via the command line, and then save out the results in a format
that can be easily read and used by python.

See help for the __main__ function below for more details about how to
run this script.

Take a look at 'analyze_cnmfe_matlab.ipynb' for an example of how
to use the results that will be generated.

Created on Nov 1 2017

@author: tamachado@stanford.edu
"""
import os
import argparse
from subprocess import call

CNMFE_ANALYSIS_FOLDER = 'source_extraction'


def get_paths(base_path, mouse_name):
    """
    Recursively search all directories within base path for files
    of type '.tif' that optionally also contain a specific substring.
    """
    continued = False
    input_folders = []
    output_folders = []
    for root, dirs, files in os.walk(base_path):

        # Check if there is any intermediate analysis in this directory
        if any(CNMFE_ANALYSIS_FOLDER in dir for dir in dirs):
            if not continued:
                print('Intermediate results exist!\n',
                      'If you do not wish to continue',
                      'press ^C to stop or enter to continue')
                print('Incomplete analysis files may cause strange behavior\n',
                      'If processing is not successful, remove old analysis\n',
                      'folders and try again.')
                _ = input()
                print('Previous results will not be recomputed',
                      'but unanalyzed data will be processed...')
                continued = True
            continue

        # Otherwise look for tiff files
        found_file = False
        for file_ in files:
            if ((mouse_name == 'none' or mouse_name in file_) and
                    ('tif' in file_)):
                    if found_file:
                        raise IOError('There can only be one tiff stack ' +
                                      'to analyze per subfolder!')
                    input_folders.append(root + os.sep + file_)
                    found_file = True
    input_folders.sort()
    for root in input_folders:
        offset = len(root.split(os.sep)[-1])
        output_folders.append(root[:-offset])

    return input_folders, output_folders


if __name__ == '__main__':

    # Parse arguments
    usage = ('Recursively analyze multipage .tif stacks or mat' +
             'files saved by CNMF-E found in subdirectories' +
             'within a provided root folder.' +
             'Then convert them to a format suitable for opening in Python.')
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('path',
                        help='Root folder containing datasets to process.')
    parser.add_argument('--name',
                        default='none',
                        help='Substring in all image stacks files to process.')
    parser.add_argument('--matlab',
                        default='~/Matlab/bin/matlab -softwareopengl',
                        help='Path to MATLAB executable')
    parser.add_argument('--cnmf',
                        default='run_cnmfe_matlab',
                        help='Name of CNMF script to run in MATLAB')
    parser.add_argument('--fixer',
                        default='fix_mat_files',
                        help='Name of postprocessing script to run in MATLAB')
    args = parser.parse_args()

    if os.path.exists(args.path):
        args.path = os.path.abspath(args.path)
    else:
        raise IOError('Path to stack or output is invalid!')

    # Run CNMF-E on each dataset
    input_folders, output_folders = get_paths(args.path, args.name)
    for idx, (inp, out) in enumerate(zip(input_folders, output_folders)):
        print('Processing ', idx)
        print('in:  ', inp)
        print('out: ', out)

        # Start MATLAB
        call(args.matlab + ' -nodesktop -nosplash -r \"' + args.cnmf +
             '(\'' + inp + '\'); exit;\"', shell=True)

    # Reformat all output files
    print('Postprocess .mat files in MATLAB...')
    call(args.matlab + ' -nodesktop -nosplash -r \"' + args.fixer +
         '(\'' + args.path + '\'); exit;\"', shell=True)

    print('All files analyzed.')
