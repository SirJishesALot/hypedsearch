'''
exec.py

Author: Zachary McGrath 
Date: 6 April 2020

Executor for the program
In charge of the flow of the program
'''
import os
from src.alignment.alignment import TIME_LOG_FILE
from src import identification
from src.postprocessing import summary, review

import multiprocessing as mp

def run(args: dict) -> None:
    '''
    Executing function for the program

    Inputs:
        args:   object arguments from main. Should be validated in main. Attributes of args:
            spectra_folder:             (str) full path the the directory containing all spectra files
            database_file:              (str) full path to the .fasta database file
            output_dir:                 (str) full path the the directory to save output to
            min_peptide_len:            (int) minimum peptide length to consider
            max_peptide_len:            (int) maximum peptide length to consider
            tolerance:                  (int) the ppm tolerance to allow in search
            precursor_tolerance:        (int) the ppm tolerance to allow when matching precursors
            peak_filter:                (int) the number of peaks to filter by 
            relative_abundance_filter:  (float) the percentage of the total abundance a peak must
                                            be to pass the filter
            digest:                     (str) the digest performed
            missed_cleavages:           (int) the number of missed cleavages allowed in digest
            verbose:                    (bool) extra printing
            cores:                      (int) the number of cores allowed to use
            n:                          (int) the number of alignments to keep per spectrum
            DEBUG:                      (bool) debuging print messages. Default=False
    Outputs:
        None
    '''
    # get all the spectra file names
    spectra_files = []
    for (root, _, filenames) in os.walk(args['spectra_folder']):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))

    # make sure cores is: 1 <= cores <= cpu cores
    cores = max(1, args['cores'])
    cores = min(cores, mp.cpu_count() - 1)

    # timing data:
    TIME_LOG_FILE = 'timelog.txt'
    o = open(TIME_LOG_FILE, 'w')
    o.write('')
    o.close()
    a = open('metadata.txt', 'w')
    a.write('')
    matched_spectra = identification.id_spectra(
        spectra_files, args['database_file'], 
        min_peptide_len=args['min_peptide_len'], 
        max_peptide_len=args['max_peptide_len'], 
        ppm_tolerance=args['tolerance'], 
        precursor_tolerance=args['precursor_tolerance'],
        peak_filter=args['peak_filter'],
        relative_abundance_filter=args['relative_abundance_filter'],
        digest=args['digest'], 
        n=args['n'] * 10,
        verbose=True, 
        DEBUG=args['DEBUG'], 
        cores=cores,
        truth_set=args['truth_set'], 
        output_dir=args['output_dir']
    )

    print('\nFinished search. Writting results to {}...'.format(args['output_dir']))

    # matched_spectra = review.tie_breaker(matched_spectra, '', args['n'])

    summary.generate(matched_spectra, args['output_dir'])
    
    # create averages for everything in the timelog file
    non_hybrid_refine_time = []
    hybrid_refine_time = []
    total_scoring_time = []
    total_alignment_time = []
    spin_time = 0
    with open(TIME_LOG_FILE, 'r') as o:
        for line in o:
            split = line.split(':')
            if "Non-hybrid refine time" == split[0]:
                value = split[1]
                non_hybrid_refine_time.append(value)
            elif "Hybrid refine time" == split[0]:
                value = split[1]
                hybrid_refine_time.append(value)
            elif "Total scoring time" == split[0]:
                value = split[1]
                total_scoring_time.append(value)
            elif "alignment time" == split[0]:
                value = split[1]
                total_alignment_time.append(value)
            elif "Time to spin up cores" == split[0]:
                value = split[1]
                spin_time = value

    avg_hybrid_refine_time = 0
    avg_non_hybrid_refine_time = 0
    avg_total_scoring_time = 0
    avg_total_alignment_time = 0
    #empty the current text file
    # o = open(TIME_LOG_FILE, 'w')
    # o.write('')
    # o.close()
    o = open(TIME_LOG_FILE,'a')
    count = 1
    for i in hybrid_refine_time:
        avg_hybrid_refine_time = avg_hybrid_refine_time + float(i)
        count = count + 1
    o.write('average hybrid refine time: ' + str(avg_hybrid_refine_time / count) + '\n')
    count = 1
    for i in non_hybrid_refine_time:
        avg_non_hybrid_refine_time = avg_non_hybrid_refine_time + float(i)
        count = count + 1
    o.write('average non hybrid refine time: ' + str(avg_non_hybrid_refine_time / count) + '\n')
    count = 1
    for i in total_scoring_time:
        avg_total_scoring_time = avg_total_scoring_time + float(i)
        count = count + 1
    o.write('average total scoring time ' + str(avg_total_scoring_time / count) + '\n')
    count = 1
    for i in total_alignment_time:
        avg_total_alignment_time = avg_total_alignment_time + float(i)
        count = count + 1
    o.write('average total alignment time ' + str(avg_total_alignment_time / count) + '\n')
    count = 1

    o.write('time to spin up cores: ' + str(spin_time) + '\n')

    o.close()
    # Averaging out metadata
    non_hybrid_b = []
    non_hybrid_y = []
    non_hybrid_total = []
    hybrid_b = []
    hybrid_y = []
    hybrid_total = []

    with open('metadata.txt', 'r') as m:
        for line in m:
            split = line.split(':')
            if (split[0] == 'average non_hybrid b score'):
                non_hybrid_b.append(split[1])
            elif (split[0] == 'average non_hybrid y score'):
                non_hybrid_y.append(split[1])
            elif (split[0] == 'average non_hybrid total score'):
                non_hybrid_total.append(split[1])
            elif (split[0] == 'average hybrid total score'):
                hybrid_total.append(split[1])
            elif (split[0] == 'average hybrid b score'):
                hybrid_b.append(split[1])
            elif (split[0] == 'average hybrid y score'):
                hybrid_y.append(split[1])

    m.close()
    file1 = open('metadata.txt', 'w')
    file1.write('')
    file1.close()
    avg_non_hybrid_b = 0
    avg_non_hybrid_y = 0
    avg_non_hybrid_tot = 0
    avg_hybrid_b = 0
    avg_hybrid_y = 0
    avg_hybrid_tot = 0
    with open('metadata.txt', 'a') as m:
        count = 1
        for b in non_hybrid_b:
            avg_non_hybrid_b = avg_non_hybrid_b + float(b)
            count = count + 1
        m.write('Average non_hybrid b score: ' + str(avg_non_hybrid_b / count) + '\n')

        count = 1
        for y in non_hybrid_y:
            avg_non_hybrid_y = avg_non_hybrid_y + float(y)
            count = count + 1
        m.write('Average non_hybrid y score: ' + str(avg_non_hybrid_y / count) + '\n')

        count = 1
        for t in non_hybrid_total:
            avg_non_hybrid_tot = avg_non_hybrid_tot + float(t)
            count = count + 1
        m.write('Average non_hybrid total score: ' + str(avg_non_hybrid_tot / count) + '\n')

        count = 1
        for b in hybrid_b:
            avg_hybrid_b = avg_hybrid_b + float(b)
            count = count + 1
        m.write('Average hybrid b score: ' + str(avg_hybrid_b / count) + '\n')

        count = 1
        for y in hybrid_y:
            avg_hybrid_y = avg_hybrid_y + float(y)
            count = count + 1
        m.write('Average hybrid y score: ' + str(avg_hybrid_y / count) + '\n')

        count = 1
        for t in hybrid_total:
            avg_hybrid_tot = avg_hybrid_tot + float(t)
            count = count + 1
        m.write('Average hybrid total score: ' + str(avg_hybrid_tot / count) + '\n')
