import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from BIOINFO_M260B.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref

def process_read_pair(front, back, ref):
    min_mismatches = len(front) + 1
    min_mismatch_location = -1

    for i in range(len(ref) - len(front)):
        mismatches = [1 if front[j] != ref[i + j] else 0 for j in range(len(front))]
        n_mismatches = sum(mismatches)

        if n_mismatches < min_mismatches:
            min_mismatches = n_mismatches
            min_mismatch_location = i

    min_tail_mismatches = len(back) + 1
    min_tail_location = -1

    for k in range(130,170):
        endpoint = min_mismatch_location + k
        if endpoint + len(back) > len(ref):
            break
    # for endpoint in range(len(ref) - len(back)):
        mismatches = [1 if back[j] != ref[endpoint + j] else 0 for j in range(len(back))]
        n_mismatches = sum(mismatches)
        
        # for num in mismatches:
        #     print(str(num), end='')
        # print('\t' + str(n_mismatches))

        if n_mismatches < min_tail_mismatches:
            min_tail_mismatches = n_mismatches
            min_tail_location = endpoint

    return min_mismatches, min_mismatch_location, min_tail_mismatches, min_tail_location

def trivial_algorithm(paired_end_reads, ref):
    """
    This is a functional aligner, but it's a huge simplification that
    generate a LOT of potential bugs.  It's also very slow.

    Read the spec carefully; consider how the paired-end reads are
    generated, and ideally, write your own algorithm
    instead of trying to tweak this one (which isn't very good).

    :param paired_end_reads: Paired-end reads generated from read_reads
    :param ref: A reference genome generated from read_reference
    :return: 2 lists:
                1) a list of alignment locations for each read (all_alignment_locations).
                    The list gives the starting position of the minimum-mismatch alignment of both reads.
                2) a list of the paired-end reads set so that both reads are in their optimal orientation
                   with respect to the reference genome.
    """
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    start = time.clock()

    debug = False
    matched = 0

    # Main loop through all reads
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        found = False

        # Print the remaining time
        if count % 10 == 0:
            time_passed = (time.clock()-start)/60
            print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print('Approximately {:.3} minutes remaining'.format(remaining_time))

        # Look at the read_pair FRONT way
        front = read_pair[0]
        back = read_pair[1][::-1]

        min_mismatches, min_mismatch_location, min_tail_mismatches, min_tail_location = process_read_pair(front, back, ref)

        if (min_mismatches + min_tail_mismatches < 10):
            matched += 1
            if (debug):
                print("FRONT WAY")
                print("front_mismatch: {}  location: {}".format(min_mismatches, min_mismatch_location))
                print("back_mismatch: {}  distance: {}".format(min_tail_mismatches, min_tail_location))

            # append this match only if threshold is met
            read_alignment_locations.append(min_mismatch_location)
            output_read_pair.append(front)
            read_alignment_locations.append(min_tail_location)
            output_read_pair.append(back)
            found = True
        
        # if front way didn't find a satisfactory match:
        # else: 
        #     # Look at read_pair the BACK way
        #     front = read_pair[1]
        #     back = read_pair[0][::-1]

        #     min_mismatches, min_mismatch_location, min_tail_mismatches, min_tail_location = process_read_pair(front, back, ref)

        #     if (min_mismatches + min_tail_mismatches < 10):
        #         matched += 1
        #         if (debug):
        #             print("BACK WAY")
        #             print("front_mismatch: {}  location: {}".format(min_mismatches, min_mismatch_location))
        #             print("back_mismatch: {}  distance: {}".format(min_tail_mismatches, min_tail_location))
                
        #         read_alignment_locations.append(min_mismatch_location)
        #         output_read_pair.append(front)
        #         read_alignment_locations.append(min_tail_location)
        #         output_read_pair.append(back)
        
        if (found):
            all_read_alignment_locations.append(read_alignment_locations)
            output_read_pairs.append(output_read_pair)

    # after main loop
    print('Matched {} reads'.format(matched))
    return all_read_alignment_locations, output_read_pairs


if __name__ == "__main__":
    data_folder = 'practice_W_1'
    input_folder = join('../data/', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    
    # This will take a while; you can use an array slice for example:
    #
    #   input_reads = reads[:300]
    #
    # to generate some data quickly.

    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    alignments, reads = trivial_algorithm(input_reads, reference)

    # print(alignments)
    # print(reads)

    # Final outputs
    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
    print("saving to " + output_fn)
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
