#!/usr/bin/env python
# Based on Greg Caporaso's fastqz manifest maker script
# Edited to handle TGen North file patterns

import glob
import os.path
import re

def check_regx_list(regx_list, input_str):
    for regx_str in regx_list:
        if re.search(regx_str, input_str):
            return True
    return False

def fq_manifestor(input_dir,
                  output_fp,
                  fq_extensions=['fastq.gz', 'fq.gz'],
                  split_pattern='_',
                  f_read_pattern='_R1_',
                  r_read_pattern='_R2_',
                  skip_amount=0, # For most Tgen files this should be changed to 1
                  close_pattern_list=['S\d+','L\d+'],
                  filter_pattern=None,
                  verbose=True):

    input_dir = os.path.abspath(input_dir)
    if verbose: print("Searching directory: %s" % input_dir)

    fq_filepaths = []
    for fq_extension in fq_extensions:
        fq_filepaths += glob.glob('%s/**/*.%s' % (input_dir, fq_extension),
                                 recursive = True)
    if filter_pattern is not None:
        fq_filepaths = [fp for fp in fq_filepaths if filter_pattern in fp]

    n_fq_filepaths = len(fq_filepaths)
    if verbose: print('Found %d fastq files.' % n_fq_filepaths)

    sids_to_fps = {}

    for fq_filepath in fq_filepaths:
        fq_filename = os.path.basename(fq_filepath)
        sid_fields = re.split(split_pattern, fq_filename)

        sid = ""
        skip_count = 0
        if len(sid_fields) == 1:
            raise ValueError('Sample ID not found in file: %s' % fq_filepath)
        else:
            for sid_f in sid_fields:
                if skip_count < skip_amount:
                    skip_count+=1
                    continue
                if check_regx_list(close_pattern_list, sid_f):
                    break
                sid = sid + "_" + sid_f
        sid = sid.strip("_") # Remove start/end '_'

        if bool(re.search(f_read_pattern, fq_filename)):
            forward = True
        elif bool(re.search(r_read_pattern, fq_filename)):
            forward = False
        else:
            raise ValueError('Forward/reverse patterns not found in file: %s' % fq_filepath)

        try:
            sid_fps = sids_to_fps[sid]
        except KeyError:
            sid_fps = [None, None]

        if forward:
            sid_fps[0] = fq_filepath
        else:
            sid_fps[1] = fq_filepath
        sids_to_fps[sid] = sid_fps

    lines = ['sample-id\tforward-absolute-filepath\treverse-absolute-filepath']
    for sid, (fwd_fq_filepath, rev_fq_filepath) in sids_to_fps.items():
        if fwd_fq_filepath is None:
            raise ValueError('Missing forward read for sample: %s' % s)

        if rev_fq_filepath is None:
            raise ValueError('Missing reverse read for sample: %s' % s)

        lines.append('%s\t%s\t%s' % (sid, fwd_fq_filepath, rev_fq_filepath))

    if (len(lines) - 1) != (n_fq_filepaths / 2) and verbose:
        print("\n** WARNING**: "
              "The number of manifest records doesn't align with the number of "
              "fastq files that were found. It's possible that the match "
              "patterns aren't working correctly for your files. These can "
              "be customized when using the API.\n")

    with open(output_fp, 'w') as of:
        of.write('\n'.join(lines))
        of.write('\n')