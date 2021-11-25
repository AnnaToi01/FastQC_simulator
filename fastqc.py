import argparse
from pathlib import Path
import sys
import logging
import time
import datetime
from calculations import SetLogger
import warnings
warnings.simplefilter('ignore')

if __name__ == '__main__':
    # create argparse logger
    SetLogger(logger_name='argparse')
    logger = logging.getLogger('argparse')
    # argparse main parameters
    parser = argparse.ArgumentParser(description='FastQC simulator', epilog='Enjoy the program! :)',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type=str, help='Path to fastq file', required=True, metavar='')
    parser.add_argument('-o', '--output', type=str, help='Path to output folder for storing results',
                        required=True, metavar='')
    parser.add_argument('-a', '--adapters', type=str, help='Path to file with adapters. Default: ./adapters.txt',
                        default='./adapters.txt', metavar='')
    args = parser.parse_args()
    if not Path(args.input).exists():
        logger.warning('File {} does not exist!'.format(args.input))
        logger.info('Abort calculations. Specify correct path to file.')
        sys.exit()
    if not Path(args.adapters).exists():
        logger.warning('File {} does not exist!'.format(args.adapters))
        logger.info('Abort calculations. Specify correct path to file.')
        sys.exit()
    logger.info('Start fastq quality control...')
    start_time = time.process_time()
    Path(args.output).mkdir(parents=True, exist_ok=True)
    # insert all functions for calculations and plotting here
    end_time = time.process_time()
    total_time = str(datetime.timedelta(seconds=end_time-start_time))
    logger.info('End fastq quality_control')
    logger.info(f'Total time for FastQC analysis: {total_time}')
