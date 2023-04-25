#!/bin/env python3
import argparse
import logging
import math
import numpy as np
import struct
import tqdm

logging.basicConfig(
    format='\x1b[1;m%(levelname)1.1s\x1b[0m %(message)s',
    level=logging.INFO,
    datefmt='%H:%M:%S',
)


def showbin(bin):
    with open(bin, 'rb') as b:
        width, = struct.unpack('f', b.read(4))
        width = int(width)
        logging.info(f'Width: {width}')
        header = struct.unpack(f'{width}f', b.read(4 * width))
        logging.info(f'X ({len(header)} pts: {header[0]} {header[1]} ... {header[-1]})')
        ys = []
        vs = []
        for i in range(10000000):
            y_bytes = b.read(4)
            if not y_bytes:
                break
            y, = struct.unpack('f', y_bytes)
            values = struct.unpack(f'{width}f', b.read(4 * width))
            ys.append(y)
            vs.append(values)
            if i < 5:
                logging.info(f'Y{i}={y}: {values[0]} {values[1]} ... {values[-1]}')
        logging.info(f'...')
        logging.info(f'Y{i-1}={y}: {values[0]} {values[1]} ... {values[-1]}')
        vs = np.array(vs)
        logging.info(f'Values vary from {vs.min()} to {vs.max()}')
        mean = np.mean(vs)
        std = np.std(vs)
        logging.info(f'Mean {mean} std_dev {std}')
        count_below_mean = np.sum(vs < mean)
        count_above_mean = np.sum(mean <= vs)
        logging.info(f'Total {count_below_mean+count_above_mean} pts'
                     f' = {count_below_mean} pts below + {count_above_mean} pts above mean')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fname', type=str, help='Path to bin file')
    args = parser.parse_args()
    showbin(args.fname)
