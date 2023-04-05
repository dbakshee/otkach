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
        w, = struct.unpack('f', b.read(4))
        w = int(w)
        logging.info(f'Width: {w}')
        h = struct.unpack(f'{w}f', b.read(4 * w))
        logging.info(f'X ({len(h)}: {h[0]} {h[1]} ... {h[-1]})')
        ys = []
        vs = []
        for i in range(10000000):
            y = b.read(4)
            if not y:
                break
            y, = struct.unpack('f', y)
            v = struct.unpack(f'{w}f', b.read(4 * w))
            ys.append(y)
            vs.append(v)
            if i < 5:
                logging.info(f'Y{i}={y}: {v[0]} {v[1]} ... {v[-1]}')
        #logging.info(f'Got {len(ys)} y-values {ys[0]} {ys[1]} ... {ys[-1]}')
        logging.info(f'Got {len(ys)} y-values {ys[0]:.5} {ys[1]:.5} ... {ys[-1]:5} step {(ys[-1]-ys[0])/(len(ys)-1)}')
        print([f'{y:.4}' for y in ys])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fname', type=str, help='Path to bin file')
    args = parser.parse_args()
    showbin(args.fname)
