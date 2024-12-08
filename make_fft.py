#!/bin/env python3
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


def make_fft(dosbin, ampbin, phabin):
    with open(dosbin, 'rb') as stream_dos, open(ampbin, 'wb') as amp, open(phabin, 'wb') as pha:
        lenob, = struct.unpack('f', stream_dos.read(4))
        lenob = int(lenob)
        obs = struct.unpack(f'{lenob}f', stream_dos.read(4 * lenob))
        es = []
        dos = []
        while True:
            e = stream_dos.read(4)
            if not e:
                break
            e, = struct.unpack('f', e)
            d = struct.unpack(f'{lenob}f', stream_dos.read(4 * lenob))
            es.append(e)
            dos.append(d)
        logging.info(f'Got {len(es)} energies {es[0]}...{es[-1]}')
        logging.info(f'Got {len(obs)} obs {obs[0]}...{obs[-1]}')
        need_header = True
        lenfft = lenob // 2 + 1
        for i, e in enumerate(tqdm.tqdm(es)):
            f = np.fft.rfft(dos[i])
            k = np.array(range(len(f))) / (obs[-1] - obs[0])
            famp = np.log(np.abs(f))
            fpha = np.mod(np.angle(f), np.pi * 2)
            part = len(f) # int(5 * len(f) / k[-1])
            f = f[:part]
            k = k[:part]
            famp = famp[:part]
            fpha = fpha[:part]
            if need_header:
                amp.write(struct.pack(f'{len(k)+1}f', len(k), *k))
                pha.write(struct.pack(f'{len(k)+1}f', len(k), *k))
                need_header = False
            amp.write(struct.pack(f'{len(famp)+1}f', e, *famp))
            pha.write(struct.pack(f'{len(fpha)+1}f', e, *fpha))
        logging.info('Done')

if __name__ == '__main__':
    wxx='w02'
    ampbin=f'{wxx}/{wxx}a_fftdos_amp.bin'
    phabin=f'{wxx}/{wxx}a_fftdos_pha.bin'
    logging.info(f'Creating {ampbin} {phabin}')
    make_fft(f'{wxx}/{wxx}a_dos.bin', ampbin, phabin)
