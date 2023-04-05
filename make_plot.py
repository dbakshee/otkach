#!/bin/env python3
import argparse
import fnmatch
import logging
import os
import struct
import sqlite3 as db
import tqdm

def config_logging():
    logging.basicConfig(
        format='\x1b[1;m%(levelname)1.1s\x1b[0m %(message)s',
        level=logging.INFO,
        datefmt='%H:%M:%S',
    )


def as_ssv(array_of_float, prec):
    ans = []
    for a in array_of_float:
        ans.append(f"{a:.{prec}g}")
    return " ".join(ans)


def main1():
    with db.connect(DB) as con, open('dos.dat', 'w') as o:
        c = con.cursor()
        #d = c.execute("SELECT * FROM MAGBREAK WHERE e == 11.0").fetchall()
        es = [rec[0]
              for rec in set(c.execute("SELECT e FROM MAGBREAK").fetchall())]
        one_ob = None
        for rowno, e in enumerate(sorted(es)):
            if e > 14:
                break
            if e == 13.25:
                continue
            logging.info(f"Processing row {e}")
            req = f"SELECT ob, dos FROM MAGBREAK WHERE e == {e}"
            d = c.execute(req).fetchall()
            ob = [rec[0] for rec in sorted(d)]
            z = [rec[1] for rec in sorted(d)]
            if one_ob is None:
                one_ob = ob
            elif one_ob != ob:
                logging.error(
                    f"For e={e} len(ob) is {len(ob)} exd {len(one_ob)}")
            if rowno == 0:
                o.write(
                    f'# 0:{one_ob[0]}, 1:{one_ob[1]}, ..., {len(one_ob)-1}:{one_ob[-1]}\n')
            o.write(f'{as_ssv(z, 8)}\n')


def main():
    with db.connect(DB) as con, open('dos_t.dat', 'w') as o:
        c = con.cursor()
        obs = [rec[0]
               for rec in set(c.execute("SELECT ob FROM MAGBREAK").fetchall())]
        for rowno, ob in enumerate(sorted(obs)):
            req = f"SELECT e, dos FROM MAGBREAK WHERE ob == {ob}"
            d = sorted(c.execute(req).fetchall())
            z = [rec[1] for rec in d]
            o.write(f'{as_ssv([ob] + z, 8)}\n')


def make_bin_map(dbfile, prop, binout):
    with db.connect(dbfile) as con, open(binout, 'wb') as o:
        c = con.cursor()
        obs = sorted([rec[0] for rec in c.execute(
            "SELECT DISTINCT ob FROM MAGBREAK").fetchall()])
        es = sorted([rec[0] for rec in c.execute(
            "SELECT DISTINCT e FROM MAGBREAK").fetchall()])
        logging.info(f'Got {len(es)} energies {es[0]}...{es[-1]}')
        logging.info(f'Got {len(obs)} obs {obs[0]}...{obs[-1]}')
        o.write(struct.pack(f'{len(obs)+1}f', len(obs), *obs))
        for e in tqdm.tqdm(es):
            req = f"SELECT ob, {prop} FROM MAGBREAK WHERE e == {e}"
            row = sorted(c.execute(req).fetchall())
            vals = list(zip(*row))[1]
            try:
                o.write(struct.pack(f'{len(obs)+1}f', e, *vals))
            except Exception as e:
                logging.critical(f"FATAL ERROR for e={e} len(vals)={len(vals)}")
                raise e


def check_db(dbfile):
    with db.connect(dbfile) as con:
        c = con.cursor()
        es = sorted(c.execute("SELECT DISTINCT e FROM MAGBREAK").fetchall())
        obs = sorted(c.execute("SELECT DISTINCT ob FROM MAGBREAK").fetchall())
        logging.info(f'Got {len(es)} energies times {len(obs)} obs')
        for e, in [(8,)]:
            req = f"SELECT ob, dos FROM MAGBREAK WHERE e == {e}"
            row = sorted(c.execute(req).fetchall())
            vals = list(zip(*row))[1]
            logging.info(f"{e}: {len(vals)} values")
            assert len(vals) == len(obs)
    exit(1)

if __name__ == '__main__':
    dbfile = 'w02/w02a.db'
    config_logging()
    for prop in 'dos rxy21 rxx rxy12 r2t22 r2t11 ryy'.split():
        binfile = f'w02/w02a_{prop}.bin'
        if os.path.exists(binfile):
            logging.info(f'Skip {binfile} - already there')
            continue
        logging.info(f'Creating {binfile}')
        make_bin_map(dbfile, prop, binfile)
