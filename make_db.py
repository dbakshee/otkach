#!/bin/env python3
import argparse
import fnmatch
import logging
import os
import sqlite3 as db


def config_logging():
    logging.basicConfig(
        format='\x1b[1;m%(levelname)1.1s\x1b[0m %(message)s',
        level=logging.INFO,
        datefmt='%H:%M:%S',
    )


def process(f, dbc):
    data = []
    for l in open(f).readlines():
        try:
            # ob, e, dos, rxy21, rxx, rxy12, r2t22, r2t11, ryy
            data.append(tuple(map(float, l.replace(',', '').split())))
        except Exception:
            pass
    dbc.executemany('INSERT INTO MAGBREAK (ob, e, dos, rxy21, rxx, rxy12, r2t22, r2t11, ryy)'
                    ' VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)', data)
    if False:
        for lno, l in enumerate(open(f).readlines()):
            if l.startswith('#'):
                continue
            try:
                d = tuple(map(float, l.replace(',', '').split()))
                dbc.execute('INSERT INTO MAGBREAK (ob, e, dos, rxy21, rxx, rxy12, r2t22, r2t11, ryy)'
                            ' VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)', d)
            except db.IntegrityError as e:
                if 'columns ob, e are not unique' in str(e):
                    pass
                else:
                    raise e 
            except Exception as e:
                logging.error(f"Error in file {f} on line {lno}: {l}")
                raise e

def main(dbfile, pat):
    if not os.path.exists(dbfile):
        logging.info(f'Create {dbfile}')
        with db.connect(dbfile) as con:
            con.execute("""
                CREATE TABLE MAGBREAK (
                    ob DOUBLE,
                    e DOUBLE,
                    dos DOUBLE,
                    rxy21 DOUBLE,
                    rxx DOUBLE,
                    rxy12 DOUBLE,
                    r2t22 DOUBLE,
                    r2t11 DOUBLE,
                    ryy DOUBLE,
                    PRIMARY KEY (ob, e) ON CONFLICT IGNORE -- WITH (IGNORE_DUP_KEY = ON)
                );""")
            con.commit()
    with db.connect(dbfile) as con:
        for root, _, files in os.walk("."):
            if pat in root:
                logging.info(f'Processing {root}')
                for f in fnmatch.filter(files, "comp_*.dat"):
                    process(os.path.join(root, f), con.cursor())
                con.commit()


if __name__ == '__main__':
    config_logging()
    main('db.db', 'broadwell.run.sh.')
