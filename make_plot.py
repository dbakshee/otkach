#!/bin/env python3
import argparse
import fnmatch
import logging
import os
import sqlite3 as db

DB = 'db.db'


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
       es = [rec[0] for rec in set(c.execute("SELECT e FROM MAGBREAK").fetchall())]
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
               logging.error(f"For e={e} len(ob) is {len(ob)} exd {len(one_ob)}")
           if rowno == 0:
               o.write(f'# 0:{one_ob[0]}, 1:{one_ob[1]}, ..., {len(one_ob)-1}:{one_ob[-1]}\n')
           o.write(f'{as_ssv(z, 8)}\n')

def main():
   with db.connect(DB) as con, open('dos_t.dat', 'w') as o:
       c = con.cursor()
       obs = [rec[0] for rec in set(c.execute("SELECT ob FROM MAGBREAK").fetchall())]
       for rowno, ob in enumerate(sorted(obs)):
           req = f"SELECT e, dos FROM MAGBREAK WHERE ob == {ob}"
           d = sorted(c.execute(req).fetchall())
           z = [rec[1] for rec in d]
           o.write(f'{as_ssv([ob] + z, 8)}\n')


if __name__ == '__main__':
    config_logging()
    main()
