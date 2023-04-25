#!/bin/env python3
from cmath import exp
from matplotlib import pyplot, backends
from types import SimpleNamespace as _SimpleNamespace
import argparse
import json
import kwant
import logging
import matplotlib as mpl
import numpy as np
import os
import re
import struct
import sys
import tqdm


class SimpleNamespace(_SimpleNamespace):
    def __str__(self):
        return super().__str__().replace('SimpleNamespace(', 'SN(')


def fname_XXX(n):
    assert 'XXX' in n
    fn = 0
    while True:
        fname = n.replace('XXX', '{:03d}').format(fn)
        if not os.path.exists(fname):
            return fname
        fn += 1


def new_file_with_header_XXX(fname, header):
    fname = fname_XXX(fname)
    open(fname, 'w').write(header + '\n')
    return fname


def get_beg_end(n_items, me, team):
    items_per_worker = n_items // team
    n_big_workers = n_items % team
    n_small_workers = team - n_big_workers
    if me < n_big_workers:
        # i am a big worker
        beg = me * (items_per_worker + 1)
        end = beg + items_per_worker + 1
    else:
        beg = (n_big_workers * (items_per_worker + 1)
               + (me - n_big_workers) * items_per_worker)
        end = beg + items_per_worker
    end = min(end, n_items)
    return beg, end

# ------------------------------------------------------------------------------


class Main(object):
    def __init__(self, args):
        self.args = args
        self.SITES = dict()
        self.POT = dict()
        self.BUMPS = self.bump_locations(args)
        self.syst = self.make_system(args)
        self.params = SimpleNamespace(
            vr=self.args.vr,
            w0=self.args.w0,
            salt=self.args.salt,
            disorder=self.args.disorder,
        )

    def plot(self, file=None):
        logging.info(f'plot: Collecting potential')
        pot = [self.syst.hamiltonian(i, i, self.params) for i, _ in enumerate(self.syst.sites)]
        logging.info(f'plot: Plotting potential')
        kwant.plot(self.syst, site_symbol='o', cmap=mpl.colormaps['seismic'], site_color=pot, site_size=0.4, hop_lw=0, file=file)

    def plot_bin(self, file='last-model.bin'):
        W, L, a = self.args.W, self.args.L, self.args.a
        logging.info(f'plot: Collecting potential')
        pot = np.ndarray((W, L), float)
        for i, s in enumerate(self.syst.sites):
            x, y = int(s.pos[0] / a), int(s.pos[1] / a)
            pot[x, y] = self.syst.hamiltonian(i, i, self.params)
        logging.info(f'plot: Plotting potential')
        with open(file, 'wb') as f:
            f.write(struct.pack(f'{1+W}f', W, *(np.arange(W) * a)))
            for x in range(L):
                f.write(struct.pack(f'{1+W}f', x * a, *pot[x,:]))

    def bump_locations(self, args):
        L, W, a, overlap = args.L, args.W, args.a, 80 * 2
        bargs = [(0, 0), 80, (0-overlap, 0-overlap),
                 (L * a+overlap, W * a+overlap), args.salt]
        return self.bump_locations_impl(*bargs)

    def bump_locations_impl(self, o, ba, c, d, salt):
        """
        compute bump locations starting at `o`
        with lattice constant `ba` in rect(`c`, `d`)
        """
        b1 = (ba, 0)
        b2 = (ba / 2.0, np.sqrt(3.0)/2.0 * ba)
        cx, cy = c
        dx, dy = d
        kmin = int((dx - cx) / 2. / ba) - 100
        kmax = int((dx - cx) / ba) + 100
        lmin = 0 - 100
        lmax = int((dy - cy) / (ba * np.sqrt(3) / 2.0)) + 100
        ans = dict()
        for k in range(kmin, kmax):
            for l in range(lmin, lmax):
                x, y = b1[0] * k + b2[0] * l, b1[1] * k + b2[1] * l
                if not cx <= x <= dx:
                    continue
                if not cy <= y <= dy:
                    continue
                ans[(x,  y,)] = kwant.digest.uniform(
                    repr((x, y,)) + salt) - 0.5
        return ans

    def get_random_shift(self, site):
        x, y = site.pos
        if (x, y) not in self.SITES:
            modulation = 0.0
            for bxy, vr in self.BUMPS.items():  # BUMPS have randomized modulation
                bx, by = bxy
                dist2 = (x - bx)**2 + (y - by)**2
                #    modulation += np.exp(-dist2/(40*40)) * vr
                #    modulation += np.exp(-dist2/(35*35)) * vr
                modulation += np.exp(-dist2/(30*30)) * vr
            self.SITES[(x, y)] = modulation
        return self.SITES[(x, y)]

    def make_system(self, args):
        a, t = args.a, args.t
        W, L, WL, Wtop = args.W, args.L, args.Wbot, args.Wtop
        # Start with an empty tight-binding system and a single square lattice.
        # ------------------
        # 3                2   (width = Wtop)
        # ---            ---
        #    |  syst    |
        # ---            ---
        # 0                1   (width = WL)
        # ------------------
        # - `a` - lattice constant
        # - `t` - tight-binding coupling (hopping amplitude)
        # - `L` - lenth of the rectangle
        # - `W` - height of the rectangle (width of the sample)
        # - `Wtop` - width of the connecting channels
        lat = kwant.lattice.square(a, norbs=1)

        def hopping(site_i, site_j, params):
            xi, yi = site_i.pos
            xj, yj = site_j.pos
            # raise ValueError('PARAMS is {}, {}, {}'.format(type(params), repr(params), str(params)))
            hmag = 6.2831853 * params.Bt / 4135.67  # Ph_0=4.13567*10^3Tesla*nm^2
            return -t * exp(-0.5j * hmag * (xi - xj) * (yi + yj))

        def potential(site, params):
            w0 = params.w0
            salt = params.salt
            vr = params.vr
            x, y = site.pos
            #    v0 = 0.78 * w0 #L=160nm
            #    L = 160
            #    v0 = 1.386 * w0 #L=120nm
            #    L = 120
            #    v0 = 2.5 * w0 #L=104nm 2cos deltaU=1meV
            v0 = 3.116 * w0  # L=80nm
            #    L = 104
            L = 80
            tpl = 2 * np.pi / L
            g1x, g2x, g1y, g2y = tpl, 0, tpl / np.sqrt(3), tpl * 2 / np.sqrt(3)
            g3x, g3y = g1x - g2x, g1y - g2y
            modulation_3c = v0 * (np.cos(g1x*x + g1y*y) +
                                  np.cos(g2x*x + g2y*y) + np.cos(g3x*x + g3y*y))

            if params.disorder == 'site':
                random_delta = (kwant.digest.uniform(
                    repr(site) + salt) - 0.5)  # disorder in sites
            if params.disorder == 'bump':
                random_delta = self.get_random_shift(site)  # disorder in bumps
            self.POT[site] = vr * random_delta + modulation_3c
            return 4 * t + self.POT[site]

        def potential_lead(site, params):
            w0 = params.w0
            return 4 * t - 1.5 * 3.116 * w0  # 3cos

        syst = kwant.Builder()

        syst[(lat(x, y) for x in range(L) for y in range(W))] = potential
        syst[lat.neighbors()] = hopping

        lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
        lead[(lat(0, j) for j in range(WL))] = potential_lead
        lead[lat.neighbors()] = hopping  # depends on B

        syst.attach_lead(lead)             # 0
        syst.attach_lead(lead.reversed())  # 1

        top_lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
        top_lead[(lat(0, W-1-j) for j in range(Wtop))] = potential_lead
        top_lead[lat.neighbors()] = hopping

        syst.attach_lead(top_lead.reversed())  # 2
        syst.attach_lead(top_lead)             # 3
        return syst.finalized()

    def hall_resistance(self, energy, params):
        smatrix = kwant.smatrix(self.syst, energy, params=dict(params=params))

        def st(i, j):
            return smatrix.transmission(i, j)
        T0 = st(0, 1) + st(0, 2) + st(0, 3)
        T1 = st(1, 0) + st(1, 2) + st(1, 3)
        T2 = st(2, 0) + st(2, 1) + st(2, 3)
        T3 = st(3, 0) + st(3, 1) + st(3, 2)
        S = st(1, 0) + st(3, 0) + st(1, 2) + st(3, 2)
        alfa22 = S*T1 - (st(1, 0) + st(1, 2)) * (st(2, 1) + st(0, 1))
        alfa11 = S*T0 - (st(0, 3) + st(0, 1)) * (st(3, 0) + st(1, 0))
        alfa21 = (st(1, 0)*st(3, 2) - st(1, 2)*st(3, 0))
        alfa12 = (st(0, 1)*st(2, 3) - st(2, 1)*st(0, 3))
        detb0 = (alfa22*alfa11 - alfa12*alfa21)/S
        # print('LOCALS ', {k:v for k, v in locals().items() if k.startswith('T') or k.startswith('det')})
        rxx = -(st(1, 3) * st(2, 0) - st(2, 3) *
                st(1, 0)) / detb0  # k=0 l=1 m=3 n=2
        ryy = -(st(3, 1) * st(2, 0) - st(2, 1) * st(3, 0)) / detb0
        # rxx2 = (st(0, 1) * st(3, 2) - st(0, 2) * st(3, 1)) / detb0 # k=0 l=3 m=1 n=2
        #    rxx = (st(1, 3) * st(2, 0) - st(2, 1) * st(3, 0)) / detb0 # up to 23h 17.06.2021
        rxy12 = -alfa12 / detb0
        rxy21 = alfa21 / detb0
        r2t22 = alfa22 / detb0
        r2t11 = alfa11 / detb0
        return (rxy21, rxx, r2t22, rxy12, r2t11, ryy)

    def comp_R_Hall_DoS(self, energy, Bts, fname):
        logging.info(
            f"comp_R_Hall_DoS(e={energy}, Bts=[{Bts[0]}...{Bts[-1]},{len(Bts)}], \"{fname}\")")
        params = self.params
        for Bt in tqdm.tqdm(Bts, disable=self.args.mpi_rank > 0):
            params.Bt = 1/Bt
            rxy21, rxx, r2t22, rxy12, r2t11, rxx2 = self.hall_resistance(
                energy, params)
            local_dos = kwant.ldos(
                self.syst, energy, params=dict(params=params))
            dos = local_dos.sum() / local_dos.size
            open(fname, 'a').write(
                f'{Bt} {energy} {dos} {rxy21} {rxx} {rxy12} {r2t22} {r2t11}, {rxx2}\n')

    def check(self):
        e = 2.824086  # 5.14759966618
        x = 0.788  # -1.964
        params = self.params
        params.Bt = x
        logging.info(f'PARAMS={params}')
        logging.info(f'SYST={self.syst}')
        sm = kwant.smatrix(self.syst, e, params=dict(params=params))

        def st(i, j):
            return sm.transmission(i, j)
        T0 = st(0, 1) + st(0, 2) + st(0, 3)
        T1 = st(1, 0) + st(1, 2) + st(1, 3)
        T2 = st(2, 0) + st(2, 1) + st(2, 3)
        T3 = st(3, 0) + st(3, 1) + st(3, 2)
        S = st(1, 0) + st(3, 0) + st(1, 2) + st(3, 2)
        alfa22 = S*T1 - (st(1, 0) + st(1, 2)) * (st(2, 1) + st(0, 1))
        alfa11 = S*T0 - (st(0, 3) + st(0, 1)) * (st(3, 0) + st(1, 0))
        alfa12 = (st(0, 1)*st(2, 3) - st(2, 1)*st(0, 3))
        alfa21 = (st(1, 0)*st(3, 2) - st(1, 2)*st(3, 0))
        detb0 = (alfa22*alfa11 - alfa12*alfa21)/S
        det1 = -T1 * (T2 * T3 - st(2, 3) * st(3, 2))
        det2 = st(1, 2) * (T3 * st(2, 1) + st(2, 3) * st(3, 1))
        det3 = st(1, 3) * (T2 * st(3, 1) + st(2, 1) * st(3, 2))
        det11 = det1 + det2 + det3
        # print('LOCALS ', {k:v for k, v in locals().items() if k.startswith('T') or k.startswith('det')})
        rxx = (st(3, 1) * st(2, 0) - st(2, 1) * st(3, 0)) / det11
        rxy = -(st(0, 3) * st(2, 1) - st(0, 1) * st(2, 3)) / det11

        r2t1 = -(st(2, 2) * st(0, 0) - st(2, 0) * st(0, 2)) / det11
        r2tb = - alfa22 / det11
        rxy12 = -alfa12 / detb0
        rxy21 = alfa21 / detb0
        r2t22 = alfa22 / detb0
        r2t11 = alfa11 / detb0
        logging.info(f'check: rxy={rxy} -det11={-det11} rxx={rxx}')
        logging.info(f'check: rxy21={rxy21} detb0={detb0} r2t22={r2t22}')
        logging.info(
            f'check: rxy12={rxy12}, detb0={detb0} r2t11={r2t11}')


def main(args):
    model = Main(args)
    model.plot_bin(file='model.bin')
    if args.mpi_rank == 0:
        model.check()
        logging.info(
            f'env={json.dumps(dict(os.environ), indent=4, sort_keys=True)}')

    Bts = np.linspace(args.I[0], args.I[1], args.I[2])
    beg, end = get_beg_end(len(Bts), args.mpi_rank, args.mpi_size)
    myBts = Bts[beg:end]
    fname = new_file_with_header_XXX(
        f"comp_R_Hall_Dos_E{args.E:.2f}_{args.mpi_rank:03d}of{args.mpi_size}.XXX.dat",
        f"# 1/B E DoS Rxy21 Rxx Rxy12 R2t22 R2t11 Ryy"
    )
    model.comp_R_Hall_DoS(args.E, myBts, fname)


def config_logging():
    fmt_atty = '\x1b[1m%(asctime)s|%(levelname)1.1s|L%(lineno)04d:\x1b[0m %(message)s'
    fmt_file = '%(asctime)s|%(levelname)1.1s|L%(lineno)04d: %(message)s'
    logging.basicConfig(
        format=fmt_file,
        datefmt='%H:%M:%S',
        level=logging.INFO,
    )


def get_rank_size():
    mpivars = [(k, v)
               for k, v in os.environ.items() if re.search(r'(PMI|MPI).*RANK', k)]
    if not mpivars:
        return 0, 1
    if 'OMPI_COMM_WORLD_RANK' in os.environ:
        return (
            int(os.environ['OMPI_COMM_WORLD_RANK']),
            int(os.environ['OMPI_COMM_WORLD_SIZE']),
        )
    if 'PMI_RANK' in os.environ:
        return (
            int(os.environ['PMI_RANK']),
            int(os.environ['PMI_SIZE']),
        )
    logging.info(f'mpivars = {mpivars}')
    raise ValueError(f'NO MPI')


if __name__ == '__main__':
    config_logging()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-E', type=float,
                        help='Energy',
                        required=True,
                        )
    parser.add_argument('-I', type=str, metavar='BEG,END,NPT',
                        help='Inverse magnetic field (e.g. `1,35,100` means 1T..1/35T on 100pts)',
                        required=True,
                        )
    parser.add_argument('-vr', type=float,
                        help='Amplitude of disorder (e.g. 0 or 2)',
                        required=True,
                        )
    parser.add_argument('-w0', type=float,
                        help='Amplitude of modulation (e.g. 0.1 to 0.5)',
                        required=True,
                        )
    parser.add_argument('--salt',
                        type=str,
                        help='Random seed, default `olya`',
                        default='olya',
                        )
    parser.add_argument('-a', type=float,
                        help='Lattice constant, nm',
                        default=8.0,
                        )
    parser.add_argument('-W', type=int,
                        help='Width of the system, sites',
                        default=500,
                        )
    parser.add_argument('-L', type=int,
                        help='Length of the system, sites',
                        default=500,
                        )
    parser.add_argument('-Wtop', type=int,
                        help='Width of top leads, sites',
                        default=70,
                        )
    parser.add_argument('-Wbot', type=int,
                        help='Width of bottom leads, sites',
                        default=70,
                        )
    parser.add_argument('-t', type=float,
                        help='Hopping magnitude',
                        default=8.8853,
                        )
    parser.add_argument('--disorder', choices=['site', 'bump'], default='site',
                        help='Disorder type',
                        )
    args = parser.parse_args()
    ibbeg, ibend, nb = args.I.split(',')
    args.I = [float(ibbeg), float(ibend), int(nb)]
    args.mpi_rank, args.mpi_size = get_rank_size()
    if hasattr(kwant.solvers, 'mumps'):
        kwant.solvers.mumps.options(nrhs=1)
    logging.info(args)
    main(args)
