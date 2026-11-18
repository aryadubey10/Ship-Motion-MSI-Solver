#!/usr/bin/env python3
"""OE3035 MSI project: reproducible six-DOF frequency-domain calculation.

Examples (run at repository root):
    python tools/get_data.py KCS_L0p950
    python src/msi.py KCS_L0p950
    python src/msi.py KCS_L0p950 --dw 0.0005  # convergence check

No artificial results are generated when the HydRA inputs are absent.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import norm

G = 9.80665
HS, T0M1, GAMMA = 4.5, 11.5, 1.5
HEADINGS = (90, 135, 180)
RATIOS = (0.90, 0.95, 1.00, 1.05, 1.10)
ROOT = Path(__file__).resolve().parents[1]
COLS_CASES = ['hull','beta_deg','fn_ratio','U_mps','sigma_a','abar','we_bar','mu_msi','msi_pct']
COLS_ROLL = ['hull','omega_n','zeta_wave','Bv44']
COLS_SENS = ['study','hull','beta_deg','fn_ratio','x_P_m','y_P_m','zeta_v',
             'sigma_a','abar','we_bar','mu_msi','msi_pct']
COLS_SPIKE = ['hull','beta_deg','fn_ratio','omega']


def spectrum(omega):
    """IACS Rec.34 JONSWAP, in wave-frequency coordinates, m^2/(rad/s)."""
    o = np.asarray(omega, dtype=float)
    tp = T0M1 / (0.7757 + 0.0965 * np.sqrt(GAMMA) - 0.0144 * GAMMA)
    wp = 2*np.pi / tp
    alpha = 0.2 / (0.1160 + 0.0594*np.sqrt(GAMMA) + 0.0246*GAMMA)
    safe = np.maximum(o, 1e-12)
    sd = np.where(o <= wp, 0.07, 0.09)
    peak = np.exp(-0.5*((safe/wp-1)/sd)**2)
    return np.where(o > 0, alpha * (5/16)*HS**2 * wp**4 * safe**-5 *
                    np.exp(-1.25*(wp/safe)**4) * GAMMA**peak, 0.0)


def msi_model(sigma_a, we_bar):
    abar = np.sqrt(2/np.pi) * sigma_a
    mu = -0.819 + 2.32*(np.log10(we_bar))**2
    msi = 100 * norm.cdf((np.log10(abar/G)-mu)/0.4)
    return float(abar), float(mu), float(msi)


def speed_from_file(d):
    return np.asarray(d['speed'], dtype=float).ravel()


def interp_matrix(x, arr, target):
    """PCHIP along wave-frequency axis; no interpolation of pre-solved RAOs."""
    return PchipInterpolator(x, arr, axis=0, extrapolate=False)(target)


@dataclass
class Hydro:
    name: str
    w: np.ndarray
    beta: np.ndarray
    speed: np.ndarray
    M: np.ndarray
    C: np.ndarray
    A: np.ndarray
    B: np.ndarray
    F: np.ndarray
    R: np.ndarray
    L: float
    breadth: float
    x: float
    y: float

    @classmethod
    def load(cls, name, data_dir, particulars):
        paths = sorted(Path(data_dir).rglob(f'{name}_hydra.json'))
        if len(paths) != 1:
            raise FileNotFoundError(f'Expected one {name}_hydra.json in {data_dir}; found {len(paths)}. '
                                    f'Download the two HydRA ZIPs as described in RUN_ME_FIRST.md.')
        d = json.loads(paths[0].read_text(encoding='utf-8'))
        if not all(k in d for k in ('w','beta','speed','M','C','AM','BD','Fk_Re','Fk_Im','Fs_Re','Fs_Im')):
            raise ValueError(f'{paths[0]} lacks expected HydRA JSON keys: {list(d)}')
        w_all = np.asarray(d['w'], dtype=float).ravel()
        indices = np.flatnonzero(w_all > 0)
        w = w_all[indices]
        if not np.all(np.diff(w)>0):
            order = np.argsort(w)
            indices, w = indices[order], w[order]
        beta = np.asarray(d['beta'], dtype=float).ravel()
        speed = speed_from_file(d)
        if not np.all(np.diff(speed) >= 0):
            raise ValueError(f'{name}: speeds are out of order')
        def get(key, tail):
            arr = np.asarray(d[key], dtype=float)
            expected = (len(w_all),) + tail
            if arr.shape != expected:
                raise ValueError(f'{name}: {key} shape is {arr.shape}, expected {expected}. '
                                 'Inspect the JSON orientation before modifying code.')
            return arr[indices]
        nbt = (6,6,len(beta),len(speed))
        nfl = (6,len(beta),len(speed))
        A = get('AM', nbt)
        B = get('BD', nbt)
        F = (get('Fk_Re',nfl)+get('Fs_Re',nfl)) + 1j*(get('Fk_Im',nfl)+get('Fs_Im',nfl))
        R = None
        if 'rao_Re' in d and 'rao_Im' in d:
            R = get('rao_Re',nfl) + 1j*get('rao_Im',nfl)
        return cls(name,w,beta,speed,np.asarray(d['M'],float),np.asarray(d['C'],float),
                   A,B,F,R,particulars['Lpp'],particulars['B'],particulars['x_P'],particulars['y_P'])

    def ib(self, beta_deg):
        j = int(np.argmin(abs(np.rad2deg(self.beta)-beta_deg)))
        if abs(np.rad2deg(self.beta[j])-beta_deg) > 0.01:
            raise ValueError(f'{self.name}: heading {beta_deg} degrees is missing: {np.rad2deg(self.beta)}')
        return j

    def roll_parameters(self, zeta_v=0.05):
        ib = self.ib(90)
        A44 = PchipInterpolator(self.w, self.A[:,3,3,ib,2], extrapolate=False)
        den = lambda w: self.M[3,3] + A44(w)
        f = lambda w: w*w*den(w)-self.C[3,3]
        grid = np.linspace(self.w[0], self.w[-1], 1000)
        pairs = [(a,b) for a,b in zip(grid[:-1],grid[1:]) if f(a)*f(b)<0]
        if not pairs:
            raise ValueError(f'{self.name}: no roll natural-frequency root in HydRA frequency range')
        # Physical roll root closest to a first estimate ignoring added mass.
        guess = np.sqrt(self.C[3,3] / max(self.M[3,3],1))
        lo, hi = min(pairs, key=lambda ab: abs(sum(ab)/2-guess))
        wn = float(brentq(f, lo, hi))
        critical = 2*np.sqrt(self.C[3,3]*den(wn))
        bwave = float(PchipInterpolator(self.w, self.B[:,3,3,ib,2])(wn))
        return dict(omega_n=wn, zeta_wave=float(bwave/critical),
                    Bv44=float(zeta_v*critical), Bcritical=float(critical))

    def solve(self, omega, beta_deg, iu, Bv44=0, validate=False):
        """Complex 6DOF RAO per m wave amplitude at arbitrary *wave* frequencies."""
        ib = self.ib(beta_deg)
        om = np.atleast_1d(np.asarray(omega, float))
        A = interp_matrix(self.w,self.A[:,:,:,ib,iu],om)
        B = interp_matrix(self.w,self.B[:,:,:,ib,iu],om)
        F = (interp_matrix(self.w,self.F[:, :,ib,iu].real,om) +
             1j*interp_matrix(self.w,self.F[:,:,ib,iu].imag,om))
        we = om - om*om*self.speed[iu]*np.cos(np.deg2rad(beta_deg))/G
        Bd = B.copy()
        Bd[:,3,3] += Bv44
        Z = -we[:,None,None]**2*(self.M[None,:,:]+A) + 1j*we[:,None,None]*Bd + self.C[None,:,:]
        rao = np.linalg.solve(Z,F[...,None])[...,0]
        if validate and self.R is not None:
            exact = self.R[:,:,ib,iu]
            err = np.linalg.norm(rao-exact)/max(np.linalg.norm(exact),1e-12)
            if err > 1e-3:
                raise ValueError(f'{self.name}: zero-viscous RAOs mismatch HydRA at beta={beta_deg}, '
                                 f'iu={iu} by {err:.2%}. Check input conventions or data orientation.')
        return rao, we

    def response(self, beta_deg, iu, zeta_v=0.05, x=None, y=None, dw=0.002):
        x = self.x if x is None else x
        y = self.y if y is None else y
        roll = self.roll_parameters(zeta_v)
        # The input only resolves wave frequencies 0.1 ... 1.5 rad/s.
        # No unsupported extrapolation into high waves.
        n = int(np.ceil((self.w[-1]-self.w[0])/dw))
        omega = np.linspace(self.w[0],self.w[-1],n+1)
        rao, we = self.solve(omega,beta_deg,iu,roll['Bv44'])
        # Body-frame positive pitch is bow-down (lecture Chapter 8).
        vertical = rao[:,2] + y*rao[:,3] - x*rao[:,4]
        # Cosine short-wave taper for head and bow-quartering seas.
        cutoff = np.sqrt(8*np.pi*G/self.L)
        if beta_deg != 90:
            first = 0.9*cutoff
            taper = np.ones_like(omega)
            on = (omega > first) & (omega < cutoff)
            taper[on] = 0.5*(1+np.cos(np.pi*(omega[on]-first)/(cutoff-first)))
            taper[omega>=cutoff] = 0
            vertical = vertical*taper
        S = np.abs(vertical)**2 * we**4 * spectrum(omega)
        m0 = float(np.trapezoid(S,omega))
        m1 = float(np.trapezoid(we*S,omega))
        if not (m0>0 and m1>0):
            raise ValueError(f'{self.name} {beta_deg}deg {iu}: non-positive spectral moment')
        sig = np.sqrt(m0)
        webar = m1/m0
        abar, mu, msi = msi_model(sig,webar)
        # One-sided acceleration spectrum on encounter-frequency axis.
        jac = 1-2*omega*self.speed[iu]*np.cos(np.deg2rad(beta_deg))/G
        Se = S/np.abs(jac)
        return dict(sigma_a=sig,abar=abar,we_bar=webar,mu_msi=mu,msi_pct=msi,
                    omega=omega,omega_e=we,S_a_wave=S,S_a_encounter=Se)


def write_rows(path, names, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf8') as f:
        w = csv.DictWriter(f,fieldnames=names)
        w.writeheader()
        for r in rows:
            w.writerow({k:(f'{v:.10g}' if isinstance(v, (float,np.floating)) else v)
                        for k,v in r.items() if k in names})


def plot_spectrum(figdir):
    w = np.linspace(0.001,6,6000)
    plt.figure(figsize=(7,4))
    plt.plot(w,spectrum(w),label='IACS JONSWAP (gamma=1.5)')
    plt.xlim(0,2); plt.xlabel('Wave frequency (rad/s)')
    plt.ylabel('Wave spectral density (m² s)')
    plt.title('Design sea: Hs=4.5 m; T0m1=11.5 s')
    plt.grid(alpha=.25); plt.tight_layout()
    plt.savefig(figdir/'sea_spectrum.png',dpi=180);plt.close()
    m0=np.trapezoid(spectrum(w),w)
    mm1=np.trapezoid(spectrum(w)/w,w)
    print(f'JONSWAP check: m0 = {m0:.5f} m^2; T0m1 = {2*np.pi*mm1/m0:.4f} s')
    if abs(m0-1.2653)>0.006 or abs(2*np.pi*mm1/m0-11.5)>0.06:
        raise AssertionError('JONSWAP verification failed')


def make_figures(hydros,cases,figdir,dw):
    colors=['#1D4ED8','#D97706']
    hnames=list(hydros)
    # 4a roll RAO, with/without extra viscous damping
    plt.figure(figsize=(8,5))
    for i,(h,hyd) in enumerate(hydros.items()):
        omega=np.linspace(hyd.w[0],hyd.w[-1],int((hyd.w[-1]-hyd.w[0])/dw)+1)
        for zeta,ls in ((0,'--'),(.05,'-')):
            Bv=hyd.roll_parameters(zeta)['Bv44']
            R,_=hyd.solve(omega,90,2,Bv)
            plt.plot(omega,np.abs(R[:,3]),ls,color=colors[i],
                     label=f'{h}: zeta_v={zeta:.2f}')
    plt.xlabel('Wave / encounter frequency (rad/s)')
    plt.ylabel('Roll RAO (rad/m)')
    plt.title('Beam seas at design Froude number')
    plt.legend(fontsize=8);plt.grid(alpha=.25);plt.tight_layout()
    plt.savefig(figdir/'roll_rao.png',dpi=180);plt.close()

    # 4b MSI vs speed, one curve per heading, and both hulls.
    fig,ax=plt.subplots(figsize=(8,5))
    for i,h in enumerate(hnames):
        for beta,ls in zip(HEADINGS,('solid','dashed','dotted')):
            ss=sorted([r for r in cases if r['hull']==h and r['beta_deg']==beta],key=lambda r:r['fn_ratio'])
            ax.plot([s['fn_ratio'] for s in ss],[s['msi_pct'] for s in ss],
                    marker='o',ms=3,color=colors[i],ls=ls,label=f'{h}, {beta}°')
    ax.set_xlabel('Fn/Fn_des'); ax.set_ylabel('MSI in 2 h (%)')
    ax.set_title('MSI versus speed at three wave headings')
    ax.grid(alpha=.25);ax.legend(fontsize=7,ncol=2);fig.tight_layout()
    fig.savefig(figdir/'msi_vs_speed.png',dpi=180);plt.close(fig)

    # 4c - acceleration spectra by heading, two hulls
    for beta in HEADINGS:
        plt.figure(figsize=(8,5))
        for i,(h,hyd) in enumerate(hydros.items()):
            res=hyd.response(beta,2,dw=dw)
            plt.plot(res['omega_e'],res['S_a_encounter'],color=colors[i],label=h)
        plt.xlabel('Encounter frequency (rad/s)')
        plt.ylabel('Acceleration spectrum (m² s⁻³)')
        plt.title(f'Cabin acceleration spectrum, beta={beta}°, Fn/Fn_des=1')
        plt.grid(alpha=.25);plt.legend(fontsize=8);plt.tight_layout()
        plt.savefig(figdir/f'accel_spectrum_{beta}.png',dpi=180);plt.close()

    # 4d MSI contours on average frequency / mean absolute acceleration plane
    ff=np.linspace(.1,2.5,220)
    aa=np.linspace(.02,5,220)
    xx,yy=np.meshgrid(ff,aa)
    mus=-.819+2.32*np.log10(xx)**2
    val=100*norm.cdf((np.log10(yy/G)-mus)/.4)
    fig,ax=plt.subplots(figsize=(8,5))
    cs=ax.contour(xx,yy,val,levels=[1,5,10,20,40,60,80,95],colors='0.5')
    ax.clabel(cs,inline=True,fontsize=8,fmt='%d%%')
    for i,h in enumerate(hnames):
        for beta,marker in zip(HEADINGS,('o','s','^')):
            rr=[r for r in cases if r['hull']==h and r['beta_deg']==beta]
            ax.scatter([r['we_bar'] for r in rr],[r['abar'] for r in rr],
                       color=colors[i],marker=marker,s=32,label=f'{h}, {beta}°')
    ax.set_xlabel('Mean encounter frequency (rad/s)')
    ax.set_ylabel('Mean absolute acceleration (m/s²)')
    ax.set_title('All 30 computed cases on equal-MSI contours')
    ax.set_xlim(.1,2.5);ax.set_ylim(0,5)
    ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.15);fig.tight_layout()
    fig.savefig(figdir/'msi_contours.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('variant',nargs='?',default='KCS_L0p950')
    parser.add_argument('--dw',type=float,default=.002,help='Wave-frequency integration step in rad/s')
    parser.add_argument('--data-dir',type=Path,default=ROOT/'data')
    args=parser.parse_args()
    if args.dw<=0 or args.dw>.02:
        parser.error('Use integration step 0 < --dw <= 0.02 rad/s')
    baseline=args.variant.split('_')[0]+'_L1p000'
    if baseline==args.variant:
        parser.error('Specify a non-baseline variant')
    members=json.loads((ROOT/'tools'/'members.json').read_text())
    if args.variant not in members or baseline not in members:
        parser.error('Unknown hull: '+args.variant)
    hydros={h:Hydro.load(h,args.data_dir,members[h]) for h in (baseline,args.variant)}
    # Fail before overwriting CSVs if either input dataset is invalid.
    for h,hyd in hydros.items():
        if len(hyd.speed)!=5 or not np.allclose(hyd.speed,members[h]['speeds'],atol=.02):
            raise ValueError(f'{h}: speed grid differs from published project members.json')
        for b in HEADINGS:
            for i in range(5):
                hyd.solve(hyd.w,b,i,validate=True)
        roll=hyd.roll_parameters()
        print(f'{h}: roll wn={roll["omega_n"]:.5f} rad/s; added Bv44={roll["Bv44"]:.5e} N m s')
    figdir=ROOT/'report'/'figures';figdir.mkdir(parents=True,exist_ok=True)
    plot_spectrum(figdir)
    roll_rows=[];case_rows=[];sensitivity=[]
    for h,hyd in hydros.items():
        roll=hyd.roll_parameters()
        roll_rows.append(dict(hull=h,**{k:roll[k] for k in ('omega_n','zeta_wave','Bv44')}))
        for beta in HEADINGS:
            for i,ratio in enumerate(RATIOS):
                out=hyd.response(beta,i,dw=args.dw)
                case_rows.append(dict(hull=h,beta_deg=beta,fn_ratio=ratio,U_mps=float(hyd.speed[i]),
                                      **{k:out[k] for k in ('sigma_a','abar','we_bar','mu_msi','msi_pct')}))
    hyd=hydros[args.variant]
    for study,beta,x,y,zeta in [
        ('cabin_midship',180,0,hyd.y,.05),
        ('cabin_centreline',90,hyd.x,0,.05),
        ('cabin_starboard',90,hyd.x,-hyd.y,.05),
        ('damping_low',90,hyd.x,hyd.y,.025),
        ('damping_high',90,hyd.x,hyd.y,.10),
    ]:
        out=hyd.response(beta,2,zeta_v=zeta,x=x,y=y,dw=args.dw)
        sensitivity.append(dict(study=study,hull=hyd.name,beta_deg=beta,fn_ratio=1.0,
                                x_P_m=x,y_P_m=y,zeta_v=zeta,
                                **{k:out[k] for k in ('sigma_a','abar','we_bar','mu_msi','msi_pct')}))
    resdir=ROOT/'results'
    write_rows(resdir/'msi_cases.csv',COLS_CASES,case_rows)
    write_rows(resdir/'roll_damping.csv',COLS_ROLL,roll_rows)
    write_rows(resdir/'sensitivity.csv',COLS_SENS,sensitivity)
    # No data points are replaced. Do not declare invented 'spikes'.
    write_rows(resdir/'spikes.csv',COLS_SPIKE,[])
    make_figures(hydros,case_rows,figdir,args.dw)
    print('Wrote 30 cases, 2 roll rows, 5 sensitivity studies, 7 figures.')
    print('Review raw coefficients for irregular-frequency spikes before submission.')
    print('Next: python src/build_report.py')


if __name__=='__main__':
    main()
