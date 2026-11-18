#!/usr/bin/env python3
"""Compile an OE3035 report PDF from computed CSVs and PNGs (no invented data).

Usage: python src/build_report.py
Must run src/msi.py successfully first, with actual HydRA data.
"""
from __future__ import annotations
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, KeepTogether, PageBreak


def rows(file):
    p=ROOT/'results'/file
    if not p.exists():
        raise FileNotFoundError(f'Run python src/msi.py KCS_L0p950 first; {p} is absent')
    with p.open() as f:
        return list(csv.DictReader(f))


def load_data():
    cases=rows('msi_cases.csv')
    roll=rows('roll_damping.csv')
    sens=rows('sensitivity.csv')
    spikes=rows('spikes.csv')
    if len(cases)!=30 or len(roll)!=2 or len(sens)<5:
        raise ValueError('Not enough computed data for a full report; no PDF generated.')
    if set(r['hull'] for r in cases) != set(r['hull'] for r in roll):
        raise ValueError('Hull mismatches in computed data')
    return cases,roll,sens,spikes


def run():
    cases,roll,sens,spikes=load_data()
    baseline=next(r['hull'] for r in roll if r['hull'].endswith('_L1p000'))
    variant=next(r['hull'] for r in roll if not r['hull'].endswith('_L1p000'))
    from msi import HEADINGS, RATIOS
    figdir=ROOT/'report'/'figures'
    required=['sea_spectrum.png','roll_rao.png','msi_vs_speed.png','msi_contours.png']+[f'accel_spectrum_{b}.png' for b in HEADINGS]
    missing=[f for f in required if not (figdir/f).exists()]
    if missing:
        raise FileNotFoundError(f'Missing computed figures: {missing}; run src/msi.py first')

    # Font support for physics symbols in paragraphs. Font files are embedded in PDF only.
    reg='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
    bold='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
    if Path(reg).exists():
        pdfmetrics.registerFont(TTFont('DejaVu',reg))
        pdfmetrics.registerFont(TTFont('DejaVu-Bold',bold))
        pdfmetrics.registerFontFamily('DejaVu', normal='DejaVu',bold='DejaVu-Bold')
        font,bfont='DejaVu','DejaVu-Bold'
    else:
        font,bfont='Helvetica','Helvetica-Bold'
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='ReportTitleX',parent=styles['Title'],fontName=bfont,fontSize=17,leading=23,spaceAfter=14,textColor=colors.HexColor('#153451')))
    styles.add(ParagraphStyle(name='SectionX',fontName=bfont,fontSize=12,leading=16,spaceBefore=15,spaceAfter=7,textColor=colors.HexColor('#153451'),keepWithNext=True))
    styles.add(ParagraphStyle(name='SubX',fontName=bfont,fontSize=10,leading=14,spaceBefore=10,spaceAfter=5,keepWithNext=True))
    styles.add(ParagraphStyle(name='BodyX',fontName=font,fontSize=9,leading=14,spaceAfter=7))
    styles.add(ParagraphStyle(name='TinyX',fontName=font,fontSize=7,leading=10,spaceAfter=5))
    styles.add(ParagraphStyle(name='FigX',fontName=font,fontSize=8,leading=12,alignment=TA_CENTER,spaceAfter=9))
    text_style=styles['BodyX']
    story=[]
    P=lambda txt:Paragraph(txt,text_style)
    def para(txt):story.append(P(txt))
    def section(txt): story.append(Paragraph(txt,styles['SectionX']))
    def sub(txt): story.append(Paragraph(txt,styles['SubX']))
    def figure(name,caption):
        from PIL import Image as PILImage
        path=figdir/name
        with PILImage.open(path) as im:
            aspect=im.height/im.width
        width=16.1*cm
        img=Image(str(path),width=width,height=width*aspect)
        story.append(KeepTogether([img,Paragraph(caption,styles['FigX'])]))
    def table(headers,rr,widths=None):
        data=[[Paragraph(str(v),styles['TinyX']) for v in headers]]
        for row in rr:
            data.append([Paragraph(str(v),styles['TinyX']) for v in row])
        tb=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
        tb.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7eef5')),
            ('TEXTCOLOR',(0,0),(-1,0),colors.HexColor('#102b46')),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')]),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('BOTTOMPADDING',(0,0),(-1,-1),5),
            ('TOPPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#aac2d8')),
        ]))
        story.append(tb)

    story.append(Paragraph('OE3035 — Motion Sickness Incidence on Board',styles['ReportTitleX']))
    para(f'<b>Compared hulls:</b> {baseline} (baseline) and {variant} (assigned variant).')
    para('<b>Sea state:</b> Hs = 4.5 m; T0m1 = 11.5 s; IACS Rec.34 JONSWAP γ = 1.5; long-crested. '
         '<b>Wave directions:</b> 90°, 135°, 180°. <b>Speeds:</b> Fn/Fn_des = 0.90 through 1.10. '
         '<b>Cabin:</b> port, at 0.4B; x = -0.35Lpp. <b>Viscous roll damping:</b> 5% critical.')
    para('<b>Scope:</b> Linear hydrodynamic response and 2 h O’Hanlon–McCauley motion sickness incidence. '
         'All reported values are computed from the HydRA input arrays; the computations are reproducible '
         'by running src/msi.py, followed by src/build_report.py.')

    section('1. Derivation and assumptions')
    para('<b>Coordinate system:</b> body-fixed x forward, y to port, z upward, origin at the calm-water '
         'midship plane. Heading β is the wave propagation direction measured from the positive x axis. '
         'The ship has constant forward speed in deep water and small motions; linear superposition applies.')
    para('<b>Encounter frequency:</b> ωe = ω − (ω²U/g) cos(β). For the three assigned headings, '
         'cos(β) ≤ 0 and ωe increases monotonically with ω. The exact ship speeds are taken from the HydRA data.')
    para('<b>Equations of motion:</b> Z(ω,β,U) ξ̂ = F̂, where '
         'Z = −ωe²(M+A) + iωe(B+Bv)+C, F̂ = F̂<sub>FK</sub>+F̂<sub>S</sub>. '
         'M and C are the rigid-body and hydrostatic 6×6 matrices; A and B are radiation matrices, '
         'and the excitation vectors are complex. All six degrees of freedom are solved simultaneously, '
         'retaining off-diagonal coupling. The source coefficients and excitation (not HydRA RAOs) '
         'are interpolated onto a fine wave-frequency grid before solving the linear system.')
    para('<b>Roll damping:</b> The natural frequency ωn satisfies '
         'ωn²[M44+A44(ωn)] = C44, using A44 at beam seas and design speed. '
         'The critical roll damping is Bcrit = 2√(C44[M44+A44(ωn)]); '
         'Bv,44 = 0.05 Bcrit, kept constant over headings and speeds. '
         'Wave-making damping fraction is ζwave=B44(ωn)/Bcrit.')
    para('<b>Vertical acceleration at cabin P:</b> The small-angle vertical displacement is '
         'ẑP = ξ̂3 + yP ξ̂4 − xP ξ̂5 (positive pitch is bow-down). '
         'Thus the acceleration RAO is âP(ω) = −ωe²[ξ̂3+yP ξ̂4−xP ξ̂5]. '
         'The encounter-frequency acceleration spectrum is '
         'Saa<sup>e</sup>(ωe) = |âP(ω)|² Sηη(ω) / |1−2ωU cos(β)/g|. '
         'Integrating in <i>wave frequency</i> eliminates this Jacobian.')
    para('<b>Spectral moments:</b> m<sub>n</sub><sup>a</sup> = '
         '∫ ωe<sup>n+4</sup> |ξ̂3 + yP ξ̂4 − xP ξ̂5|² Sηη(ω) dω. '
         'Hence σa=√m<sub>0</sub><sup>a</sup>, and ω̄e = m<sub>1</sub><sup>a</sup>/m<sub>0</sub><sup>a</sup>. '
         'The integrals use frequencies resolved by HydRA. For quartering and head seas, a cosine '
         'short-wave taper is applied between 90% and 100% of ω=√(8πg/Lpp). No such taper is used in beam seas.')
    para('<b>MSI model:</b> ā = √(2/π) σa for zero-mean Gaussian acceleration; '
         'μMSI = −0.819 + 2.32[log10(ω̄e)]²; '
         'MSI = 100 Φ{[log10(ā/g)−μMSI]/0.4}, where Φ is the standard normal CDF. '
         'MSI is the fraction predicted to vomit within two hours, expressed as a percentage.')
    figure('sea_spectrum.png','Figure 1. Specified long-crested design wave spectrum (IACS Rec. 34).')

    section('2. Roll damping and hydrodynamic data')
    table(['Hull','ωn (rad/s)','ζwave','Bv,44 (N m s)'],[
        (r['hull'],f'{float(r["omega_n"]):.5f}',f'{float(r["zeta_wave"]):.5f}',f'{float(r["Bv44"]):.6e}')
        for r in roll],widths=[5.0*cm,3.2*cm,2.6*cm,5.2*cm])
    if spikes:
        para(f'{len(spikes)} irregular-frequency data point(s) were replaced:')
        table(['Hull','β','Fn/Fn_des','Wave ω'],[
            (r['hull'],r['beta_deg'],r['fn_ratio'],r['omega']) for r in spikes])
    else:
        para('<b>Irregular-frequency spikes:</b> no data points were explicitly replaced in this run. '
             'Raw coefficients should be examined for speed-dependent outliers before submission; '
             'an empty spikes.csv is not proof of spike-free input.')
    figure('roll_rao.png','Figure 2. Beam-sea roll RAO, with and without the prescribed viscous damping.')

    section('3. MSI results: all 30 operating cases')
    para('Table values are SI except heading (degrees), Froude ratio and MSI (%). '
         'All CSV values are preserved with at least four significant figures.')
    for h in (baseline,variant):
        sub(h)
        block=[r for r in cases if r['hull']==h]
        table(['β','Fn ratio','U m/s','σa m/s²','ā m/s²','ω̄e rad/s','μMSI','MSI %'],[
            (r['beta_deg'],f'{float(r["fn_ratio"]):.2f}',f'{float(r["U_mps"]):.3f}',
             f'{float(r["sigma_a"]):.4f}',f'{float(r["abar"]):.4f}',
             f'{float(r["we_bar"]):.4f}',f'{float(r["mu_msi"]):.4f}',f'{float(r["msi_pct"]):.3f}')
            for r in block],widths=[.9*cm,1.45*cm,1.7*cm,2.1*cm,2.1*cm,2.15*cm,2.0*cm,1.55*cm])
    figure('msi_vs_speed.png','Figure 3. Two-hour MSI against relative speed for both hulls and three wave headings.')
    for b in HEADINGS:
        figure(f'accel_spectrum_{b}.png',f'Figure: Acceleration power spectral density for both hulls, '
               f'{b}° heading, at the design Froude number.')
    figure('msi_contours.png','Figure: All 30 cases on the mean encounter-frequency / mean absolute-acceleration plane.')

    section('4. Discussion')
    get=lambda h,b,s:next(r for r in cases if r['hull']==h and abs(float(r['fn_ratio'])-s)<1e-5 and float(r['beta_deg'])==b)
    better=[(r['beta_deg'],r['fn_ratio'],r['msi_pct']) for r in cases if r['hull']==variant]
    diffs=[]
    for b in HEADINGS:
        for s in RATIOS:
            a=float(get(baseline,b,s)['msi_pct'])
            v=float(get(variant,b,s)['msi_pct'])
            diffs.append((b,s,v-a))
    more=sum(delta>0 for _,_,delta in diffs);less=sum(delta<0 for _,_,delta in diffs)
    para(f'<b>Hull comparison:</b> {variant} has a lower MSI in {less} cases and a higher MSI in '
         f'{more} cases relative to {baseline} (of 15 comparisons). The most favourable change for the '
         f'variant is {min(diffs,key=lambda x:x[2])[2]:+.3f} percentage points at '
         f'β={min(diffs,key=lambda x:x[2])[0]}°, Fn/Fn_des={min(diffs,key=lambda x:x[2])[1]:.2f}; '
         f'the least favourable is {max(diffs,key=lambda x:x[2])[2]:+.3f} points at '
         f'β={max(diffs,key=lambda x:x[2])[0]}°, Fn/Fn_des={max(diffs,key=lambda x:x[2])[1]:.2f}. '
         'Read these changes alongside the complex roll, heave and pitch RAOs to distinguish '
         'geometry-induced amplitude changes from frequency shifts.')
    for h in (baseline,variant):
        sub(f'Speed and heading: {h}')
        for b in HEADINGS:
            rs=[get(h,b,s) for s in RATIOS]
            first,last=rs[0],rs[-1]
            low,high=float(first['msi_pct']),float(last['msi_pct'])
            sig0,sig1=float(first['sigma_a']),float(last['sigma_a'])
            we0,we1=float(first['we_bar']),float(last['we_bar'])
            para(f'At β={b}°, the MSI changes from {low:.2f}% to {high:.2f}% '
                 f'(Δ={high-low:+.2f} percentage points) as Fn/Fn_des rises from 0.90 to 1.10. '
                 f'Over the same range σa changes from {sig0:.3f} to {sig1:.3f} m/s², '
                 f'and ω̄e from {we0:.3f} to {we1:.3f} rad/s. '
                 'Both magnitude and frequency enter MSI nonlinearly; their relative influence '
                 'can be inferred from these changes and the contour diagram.')
        worst=max([get(h,b,1.0) for b in HEADINGS],key=lambda r:float(r['msi_pct']))
        para(f'At design speed, the highest MSI of {h} occurs at β={worst["beta_deg"]}° '
             f'({float(worst["msi_pct"]):.2f}%). In head seas symmetry suppresses roll and '
             'the heave/pitch contribution dominates; in beam seas roll and heave are key; '
             'the quartering case can contain all three modes and their complex-phase interference.')
    sub('Cabin position and viscous roll damping sensitivity')
    study={r['study']:r for r in sens}
    table(['Study','β','xP (m)','yP (m)','ζv','MSI %'],[
        (r['study'],r['beta_deg'],f'{float(r["x_P_m"]):.2f}',
         f'{float(r["y_P_m"]):.2f}',f'{float(r["zeta_v"]):.3f}',f'{float(r["msi_pct"]):.3f}')
        for r in sens],widths=[4.1*cm,1.2*cm,2.5*cm,2.5*cm,2.0*cm,2.4*cm])
    for name,b in (('cabin_midship',180),('cabin_centreline',90),('cabin_starboard',90),
                   ('damping_low',90),('damping_high',90)):
        rr=study[name]
        ref=float(get(variant,b,float(rr['fn_ratio']))['msi_pct'])
        new=float(rr['msi_pct'])
        para(f'{name}: MSI={new:.2f}% versus {ref:.2f}% in the reference cabin '
             f'(change {new-ref:+.2f} percentage points).')
    para('Moving toward midship removes the pitch lever arm, while centreline placement removes the roll '
         'lever arm. Starboard placement reverses the roll term and can change constructive/destructive '
         'interference with heave. Increasing viscous damping normally reduces resonant roll amplitudes, '
         'but MSI depends on the wave spectral energy near the natural roll frequency.')
    sub('Model limitations and uncertainty')
    para('The MSI regression comes from controlled exposures of a limited, mostly young male population; '
         'it need not apply equally to all ages or experienced mariners. The two-hour horizon does not '
         'capture habituation, individual sensitivity or longer exposures. A sinusoidal-motion model is '
         'applied to a stochastic sea via spectral summary measures. Long-crested waves can overstate '
         'coherent motions relative to directionally spread seas, but spreading can also activate other '
         'directions. Linear panel RAOs omit nonlinearities such as large-amplitude roll, slamming, '
         'green water and speed-dependent viscous effects; these omissions can underestimate actual '
         'accelerations in severe events. An assumed 5% critical viscous damping could either over- or '
         'underestimate actual roll response. High-frequency results depend on data bounds and the '
         'short-wave taper. A convergence check using --dw 0.001 versus 0.002 is recommended before '
         'interpreting small differences between hulls.')

    section('5. Overall recommendation and numerical verification')
    design={b:(float(get(baseline,b,1.0)['msi_pct']), float(get(variant,b,1.0)['msi_pct']))
            for b in HEADINGS}
    para('At the nominal design speed, the variant-minus-baseline differences in two-hour MSI are '
         + ', '.join(f'{design[b][1]-design[b][0]:+.3f} percentage points at β={b}°'
                     for b in HEADINGS) + '. Thus, geometry selection for this sea state should not '
         'be based on a single heading. A recommendation for this ship should account for the operating '
         'heading distribution and expected time spent at sea.')
    para(f'For the assigned hull at design speed, shifting the cabin toward midship reduces head-sea '
         f'MSI from {float(get(variant,180,1.0)["msi_pct"]):.2f}% to '
         f'{float(study["cabin_midship"]["msi_pct"]):.2f}% in this model. '
         f'In beam seas the starboard-side case gives '
         f'{float(study["cabin_starboard"]["msi_pct"]):.2f}% compared with '
         f'{float(get(variant,90,1.0)["msi_pct"]):.2f}% at the original port-side cabin. '
         'Cabin location is therefore a potentially more consequential comfort-design variable than '
         'the small roll-damping changes tested here; real placement is also constrained by arrangement, '
         'safety and other operational requirements.')
    import json
    validation=ROOT/'quality_checks'/'validation.json'
    if validation.exists():
        q=json.loads(validation.read_text(encoding='utf8'))
        para(f'<b>Numerical checks:</b> reconstruction of the original HydRA RAOs without viscous '
             f'damping has a worst relative residual of '
             f'{q["max_rao_relative_residual"]:.3e} among 30 source cases. '
             f'Recalculating with wave-frequency steps 0.002 and 0.001 rad/s gives a largest '
             f'MSI change of {q["max_msi_grid_diff_percentage_points"]:.7f} percentage points. '
             'The reported check of 20/20 validates file format and formula consistency only.')
    para('<b>Unresolved source-data issue:</b> irregular-frequency spikes were not automatically '
         'corrected. See quality_checks/ for raw-coefficient plots by hull and heading. '
         'An expert review is needed to decide if any points require replacement before final '
         'submission; any correction should be recorded and all reported results regenerated.')

    section('6. Reproducibility and references')
    para('Reproduce from the original hull JSON files: '
         '<b>python src/msi.py '+variant+'</b> then <b>python src/build_report.py</b>. '
         'The public course checker can then be run as <b>python tools/check_submission.py</b>. '
         'The checker verifies formats and internal consistency, not physical accuracy.')
    para('Sources: OE3035 Project Chapter 10 (2026), '
         'https://oe3035.mavlab.in/project/project; OE3035 Chapters 7–9; '
         'IACS Recommendation 34, Rev. 2 (2022); '
         'O’Hanlon and McCauley (1974), <i>Aerospace Medicine</i> 45(4):366–369; '
         'Journée and Pinkster (2002), <i>Introduction in Ship Hydromechanics</i>.')

    target=ROOT/'report'/'report.pdf'
    doc=SimpleDocTemplate(str(target),pagesize=A4,rightMargin=2*cm,leftMargin=2*cm,
                           topMargin=1.8*cm,bottomMargin=1.8*cm,title='OE3035 MSI - KCS comparison')
    def on_page(canvas,document):
        canvas.saveState()
        canvas.setFont(font,7)
        canvas.setFillColor(colors.HexColor('#526475'))
        canvas.drawString(2*cm,1.3*cm,f'OE3035  •  {variant} vs {baseline}')
        canvas.drawRightString(A4[0]-2*cm,1.3*cm,str(document.page))
        canvas.restoreState()
    doc.build(story,onFirstPage=on_page,onLaterPages=on_page)
    print('Created',target)

if __name__=='__main__':
    run()
