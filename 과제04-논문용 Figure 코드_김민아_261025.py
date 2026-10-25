"""Reproduce the final A–H fine-mapping figure directly from the supplied data.

Inputs beside this script:
    source_data.xlsx                      Original individual-plant measurements
    family_genotype_configuration.csv     Founding-family marker configurations
                                          (optional: if absent, it is written from CFG_SPEC below)

Run:
    python build_final_figure.py

Outputs:
    outputs/Integrated_mapping_additive_trait_figure.png (600 dpi)
    outputs/Integrated_mapping_additive_trait_figure.pdf
    outputs/Integrated_mapping_additive_trait_figure.svg
    outputs/*csv                          Recalculated statistics and sample lists

Genotypes: F3/F7 use 1.28 Mb; F14 uses 1.67 Mb.
F–H: HD x axis, CL y axis, genotype colors, panicle-length point shapes.
Ellipses summarize genotype-specific HD–CL covariance (68% contours).
A candidate interval of 1.16–1.41 Mb is displayed as specified by the researcher;
this graphical annotation is not an automatically estimated breakpoint interval.
All statistics are calculated from workbook records; no plotted effect or P value
is taken from the reference image or hard-coded as a numerical plot input.
"""

from pathlib import Path
import argparse, itertools
import numpy as np
import pandas as pd
import openpyxl
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.lines import Line2D
BASE = Path(__file__).resolve().parent
ROOT = BASE / 'outputs'
ROOT.mkdir(parents=True, exist_ok=True)
ap = argparse.ArgumentParser()
ap.add_argument('--input', type=Path, default=BASE / 'source_data.xlsx')
opt = ap.parse_args()
# 1. 원자료 읽기와 결측치 제외
rows = list(openpyxl.load_workbook(opt.input, data_only=True).active.values)
raw = pd.DataFrame(rows[1:], columns=rows[0]).apply(lambda c: c.map(lambda v: np.nan if v in ('NA', '', None) else v))
raw['family'] = raw['sample'].str.split('-').str[0].astype(int)
traits = ['HD', 'CL', 'PL', 'PN']
genos = ['A', 'H', 'B']
chosen = [3, 7, 14]
marker = 'genotype'
family_markers = {3: '1.28Mb', 7: '1.28Mb', 14: '1.67Mb'}
available_markers = ['0.96Mb', '1.16Mb', '1.28Mb', '1.55Mb', '1.67Mb']
for t in traits:
    raw[t] = pd.to_numeric(raw[t], errors='raise')
assert len(raw) == 1344 and raw['sample'].is_unique
parts = []
for f in chosen:
    a = raw[raw.family.eq(f)][['sample', 'family'] + traits].copy()
    a['genotype'] = raw.loc[a.index, family_markers[f]]
    a['analysis_marker'] = family_markers[f]
    parts.append(a)
subset = pd.concat(parts)
missing = subset[traits + [marker]].isna().any(axis=1)
excluded = subset[missing].copy()
excluded['reason'] = excluded.apply(lambda r: '; '.join((k + ' missing' for k in traits + [marker] if pd.isna(r[k]))), axis=1)
complete = subset[~missing].copy()
assert not complete.isna().any().any() and complete.groupby('family').size().tolist() == [95, 96, 96]
complete.to_csv(ROOT / 'complete_cases_for_all_trait_panels.csv', index=False)
excluded.to_csv(ROOT / 'excluded_selected_plants.csv', index=False)
# 2. 패밀리별 출수기 ANOVA와 BH 보정
primary = {}
screen_rows = []
screening_cases = []
all_exclusions = []

def bh(p):
    p = np.asarray(p, float)
    ix = np.argsort(p)
    q = np.empty(len(p))
    q[ix] = np.minimum(1, np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q
for f in range(1, 15):
    sub = raw[raw.family.eq(f)]
    avail = [m for m in available_markers if sub[m].notna().any()]
    m = family_markers[f] if f in chosen else '1.28Mb' if '1.28Mb' in avail else avail[0]
    primary[f] = m
    if f in chosen:
        a = complete[complete.family.eq(f)].copy()
    else:
        a = sub.dropna(subset=['HD', m]).copy()
        a['genotype'] = a[m]
    g = [a.loc[a[marker].eq(z), 'HD'].to_numpy(float) for z in genos]
    test = stats.f_oneway(*g)
    screen_rows.append(dict(family=f, line=f'J{26336 + f}', marker=m, n=len(a), n_A=len(g[0]), n_H=len(g[1]), n_B=len(g[2]), mean_HD=a.HD.mean(), HD_SE=stats.sem(a.HD), F=test.statistic, p=test.pvalue, neg_log10_p=-np.log10(test.pvalue), additive=(g[0].mean() - g[2].mean()) / 2))
    b = a[['sample', 'family', 'HD', 'genotype']].copy()
    b['marker'] = m
    assert not b.isna().any().any()
    screening_cases.append(b)
    ids = set(a['sample'])
    drop = sub[~sub['sample'].isin(ids)].copy()
    for r in drop.itertuples():
        all_exclusions.append(dict(sample=r.sample, family=f, reason=excluded.set_index('sample').loc[r.sample, 'reason'] if f in chosen else 'HD or representative marker missing'))
screen = pd.DataFrame(screen_rows)
screen['q_BH_14'] = bh(screen.p)
screen.to_csv(ROOT / 'recalculated_family_heading_associations.csv', index=False)
pd.concat(screening_cases).to_csv(ROOT / 'complete_cases_for_family_heading_tests.csv', index=False)
pd.DataFrame(all_exclusions).to_csv(ROOT / 'excluded_from_family_heading_tests.csv', index=False)
# 3. 창시 F₂ 유전자형 지도 읽기
CFG_PATH = BASE / 'family_genotype_configuration.csv'
CFG_POSITIONS = [0.52, 0.96, 1.12, 1.16, 1.20, 1.28, 1.31, 1.41, 1.44, 1.55, 1.57, 1.60, 1.62, 1.66, 1.67, 2.00, 2.83]
# 가계별 구성: (왼쪽 유전자형, 왼쪽 구간의 마지막 마커, 오른쪽 유전자형)
CFG_SPEC = {1: ('A', 1.62, 'H'), 2: ('A', 1.55, 'H'), 3: ('H', 1.62, 'A'), 4: ('H', 1.41, 'A'),
            5: ('H', 1.31, 'A'), 6: ('B', 1.44, 'H'), 7: ('H', 1.28, 'B'), 8: ('H', 1.66, 'B'),
            9: ('H', 1.60, 'B'), 10: ('H', 1.67, 'A'), 11: ('H', 1.16, 'A'), 12: ('H', 0.96, 'A'),
            13: ('B', 1.16, 'H'), 14: ('H', 1.67, 'B')}
if not CFG_PATH.exists():
    cfg_rows = []
    for f, (left, last, right) in CFG_SPEC.items():
        r = {'family': f, 'line': f'J{26336 + f}'}
        r.update({f'{p:.2f}Mb': left if p <= last else right for p in CFG_POSITIONS})
        cfg_rows.append(r)
    pd.DataFrame(cfg_rows).to_csv(CFG_PATH, index=False)
    print(f'NOTE: {CFG_PATH.name} not found -> written from the embedded CFG_SPEC table')
cfg = pd.read_csv(CFG_PATH)
positions = [float(c[:-2]) for c in cfg.columns[2:]]
cfg = cfg.set_index('family')
assert len(positions) == 17
for f in range(1, 15):
    assert cfg.loc[f, f'{float(primary[f][:-2]):.2f}Mb'] == 'H'
colors = {'A': '#6D91D2', 'H': '#82B77E', 'B': '#D89163'}
dark = '#304858'
sigcolor = '#347F7A'
days = np.arange(70, 87)
idx = pd.MultiIndex.from_product([chosen, days, genos], names=['family', 'HD', 'genotype'])
counts = complete.rename(columns={marker: 'genotype'}).groupby(['family', 'HD', 'genotype']).size().rename('n').reindex(idx, fill_value=0).reset_index()
counts.to_csv(ROOT / 'complete_case_daily_counts.csv', index=False)
cor = []
for f, t in itertools.product(chosen, ['CL', 'PL', 'PN']):
    a = complete[complete.family.eq(f)]
    pr = stats.pearsonr(a.HD, a[t])
    cor.append(dict(family=f, trait=t, n=len(a), r=pr.statistic, p=pr.pvalue))
cor = pd.DataFrame(cor)
cor['q_BH_9'] = bh(cor.p)
cor.to_csv(ROOT / 'complete_case_heading_trait_correlations.csv', index=False)
statsrows = []
for f in chosen:
    a = complete[complete.family.eq(f)]
    for t in traits:
        arrays = [a.loc[a[marker].eq(g), t].to_numpy(float) for g in genos]
        tt = stats.f_oneway(*arrays)
        statsrows.append(dict(family=f, marker=family_markers[f], trait=t, n=len(a), p=tt.pvalue, A_mean=arrays[0].mean(), H_mean=arrays[1].mean(), B_mean=arrays[2].mean()))
pd.DataFrame(statsrows).to_csv(ROOT / 'complete_case_genotype_trait_tests.csv', index=False)
# 4. 최종 피겨의 색상과 글꼴 설정
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.65, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

def pstr(p):
    return f'{p:.2e}' if p < 0.01 else f'{p:.3f}'
PL_CLASSES = [('≤20', 0, 20, 'o'), ('20–22', 20, 22, '^'), ('22–24', 22, 24, 's'), ('>24', 24, np.inf, 'D')]
shape_counts = []
for f in chosen:
    for label_pl, lo, hi, shape in PL_CLASSES:
        a = complete[complete.family.eq(f)]
        shape_counts.append(dict(family=f, panicle_length_class_cm=label_pl, shape=shape, n=int(((a.PL > lo) & (a.PL <= hi)).sum())))
pd.DataFrame(shape_counts).to_csv(ROOT / 'panicle_length_shape_class_counts.csv', index=False)

def label(ax, l):
    ax.text(-0.2, 1.1, l, transform=ax.transAxes, fontsize=11, fontweight='bold', va='top')

# 5. 패널 A: 유전자형 지도와 평균 출수일수
def map_and_significance(fig, outer):
    sub = outer.subgridspec(1, 2, width_ratios=[2.1, 1], wspace=0.1)
    ax = fig.add_subplot(sub[0])
    bar = fig.add_subplot(sub[1], sharey=ax)
    for row, f in enumerate(range(1, 15)):
        row_colors = colors
        calls = [cfg.loc[f, f'{p:.2f}Mb'] for p in positions]
        transitions = [j for j in range(16) if calls[j] != calls[j + 1]]
        assert len(transitions) == 1
        j = transitions[0]
        ax.add_patch(Rectangle((-0.4, row - 0.34), j + 0.4, 0.68, fc=row_colors[calls[0]], ec='none'))
        ax.add_patch(Rectangle((j + 1, row - 0.34), 16.4 - (j + 1), 0.68, fc=row_colors[calls[-1]], ec='none'))
        ax.add_patch(Rectangle((j, row - 0.34), 0.5, 0.68, fc=row_colors[calls[0]], ec='none'))
        ax.add_patch(Rectangle((j + 0.5, row - 0.34), 0.5, 0.68, fc=row_colors[calls[-1]], ec='none'))
        ax.add_patch(Rectangle((-0.4, row - 0.34), 16.8, 0.68, fill=False, ec=dark, lw=0.4))
        ax.text((j - 0.4) / 2, row, calls[0], ha='center', va='center', fontsize=5.9, color='#1C3344')
        ax.text((j + 1 + 16.4) / 2, row, calls[-1], ha='center', va='center', fontsize=5.9, color='#1C3344')
        r = screen[screen.family.eq(f)].iloc[0]
        col = sigcolor if r.q_BH_14 < 0.05 else '#BBC3C9'
        bar.barh(row, r.mean_HD, xerr=r.HD_SE, height=0.68, color=col, edgecolor=dark, linewidth=1.0 if f in chosen else 0.35, error_kw=dict(ecolor=dark, lw=0.6, capsize=1.5))
        bar.text(r.mean_HD - 2, row, f'{r.mean_HD:.2f} days | P = ' + pstr(r.p) + ' ' + stars(r.p), ha='right', va='center', fontsize=5.4, color='white' if r.q_BH_14 < 0.05 else '#304858', linespacing=1.2)
    for x in [3, 7]:
        ax.axvline(x, color='#A9513B', ls=(0, (3, 2)), lw=0.85)
    ax.add_patch(Rectangle((3, 14.0), 4, 0.4, fc='#A9513B', ec='#A9513B', lw=0.65, clip_on=False))
    ax.text(7.3, 14.2, 'Candidate interval: 1.16–1.41 Mb (~250 kb)', ha='left', va='center', fontsize=5.8, color='#A9513B')
    ax.set_xlim(-0.45, 16.45)
    ax.set_ylim(14.7, -0.6)
    ax.set_yticks(range(14), [f'Family {f}' for f in range(1, 15)], fontsize=5.9)
    for obj, f in zip(ax.get_yticklabels(), range(1, 15)):
        obj.set_color(sigcolor if screen.loc[screen.family.eq(f), 'q_BH_14'].iloc[0] < 0.05 else dark)
        obj.set_fontweight('bold' if f in chosen else 'normal')
    ax.set_xticks(range(17), [f'{p:.2f}' for p in positions], rotation=60, ha='left')
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=3)
    ax.spines[:].set_visible(False)
    for i, obj in enumerate(ax.get_xticklabels()):
        if i in [3, 7]:
            obj.set_color('#A9513B')
            obj.set_fontweight('bold')
    ax.text(-0.14, 1.25, 'A', transform=ax.transAxes, fontweight='bold', fontsize=12, va='top')
    ax.text(0.5, -0.095, 'Marker position (Mb)', transform=ax.transAxes, ha='center', va='top', fontsize=7, color=dark)
    bar.set_xlim(0, 90)
    bar.set_xticks(np.arange(0, 91, 10))
    bar.tick_params(axis='y', left=False, labelleft=False)
    bar.set_xlabel('Mean days to heading (days)', fontsize=7)
    bar.grid(axis='x', color='#E8ECEF', lw=0.4)
    bar.set_axisbelow(True)
    bar.spines[['left', 'right', 'top']].set_visible(False)
    bar.text(0.5, 1.09, 'Mean heading date', transform=bar.transAxes, fontweight='bold', fontsize=7, ha='center', va='bottom')
    return (ax, bar)

# 6. 패널 C–E: 실제 개체 수에 따른 출수일수 분포
def draw_hist(fig, outer, f, letter):
    ax = fig.add_subplot(outer)
    a = complete[complete.family.eq(f)]
    bottom = np.zeros(len(days))
    r = screen[screen.family.eq(f)].iloc[0]
    for g in genos:
        n = counts[counts.family.eq(f) & counts.genotype.eq(g)].set_index('HD').n.reindex(days, fill_value=0).to_numpy()
        ax.bar(days, n, bottom=bottom, width=0.84, color=colors[g], edgecolor='white', lw=0.35)
        bottom += n
    ax.set_ylim(0, 23)
    ax.set_xlim(69.4, 86.6)
    ax.set_xticks([70, 75, 80, 85])
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_xlabel('Days to heading (days)')
    ax.set_ylabel('Number of plants')
    ax.grid(axis='y', color='#E6EBEF', lw=0.4, ls='--')
    ax.set_axisbelow(True)
    ax.set_title(f'Family {f} | {family_markers[f][:-2]} Mb | n = {len(a)}', loc='left', fontweight='bold', pad=16)
    label(ax, letter)
    return ax

# 7. 패널 F–H: 출수일수·간장·수장
def draw_joint(fig, outer, f, letter):
    sub = outer.subgridspec(2, 2, height_ratios=[0.18, 1], width_ratios=[1, 0.14], hspace=0.02, wspace=0.025)
    top = fig.add_subplot(sub[0, 0])
    ax = fig.add_subplot(sub[1, 0])
    right = fig.add_subplot(sub[1, 1], sharey=ax)
    a = complete[complete.family.eq(f)]
    rng = np.random.default_rng(850 + f)
    for g in genos:
        b = a[a[marker].eq(g)]
        for pl_label, lo, hi, shape in PL_CLASSES:
            points = b[b.PL.gt(lo) & b.PL.le(hi)]
            ax.scatter(points.HD + rng.uniform(-0.07, 0.07, len(points)), points.CL, s=25, marker=shape, c=colors[g], ec=dark, linewidths=0.4, alpha=0.8, zorder=3)
        xx = np.linspace(69.4, 86.6, 200)
        yy = np.linspace(47, 83, 200)
        dx = stats.gaussian_kde(b.HD)(xx)
        dy = stats.gaussian_kde(b.CL)(yy)
        top.fill_between(xx, dx, color=colors[g], alpha=0.2, lw=0)
        top.plot(xx, dx, color=colors[g], lw=0.7)
        right.fill_betweenx(yy, dy, color=colors[g], alpha=0.2, lw=0)
        right.plot(dy, yy, color=colors[g], lw=0.7)
        cov = np.cov(b[['HD', 'CL']], rowvar=False)
        vals, vec = np.linalg.eigh(cov)
        order = np.argsort(vals)[::-1]
        vals = vals[order]
        vec = vec[:, order]
        ang = np.degrees(np.arctan2(vec[1, 0], vec[0, 0]))
        wh = 2 * np.sqrt(stats.chi2.ppf(0.68, 2)) * np.sqrt(vals)
        ax.add_patch(Ellipse((b.HD.mean(), b.CL.mean()), *wh, angle=ang, fc=colors[g], ec='none', alpha=0.1, zorder=1))
        ax.add_patch(Ellipse((b.HD.mean(), b.CL.mean()), *wh, angle=ang, fill=False, ec=colors[g], lw=0.7, alpha=0.8, zorder=2))
    x = a.HD.to_numpy(float)
    y = a.CL.to_numpy(float)
    xx = np.linspace(x.min(), x.max(), 200)
    reg = stats.linregress(x, y)
    pred = reg.intercept + reg.slope * xx
    res = y - (reg.intercept + reg.slope * x)
    mse = res @ res / (len(a) - 2)
    se = np.sqrt(mse * (1 / len(a) + (xx - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2)))
    tc = stats.t.ppf(0.975, len(a) - 2)
    ax.plot(xx, pred, color=dark, lw=1)
    ax.fill_between(xx, pred - tc * se, pred + tc * se, color=dark, alpha=0.12, lw=0, zorder=0)
    cr = cor[cor.family.eq(f) & cor.trait.eq('CL')].iloc[0]
    ax.text(0.03, 0.97, f'r = {cr.r:.2f}\nP = {pstr(cr.p)}', ha='left', va='top', transform=ax.transAxes, fontsize=5.8, bbox=dict(fc='white', ec='none', alpha=0.8, pad=1))
    ax.set_xlim(69.4, 86.6)
    ax.set_ylim(47, 83)
    ax.set_xticks([70, 75, 80, 85])
    ax.set_yticks([50, 60, 70, 80])
    ax.set_xlabel('Days to heading (days)')
    ax.set_ylabel('Culm length (cm)')
    ax.grid(color='#E6EBEF', lw=0.4)
    ax.set_axisbelow(True)
    top.set_xlim(69.4, 86.6)
    top.axis('off')
    right.axis('off')
    top.text(-0.2, 1.12, letter, transform=top.transAxes, fontweight='bold', fontsize=11, va='top')

def save(fig, name):
    for ext in ['png', 'pdf', 'svg']:
        temp = ROOT / (name + '.tmp.' + ext)
        fig.savefig(temp, dpi=600)
        temp.replace(ROOT / (name + '.' + ext))
    fig.savefig(ROOT / (name + '_preview.png'), dpi=150)
    plt.close(fig)
# 8. 패널 B에 사용할 상가효과와 95% 신뢰구간 계산
effect_rows = []
for f in range(1, 15):
    a = pd.concat(screening_cases).query('family == @f')
    y = a.HD.to_numpy(float)
    g = a.genotype.to_numpy()
    X = np.column_stack([np.ones(len(a)), np.select([g == 'A', g == 'B'], [1.0, -1.0], default=0.0), (g == 'H').astype(float)])
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ coef
    df = len(a) - 3
    mse = res @ res / df
    se = np.sqrt((mse * np.linalg.inv(X.T @ X))[1, 1])
    effect = coef[1]
    tc = stats.t.ppf(0.975, df)
    assert np.isclose(effect, (y[g == 'A'].mean() - y[g == 'B'].mean()) / 2)
    effect_rows.append(dict(family=f, marker=primary[f], n=len(a), additive_days=effect, additive_SE=se, CI95_low=effect - tc * se, CI95_high=effect + tc * se, p_additive=2 * stats.t.sf(abs(effect / se), df)))
effects = pd.DataFrame(effect_rows)
effects['q_additive_BH14'] = bh(effects.p_additive)
effects.to_csv(ROOT / 'additive_effect_estimates.csv', index=False)
assert effects.loc[effects.q_additive_BH14.lt(0.05), 'family'].tolist() == [3, 4, 5, 7, 8, 9, 10, 13, 14]

# 9. 패널 B: 상가효과 포리스트 플롯
def draw_effects(fig, outer, map_ax):
    sub = outer.subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.12)
    ax = fig.add_subplot(sub[0], sharey=map_ax)
    tab = fig.add_subplot(sub[1], sharey=map_ax)
    for i, r in effects.iterrows():
        col = '#347F7A' if r.q_additive_BH14 < 0.05 else '#BBC3C9'
        ax.errorbar(r.additive_days, i, xerr=[[r.additive_days - r.CI95_low], [r.CI95_high - r.additive_days]], fmt='o', ms=4.2, color=col, ecolor=col, elinewidth=1.1, capsize=0, zorder=3)
        tab.text(0, i, f'{r.additive_days:+.2f} [{r.CI95_low:+.2f}, {r.CI95_high:+.2f}]', ha='left', va='center', fontsize=5.8, color=dark, fontweight='bold' if r.family in chosen else 'normal')
    ax.axvline(0, ls=(0, (3, 2)), lw=0.8, color='#747D84')
    ax.set_xlim(-1, 3)
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.spines[['left', 'top', 'right']].set_visible(False)
    ax.set_xlabel('Additive effect (days)', fontsize=7)
    tab.axis('off')
    tab.set_xlim(0, 1)
    ax.text(-0.12, 1.25, 'B', transform=ax.transAxes, fontweight='bold', fontsize=12, va='top')
    tab.text(0, 1.09, 'Effect [95% CI]', transform=tab.transAxes, fontsize=6.3, fontweight='bold', va='bottom')
    tab.legend(handles=[Line2D([], [], marker='o', ms=4.5, ls='', color='#347F7A', label='q < 0.05'), Line2D([], [], marker='o', ms=4.5, ls='', color='#BBC3C9', label='q ≥ 0.05')], title='BH-adjusted additive-effect test', title_fontsize=6.5, ncol=2, frameon=False, fontsize=6.7, loc='upper center', bbox_to_anchor=(0.5, 0.01), borderaxespad=0, columnspacing=1.3, handletextpad=0.4)
    return ax

def legends(fig, component=False):
    legend_y = 0.065 if component else 0.045
    fig.legend(handles=[Line2D([], [], marker='s', ms=5, ls='', mfc=colors[g], mec='none', label=g) for g in genos], title='Genotype', title_fontsize=7, ncol=3, frameon=False, fontsize=7, loc='lower center', bbox_to_anchor=(0.2, legend_y), columnspacing=2)
    if not component:
        fig.legend(handles=[Line2D([], [], marker=shape, ms=5, ls='', mfc='#BCC5CC', mec=dark, mew=0.4, label=lab) for lab, lo, hi, shape in PL_CLASSES], title='Panicle length (cm)', title_fontsize=7, ncol=4, frameon=False, fontsize=7, loc='lower center', bbox_to_anchor=(0.49, legend_y), columnspacing=1.5)
# 10. 전체 A–H 배치 및 저장
fig = plt.figure(figsize=(13.8, 10.6))
grid = fig.add_gridspec(3, 3, left=0.075, right=0.985, top=0.87, bottom=0.15, height_ratios=[1.8, 0.86, 0.98], hspace=0.73, wspace=0.4)
top = grid[0, :].subgridspec(1, 2, width_ratios=[2.1, 1], wspace=0.18)
map_ax, bar_ax = map_and_significance(fig, top[0])
draw_effects(fig, top[1], map_ax)
for i, f in enumerate(chosen):
    draw_hist(fig, grid[1, i], f, chr(67 + i))
    draw_joint(fig, grid[2, i], f, chr(70 + i))
fig.suptitle('Fine-mapping configurations, additive effects, and heading-date trait relationships', x=0.53, y=0.985, fontweight='bold', fontsize=13)
legends(fig)
save(fig, 'Integrated_mapping_additive_trait_figure')
print("완료: outputs 폴더에 최종 A–H 피겨와 통계표를 저장했습니다.")