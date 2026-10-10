"""Tables and a convergence figure for the TMOP shortlist, in the Paper A style.

Reads the CSV files written by `scripts/plot_benchmarks.jl` and writes LaTeX
tables and one PDF into a paper's `generated/` directory. The mesh and
element-metric figures are produced by the Julia script itself.

    python scripts/export_shortlist.py --paper /path/to/PaperA
"""
import argparse, csv, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'

# Paper A convention.
C_UNIFORM, C_TMOP, C_GCEE = '#444444', '#24628b', '#a8372c'

CASES = [
    ('convergence_poisson',       r'Circular front, $\beta=20$'),
    ('convergence_poisson_sharp', r'Circular front, $\beta=200$'),
    ('convergence_shock',         r'Oblique shock, $\mu_{107}$'),
]

# Published target, metric and admitted size ratio for each geometric case.
QUALITY_SOURCE = {
    'ring-paper':    ('Annulus, hr \\S4.2',      '$10\\!:\\!1$'),
    'ring-tid9':     ('Annulus, \\texttt{tid 9}', '$10\\!:\\!1$'),
    'ring-tid10':    ('Annulus, \\texttt{tid 10}', '$10\\!:\\!1$'),
    'band-size':     ('Band $W_1$, IMR \\S6.1',  '$7\\!:\\!1$'),
    'band-shape':    ('Band $W_2$, IMR \\S6.1',  '$3\\!:\\!1$'),
    'band-both':     ('Band $W_3$, IMR \\S6.1',  '$3\\!:\\!1$'),
    'material':      ('Interface, \\texttt{tid 6}', '$9\\!:\\!1$'),
    'discrete-size': ('Size, \\texttt{tid 5}',   '$9\\!:\\!1$'),
}
UNIFORM_AREA = 1.0 / (12 * 12)   # ngeo = 12 in plot_benchmarks.jl


def read(stem):
    with (DATA / (stem + '.csv')).open() as fh:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def read_cases():
    with (DATA / 'cases_gcee.csv').open() as fh:
        return list(csv.DictReader(fh))


def convergence_table(out):
    lines = [r'\begin{tabular}{llrrrrrr}', r'\toprule',
             r'& & \multicolumn{3}{c}{$\|u-u_h\|_{L^2}$}'
             r' & \multicolumn{3}{c}{$|u-u_h|_{H^1}$}\\',
             r'\cmidrule(lr){3-5}\cmidrule(lr){6-8}',
             r'Case & $n$ & uniform & TMOP & GCEE & uniform & TMOP & GCEE\\',
             r'\midrule']
    for k, (stem, label) in enumerate(CASES):
        rows = read(stem)
        for j, r in enumerate(rows):
            head = label if j == 0 else ''
            lines.append(
                f"{head} & {int(r['n'])} & {r['L2_uniform']:.4f} & {r['L2_TMOP']:.4f} "
                f"& {r['L2_GCEE']:.4f} & {r['energy_uniform']:.3f} & {r['energy_TMOP']:.3f} "
                f"& {r['energy_GCEE']:.3f}\\\\")
        if k < len(CASES) - 1:
            lines.append(r'\midrule')
    lines += [r'\bottomrule', r'\end{tabular}']
    (out / 'shortlist_table.tex').write_text('\n'.join(lines) + '\n')


def quality_table(out):
    rows = read_cases()
    lines = [r'\begin{tabular}{llcrrrrrr}', r'\toprule',
             r'& & & \multicolumn{2}{c}{$\min_K |K|\ (\times 10^{-3})$}'
             r' & \multicolumn{2}{c}{$|K_{\mathrm{unif}}|/\min_K|K|$}'
             r' & \multicolumn{2}{c}{$\max_K \mu$}\\',
             r'\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}',
             r'Target & $\mu$ & ratio & TMOP & GCEE & TMOP & GCEE & TMOP & GCEE\\',
             r'\midrule']
    for r in rows:
        label, ratio = QUALITY_SOURCE[r['case']]
        at, ag = float(r['min_area_TMOP']), float(r['min_area_GCEE'])
        mt, mg = float(r['max_mu_TMOP']), float(r['max_mu_GCEE'])
        fmt = lambda v: f"{v:.0f}" if v >= 100 else f"{v:.2f}"
        lines.append(
            f"{label} & $\\mu_{{{int(r['metric'])}}}$ & {ratio} & "
            f"{1e3*at:.2f} & {1e3*ag:.2f} & "
            f"{UNIFORM_AREA/at:.1f} & {UNIFORM_AREA/ag:.1f} & "
            f"{fmt(mt)} & {fmt(mg)}\\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (out / 'shortlist_quality.tex').write_text('\n'.join(lines) + '\n')


def convergence_figure(out):
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(9.4, 5.2), layout='constrained')
    for col, (stem, label) in enumerate(CASES):
        rows = read(stem)
        h = [r['h'] for r in rows]
        series = [
            ('L2_uniform', 'L2_TMOP', 'L2_GCEE', 2, r'$\|u-u_h\|_{L^2}$'),
            ('energy_uniform', 'energy_TMOP', 'energy_GCEE', 1, r'$|u-u_h|_{H^1}$'),
        ]
        for row, (ku, kt, kg, p, ylab) in enumerate(series):
            ax = axes[row, col]
            for key, color, marker, name in (
                    (ku, C_UNIFORM, 'o', 'uniform'),
                    (kt, C_TMOP, 's', 'TMOP target'),
                    (kg, C_GCEE, 'D', 'GCEE target')):
                ax.loglog(h, [r[key] for r in rows], '-', marker=marker, ms=3.4,
                          color=color, lw=1.2, label=name)
            ref = [rows[0][ku] * (x / h[0]) ** p for x in h]
            ax.loglog(h, ref, ':', color='#999999', lw=1.0,
                      label=rf'$O(h^{p})$' if p > 1 else r'$O(h)$')
            ax.set_xlabel('$h$')
            if col == 0:
                ax.set_ylabel(ylab)
            if row == 0:
                ax.set_title(label, fontsize=9)
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.savefig(out / 'shortlist_convergence.pdf', bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paper', required=True)
    out = pathlib.Path(ap.parse_args().paper) / 'generated'
    out.mkdir(parents=True, exist_ok=True)
    convergence_table(out)
    quality_table(out)
    convergence_figure(out)
    print('wrote shortlist_table.tex, shortlist_quality.tex, shortlist_convergence.pdf ->', out)


if __name__ == '__main__':
    main()
