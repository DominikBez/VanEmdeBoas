# Copyright Dominik Bez 2022.

import math
from itertools import product
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

FOLDER = '../presentation/assets/'

# Helper

def create_marker():
    markers = [str(item[0]) + '-' for item in matplotlib.markers.MarkerStyle.markers.items()]
    del markers[1]
    cur_marker = 0

    def closure_marker():
        nonlocal cur_marker
        mark = markers[cur_marker]
        cur_marker += 1
        return mark

    def closure_reset():
        nonlocal cur_marker
        cur_marker = 0

    return closure_marker, closure_reset

marker, reset_marker = create_marker()

# End Helper

data = pd.read_csv('../results/veb.csv')
# allTypes = data['tree'].unique()
distributions = data['distribution'].unique()
sizes = data['size'].unique()
xticks = [sizes[i] for i in range(0, len(sizes), 2)]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
typesList = [['8std::set', '16std::set', '32std::set', '64std::set', '8uVEB', '16uVEB', '32uVEB', '64uVEB'],
    ['6uVEB', '8uVEB', '12uVEB', '13uVEB', '14uVEB', '16uVEB'],
    ['64uVEB', '64VEB', '127uVEB', '128uVEB', '128VEB'],
    ['32std::set', '32uVEB', 'uVEB32'],
    ['32uVEB', '32uVEBL', 'uVEB32', 'uVEB32L', 'uVEB32LT', 'uVEB32LFG', 'uVEB32LL']]

for op in ['insert', 'lookup', 'remove']:
    pdf = PdfPages(FOLDER + op + 'VEB.pdf')
    plt.figure(figsize=(9, 5.0625))
    plt.title("{K}std::set = std::set<uint{K}_t>, {K}uVEB = VanEmdeBoas<{K}>, {K}VEB = VanEmdeBoas<{K}, int{K}_t>,"
        + " uVEB32 = VanEmdeBoas32<>, 32uVEBL = VanEmdeBoasLocked<32>, uVEB32L = VanEmdeBoas32Locked<>, "
        + " uVEB32LT = VanEmdeBoas32LockedTop<>, uVEB32LFG = VanEmdeBoas32LockedFineGrained<>,"
        + " uVEB32LL = VanEmdeBoas32Lockless.\n No #defines => sf::contention_free_shared_mutex is used often;"
        + " also bytell_hash_map by Malte Skarupke is used for VanEmdeBoas and VanEmdeBoasLocked"
        + " (not VanEmdeBoas32 and its parallel variants)"
        + "\n Random distributions: uniform, cluster = random placed clusters with 1000 succeeding elements,"
        + " normal = normal distribution with mean ~0/2^31 for signed/unsigned and std (2^31)/10,"
        + " incProb = linear increasing probability where the smallest value has probability 0,"
        + " decProb = linear decreasing probability where the largest value has probability 0"
        + "\n There are ten iterations for each data point.\n Hardware: i7-7700HQ, 16GB DDR4 Windows Laptop",
        wrap=True, horizontalalignment="center", pad=1,
        verticalalignment="center")
    fig = plt.gca()
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    fig.axis("off")
    plt.tight_layout(pad=10)
    pdf.savefig()
    plt.clf()

    for dist, types in product(distributions, typesList):
        typesFrames = list(filter(lambda tf: not tf[1].empty,
            map(lambda tf: (tf[0], tf[1][tf[1]['distribution'] == dist]),
            map(lambda t: (t, data[data['tree'] == t]), types))))
        if len(typesFrames) < 3:
            continue
        typesFiltered, *_ = zip(*typesFrames)
        fig, ax = plt.subplots()
        ax.set_prop_cycle('color', (c for c in colors for _ in range(2)))
        for t, filtered in typesFrames:
            std = filtered.groupby('size').std()
            mean = filtered.groupby('size').mean().reset_index()
            m = marker()
            mean.plot(ax=ax, x='size', y=op, figsize=(9, 5.0625), yerr=std, legend=False, style=m)
            mean.plot(ax=ax, x='size', y=op, figsize=(9, 5.0625), logx=True, style=m)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('log2(Size)')
        ax.set_xticks(xticks)
        ax.set_xticklabels([int(math.log(x, 2)) for x in xticks])
        ax.legend(typesFiltered)
        opCapitalized = op.capitalize()
        ax.set_ylabel(f'Time per {opCapitalized} [ns]')
        if op == 'lookup':
            plt.suptitle(f"Time of 10000 Lookups in a Tree with 'Size' Elements ({dist} distribution)")
        else:
            plt.suptitle(f"Time to {opCapitalized} 'Size' Elements ({dist} distribution)")
        plt.tight_layout(pad=1)
        reset_marker()
        pdf.savefig(fig, transparent=True)
        plt.close(fig)

        # Zoomed in
        fig, ax = plt.subplots()
        ax.set_prop_cycle('color', (c for c in colors for _ in range(2)))
        for t, filtered in typesFrames:
            std = filtered.groupby('size').std()
            mean = filtered.groupby('size').mean().reset_index()
            m = marker()
            mean.plot(ax=ax, x='size', y=op, figsize=(9, 5.0625), yerr=std, legend=False, style=m)
            mean.plot(ax=ax, x='size', y=op, figsize=(9, 5.0625), logx=True, style=m)
        quantileData = data[data['distribution'] == dist]
        quantileData = quantileData.loc[quantileData['tree'].isin(typesFiltered)][op]
        ax.set_ylim(bottom=quantileData.quantile(0.05), top=quantileData.quantile(0.80))
        ax.set_xlabel('log2(Size)')
        ax.set_xticks(xticks)
        ax.set_xticklabels([int(math.log(x, 2)) for x in xticks])
        ax.legend(typesFiltered)
        opCapitalized = op.capitalize()
        ax.set_ylabel(f'Time per {opCapitalized} [ns]')
        if op == 'lookup':
            plt.suptitle(f"Time of 10000 Lookups in a Tree with 'Size' Elements (Zoomed in; {dist} distribution)")
        else:
            plt.suptitle(f"Time to {opCapitalized} 'Size' Elements (Zoomed in; {dist} distribution)")
        plt.tight_layout(pad=1)
        reset_marker()
        pdf.savefig(fig, transparent=True)
        plt.close(fig)
    pdf.close()
