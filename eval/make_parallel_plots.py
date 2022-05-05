# Copyright Dominik Bez 2022.

import math
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

data = pd.read_csv('../results/parallelVeb.csv')
types = data['tree'].unique()
distributions = data['distribution'].unique()
sizes = data['size'].unique()
maxSize = max(sizes)
seqData = pd.read_csv('../results/veb.csv')
seqData = seqData[seqData['tree'] == 'uVEB32']
seqData = seqData[seqData['size'] == maxSize]
threads = data['threads'].unique()
maxThreads = max(threads)
xticks = [sizes[i] for i in range(0, len(sizes), 2)]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for op in ['insert', 'lookup', 'remove']:
    opCapitalized = op.capitalize()
    pdf = PdfPages(FOLDER + 'parallel' + opCapitalized + 'VEB.pdf')
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
        + "\n The operation is parallelized with OpenMP (static scheduling) and the time until completion"
        + "is measured. Unlike for the sequential plots, the inserted elements are shuffled before insertion."
        + " This is relevant for the cluster distribution since it increases the probability that more than"
        + " one thread tries to access the same bottom data structure at the same time because without shuffling"
        + " most clusters are inserted by a single thread since the scheduling is static."
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

    for dist in distributions:
        fig, ax = plt.subplots()
        ax.set_prop_cycle('color', (c for c in colors for _ in range(2)))
        maxThreadsData = data[data['threads'] == maxThreads]
        maxThreadsData = maxThreadsData[maxThreadsData['distribution'] == dist]
        for t in types:
            filtered = maxThreadsData[maxThreadsData['tree'] == t]
            std = filtered.groupby('size').std()
            mean = filtered.groupby('size').mean().reset_index()
            m = marker()
            mean.plot(ax=ax, x='size', y=op, figsize=(9, 5.0625), yerr=std, legend=False, style=m)
            mean.plot(ax=ax, x='size', y=op, figsize=(9, 5.0625), logx=True, style=m)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('log2(Size)')
        ax.set_xticks(xticks)
        ax.set_xticklabels([int(math.log(x, 2)) for x in xticks])
        ax.legend(types)
        ax.set_ylabel(f'Time per {opCapitalized} [ns]')
        if op == 'lookup':
            plt.suptitle(f"Time of 10000 Parallel Lookups in a Tree with 'Size' Elements ({dist} distribution; {maxThreads} Threads)")
        else:
            plt.suptitle(f"Time to {opCapitalized} 'Size' Elements Parallel ({dist} distribution; {maxThreads} Threads)")
        plt.tight_layout(pad=1)
        reset_marker()
        pdf.savefig(fig, transparent=True)

        ax.set_ylim(bottom=maxThreadsData[op].quantile(0.005), top=maxThreadsData[op].quantile(0.80))
        if op == 'lookup':
            plt.suptitle(f"Time of 10000 Parallel Lookups in a Tree with 'Size' Elements (Zoomed in; {dist} distribution; {maxThreads} Threads)")
        else:
            plt.suptitle(f"Time to {opCapitalized} 'Size' Elements Parallel (Zoomed in; {dist} distribution; {maxThreads} Threads)")
        pdf.savefig(fig, transparent=True)
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.set_prop_cycle('color', (c for c in colors for _ in range(2)))
        maxSizeData = data[data['size'] == maxSize]
        maxSizeData = maxSizeData[maxSizeData['distribution'] == dist]
        for t in types:
            filtered = maxSizeData[maxSizeData['tree'] == t]
            std = filtered.groupby('threads').std()
            mean = filtered.groupby('threads').mean().reset_index()
            m = marker()
            mean.plot(ax=ax, x='threads', y=op, figsize=(9, 5.0625), yerr=std, legend=False, style=m)
            mean.plot(ax=ax, x='threads', y=op, figsize=(9, 5.0625), style=m)
        seqTime = seqData[seqData['distribution'] == dist][op].mean()
        ax.axhline(y=seqTime, label='uVEB32', linestyle='--')
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Number of Threads')
        ax.legend(list(types) + ['uVEB32'])
        ax.set_ylabel(f'Time per {opCapitalized} [ns]')
        if op == 'lookup':
            plt.suptitle(f"Time of 10000 Parallel Lookups in a Tree with 'Size' Elements ({dist} distribution; {maxSize} Elements)")
        else:
            plt.suptitle(f"Time to {opCapitalized} 'Size' Elements Parallel ({dist} distribution; {maxSize} Elements)")
        plt.tight_layout(pad=1)
        reset_marker()
        pdf.savefig(fig, transparent=True)

        ax.set_ylim(bottom=0, top=maxSizeData[op].quantile(0.80))
        if op == 'lookup':
            plt.suptitle(f"Time of 10000 Parallel Lookups in a Tree with 'Size' Elements (Zoomed in; {dist} distribution; {maxSize} Elements)", fontsize=11)
        else:
            plt.suptitle(f"Time to {opCapitalized} 'Size' Elements Parallel (Zoomed in; {dist} distribution; {maxSize} Elements)")
        pdf.savefig(fig, transparent=True)
        plt.close(fig)

        fig, ax = plt.subplots()
        maxSizeData.loc[:, op] = seqTime / maxSizeData[op]
        for t in types:
            filtered = maxSizeData[maxSizeData['tree'] == t]
            filtered = filtered.groupby('threads').mean().reset_index()
            filtered.plot(ax=ax, x='threads', y=op, figsize=(9, 5.0625), style=marker())
        plt.xlabel('Number of Threads')
        plt.ylabel('Absolute Speedup')
        ax.axhline(y=1, color='grey', linestyle='--')
        ax.legend(types)
        if op == 'lookup':
            plt.suptitle(f"Speedup over uVEB32 of 10000 Parallel Lookups in a Tree with 'Size' Elements ({dist} distribution; {maxSize} Elements)", fontsize=10)
        else:
            plt.suptitle(f"Speedup over uVEB32 to {opCapitalized} 'Size' Elements Parallel ({dist} distribution; {maxSize} Elements)")
        plt.ylim(bottom=0)
        plt.tight_layout(pad=1)
        reset_marker()
        pdf.savefig(fig, transparent=True)
        plt.close(fig)
    pdf.close()
