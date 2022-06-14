import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import seaborn as sns
sns.set_style("darkgrid")

from sklearn import metrics


sstar_1src_accuracy = pd.read_csv(snakemake.input.sstar_1src_accuracy, sep="\t").dropna()
sprime_1src_accuracy = pd.read_csv(snakemake.input.sprime_1src_accuracy, sep="\t").dropna()
skovhmm_1src_accuracy = pd.read_csv(snakemake.input.skovhmm_1src_accuracy, sep="\t").dropna()

sstar_1src_accuracy_grouped = sstar_1src_accuracy.groupby(['demography', 'scenario', 'sample', 'cutoff'], as_index=False)
sprime_1src_accuracy_grouped = sprime_1src_accuracy.groupby(['demography', 'sample', 'cutoff'], as_index=False)
skovhmm_1src_accuracy_grouped = skovhmm_1src_accuracy.groupby(['demography', 'sample', 'cutoff'], as_index=False)

sstar_1src_accuracy_mean = sstar_1src_accuracy_grouped.mean()
sprime_1src_accuracy_mean = sprime_1src_accuracy_grouped.mean()
skovhmm_1src_accuracy_mean = skovhmm_1src_accuracy_grouped.mean()

sstar_1src_accuracy_mean.to_csv(snakemake.output.sstar_1src_accuracy_mean, sep="\t", index=False)
sprime_1src_accuracy_mean.to_csv(snakemake.output.sprime_1src_accuracy_mean, sep="\t", index=False)
skovhmm_1src_accuracy_mean.to_csv(snakemake.output.skovhmm_1src_accuracy_mean, sep="\t", index=False)

sprime_1src_accuracy_mean['scenario'] = ['true'] * len(sprime_1src_accuracy_mean)
skovhmm_1src_accuracy_mean['scenario'] = ['true'] * len(skovhmm_1src_accuracy_mean)

sstar_2src_accuracy = pd.read_csv(snakemake.input.sstar_2src_accuracy, sep="\t").dropna()
sprime_2src_accuracy = pd.read_csv(snakemake.input.sprime_2src_accuracy, sep="\t").dropna()
archaicseeker2_2src_accuracy = pd.read_csv(snakemake.input.archaicseeker2_2src_accuracy, sep="\t").dropna()

sstar_2src_accuracy_grouped = sstar_2src_accuracy.groupby(['demography', 'sample', 'cutoff', 'src'], as_index=False)
sprime_2src_accuracy_grouped = sprime_2src_accuracy.groupby(['demography', 'sample', 'cutoff', 'src'], as_index=False)
archaicseeker2_2src_accuracy_grouped = archaicseeker2_2src_accuracy.groupby(['demography', 'sample', 'cutoff', 'src'], as_index=False)

sstar_2src_accuracy_mean = sstar_2src_accuracy_grouped.mean()
sprime_2src_accuracy_mean = sprime_2src_accuracy_grouped.mean()
archaicseeker2_2src_accuracy_mean = archaicseeker2_2src_accuracy_grouped.mean()

sstar_2src_accuracy_mean.to_csv(snakemake.output.sstar_2src_accuracy_mean, sep="\t", index=False)
sprime_2src_accuracy_mean.to_csv(snakemake.output.sprime_2src_accuracy_mean, sep="\t", index=False)
archaicseeker2_2src_accuracy_mean.to_csv(snakemake.output.archaicseeker2_2src_accuracy_mean, sep="\t", index=False)

methods1 = ['sstar', 'sprime', 'skovhmm']
demography1 = ['HumanNeanderthal', 'BonoboGhost']
samples = ['nref_10_ntgt_1', 'nref_50_ntgt_1', 'nref_10_ntgt_10', 'nref_50_ntgt_10']
scenarios = ['true', 'const', 'ref_tgt_only']
accuracy1 = {
    'sstar': sstar_1src_accuracy_mean,
    'sprime': sprime_1src_accuracy_mean,
    'skovhmm': skovhmm_1src_accuracy_mean,
}
methods2 = [
    'sstar', 
    'sprime', 
    'archaicseeker2'
]
demography2 = ['HumanNeanderthalDenisovan', 'ChimpBonoboGhost']
accuracy2 = {
    'sstar': sstar_2src_accuracy_mean,
    'sprime': sprime_2src_accuracy_mean,
    'archaicseeker2': archaicseeker2_2src_accuracy_mean,
}

fig, axs = plt.subplots(nrows=4, ncols=3, constrained_layout=True, figsize=(7.5,7.5), dpi=300)
gridspec = axs[0, 0].get_subplotspec().get_gridspec()
for a in axs[:,2]:
    a.remove()

markers = {
    'nref_10_ntgt_1': {'symbol':'.', 'size': 6},
    'nref_50_ntgt_1': {'symbol':'*', 'size': 6},
    'nref_10_ntgt_10': {'symbol':'p', 'size': 4},
    'nref_50_ntgt_10': {'symbol':'d', 'size': 4},
}

colors = {
    'sstar': {'true': 'blue', 'const': 'cyan', 'ref_tgt_only': 'purple'},
    'skovhmm': 'green',
    'sprime': 'orange',
    'archaicseeker2': 'magenta',
}

linestyles = {
    'const': 'dotted',
    'true': 'solid',
    'ref_tgt_only': (0, (3, 1, 1, 1, 1, 1)),
}

titles = {
    'HumanNeanderthal': 'Human-Neanderthal model',
    'BonoboGhost': 'Bonobo-Ghost model',
    'HumanNeanderthalDenisovan': 'Human-Neanderthal-Denisovan model',
    'ChimpBonoboGhost': 'Chimpanzee-Ghost-Bonobo model',
}

auc = open(snakemake.output.auc1, 'w')
auc.write('method\tdemography\tscenario\tsample\tAUC\n')

j = 0
for d in demography1:
    for s in samples:
        for sc in scenarios:
            for m in methods1:
                if (m == 'sstar') or (m == 'skovhmm'):
                    if (s == 'nref_10_ntgt_10') or (s == 'nref_50_ntgt_10'): continue
                if m == 'sstar': color = colors[m][sc]
                else: color = colors[m]
                df = accuracy1[m][
                        (accuracy1[m]['demography'] == d) &
                        (accuracy1[m]['sample'] == s) &
                        (accuracy1[m]['scenario'] == sc)
                    ].sort_values(by='recall', ascending=False)
                recall = df['recall']
                precision = df['precision']
                if (m == 'sprime') or (m == 'skovhmm'):
                    if sc != 'true': continue
                auc_score = metrics.auc(recall/100, precision/100)
                auc.write(f'{m}\t{d}\t{sc}\t{s}\t{auc_score}\n')
                axs[j,0].plot(recall, precision,
                    marker=markers[s]['symbol'], ms=markers[s]['size'],
                    c=color)

    axs[j,0].set_xlabel('Recall (%)', fontsize=10)
    axs[j,0].set_ylabel('Precision (%)', fontsize=10)
    axs[j,0].set_xlim([-5, 105])
    axs[j,0].set_ylim([-5, 105])
    axs[j,0].set_title(titles[d], fontsize=8, weight='bold')
    if j == 0: axs[j,0].text(-35, 110, 'B', fontsize=10, weight='bold')
    if j == 1: axs[j,0].text(-35, 110, 'D', fontsize=10, weight='bold')
    j += 1

auc.close()

auc = open(snakemake.output.auc2, 'w')
auc.write('method\tdemography\tsample\tsrc\tAUC\n')

for d in demography2:
    for s in samples:
        for m in methods2:
            if (m == 'sstar') or (m == 'skovhmm'):
                if (s == 'nref_10_ntgt_10') or (s == 'nref_50_ntgt_10'): continue
            if m == 'sstar': color = colors[m]['true']
            else: color = colors[m]
            src1_df = accuracy2[m][
                         (accuracy2[m]['demography'] == d) &
                         (accuracy2[m]['sample'] == s) &
                         (accuracy2[m]['src'] == 'src1')
                     ].sort_values(by='recall', ascending=False)
            src2_df = accuracy2[m][
                         (accuracy2[m]['demography'] == d) &
                         (accuracy2[m]['sample'] == s) &
                         (accuracy2[m]['src'] == 'src2')
                     ].sort_values(by='recall', ascending=False)
            src1_recall = src1_df['recall']
            src1_precision = src1_df['precision']
            src2_recall = src2_df['recall']
            src2_precision = src2_df['precision']
            if m != 'archaicseeker2':
                src1_auc_score = metrics.auc(src1_recall/100, src1_precision/100)
                src2_auc_score = metrics.auc(src2_recall/100, src2_precision/100)
                auc.write(f'{m}\t{d}\t{s}\tsrc1\t{src1_auc_score}\n')
                auc.write(f'{m}\t{d}\t{s}\tsrc2\t{src2_auc_score}\n')
            else:
                auc.write(f'{m}\t{d}\t{s}\tsrc1\tNA\n')
                auc.write(f'{m}\t{d}\t{s}\tsrc2\tNA\n')
            axs[j,0].plot(src1_recall, src1_precision,
                     marker=markers[s]['symbol'], ms=markers[s]['size'],
                     c=color, markerfacecolor='white')
            axs[j,0].plot(src2_recall, src2_precision,
                     marker=markers[s]['symbol'], ms=markers[s]['size'],
                     c=color, linestyle='dashdot')

    axs[j,0].set_xlabel('Recall (%)', fontsize=10)
    axs[j,0].set_ylabel('Precision (%)', fontsize=10)
    axs[j,0].set_xlim([-5, 105])
    axs[j,0].set_ylim([-5, 105])
    axs[j,0].set_title(titles[d], fontsize=8, weight='bold')
    if j == 2: axs[j,0].text(-35, 110, 'F', fontsize=10, weight='bold')
    if j == 3: axs[j,0].text(-35, 110, 'H', fontsize=10, weight='bold')
    j += 1

auc.close()

auc_1src = pd.read_csv(snakemake.output.auc1, sep="\t").dropna()
auc_2src = pd.read_csv(snakemake.output.auc2, sep="\t").dropna()

i = 0
ytick_labels = []
yticks = []
for d in demography1:
    k = 0
    for s in samples:
        for sc in scenarios:
            for m in methods1:
                if (m == 'sstar') or (m == 'skovhmm'):
                    if (s == 'nref_10_ntgt_10') or (s == 'nref_50_ntgt_10'): continue
                if m == 'sstar': color = colors[m][sc]
                else: color = colors[m]
                if (m == 'sprime') or (m == 'skovhmm'):
                    if sc != 'true': continue
                auc = auc_1src[
                    (auc_1src['method'] == m) &
                    (auc_1src['demography'] == d) &
                    (auc_1src['scenario'] == sc) &
                    (auc_1src['sample'] == s)
                ]['AUC']
                bars = axs[i,1].barh(k, auc, color=color)
                axs[i,1].bar_label(bars, fontsize=5)
                if s == 'nref_10_ntgt_1': ytick_labels.append('Nref=10,Ntgt=1')
                elif s == 'nref_10_ntgt_10': ytick_labels.append('Nref=10,Ntgt=10')
                elif s == 'nref_50_ntgt_1': ytick_labels.append('Nref=50,Ntgt=1')
                elif s == 'nref_50_ntgt_10': ytick_labels.append('Nref=50,Ntgt=10')
                yticks.append(k)
                k += 1

    axs[i,1].set_xlabel('AUC', fontsize=10)
    axs[i,1].set_ylabel('Sample size', fontsize=10)
    axs[i,1].set_yticks(yticks)
    axs[i,1].set_yticklabels(ytick_labels, fontsize=5)
    axs[i,1].set_title(titles[d], fontsize=8, weight='bold')
    if i == 0: axs[i,1].text(-0.25, 12.5, 'C', fontsize=10, weight='bold')
    if i == 1: axs[i,1].text(-0.25, 12.5, 'E', fontsize=10, weight='bold')
    axs[i,1].set_xlim([0,1])
    i += 1

ytick_labels = []
yticks = []
for d in demography2:
    k = 0
    for s in samples:
        for m in methods2:
            if (m == 'sstar') or (m == 'skovhmm'):
                if (s == 'nref_10_ntgt_10') or (s == 'nref_50_ntgt_10'): continue
            if m == 'sstar': color = colors[m]['true']
            else: color = colors[m]
            if m == 'archaicseeker2': continue
            auc_src1 = auc_2src[
                    (auc_2src['method'] == m) &
                    (auc_2src['demography'] == d) &
                    (auc_2src['src'] == 'src1') &
                    (auc_2src['sample'] == s)
                ]['AUC']
            auc_src2 = auc_2src[
                    (auc_2src['method'] == m) &
                    (auc_2src['demography'] == d) &
                    (auc_2src['src'] == 'src2') &
                    (auc_2src['sample'] == s)
                ]['AUC']
            bars = axs[i,1].barh(k, auc_src1, color=color)
            axs[i,1].bar_label(bars, fontsize=5)
            if s == 'nref_10_ntgt_1': ytick_labels.append('Nref=10,Ntgt=1(src1)')
            elif s == 'nref_10_ntgt_10': ytick_labels.append('Nref=10,Ntgt=10(src1)')
            elif s == 'nref_50_ntgt_1': ytick_labels.append('Nref=50,Ntgt=1(src1)')
            elif s == 'nref_50_ntgt_10': ytick_labels.append('Nref=50,Ntgt=10(src1)')
            yticks.append(k)
            k += 1
            bars = axs[i,1].barh(k, auc_src2, color=color)
            axs[i,1].bar_label(bars, fontsize=5)
            if s == 'nref_10_ntgt_1': ytick_labels.append('Nref=10,Ntgt=1(src2)')
            elif s == 'nref_10_ntgt_10': ytick_labels.append('Nref=10,Ntgt=10(src2)')
            elif s == 'nref_50_ntgt_1': ytick_labels.append('Nref=50,Ntgt=1(src2)')
            elif s == 'nref_50_ntgt_10': ytick_labels.append('Nref=50,Ntgt=10(src2)')
            yticks.append(k)
            k += 1

    axs[i,1].set_xlabel('AUC', fontsize=10)
    axs[i,1].set_ylabel('Sample size', fontsize=10)
    axs[i,1].set_yticks(yticks)
    axs[i,1].set_yticklabels(ytick_labels, fontsize=5)
    axs[i,1].set_title(titles[d], fontsize=8, weight='bold')
    if i == 2: axs[i,1].text(-0.25, 12.5, 'G', fontsize=10, weight='bold')
    if i == 3: axs[i,1].text(-0.25, 12.5, 'I', fontsize=10, weight='bold')
    axs[i,1].set_xlim([0,1])
    i += 1

# legend
subfig = fig.add_subfigure(gridspec[:,2])
handles, labels = subfig.gca().get_legend_handles_labels()
sstar_line = plt.Line2D([0], [0], label='sstar (full)', color=colors['sstar']['true'])
sstar_line2 = plt.Line2D([0], [0], label='sstar (constant)', color=colors['sstar']['const'])
sstar_line3 = plt.Line2D([0], [0], label='sstar (only ref & tgt)', color=colors['sstar']['ref_tgt_only'])
skovhmm_line = plt.Line2D([0], [0], label='SkovHMM', color=colors['skovhmm'])
sprime_line = plt.Line2D([0], [0], label='SPrime', color=colors['sprime'])
archaicseeker2_line = plt.Line2D([0], [0], label='ArchaicSeeker 2.0', color=colors['archaicseeker2'])
nref_10_ntgt_1 = plt.Line2D([0], [0], marker=markers['nref_10_ntgt_1']['symbol'],
                            ms=5, label='Nref=10, Ntgt=1', color='black', linewidth=0)
nref_50_ntgt_1 = plt.Line2D([0], [0], marker=markers['nref_50_ntgt_1']['symbol'],
                            ms=5, label='Nref=50, Ntgt=1', color='black', linewidth=0)
nref_10_ntgt_10 = plt.Line2D([0], [0], marker=markers['nref_10_ntgt_10']['symbol'],
                             ms=5, label='Nref=10, Ntgt=10', color='black', linewidth=0)
nref_50_ntgt_10 = plt.Line2D([0], [0], marker=markers['nref_50_ntgt_10']['symbol'],
                             ms=4, label='Nref=50, Ntgt=10', color='black', linewidth=0)
src1 = plt.Line2D([0], [0], label='src1', color='black', marker='o', ms=4, markerfacecolor='white')
src2 = plt.Line2D([0], [0], label='src2', color='black', marker='o', ms=4, linestyle='dotted')

handles.extend([sstar_line, sstar_line2, sstar_line3, skovhmm_line, sprime_line, archaicseeker2_line,
                nref_10_ntgt_1, nref_50_ntgt_1, nref_10_ntgt_10, nref_50_ntgt_10, src1, src2])
subfig.legend(handles=handles, fontsize=8, handlelength=1.5)

fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0, wspace=0.1)
plt.savefig(snakemake.output.accuracy, bbox_inches='tight')
