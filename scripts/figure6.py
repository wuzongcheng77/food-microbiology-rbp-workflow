"""Figure 6: static publication exports from audited CSVs; no statistical computation."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):
    with p.open(encoding='utf-8',newline='') as f: return list(csv.DictReader(f))

def label_panel_c(ax, bars, percentages):
    """Place high-frequency labels inside; retain outside labels for shorter bars."""
    return [ax.annotate(f'{value:.2f}%',
                xy=(bar.get_x()+bar.get_width(),bar.get_y()+bar.get_height()/2),
                xytext=(-3 if value >= 95 else 3,0),textcoords='offset points',
                ha='right' if value >= 95 else 'left',va='center',fontsize=9,
                color='white' if value >= 95 else 'black')
            for bar,value in zip(bars,percentages)]

def make(source,out):
    out.mkdir(parents=True,exist_ok=False)
    rows=read(source/'Supplementary_Table_S4.csv'); by={r['analysis']:r for r in rows}
    freq=read(source/'candidate_inclusion_frequency.csv')
    nonzero=[r for r in freq if float(r['top5_frequency'])>0]
    nonzero.sort(key=lambda r:(-float(r['top5_frequency']),int(r['frozen_rank'])))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'svg.hashsalt':'food-microbiology-rbp-v1','axes.spines.top':False,
                         'axes.spines.right':False,'axes.titlesize':12})
    fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    names=['target_context','docking_support','solubility_no_chaperone','E7_residual']
    labels=['Target context','Docking support','Solubility / no chaperone','E7 residual']
    def bars(ax,title,labs,vals,xlabel,limit=None):
        b=ax.barh(labs,vals,color='#4477AA',edgecolor='#243746',linewidth=.6)
        ax.invert_yaxis(); ax.set_title(title,loc='left',fontweight='bold',pad=12)
        ax.set_xlabel(xlabel); ax.set_xlim(0,limit or max(vals)*1.22)
        ax.bar_label(b,labels=[str(v) for v in vals],padding=4,fontsize=9)
        ax.set_axisbelow(True); ax.grid(axis='x',color='#E3E3E3',linewidth=.6)
    bars(axs[0,0],'A  Leave-one-component-out',labels,
         [int(by['without_'+n]['target_rank']) for n in names],'Q9T0Z9 rank (lower is better)')
    bars(axs[0,1],'B  Single-component scoring baselines',labels,
         [int(by['only_'+n]['target_rank']) for n in names],'Q9T0Z9 rank (lower is better)')
    vals=[100*float(r['top5_frequency']) for r in nonzero]
    ax=axs[1,0]; b=ax.barh([r['accession'] for r in nonzero],vals,color='#4477AA',edgecolor='#243746',linewidth=.6)
    ax.invert_yaxis(); ax.set_xlim(0,120); ax.set_xticks([0,25,50,75,100])
    ax.set_title('C  Top5 inclusion under weight perturbation',loc='left',fontweight='bold',pad=12)
    ax.set_xlabel('Inclusion frequency (%) across 10,000 trials')
    label_panel_c(ax,b,vals)
    ax.set_axisbelow(True); ax.grid(axis='x',color='#E3E3E3',linewidth=.6)
    bars(axs[1,1],'D  Two-layer masking audit',['Full','Identifier masked','Strict context blind'],
        [int(by[n]['target_rank']) for n in ['full','identifier_masked','strict_context_blind']],
        'Q9T0Z9 rank (lower is better)')
    n=rows[0]['candidates']
    fig.suptitle(f'Figure 6 | Frozen G_pre1 robustness across all {n} candidates',fontsize=16,fontweight='bold')
    fig.supxlabel('Scoring-source comparisons; no classifier accuracy is estimated. Panel C shows all nonzero frequencies.',fontsize=10)
    fig.savefig(out/'Figure6.png',dpi=600)
    fig.savefig(out/'Figure6.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    fig.savefig(out/'Figure6.svg',metadata={'Date':None})
    fig.savefig(out/'Figure6_preview.png',dpi=120)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();make(a.source,a.out)
