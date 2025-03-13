from caltable import DataUnit
from caltable import Engines
from caltable import FileUnit
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.lines import Line2D
from itertools import chain
import pandas as pd
import numpy as np
import io
import base64
import math
import json

@DataUnit.register(['sasa', 'corex', 'bfactor', 'sequence_entropy', 'apl-aggregate',
                    'apl-residue-likelihood', 'apl-peptide-likelihood'])
class ProteinValuesTypeEngine(Engines.NumArrayTypeEngine):
    def __init__(self, value, iotype):
        super().__init__(value, iotype)
    
    def _plot(self, value, title=None, x_label=None, y_label=None, grid=False):
        matplotlib.use('agg')
        fig = plt.figure()
        ax = fig.subplots()
        ax.plot(list(range(len(value))), value)
        ax.set_title(self.iotype.name if title is None else title)
        if x_label is not None: ax.set_xlabel(x_label)
        if y_label is not None: ax.set_xlabel(y_label)
        if grid: ax.grid()
        return fig
        
    def _repr_markdown_(self):
        return self.view_html()
    
    def view_html(self, title=None, x_label=None, y_label=None, grid=False):
        fig = self._plot(self.value, title=title, x_label=x_label, y_label=y_label, grid=grid)
        img_stream = io.BytesIO()
        fig.savefig(img_stream, format='jpg', bbox_inches='tight')
        img_stream.seek(0)
        img_base64 = base64.b64encode(img_stream.read()).decode()
        _image_base64 = f'<img src="data:image/jpg;base64,{img_base64}">'
        fig.clear()
        plt.close(fig)
        return _image_base64
    
    def file(self, name=None, ext='csv'):
        return FileUnit(pd.DataFrame({name:self.value}).to_csv(lineterminator='\n'), name=name, ext=ext)
    
@DataUnit.register(['mhcii', 'apl-mhc-combined'])
class ProteinMHCIITableTypeEngine(Engines.StringJSONTypeEngine):
    
    def _plot(self, value, title=None, x_label=None, y_label=None, grid=False):
        matplotlib.use('agg')
        fig = plt.figure()
        ax = fig.subplots()
        _df = pd.read_json(io.StringIO(value))
        _df.plot(ax=ax)
        ax.set_title(self.iotype.name if title is None else title)
        if x_label is not None: ax.set_xlabel(x_label)
        if y_label is not None: ax.set_xlabel(y_label)
        if grid: ax.grid()
        return fig
        
    def _repr_markdown_(self):
        return self.view_html()
    
    def view_html(self, title=None, x_label=None, y_label=None, grid=False):
        fig = self._plot(self.value, title=title, x_label=x_label, y_label=y_label, grid=grid)
        img_stream = io.BytesIO()
        fig.savefig(img_stream, format='jpg', bbox_inches='tight')
        img_stream.seek(0)
        img_base64 = base64.b64encode(img_stream.read()).decode()
        _image_base64 = f'<img src="data:image/jpg;base64,{img_base64}">'
        fig.clear()
        plt.close(fig)
        return _image_base64
    
    def file(self, name=None, ext='csv'):
        _df = pd.read_json(io.StringIO(self.value))
        return FileUnit(_df.to_csv(lineterminator='\n'), name=name, ext=ext)
    
    
@DataUnit.register(['apl-table'])
class ProteinAPLTableTypeEngine(Engines.StringJSONTypeEngine):
    
    def _build_df(self, value):
        labels = {
            'APL': 'Likelihood',
            'Aggregate': 'Aggregate',
            'COREX': 'COREX',
            'SASA': 'ASA',
            'B-Factor': 'B-factor',
            'Sequence Entropy': 'Conservation',
        }
        _dict = json.loads(value)
        _max_len = max([len(_dict[key]) for key in labels])
        _dict = {label:_dict[key]+([math.nan]*(_max_len-len(_dict[key]))) for key,label in labels.items()}
        return pd.DataFrame(_dict)
    
    def _plot(self, value, title=None, x_label=None, y_label=None, grid=False):
        labels = {
            'APL': ('Likelihood', 'm-', 2),
            'Aggregate': ('Aggregate', 'b-', 2),
            'COREX': ('COREX', 'r-', 1),
            'SASA': ('ASA', 'b-', 1),
            'B-Factor': ('B-factor', 'c-', 1),
            'Sequence Entropy': ('Seq Entropy', 'y-', 1),
        }
        matplotlib.use('agg')
        data = json.loads(self.value)
        antigen_tag = data['Antigen']
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.4, 6.8), gridspec_kw={'height_ratios': [3, 1]})
        fig.subplots_adjust(hspace=0.4)
        ax.set_xlim(0, len(data['Aggregate']) + 1)
        _flatten_values = list(chain(*[data[key] for key in labels]))
        ax.set_ylim(min(_flatten_values) - 0.1, max(_flatten_values) + 0.1)
        ax.set_xlabel('Primary Sequence')
        ax.set_ylabel('Score')
        ax.set_title(f"{antigen_tag} Residue APL")
        if (len(data['APL']) > 500):
            ax.xaxis.set_major_locator(MultipleLocator(100))
        else:
            ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        
        for key, (label, linestyle, linewidth) in labels.items():
            ax.plot(range(1, len(data['APL'])+1), data[key][:len(data['APL'])], linestyle, label=label, linewidth=linewidth)
        
        ax.axhline(0, c='black', linestyle = '--', linewidth=1)
        ax.legend(loc='lower right', bbox_to_anchor=(1, 0.0), prop={'size':8})
        
        def plot_one_line(ax_, plot_score, _labels, title='', percentile=0.5):
            threshold = np.unique(plot_score)[-int(percentile*len(plot_score))]
            if threshold > 0:
                ax_.axhline(threshold, c='m', linestyle = '--', linewidth=1)
            for x, label in zip(range(len(plot_score)), _labels):
                if plot_score[x] >= threshold and label:
                    ax_.axvline(x+1, ymin=0, ymax=1.0, c='green', linewidth=2)
                elif plot_score[x] < threshold and label:
                    ax_.axvline(x+1, ymin=0, ymax=1.0, c='blue', linewidth=2)
                    
            ax_.plot(range(1,len(plot_score)+1),plot_score, linestyle = '-', linewidth=2, color='0.5')
            ax_.set_title(title)
            ax_.set_adjustable(adjustable='datalim')
            ax_.set_ylim(0, 1.1)
            label_lines = [Line2D([0], [0], color="green", lw=2),
                           Line2D([0], [0], color="blue", lw=2),]
            ax_.legend(label_lines, ['TP', 'FN'], loc=4)
        plot_score = normalize_score(data['Peptide-Likelihood'])
        percentile = float(data['APL-Threshold'])
        _labels = [i['label'] for i in json.loads(data['Regular-Mers']).values()]
        plot_one_line(bx, plot_score, _labels=_labels, title=f"{antigen_tag} Peptide APL", percentile=percentile)
        
        return fig
        
    def _repr_markdown_(self):
        return self.view_html()
    
    def view_html(self, title=None, x_label=None, y_label=None, grid=False):
        fig = self._plot(self.value, title=title, x_label=x_label, y_label=y_label, grid=grid)
        img_stream = io.BytesIO()
        fig.savefig(img_stream, format='jpg', bbox_inches='tight')
        img_stream.seek(0)
        img_base64 = base64.b64encode(img_stream.read()).decode()
        _image_base64 = f'<img src="data:image/jpg;base64,{img_base64}">'
        fig.clear()
        plt.close(fig)
        return _image_base64
    
    def file(self, name=None, ext='csv'):
        _df = self._build_df(self.value)
        return FileUnit(_df.to_csv(lineterminator='\n'), name=name, ext=ext)
    
def normalize_score(L):
    max_L = max(L)
    min_L = min(L)
    norm_L = [(x-min_L)/(max_L-min_L) for x in L]
    return norm_L

@DataUnit.register(['aplmhc-table'])
class ProteinAPLMHCTableTypeEngine(Engines.StringJSONTypeEngine):
    
    def _build_df(self, value):
        data = json.loads(value)
        mhc_df = pd.read_json(io.StringIO(data['MHC'])).set_index('Peptide')
        combine_df = pd.read_json(io.StringIO(data['APL-MHC'])).set_index('Peptide')
        combine_df.columns = [f'combined-{col}' for col in combine_df.columns]
        df = pd.merge(left=mhc_df, right=combine_df, left_index=True, right_index=True)
        df['likelihood'] = data['Peptide-Likelihood']
        return df
    
    def _plot(self, value, title=None, x_label=None, y_label=None, grid=False):
        matplotlib.use('agg')
        data = json.loads(self.value)
        mhc_df = pd.read_json(io.StringIO(data['MHC'])).set_index('Peptide')
        combine_df = pd.read_json(io.StringIO(data['APL-MHC'])).set_index('Peptide')
        fig, (axes) = plt.subplots(len(mhc_df.columns) + len(combine_df.columns) + 1, 1, sharex=True)
        fig.subplots_adjust(hspace=0.9)
        percentile = float(data['APL-Threshold'])
        antigen_tag = data['Antigen']
        labels = data['Labels']

        def plot_one_line(ax, plot_score, labels, title='', percentile=0.5):
            threshold = np.unique(plot_score)[-int(percentile*len(plot_score))]
            if threshold > 0:
                ax.axhline(threshold, c='m', linestyle = '--', linewidth=1)
            for x, label in zip(range(len(plot_score)), labels):
                if plot_score[x] >= threshold and label:
                    ax.axvline(x+1, ymin=0, ymax=1.0, c='green', linewidth=2)
                elif plot_score[x] < threshold and label:
                    ax.axvline(x+1, ymin=0, ymax=1.0, c='blue', linewidth=2)
                    
            ax.plot(range(1,len(plot_score)+1),plot_score, linestyle = '-', linewidth=2, color='0.5')
            ax.set_title(title)
            ax.set_adjustable(adjustable='datalim')
            ax.set_ylim(0, 1.1)
            label_lines = [Line2D([0], [0], color="green", lw=2),
                           Line2D([0], [0], color="blue", lw=2),]
            ax.legend(label_lines, ['TP', 'FN'], loc=4, bbox_to_anchor=(1.2,-0.2))
        plot_score = normalize_score(data['Peptide-Likelihood'])
        plot_one_line(axes[0], plot_score, labels=labels, title=f"{antigen_tag} Peptide APL", percentile=percentile)
        for i in range(1, 1+len(mhc_df.columns)):
            plot_score = normalize_score(mhc_df.iloc[:,i-1].values)
            plot_one_line(axes[i], plot_score, labels=labels, title=f"{antigen_tag} MHC-II ({mhc_df.columns[i-1]})", percentile=percentile)
        for i in range(1, 1+len(combine_df.columns)):
            plot_score = normalize_score(combine_df.iloc[:,i-1].values)
            plot_one_line(axes[i+len(mhc_df.columns)], plot_score, labels=labels, title=f"{antigen_tag} APL & MHC-II ({combine_df.columns[i-1]}) Combined", percentile=percentile)
        return fig
        
        
    def _repr_markdown_(self):
        return self.view_html()
    
    def view_html(self, title=None, x_label=None, y_label=None, grid=False):
        fig = self._plot(self.value, title=title, x_label=x_label, y_label=y_label, grid=grid)
        img_stream = io.BytesIO()
        fig.savefig(img_stream, format='jpg', bbox_inches='tight')
        img_stream.seek(0)
        img_base64 = base64.b64encode(img_stream.read()).decode()
        _image_base64 = f'<img src="data:image/jpg;base64,{img_base64}">'
        fig.clear()
        plt.close(fig)
        return _image_base64
    
    def file(self, name=None, ext='csv'):
        _df = self._build_df(self.value)
        return FileUnit(_df.to_csv(lineterminator='\n'), name=name, ext=ext)