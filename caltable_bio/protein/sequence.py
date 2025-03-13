import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from caltable import DataUnit
from caltable import Engines
from caltable import FileUnit

@DataUnit.register(['fasta'])
class SequenceAlignmentTypeEngine(Engines.StringTypeEngine):
    def __init__(self, value, iotype):
        super().__init__(value, iotype)
    
    def _fetch_sequences(self, _alignment):
        return dict(re.findall(r">(.*?)\n([^>]+)\n", _alignment))
    def _generate_heatmap_data(self, sequences):
        """Converts sequences into numerical matrix for visualization."""
        residues = sorted(set("".join(sequences.values())))  # Unique residues
        residue_to_int = {res: i for i, res in enumerate(residues)}  # Map residue to int
        _max_len = max([len(sequence) for sequence in sequences.values()])
        # Convert sequences to a numerical matrix
        alignment_matrix = [
            [residue_to_int[res] for res in seq] + [-1] *(_max_len - len(seq))
            for seq in sequences.values()
        ]
        residues.insert(0, "*")
        return np.array(alignment_matrix), list(sequences.keys()), residues
    def _plot_sequences(self, sequences):
        alignment_matrix, sequence_names, residues = self._generate_heatmap_data(sequences)
        fig = go.Figure(
            data=go.Heatmap(
                z=alignment_matrix, 
                x=list(range(alignment_matrix.shape[1])),  # Columns represent positions
                y=sequence_names,  # Rows are sequence names
                colorscale="Viridis",  # Color scheme
                colorbar=dict(
                    title="Residues",
                    tickvals=list(range(len(residues))),
                    ticktext=residues,
                )
            )
        )
        # Update layout for better readability
        fig.update_layout(
            title="FASTA Alignment Visualization",
            xaxis_title="Position",
            yaxis_title="Sequences",
            yaxis=dict(showgrid=False, automargin=True),
            xaxis=dict(
                title="Position",
                tickmode="array",
                showgrid=False,
                scaleanchor="y",  # Lock x and y axis scales to make squares
                scaleratio=1
            ),
        )

        return fig
    def _plot(self, value):
        _sequences = self._fetch_sequences(value)
        fig = self._plot_sequences(_sequences)
        return fig
        
    def _repr_markdown_(self):
        self._plot(self.value).show()
        return f'{pd.DataFrame(list(self._fetch_sequences(self.value).items()), columns=["ID", "Sequence"])._repr_html_()}'
    
    def view_html(self, **kwargs):
        return self._plot(self.value).to_html()
    
    def file(self, name=None, ext='fasta'):
        return super().file(name, ext)
    
    
@DataUnit.register(['protein-seq'])
class SequenceTypeEngine(Engines.StringTypeEngine):
    def __init__(self, value, iotype):
        super().__init__(value, iotype)
    
    def _render_sequence(self, sequence, fragment_length=10, column_num=6):
        fragments = [sequence[i:i+fragment_length] for i in range(0, len(sequence), fragment_length)]
        index_content, fragment_content = [], []
        for index, fragment in enumerate(fragments):
            _index = str(index)
            _index = f'<b>{_index}</b>' + ''.join(['&nbsp;']*(fragment_length - len(_index)))
            index_content.append(_index)
            fragment_content.append(f'<span style="background-color:grey;"><b>{fragment[0]}</b></span>{fragment[1:]}')
        index_content = [index_content[i:i+column_num] for i in range(0, len(index_content), column_num)]
        fragment_content = [fragment_content[i:i+column_num] for i in range(0, len(fragment_content), column_num)]
        _content = '<br><br>'.join(['&nbsp;'.join(ind) + '<br>' + '&nbsp;'.join(frag) for (ind, frag) in zip(index_content, fragment_content)])
        _content = f'<div style="overflow:scroll; font-family: Courier,monospace;">{_content}</div>'
        return _content
        
    def _repr_markdown_(self):
        return self._render_sequence(self.value, fragment_length=10, column_num=6)
    
    def view_html(self, **kwargs):
        return self._render_sequence(self.value, fragment_length=10, column_num=6)
    
    def file(self, name=None, ext='fasta'):
        _data = f'>{name}\n{self.value}'
        return FileUnit(_data, name=name, ext=ext)
    
@DataUnit.register(['protein-regulared-peptide'])
class ProteinRegularMersTypeEngine(Engines.StringJSONTypeEngine):
    pass

@DataUnit.register(['protein-peptides'])
class ProteinMersTypeEngine(Engines.StringTypeEngine):
    pass