from caltable import DataUnit
from caltable import Engines
import pandas as pd
import py3Dmol


@DataUnit.register('pdb')
class ProteinPDBTypeEngine(Engines.StringTableTypeEngine):
    def __init__(self, value, iotype):
        super().__init__(value, iotype, sep=' +', end='\n')
    @property
    def preview(self): return f'PDB:{len(self)} lines'
    
    def _plot(self, value, width=400, height=300):
        view = py3Dmol.view(width=width, height=height)
        view.addModelsAsFrames(value)
        view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
        view.zoomTo()
        return view
        
    def _repr_markdown_(self):
        self._plot(self.value, width=400, height=300).show()
        return f'{pd.DataFrame(self.table_value)._repr_html_()}'
    
    def view_html(self, width=400, height=300, **kwargs):
        return self._plot(self.value, width=width, height=height).write_html()
    
    def file(self, name=None, ext='pdb'):
        return super().file(name, ext)
    