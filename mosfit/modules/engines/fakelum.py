"""Definitions for the `FakeLum` class."""

# Important: Only define one ``Module`` class per file.
from mosfit.modules.engines.engine import Engine


class FakeLum(Engine):
    """
    Honestly I'm out of ideas on how to make Kasen SEDs work without
    luminosities so I'm just going to make some fake ones and see if
    that helps :(
    """

    def __init__(self, **kwargs):
        """Initialize module."""
        super(Engine, self).__init__(**kwargs)
        self._wants_dense = True # This function wants a dense times array
    

    def process(self, **kwargs):

        self._times = kwargs[self.key('dense_times')]
        luminosities = [ 0.0 for t in self._times]
        
        return {self.dense_key('luminosities'): luminosities}
