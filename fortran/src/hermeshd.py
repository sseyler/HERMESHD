import _hermeshd
import f90wrap.runtime
import logging

class Hermeshd(f90wrap.runtime.FortranModule):
    """
    Module hermeshd
    
    
    Defined at hermeshd.f90 lines 2-269
    
    """
    @staticmethod
    def main(comm):
        """
        main(comm)
        
        
        Defined at hermeshd.f90 lines 27-48
        
        Parameters
        ----------
        comm : int
        
        -----------------------------
        """
        _hermeshd.f90wrap_main(comm=comm)
    
    _dt_array_initialisers = []
    

hermeshd = Hermeshd()

