import _hermeshd
import f90wrap.runtime
import logging

class Hermeshd(f90wrap.runtime.FortranModule):
    """
    Module hermeshd
    
    
    Defined at hermeshd.f90 lines 2-263
    
    """
    @staticmethod
    def main():
        """
        main()
        
        
        Defined at hermeshd.f90 lines 27-45
        
        
        -----------------------------
        """
        _hermeshd.f90wrap_main()
    
    _dt_array_initialisers = []
    

hermeshd = Hermeshd()

