import _hermeshd
import f90wrap.runtime
import logging

class Hermeshd(f90wrap.runtime.FortranModule):
    """
    Module hermeshd
    
    
    Defined at hermeshd.f90 lines 2-217
    
    """
    @staticmethod
    def main(comm):
        """
        main(comm)
        
        
        Defined at hermeshd.f90 lines 29-51
        
        Parameters
        ----------
        comm : int
        
        -----------------------------
        """
        _hermeshd.f90wrap_main(comm=comm)
    
    @staticmethod
    def step(q_io, q_1, q_2, t, dt):
        """
        step(q_io, q_1, q_2, t, dt)
        
        
        Defined at hermeshd.f90 lines 57-65
        
        Parameters
        ----------
        q_io : float array
        q_1 : float array
        q_2 : float array
        t : float
        dt : float
        
        """
        _hermeshd.f90wrap_step(q_io=q_io, q_1=q_1, q_2=q_2, t=t, dt=dt)
    
    @staticmethod
    def setup(q_io, t, dt, t1, t_start, dtout, nout, comm):
        """
        setup(q_io, t, dt, t1, t_start, dtout, nout, comm)
        
        
        Defined at hermeshd.f90 lines 71-124
        
        Parameters
        ----------
        q_io : float array
        t : float
        dt : float
        t1 : float
        t_start : float
        dtout : float
        nout : int
        comm : int
        
        -------------------------------------------------
         1. Initialize general simulation variables
        -------------------------------------------------
        """
        _hermeshd.f90wrap_setup(q_io=q_io, t=t, dt=dt, t1=t1, t_start=t_start, \
            dtout=dtout, nout=nout, comm=comm)
    
    @staticmethod
    def cleanup(t_start):
        """
        cleanup(t_start)
        
        
        Defined at hermeshd.f90 lines 130-148
        
        Parameters
        ----------
        t_start : float
        
        -------------------------------------------------
         1. De-allocate system resources for RNG
        -------------------------------------------------
        """
        _hermeshd.f90wrap_cleanup(t_start=t_start)
    
    @staticmethod
    def generate_output(q_r, t, dt, t1, dtout, nout):
        """
        generate_output(q_r, t, dt, t1, dtout, nout)
        
        
        Defined at hermeshd.f90 lines 156-215
        
        Parameters
        ----------
        q_r : float array
        t : float
        dt : float
        t1 : float
        dtout : float
        nout : int
        
        """
        _hermeshd.f90wrap_generate_output(q_r=q_r, t=t, dt=dt, t1=t1, dtout=dtout, \
            nout=nout)
    
    _dt_array_initialisers = []
    

hermeshd = Hermeshd()

