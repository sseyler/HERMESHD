import _hermeshd
import f90wrap.runtime
import logging

class HermesHD(f90wrap.runtime.FortranModule):

    def __init__(self, comm=None):
        self.comm = comm


    @staticmethod
    def setup(self):


    @staticmethod
    def run(self):


    @staticmethod
    def cleanup(self):



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


hermeshd = HermesHD()
