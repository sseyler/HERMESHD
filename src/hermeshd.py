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

class Spatial(f90wrap.runtime.FortranModule):
    """
    Module spatial
    
    
    Defined at spatial.f90 lines 1-119
    
    """
    @staticmethod
    def xc(i):
        """
        xc = xc(i)
        
        
        Defined at spatial.f90 lines 90-93
        
        Parameters
        ----------
        i : int
        
        Returns
        -------
        xc : float
        
        """
        xc = _hermeshd.f90wrap_xc(i=i)
        return xc
    
    @staticmethod
    def yc(j):
        """
        yc = yc(j)
        
        
        Defined at spatial.f90 lines 95-98
        
        Parameters
        ----------
        j : int
        
        Returns
        -------
        yc : float
        
        """
        yc = _hermeshd.f90wrap_yc(j=j)
        return yc
    
    @staticmethod
    def zc(k):
        """
        zc = zc(k)
        
        
        Defined at spatial.f90 lines 100-103
        
        Parameters
        ----------
        k : int
        
        Returns
        -------
        zc : float
        
        """
        zc = _hermeshd.f90wrap_zc(k=k)
        return zc
    
    @property
    def q_r0(self):
        """
        Element q_r0 ftype=real pytype=float
        
        
        Defined at spatial.f90 line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hermeshd.f90wrap_spatial__array__q_r0(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            q_r0 = self._arrays[array_handle]
        else:
            q_r0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hermeshd.f90wrap_spatial__array__q_r0)
            self._arrays[array_handle] = q_r0
        return q_r0
    
    @q_r0.setter
    def q_r0(self, q_r0):
        self.q_r0[...] = q_r0
    
    @property
    def q_r1(self):
        """
        Element q_r1 ftype=real pytype=float
        
        
        Defined at spatial.f90 line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hermeshd.f90wrap_spatial__array__q_r1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            q_r1 = self._arrays[array_handle]
        else:
            q_r1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hermeshd.f90wrap_spatial__array__q_r1)
            self._arrays[array_handle] = q_r1
        return q_r1
    
    @q_r1.setter
    def q_r1(self, q_r1):
        self.q_r1[...] = q_r1
    
    @property
    def q_r2(self):
        """
        Element q_r2 ftype=real pytype=float
        
        
        Defined at spatial.f90 line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hermeshd.f90wrap_spatial__array__q_r2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            q_r2 = self._arrays[array_handle]
        else:
            q_r2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hermeshd.f90wrap_spatial__array__q_r2)
            self._arrays[array_handle] = q_r2
        return q_r2
    
    @q_r2.setter
    def q_r2(self, q_r2):
        self.q_r2[...] = q_r2
    
    @property
    def q_r3(self):
        """
        Element q_r3 ftype=real pytype=float
        
        
        Defined at spatial.f90 line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hermeshd.f90wrap_spatial__array__q_r3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            q_r3 = self._arrays[array_handle]
        else:
            q_r3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hermeshd.f90wrap_spatial__array__q_r3)
            self._arrays[array_handle] = q_r3
        return q_r3
    
    @q_r3.setter
    def q_r3(self, q_r3):
        self.q_r3[...] = q_r3
    
    @property
    def dz(self):
        """
        Element dz ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dz()
    
    @dz.setter
    def dz(self, dz):
        _hermeshd.f90wrap_spatial__set__dz(dz)
    
    @property
    def dy(self):
        """
        Element dy ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dy()
    
    @dy.setter
    def dy(self, dy):
        _hermeshd.f90wrap_spatial__set__dy(dy)
    
    @property
    def dx(self):
        """
        Element dx ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dx()
    
    @dx.setter
    def dx(self, dx):
        _hermeshd.f90wrap_spatial__set__dx(dx)
    
    @property
    def dxi(self):
        """
        Element dxi ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dxi()
    
    @dxi.setter
    def dxi(self, dxi):
        _hermeshd.f90wrap_spatial__set__dxi(dxi)
    
    @property
    def dyi(self):
        """
        Element dyi ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dyi()
    
    @dyi.setter
    def dyi(self, dyi):
        _hermeshd.f90wrap_spatial__set__dyi(dyi)
    
    @property
    def dzi(self):
        """
        Element dzi ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dzi()
    
    @dzi.setter
    def dzi(self, dzi):
        _hermeshd.f90wrap_spatial__set__dzi(dzi)
    
    @property
    def dvi(self):
        """
        Element dvi ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 17
        
        """
        return _hermeshd.f90wrap_spatial__get__dvi()
    
    @dvi.setter
    def dvi(self, dvi):
        _hermeshd.f90wrap_spatial__set__dvi(dvi)
    
    @property
    def loc_lxd(self):
        """
        Element loc_lxd ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 20
        
        """
        return _hermeshd.f90wrap_spatial__get__loc_lxd()
    
    @loc_lxd.setter
    def loc_lxd(self, loc_lxd):
        _hermeshd.f90wrap_spatial__set__loc_lxd(loc_lxd)
    
    @property
    def loc_lyd(self):
        """
        Element loc_lyd ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 20
        
        """
        return _hermeshd.f90wrap_spatial__get__loc_lyd()
    
    @loc_lyd.setter
    def loc_lyd(self, loc_lyd):
        _hermeshd.f90wrap_spatial__set__loc_lyd(loc_lyd)
    
    @property
    def loc_lzd(self):
        """
        Element loc_lzd ftype=real  pytype=float
        
        
        Defined at spatial.f90 line 20
        
        """
        return _hermeshd.f90wrap_spatial__get__loc_lzd()
    
    @loc_lzd.setter
    def loc_lzd(self, loc_lzd):
        _hermeshd.f90wrap_spatial__set__loc_lzd(loc_lzd)
    
    def __str__(self):
        ret = ['<spatial>{\n']
        ret.append('    q_r0 : ')
        ret.append(repr(self.q_r0))
        ret.append(',\n    q_r1 : ')
        ret.append(repr(self.q_r1))
        ret.append(',\n    q_r2 : ')
        ret.append(repr(self.q_r2))
        ret.append(',\n    q_r3 : ')
        ret.append(repr(self.q_r3))
        ret.append(',\n    dz : ')
        ret.append(repr(self.dz))
        ret.append(',\n    dy : ')
        ret.append(repr(self.dy))
        ret.append(',\n    dx : ')
        ret.append(repr(self.dx))
        ret.append(',\n    dxi : ')
        ret.append(repr(self.dxi))
        ret.append(',\n    dyi : ')
        ret.append(repr(self.dyi))
        ret.append(',\n    dzi : ')
        ret.append(repr(self.dzi))
        ret.append(',\n    dvi : ')
        ret.append(repr(self.dvi))
        ret.append(',\n    loc_lxd : ')
        ret.append(repr(self.loc_lxd))
        ret.append(',\n    loc_lyd : ')
        ret.append(repr(self.loc_lyd))
        ret.append(',\n    loc_lzd : ')
        ret.append(repr(self.loc_lzd))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

spatial = Spatial()

