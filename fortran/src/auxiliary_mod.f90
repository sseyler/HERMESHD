module auxiliary_mod

use parameters_mod

contains


!-----------------------------------------------------------

    real function xc(i)
        integer i
        xc = loc_lxd + (i - 0.5)*dx
    end function xc

!-----------------------------------------------------------

    real function yc(j)
        integer j
        yc = loc_lyd + (j - 0.5)*dy
    end function yc

!-----------------------------------------------------------

    real function zc(k)
        integer k
        zc = loc_lzd + (k - 0.5)*dz
    end function zc

!-----------------------------------------------------------

    real function rz(i,j)
        integer i,j
        rz = sqrt(yc(j)**2)
    end function rz

!-----------------------------------------------------------

    real function r(i,j)
        integer i,j,k
        r = sqrt(xc(i)**2 + yc(j)**2)
    end function r

!-----------------------------------------------------------

    real function theta(i,j)
        integer i,j
        theta = atan2(yc(j),xc(i))
    end function theta

!-----------------------------------------------------------

    real function xvtk(i)
        integer i
        xvtk = loc_lxd + (i - 0.5)*dxvtk
    end function xvtk

!-----------------------------------------------------------

    real function yvtk(j)
        integer j
        yvtk = loc_lyd + (j - 0.5)*dyvtk
    end function yvtk

!-----------------------------------------------------------

    real function zvtk(k)
        integer k
        zvtk = loc_lzd + (k - 0.5)*dzvtk
    end function zvtk

!-----------------------------------------------------------

    real function rvtk(i,j)
        integer i,j,k
        rvtk = sqrt(xvtk(i)**2 + yvtk(j)**2)
    end function rvtk

!-----------------------------------------------------------

    real function thetavtk(i,j)
        integer i,j
        thetavtk = atan2(yvtk(j),xvtk(i))
    end function thetavtk


!------------More realistic COBRA  current driver -----------

    real function E_z(t)
        real t
        if(t.le.tr) E_z = Ez0*sin(0.5*pi*t/tr)*cos(0.5*pi*t/tr)
        if(t.ge.tr) E_z = 0.
    end function E_z

!------------More realistic COBRA  current driver -----------


    real function Icur(t)
        real t
        if(t.le.tr) Icur = Ipeak*sin(0.5*pi*t/tr)**2
        if(t.ge.tr) Icur = Ipeak
!        if(t.ge.tr.and.t.le.2*tr)Icur = Ipeak
!        if(t.ge.2*tr.and.t.le.3*tr)Icur = Ipeak*sin(0.5*pi*(t-tr)/tr)
!        if(t.ge.3*tr)Icur = 0.
    end function Icur

!-----------------------------------------------------------


end module auxiliary_mod
