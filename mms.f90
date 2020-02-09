module mms
  use prec_def
  implicit none
  real(kind = dp), parameter :: pi = acos(-1.d0)
  real(kind = dp), parameter :: lam_0 = 1.d0
  real(kind = dp), parameter :: mu_0  = 1.d0
  real(kind = dp), parameter :: omega_0  = sqrt(2.d0)
  real(kind = dp), parameter :: kx_0  = 1.d0
  real(kind = dp), parameter :: ky_0  = 1.d0
  real(kind = dp), parameter :: x_off_0  =  0.d0
  real(kind = dp), parameter :: y_off_0  =  0.d0
  real(kind = dp), parameter :: amp_u1_0 =  1.d0
  real(kind = dp), parameter :: amp_u2_0 = -1.d0

contains

  real(kind = dp)  function u_mms(x,y,t,ivar,idt)
    implicit none
    integer  :: ivar,idt,q
    real(kind = dp) :: x,y,t
    if(ivar.eq.1) then
     u_mms = amp_u1_0*sin(kx_0*x+x_off_0)*sin(ky_0*y+y_off_0)
    else
     u_mms = amp_u2_0*sin(kx_0*x+x_off_0)*sin(ky_0*y+y_off_0)
    end if
    q = MODULO(idt,4)
    SELECT CASE (q)
    CASE (3)
     u_mms = sin(omega_0*t)*u_mms
    CASE (0)
     u_mms = cos(omega_0*t)*u_mms
    CASE (1)
     u_mms = -sin(omega_0*t)*u_mms
    CASE (2)
     u_mms = -cos(omega_0*t)*u_mms
    END SELECT
    if (idt.gt.0) u_mms = u_mms*omega_0**idt
    return
  end function u_mms

end module mms
