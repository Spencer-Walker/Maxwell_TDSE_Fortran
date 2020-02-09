module stuff_module
  use prec_def
  use mms
  integer, parameter :: n_plot = 10
  integer, parameter :: nxl = 160+1
  integer, parameter :: nyl = 160+1
  integer, parameter :: nvar = 3
  integer, parameter :: nvar_p = 2
  integer, parameter :: nvar_j = 2
  real(dp), parameter :: xmin = -1.5d0*pi
  real(dp), parameter :: xmax = 1.5d0*pi
  real(dp), parameter :: ymin = -1.5d0*pi
  real(dp), parameter :: ymax = 1.5d0*pi
  real(dp), parameter :: CFL = 0.3d0
  real(dp), parameter :: Tfinal = 10.d0
  real(dp), parameter :: tau = 1.0d0
  ! SDC Stuff 
  integer, parameter :: nsdc = 2
  integer, parameter :: ncor = 2*nsdc
  ! Realistic material parameters ...
  real(dp), parameter :: mu0 = 1.d0
  real(dp), parameter :: ep0 = 1.d0
  real(dp), parameter :: epinf = 1.d0
  real(dp), parameter :: tau_inv = 0.d0
  real(dp), parameter :: omgega_0_square = 1.d0
  real(dp), parameter :: omgega_p_square = 1000.d0
  real(dp), parameter :: a2 = 1000.d0
  real(dp), parameter :: a4 = 1000.d0
  !
  real(dp) :: D1x_1(6,9)
  real(dp) :: D1x_N(6,9)
  real(dp) :: D1x_mid(-3:3)
  real(dp) :: Hx_1(6),HIx_1(6)
  real(dp) :: Hx_N(6),HIx_N(6)
  real(dp) :: D1y_1(6,9)
  real(dp) :: D1y_N(6,9)
  real(dp) :: D1y_mid(-3:3)
  real(dp) :: Hy_1(6),HIy_1(6)
  real(dp) :: Hy_N(6),HIy_N(6)
  !
  real(dp) :: IntMatrix(nsdc+1,nsdc)
  real(dp) :: dtm(nsdc+1)

end module stuff_module

program main
  !
  !
  ! This program solves
  !              \mu_0 (H_z)_t + (E_y)_x - (E_x)_y = 0,
  ! \eps_0 \eps_\infty (E_x)_t - (H_z)_y = -\eps_0 (J_x),
  ! \eps_0 \eps_\infty (E_y)_t + (H_z)_x = -\eps_0 (J_y),
  !                              (P_x)_t = (J_x),
  !                              (P_y)_t = (J_y),
  ! (J_x)_t + \tau^{-1} (J_x) + \omega_0^2 (P_x) F(P) = \omega_p^2 (E_x),
  ! (J_y)_t + \tau^{-1} (J_y) + \omega_0^2 (P_y) F(P) = \omega_p^2 (E_y).
  !
  ! Here F(P) is a scalar function of P = \sqrt((P_x)^2+(P_y)^2).
  !
  ! We use a PEC and set homogenuous Dirichlet bc for
  ! E_x on top and bottom and for E_y to the left and right.
  !
  ! Spatial discretization is done with a 6th order accurate
  ! SBP-SAT method.
  ! Temporal discretization is by leap-frog
  !
  !
  use prec_def
  use MMS
  use stuff_module
  implicit none
  !
  real(dp) :: x(nxl),y(nyl)
  real(dp) :: u(nxl,nyl,nvar)
  real(dp) :: p(nxl,nyl,nvar_p)
  real(dp) :: jf(nxl,nyl,nvar_j)
  real(dp) :: hx,hy,dt,t
  real(dp) :: w(nxl,nyl)
  real(dp) :: energy
  character(7) :: charit
  integer :: i,j,it,nt,ivar,row,col,m,k
  real(dp) :: compute_l2_sq,time1,time2

  real(dp), allocatable, dimension(:,:,:,:) :: FEu_sdc,UK_sdc
  real(dp), allocatable, dimension(:,:,:,:) :: FEp_sdc,PK_sdc
  real(dp), allocatable, dimension(:,:,:,:) :: FEj_sdc,JK_sdc
  real(dp), allocatable, dimension(:,:,:)   :: u_sdc,du_sdc,dfu_sdc
  real(dp), allocatable, dimension(:,:,:)   :: p_sdc,dp_sdc,dfp_sdc
  real(dp), allocatable, dimension(:,:,:)   :: j_sdc,dj_sdc,dfj_sdc

  allocate(FEu_sdc(nxl,nyl,nvar,nsdc),FEp_sdc(nxl,nyl,nvar_p,nsdc), &
    FEj_sdc(nxl,nyl,nvar_j,nsdc),UK_sdc(nxl,nyl,nvar,nsdc),&
    PK_sdc(nxl,nyl,nvar_p,nsdc),JK_sdc(nxl,nyl,nvar_j,nsdc))
  allocate(u_sdc(nxl,nyl,nvar),du_sdc(nxl,nyl,nvar),dfu_sdc(nxl,nyl,nvar))
  allocate(p_sdc(nxl,nyl,nvar_p),dp_sdc(nxl,nyl,nvar_p),dfp_sdc(nxl,nyl,nvar_p))
  allocate(j_sdc(nxl,nyl,nvar_j),dj_sdc(nxl,nyl,nvar_j),dfj_sdc(nxl,nyl,nvar_j))

  call kalle_anka(nsdc,dtm,IntMatrix)
  !
  hx = (xmax-xmin)/real(nxl-1,dp)
  hy = (ymax-ymin)/real(nyl-1,dp)

  call SBP_6_order(Hx_1,Hx_N,Hix_1,Hix_N,D1x_1,D1x_N,D1x_mid,hx)
  call SBP_6_order(Hy_1,Hy_N,Hiy_1,Hiy_N,D1y_1,D1y_N,D1y_mid,hy)
  do i = 1,nxl
   x(i) = xmin + real(i-1,dp)*hx
  end do
  do j = 1,nyl
   y(j) = ymin + real(j-1,dp)*hy
  end do
  call printdble1d(x,nxl,"x.txt")
  call printdble1d(y,nyl,"y.txt")
  t = 0.d0
  u = 0.d0
  p = 0.d0
  jf = 0.d0
  ivar = 1
  call random_number(u(:,:,1))
  u(:,:,1) = 0.1d0*(u(:,:,1) - 0.5d0)
  do j = 1,nyl
   do i = 1,nxl
    u(i,j,ivar) = u(i,j,ivar) + u_mms(x(i),y(j),t,ivar,0) + 0.d0*exp(-10.d0*(x(i)**2+y(j)**2)) !
   end do
  end do
  
  
  WRITE(charit,"(I7.7)") nxl
  open(21,file = 'result'//charit//'.txt',status='unknown')
  
  ! call readdble2d(u(:,:,1),1,nxl,1,nyl,'h.txt')

  energy = 0.d0
  w = u(:,:,1)**2
  energy = energy + mu0*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  do ivar = 2,nvar
   w = u(:,:,ivar)**2
   energy = energy + ep0*epinf*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  end do
  do ivar = 1,nvar_j
   w = jf(:,:,ivar)**2
   energy = energy + (ep0/omgega_p_square)*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  end do
  do ivar = 1,nvar_p
   w = p(:,:,ivar)**2
   energy = energy + &
     (ep0*omgega_0_square/omgega_p_square)*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  end do
  if (a2 .gt. 0.d0 ) then
   w = (p(:,:,1)**2+p(:,:,2)**2)**2
   energy = energy + &
     (ep0*omgega_0_square/omgega_p_square)*a2*0.5d0*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  end if
  if (a4 .gt. 0.d0 ) then
   w = (p(:,:,1)**2+p(:,:,2)**2)**3
   energy = energy + &
     (ep0*omgega_0_square/omgega_p_square)*(a4/3.d0)&
     *compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  end if

  write(*,*) energy

  write(21,fmt='(4(E24.16))') t+dt,maxval(abs(u(:,:,ivar)-w)),&
        sqrt(hx*hy*sum((u(2:nxl,2:nyl,ivar)-w(2:nxl,2:nyl))**2)),energy 
  write(21,fmt='(4(E24.16))') 0.d0,u((nxl+1)/2,(nyl+1)/2,1),p((nxl+1)/2,(nyl+1)/2,1),energy 
  write(21,fmt='(5(E24.16))') 0.d0,u((nxl-1)/4+1,(nxl-1)/4+1,1),p((nxl-1)/4+1,(nxl-1)/4+1,1),jf((nxl-1)/4+1,(nxl-1)/4+1,1),energy 

  dt = CFL*min(hx,hy) ! DEAA fix up with speed of light...
  nt = ceiling(Tfinal/dt)
  dt = Tfinal / real(nt,dp)
  ! Scale integration matrix and substeps by dt.
  do row = 1,nsdc+1
   dtm(row) = dtm(row)*dt
   do col = 1,nsdc
    IntMatrix(row,col) = IntMatrix(row,col)*dt
   end do
  end do

  call cpu_time(time1)
  
  do it = 1,nt

   t = real(it-1,dp)*dt
   write(*,*) t
   ! Special first step to get to first Gauss point
   call compute_time_derivatives(FEu_sdc(:,:,:,1),FEp_sdc(:,:,:,1),FEj_sdc(:,:,:,1),&
     u,p,jf)
   Uk_sdc(:,:,:,1) = u  + dtm(1)*FEu_sdc(:,:,:,1)
   Pk_sdc(:,:,:,1) = p  + dtm(1)*FEp_sdc(:,:,:,1)
   Jk_sdc(:,:,:,1) = jf + dtm(1)*FEj_sdc(:,:,:,1)

   call compute_time_derivatives(FEu_sdc(:,:,:,1),FEp_sdc(:,:,:,1),FEj_sdc(:,:,:,1),&
     Uk_sdc(:,:,:,1),Pk_sdc(:,:,:,1),Jk_sdc(:,:,:,1))
   ! Prediction phase. Here we use first order forward Euler to generate
   ! solution on the Gauss points. This first-order accurate
   ! solution will then be improved during the upcoming correction phase.
   do m = 1,nsdc-1
    Uk_sdc(:,:,:,m+1) = Uk_sdc(:,:,:,m)+dtm(m+1)*FEu_sdc(:,:,:,m)
    pk_sdc(:,:,:,m+1) = pk_sdc(:,:,:,m)+dtm(m+1)*FEp_sdc(:,:,:,m)
    jk_sdc(:,:,:,m+1) = jk_sdc(:,:,:,m)+dtm(m+1)*FEj_sdc(:,:,:,m)
    call compute_time_derivatives(FEu_sdc(:,:,:,m+1),FEp_sdc(:,:,:,m+1),&
      FEj_sdc(:,:,:,m+1),Uk_sdc(:,:,:,m+1),Pk_sdc(:,:,:,m+1),Jk_sdc(:,:,:,m+1))
   end do
   !  End predicition phase.

   !  Start correction phase. To get the max possible
   !  accuracy 2nsdc-2, we choose ncor = 2nsdc-3 above so there will be ncor
   !  (0 through ncor-1) sweeps through this phase. After prediction phase
   !  we have a first order accurate solution. One correction should give
   !  second order accuracy. Two should give third order, etc. So 2nsdc-3
   !  is order 2nsdc-2. Can choose ncor smaller, but fewer corrections
   !  means correspondingly less accuracy.
   !
   do k = 1,ncor
    ! Again, special first step
    DFu_sdc = FEu_sdc(:,:,:,1)
    DFp_sdc = FEp_sdc(:,:,:,1)
    DFj_sdc = FEj_sdc(:,:,:,1)
    !
    Uk_sdc(:,:,:,1) = u  + IntMatrix(1,1)*FEu_sdc(:,:,:,1)
    pk_sdc(:,:,:,1) = p  + IntMatrix(1,1)*FEp_sdc(:,:,:,1)
    jk_sdc(:,:,:,1) = jf + IntMatrix(1,1)*FEj_sdc(:,:,:,1)
    do col = 2,nsdc
     Uk_sdc(:,:,:,1) = Uk_sdc(:,:,:,1) +IntMatrix(1,col)*FEu_sdc(:,:,:,col)
     Pk_sdc(:,:,:,1) = Pk_sdc(:,:,:,1) +IntMatrix(1,col)*FEp_sdc(:,:,:,col)
     Jk_sdc(:,:,:,1) = Jk_sdc(:,:,:,1) +IntMatrix(1,col)*FEj_sdc(:,:,:,col)
    end do
    call compute_time_derivatives(FEu_sdc(:,:,:,1),FEp_sdc(:,:,:,1),FEj_sdc(:,:,:,1),&
      Uk_sdc(:,:,:,1),Pk_sdc(:,:,:,1),Jk_sdc(:,:,:,1))
    DFu_sdc = FEu_sdc(:,:,:,1) - DFu_sdc
    DFp_sdc = FEp_sdc(:,:,:,1) - DFp_sdc
    DFj_sdc = FEj_sdc(:,:,:,1) - DFj_sdc
    
    do m=1,nsdc-1
     Uk_sdc(:,:,:,m+1) = IntMatrix(m+1,1)*FEu_sdc(:,:,:,1)
     pk_sdc(:,:,:,m+1) = IntMatrix(m+1,1)*FEp_sdc(:,:,:,1)
     jk_sdc(:,:,:,m+1) = IntMatrix(m+1,1)*FEj_sdc(:,:,:,1)
     
     do col = 2,nsdc
      Uk_sdc(:,:,:,m+1) = Uk_sdc(:,:,:,m+1)+IntMatrix(m+1,col)*FEu_sdc(:,:,:,col)
      pk_sdc(:,:,:,m+1) = pk_sdc(:,:,:,m+1)+IntMatrix(m+1,col)*FEp_sdc(:,:,:,col)
      jk_sdc(:,:,:,m+1) = jk_sdc(:,:,:,m+1)+IntMatrix(m+1,col)*FEj_sdc(:,:,:,col)
     end do
     Uk_sdc(:,:,:,m+1) = Uk_sdc(:,:,:,m) + dtm(m+1)*DFu_sdc+Uk_sdc(:,:,:,m+1)
     pk_sdc(:,:,:,m+1) = pk_sdc(:,:,:,m) + dtm(m+1)*DFp_sdc+pk_sdc(:,:,:,m+1)
     jk_sdc(:,:,:,m+1) = jk_sdc(:,:,:,m) + dtm(m+1)*DFj_sdc+jk_sdc(:,:,:,m+1)
     DFu_sdc = FEu_sdc(:,:,:,m+1)
     DFp_sdc = FEp_sdc(:,:,:,m+1)
     DFj_sdc = FEj_sdc(:,:,:,m+1)
     call compute_time_derivatives(FEu_sdc(:,:,:,m+1),FEp_sdc(:,:,:,m+1),&
       FEj_sdc(:,:,:,m+1),Uk_sdc(:,:,:,m+1),Pk_sdc(:,:,:,m+1),Jk_sdc(:,:,:,m+1))
     DFu_sdc = FEu_sdc(:,:,:,m+1) - DFu_sdc
     DFp_sdc = FEp_sdc(:,:,:,m+1) - DFp_sdc
     DFj_sdc = FEj_sdc(:,:,:,m+1) - DFj_sdc
    end do
   end do

   ! The update.
   do row = 1,nsdc+1
    u  = u  + IntMatrix(row,1)*FEu_sdc(:,:,:,1)
    p  = p  + IntMatrix(row,1)*FEp_sdc(:,:,:,1)
    jf = jf + IntMatrix(row,1)*FEj_sdc(:,:,:,1)
    do col = 2,nsdc
     u  = u  + IntMatrix(row,col)*FEu_sdc(:,:,:,col)
     p  = p  + IntMatrix(row,col)*FEp_sdc(:,:,:,col)
     jf = jf + IntMatrix(row,col)*FEj_sdc(:,:,:,col)
    end do
   end do
   
   ! Compute energy
   energy = 0.d0
   w = u(:,:,1)**2
   energy = energy + mu0*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
   do ivar = 2,nvar
    w = u(:,:,ivar)**2
    energy = energy + ep0*epinf*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
   end do
   do ivar = 1,nvar_j
    w = jf(:,:,ivar)**2
    energy = energy + (ep0/omgega_p_square)*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
   end do
   do ivar = 1,nvar_p
    w = p(:,:,ivar)**2
    energy = energy + &
      (ep0*omgega_0_square/omgega_p_square)*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
   end do
   if (a2 .gt. 0.d0 ) then
    w = (p(:,:,1)**2+p(:,:,2)**2)**2
    energy = energy + &
      (ep0*omgega_0_square/omgega_p_square)*a2*0.5d0*compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
   end if
   if (a4 .gt. 0.d0 ) then
    w = (p(:,:,1)**2+p(:,:,2)**2)**3
    energy = energy + &
      (ep0*omgega_0_square/omgega_p_square)*(a4/3.d0)&
      *compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
   end if
!   write(*,*) energy

   ivar = 1
   do j = 1,nyl
    do i = 1,nxl
     w(i,j) = u_mms(x(i),y(j),t+dt,ivar,0) 
    end do
   end do
   
!   write(21,fmt='(4(E24.16))') t+dt,maxval(abs(u(:,:,ivar)-w)),&
!     sqrt(hx*hy*sum((u(2:nxl,2:nyl,ivar)-w(2:nxl,2:nyl))**2)),energy 

   write(21,fmt='(5(E24.16))') t+dt,u((nxl-1)/4+1,(nxl-1)/4+1,1),p((nxl-1)/4+1,(nxl-1)/4+1,1),jf((nxl-1)/4+1,(nxl-1)/4+1,1),energy 
   
   if ((mod(it,n_plot) .eq.0) .or. (it.eq.nt)) then
    WRITE(charit,"(I7.7)") it
    call printdble2d(u(:,:,1),1,nxl,1,nyl,"Hz_" // charit // ".txt")
    call printdble2d(u(:,:,2),1,nxl,1,nyl,"Ex_" // charit // ".txt")
    call printdble2d(u(:,:,2),1,nxl,1,nyl,"Ey_" // charit // ".txt")
   end if
  end do

  call cpu_time(time2)
  write(*,*) (time2-time1),nt,(nxl+1)*(nyl+1)*3,(time2-time1)/dble(nt*(nxl+1)*(nyl+1)*3)
  close(21)
  deallocate(FEu_sdc,FEp_sdc,FEj_sdc,UK_sdc,PK_sdc,JK_sdc)
  deallocate(u_sdc,du_sdc,dfu_sdc)
  deallocate(p_sdc,dp_sdc,dfp_sdc)
  deallocate(j_sdc,dj_sdc,dfj_sdc)


end program main


subroutine compute_time_derivatives(ut,pt,jft,u,p,jf)
  use prec_def
  use stuff_module
  implicit none
  real(dp) :: u(nxl,nyl,nvar)
  real(dp) :: ut(nxl,nyl,nvar)
  real(dp) :: p(nxl,nyl,nvar_p)
  real(dp) :: pt(nxl,nyl,nvar_p)
  real(dp) :: jf(nxl,nyl,nvar_j)
  real(dp) :: jft(nxl,nyl,nvar_j)
  real(dp) :: ux(nxl,nyl,nvar),uy(nxl,nyl,nvar),w(nxl,nyl)
  real(dp) :: ux_tmp,uy_tmp
  integer :: ivar,i,j,k
  ! compute x-derivatives of H_z and E_y (i.e. variable 1 and 3)
  do ivar  = 1,nvar,2
   do j = 1,nyl
    ! Left part
    ux(1:6,j,ivar) = matmul(D1x_1,u(1:9,j,ivar))
    ! Middle part
    do i = 7,nxl-6
     ux_tmp = 0.d0
     do k = -3,3
      ux_tmp = ux_tmp + D1x_mid(k)*u(i+k,j,ivar)
     end do
     ux(i,j,ivar) = ux_tmp
    end do
    ! Right part
    ux(nxl-5:nxl,j,ivar) = matmul(D1x_N,u(nxl-8:nxl,j,ivar))
   end do
  end do
  ! compute y-derivatives of H_z and E_x (i.e. variable 1 and 2)
  do ivar = 1,2,1
   do i = 1,nxl
    uy(i,1:6,ivar) = matmul(D1y_1,u(i,1:9,ivar))
    do j = 7,nyl-6
     uy_tmp = 0.d0
     do k = -3,3
      uy_tmp = uy_tmp + D1y_mid(k)*u(i,j+k,ivar)
     end do
     uy(i,j,ivar) = uy_tmp
    end do
    uy(i,nyl-5:nyl,ivar) = matmul(D1y_N,u(i,nyl-8:nyl,ivar))
   end do
  end do
  ! Compute time derivatives
  ut(:,:,1) = -ux(:,:,3) + uy(:,:,2)
  ut(:,:,2) =  uy(:,:,1) - ep0*jf(:,:,1)
  ut(:,:,3) = -ux(:,:,1) - ep0*jf(:,:,2)

  pt = jf
  w = p(:,:,1)**2 + p(:,:,2)**2
  w = (1.d0 + a2*w + a4*w**2)
  jft(:,:,1) = -tau_inv*jf(:,:,1) - omgega_0_square*p(:,:,1)*w &
    + omgega_p_square*u(:,:,2)
  jft(:,:,2) = -tau_inv*jf(:,:,2) - omgega_0_square*p(:,:,2)*w &
    + omgega_p_square*u(:,:,3)

  ! Impose SAT bc (weak / penalty)
  do j = 1,nyl
   ut(1,j,1)   = ut(1,j,1)   - hix_1(1)*tau*(u(1,j,3)   - 0.d0)
   ut(nxl,j,1) = ut(nxl,j,1) + hix_n(6)*tau*(u(nxl,j,3) - 0.d0)
  end do
  do i = 1,nxl
   ut(i,1,1)   = ut(i,1,1)   + hiy_1(1)*tau*(u(i,1,2)   - 0.d0)
   ut(i,nyl,1) = ut(i,nyl,1) - hiy_n(6)*tau*(u(i,nyl,2) - 0.d0)
  end do
  ut(:,:,1)   = ut(:,:,1)*(1.d0/mu0)
  ut(:,:,2:3) = ut(:,:,2:3)*(1.d0/ep0/epinf)

end subroutine compute_time_derivatives


subroutine printdble1d(u,nx,str)
  use prec_def
  implicit none
  integer, intent(in) :: nx
  real(dp), intent(in) :: u(nx)
  character(len=*), intent(in) :: str
  integer :: i
  open(2,file=trim(str),status='unknown')
  do i=1,nx,1
   write(2,fmt='(E24.16)',advance='no') u(i)
  end do
  close(2)
end subroutine printdble1d

subroutine printdble2d(u,nx1,nx2,ny1,ny2,str)
  use prec_def
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  REAL(DP), intent(in) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=ny1,ny2,1
   do i=nx1,nx2,1
    if(abs(u(i,j)) .lt. 1e-40) then
     write(2,fmt='(E24.16)',advance='no') 0.d0
    else
     write(2,fmt='(E24.16)',advance='no') u(i,j)
    end if
   end do
   write(2,'()')
  end do
  close(2)
end subroutine printdble2d

subroutine readdble2d(u,nx1,nx2,ny1,ny2,str)
  use prec_def
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  REAL(DP), intent(inout) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: j
  open(2,file=trim(str),status='old',action='read')
  do j=ny1,ny2,1
     read(2,*) u(:,j)
  end do
  close(2)
end subroutine readdble2d



real(kind = dp)  function compute_l2_sq(w,nxl,nyl,hx,hy,hx_1,hy_1,hx_n,hy_n)
  use prec_def
  implicit none
  integer  :: j,nxl,nyl
  real(dp) :: w(nxl,nyl),hx,hy
  real(dp) :: Hx_1(6)
  real(dp) :: Hx_N(6)
  real(dp) :: Hy_1(6)
  real(dp) :: Hy_N(6)
  real(dp) :: energy_integrand(nyl)

  energy_integrand = 0.d0
  do j = 1,nyl
   ! Left part
   energy_integrand(j) = sum(hx_1*w(1:6,j)) &
     + hx*sum(w(7:nxl-6,j)) + sum(hx_n*w(nxl-5:nxl,j))
  end do
  compute_l2_sq = sum(hy_1*energy_integrand(1:6)) &
    + hy*sum(energy_integrand(7:nyl-6)) &
    + sum(hy_n*energy_integrand(nyl-5:nyl))

end function compute_l2_sq

subroutine kalle_anka(cse,dtm,intmatrix)
  use prec_def
  IMPLICIT NONE
  integer :: cse
  REAL(DP) :: A(6,6), B(6)
  real(dp) :: IntMatrix(cse+1,cse),dtm(cse+1)

  select case (cse)
  case (1)
   A(1,1) = 5.000000000000000d-01
   A(2,1) = 5.000000000000000d-01
   B(1:2) =(/ 5.000000000000000d-01,  5.000000000000000d-01/)
  case (2)
   A(1,1:2) = (/ 2.499999999999999d-01, -3.867513459481287d-02/)
   A(2,1:2) = (/ 2.886751345948130d-01,  2.886751345948130d-01/)
   A(3,1:2) = (/-3.867513459481287d-02,  2.499999999999999d-01/)
   B(1:3)   = (/ 2.113248654051870d-01,  5.773502691896260d-01, 2.113248654051870d-01/)
  case (3)
   A(1,1:3) = (/ 1.388888888888884d-01, -3.597666752493867d-02, 9.789444015308249d-03/)
   A(2,1:3) = (/ 1.613743060919758d-01,  2.581988897471613d-01,-3.227486121839515d-02/)
   A(3,1:3) = (/-3.227486121839518d-02,  2.581988897471614d-01, 1.613743060919758d-01/)
   A(4,1:3) = (/ 9.789444015308235d-03, -3.597666752493861d-02, 1.388888888888883d-01/)
   B(1:4)   = (/ 1.127016653792580d-01,  3.872983346207419d-01, 3.872983346207420d-01, &
     1.127016653792580d-01/)
  case (4)
   A(1,1:4) = (/ 8.696371128436355d-02, -2.660418008499882d-02, 1.262746268940475d-02,&
     -3.555149685795708d-03/)
   A(2,1:4) = (/ 1.011544062155046d-01,  1.896404688006354d-01,-4.050789129187568d-02, &
     1.029065028033388d-02/)
   A(3,1:4) = (/-2.092619552567926d-02,  1.909167173181072d-01, 1.909167173181073d-01,&
     -2.092619552567929d-02/)
   A(4,1:4) = (/ 1.029065028033382d-02, -4.050789129187557d-02, 1.896404688006353d-01,&
     1.011544062155046d-01/)
   A(5,1:4) = (/-3.555149685795718d-03,  1.262746268940488d-02,-2.660418008499903d-02,&
     8.696371128436392d-02/)
   B(1:5)   = (/ 6.943184420297377d-02,  2.605776340045982d-01, 3.399810435848560d-01, &
     2.605776340045981d-01, 6.943184420297405d-02/)
  case (5)
   A(1,1:5) = (/ 5.923172126404737d-02, -1.957036435907604d-02, 1.125440081864303d-02,&
     -5.593793660812233d-03, 1.588112967866009d-03/)
   A(2,1:5) = (/ 6.891928440599805d-02,  1.392275319839184d-01,-3.584651543828557d-02,&
     1.591207433149566d-02,-4.357107366635624d-03/)
   A(3,1:5) = (/-1.437471766582066d-02,  1.403474840557995d-01, 1.668143368418641d-01,&
     -3.100859710164155d-02, 7.456148922639473d-03/)
   A(4,1:5) = (/ 7.456148922639435d-03, -3.100859710164146d-02, 1.668143368418642d-01,&
     1.403474840557994d-01,-1.437471766582052d-02/)
   A(5,1:5) = (/-4.357107366635649d-03,  1.591207433149570d-02,-3.584651543828551d-02,&
     1.392275319839183d-01, 6.891928440599815d-02/)
   A(6,1:5) = (/ 1.588112967866014d-03, -5.593793660812249d-03, 1.125440081864305d-02,&
     -1.957036435907612d-02, 5.923172126404717d-02/)
   B(1:6)   = (/ 4.691007703066813d-02,  1.838552679164908d-01, 2.692346550528408d-01,&
     2.692346550528411d-01, 1.838552679164910d-01, 4.691007703066788d-02/)
  end select

  IntMatrix(1:cse+1,1:cse) = A(1:cse+1,1:cse)
  dtm(1:cse+1) = B(1:cse+1)

  
end subroutine kalle_anka
