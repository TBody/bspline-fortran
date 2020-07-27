!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!
!  Multidimensional (1D-6D) B-spline interpolation of data on a regular grid.
!  Basic pure subroutine interface.
!
!### Notes
!
!  This module is based on the B-spline and spline routines from [1].
!  The original Fortran 77 routines were converted to free-form source.
!  Some of them are relatively unchanged from the originals, but some have
!  been extensively refactored. In addition, new routines for
!  1d, 4d, 5d, and 6d interpolation were also created (these are simply
!  extensions of the same algorithm into higher dimensions).
!
!### See also
!  * An object-oriented interface can be found in [[bspline_oo_module]].
!
!### References
!
!  1. DBSPLIN and DTENSBS from the
!     [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).
!     Original code is public domain.
!  2. Carl de Boor, "A Practical Guide to Splines",
!     Springer-Verlag, New York, 1978.
!  3. Carl de Boor, [Efficient Computer Manipulation of Tensor
!     Products](http://dl.acm.org/citation.cfm?id=355831),
!     ACM Transactions on Mathematical Software,
!     Vol. 5 (1979), p. 173-182.
!  4. D.E. Amos, "Computation with Splines and B-Splines",
!     SAND78-1968, Sandia Laboratories, March, 1979.
!  5. Carl de Boor,
!     [Package for calculating with B-splines](http://epubs.siam.org/doi/abs/10.1137/0714026),
!     SIAM Journal on Numerical Analysis 14, 3 (June 1977), p. 441-472.
!  6. D.E. Amos, "Quadrature subroutines for splines and B-splines",
!     Report SAND79-1825, Sandia Laboratories, December 1979.

    module bspline_sub_module

    use bspline_kinds_module, only: wp, ip
    use,intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    !Spline function order (order = polynomial degree + 1)
    integer(ip),parameter,public :: bspline_order_quadratic = 3_ip
    integer(ip),parameter,public :: bspline_order_cubic     = 4_ip
    integer(ip),parameter,public :: bspline_order_quartic   = 5_ip
    integer(ip),parameter,public :: bspline_order_quintic   = 6_ip

    !main routines:
    public :: db2ink, db2val

    public :: get_status_message

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
!  the two-dimensional gridded data
!  $$ [x(i),y(j),\mathrm{fcn}(i,j)] ~\mathrm{for}~ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y $$
!  The interpolating function and its derivatives may
!  subsequently be evaluated by the function [[db2val]].
!
!  The interpolating function is a piecewise polynomial function
!  represented as a tensor product of one-dimensional b-splines. the
!  form of this function is
!
!  $$ s(x,y) = \sum_{i=1}^{n_x} \sum_{j=1}^{n_y} a_{ij} u_i(x) v_j(y) $$
!
!  where the functions \(u_i\) and \(v_j\) are one-dimensional b-spline
!  basis functions. the coefficients \( a_{ij} \) are chosen so that
!
!  $$ s(x(i),y(j)) = \mathrm{fcn}(i,j) ~\mathrm{for}~ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y $$
!
!  Note that for each fixed value of \(y\), \( s(x,y) \) is a piecewise
!  polynomial function of \(x\) alone, and for each fixed value of \(x\), \( s(x,y) \)
!  is a piecewise polynomial function of \(y\) alone. in one dimension
!  a piecewise polynomial may be created by partitioning a given
!  interval into subintervals and defining a distinct polynomial piece
!  on each one. the points where adjacent subintervals meet are called
!  knots. each of the functions \(u_i\) and \(v_j\) above is a piecewise
!  polynomial.
!
!  Users of [[db2ink]] choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the \(x\) and
!  \(y\) directions (`kx` and `ky`). users also may define their own knot
!  sequence in \(x\) and \(y\) separately (`tx` and `ty`). if `iflag=0`, however,
!  [[db2ink]] will choose sequences of knots that result in a piecewise
!  polynomial interpolant with `kx-2` continuous partial derivatives in
!  \(x\) and `ky-2` continuous partial derivatives in \(y\). (`kx` knots are taken
!  near each endpoint in the \(x\) direction, not-a-knot end conditions
!  are used, and the remaining knots are placed at data points if `kx`
!  is even or at midpoints between data points if `kx` is odd. the \(y\)
!  direction is treated similarly.)
!
!  After a call to [[db2ink]], all information necessary to define the
!  interpolating function are contained in the parameters `nx`, `ny`, `kx`,
!  `ky`, `tx`, `ty`, and `bcoef`. These quantities should not be altered until
!  after the last call of the evaluation routine [[db2val]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db2ink(x,nx,y,ny,fcn,kx,ky,iknot,tx,ty,bcoef,iflag)

    implicit none

    integer(ip),intent(in)                  :: nx     !! Number of \(x\) abcissae
    integer(ip),intent(in)                  :: ny     !! Number of \(y\) abcissae
    integer(ip),intent(in)                  :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                  :: ky     !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)        :: y      !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in)      :: fcn    !! `(nx,ny)` matrix of function values to interpolate.
                                                      !! `fcn(i,j)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)                  :: iknot  !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db1ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)     :: tx     !! The `(nx+kx)` knots in the \(x\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)     :: ty     !! The `(ny+ky)` knots in the \(y\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:,:),intent(out)     :: bcoef  !! `(nx,ny)` matrix of coefficients of the b-spline interpolant.
    integer(ip),intent(out)                 :: iflag  !! *  0 = successful execution.
                                                      !! *  2 = `iknot` out of range.
                                                      !! *  3 = `nx` out of range.
                                                      !! *  4 = `kx` out of range.
                                                      !! *  5 = `x` not strictly increasing.
                                                      !! *  6 = `tx` not non-decreasing.
                                                      !! *  7 = `ny` out of range.
                                                      !! *  8 = `ky` out of range.
                                                      !! *  9 = `y` not strictly increasing.
                                                      !! * 10 = `ty` not non-decreasing.
                                                      !! * 700 = `size(x)`  \( \ne \) `size(fcn,1)`
                                                      !! * 701 = `size(y)`  \( \ne \) `size(fcn,2)`
                                                      !! * 706 = `size(x)`  \( \ne \) `nx`
                                                      !! * 707 = `size(y)`  \( \ne \) `ny`
                                                      !! * 712 = `size(tx)` \( \ne \) `nx+kx`
                                                      !! * 713 = `size(ty)` \( \ne \) `ny+ky`
                                                      !! * 800 = `size(x)`  \( \ne \) `size(bcoef,1)`
                                                      !! * 801 = `size(y)`  \( \ne \) `size(bcoef,2)`

    logical :: status_ok
    real(wp),dimension(:),allocatable :: temp !! work array of length `nx*ny`
    real(wp),dimension(:),allocatable :: work !! work array of length `max(2*kx*(nx+1),2*ky*(ny+1))`

    !check validity of inputs

    call check_inputs(  iknot,&
                        iflag,&
                        nx=nx,ny=ny,&
                        kx=kx,ky=ky,&
                        x=x,y=y,&
                        tx=tx,ty=ty,&
                        f2=fcn,&
                        bcoef2=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        !choose knots
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
        end if

        allocate(temp(nx*ny))
        allocate(work(max(2_ip*kx*(nx+1_ip),2_ip*ky*(ny+1_ip))))

        !construct b-spline coefficients
                         call dbtpcf(x,nx,fcn, nx,ny,tx,kx,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(y,ny,temp,ny,nx,ty,ky,bcoef,work,iflag)

        deallocate(temp)
        deallocate(work)

    end if

    end subroutine db2ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db2ink]] or one of its
!  derivatives at the point (`xval`,`yval`).
!
!  To evaluate the interpolant
!  itself, set `idx=idy=0`, to evaluate the first partial with respect
!  to `x`, set `idx=1,idy=0`, and so on.
!
!  [[db2val]] returns 0.0 if `(xval,yval)` is out of range. that is, if
!```fortran
!   xval < tx(1) .or. xval > tx(nx+kx) .or.
!   yval < ty(1) .or. yval > ty(ny+ky)
!```
!  if the knots tx and ty were chosen by [[db2ink]], then this is equivalent to:
!```fortran
!   xval < x(1) .or. xval > x(nx)+epsx .or.
!   yval < y(1) .or. yval > y(ny)+epsy
!```
!  where
!```fortran
!   epsx = 0.1*(x(nx)-x(nx-1))
!   epsy = 0.1*(y(ny)-y(ny-1))
!```
!
!  The input quantities `tx`, `ty`, `nx`, `ny`, `kx`, `ky`, and `bcoef` should be
!  unchanged since the last call of [[db2ink]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag,inbvx,inbvy,iloy,w1,w0,extrap)

    implicit none

    integer(ip),intent(in)               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer(ip),intent(in)               :: ny       !! the number of interpolation points in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    integer(ip),intent(in)               :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer(ip),intent(in)               :: ky       !! order of polynomial pieces in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    real(wp),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                  :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction.
                                                     !! (same as in last call to [[db2ink]])
    real(wp),dimension(ny+ky),intent(in) :: ty       !! sequence of knots defining the piecewise
                                                     !! polynomial in the \(y\) direction.
                                                     !! (same as in last call to [[db2ink]])
    real(wp),dimension(nx,ny),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db2ink]].
    real(wp),intent(out)                 :: f        !! interpolated value
    integer(ip),intent(out)              :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer(ip),intent(inout)            :: inbvx    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer(ip),intent(inout)            :: inbvy    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer(ip),intent(inout)            :: iloy     !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    real(wp),dimension(ky),intent(inout)              :: w1 !! work array
    real(wp),dimension(3_ip*max(kx,ky)),intent(inout) :: w0 !! work array
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    integer(ip) :: k, lefty, kcol

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(yval,ty,2_ip,extrap); if (iflag/=0_ip) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,iflag,extrap); if (iflag/=0_ip) return

    kcol = lefty - ky
    do k=1_ip,ky
        kcol = kcol + 1_ip
        call dbvalu(tx,bcoef(:,kcol),nx,kx,idx,xval,inbvx,w0,iflag,w1(k),extrap)
        if (iflag/=0_ip) return !error
    end do

    kcol = lefty - ky + 1_ip
    call dbvalu(ty(kcol:),w1,ky,ky,idy,yval,inbvy,w0,iflag,f,extrap)

    end subroutine db2val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Checks if the value is withing the range of the knot vectors.
!  This is called by the various `db*val` routines.

    pure function check_value(x,t,i,extrap) result(iflag)

    implicit none

    integer(ip)                      :: iflag   !! returns 0 if value is OK, otherwise returns `600+i`
    real(wp),intent(in)              :: x       !! the value to check
    integer(ip),intent(in)           :: i       !! 1=x, 2=y, 3=z, 4=q, 5=r, 6=s
    real(wp),dimension(:),intent(in) :: t       !! the knot vector
    logical,intent(in),optional      :: extrap  !! if extrapolation is allowed
                                                !! (if not present, default is False)

    logical :: allow_extrapolation  !! if extrapolation is allowed

    if (present(extrap)) then
        allow_extrapolation = extrap
    else
        allow_extrapolation = .false.
    end if

    if (allow_extrapolation) then
        ! in this case all values are OK
        iflag = 0_ip
    else
        if (x<t(1_ip) .or. x>t(size(t,kind=ip))) then
            iflag = 600_ip + i  ! value out of bounds (601, 602, etc.)
        else
            iflag = 0_ip
        end if
    end if

    end function check_value
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the validity of the inputs to the `db*ink` routines.
!  Prints warning message if there is an error,
!  and also sets iflag and status_ok.
!
!  Supports up to 6D: `x`,`y`,`z`,`q`,`r`,`s`
!
!### Notes
!
!  The code is new, but the logic is based on the original
!  logic in the CMLIB routines `db2ink` and `db3ink`.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine check_inputs(iknot,&
                                 iflag,&
                                 nx,ny,nz,nq,nr,ns,&
                                 kx,ky,kz,kq,kr,ks,&
                                 x,y,z,q,r,s,&
                                 tx,ty,tz,tq,tr,ts,&
                                 f1,f2,f3,f4,f5,f6,&
                                 bcoef1,bcoef2,bcoef3,bcoef4,bcoef5,bcoef6,&
                                 status_ok)

    implicit none

    integer(ip),intent(in)                              :: iknot !! = 0 if the `INK` routine is computing the knots.
    integer(ip),intent(out)                             :: iflag
    integer(ip),intent(in),optional                     :: nx,ny,nz,nq,nr,ns
    integer(ip),intent(in),optional                     :: kx,ky,kz,kq,kr,ks
    real(wp),dimension(:),intent(in),optional           :: x,y,z,q,r,s
    real(wp),dimension(:),intent(in),optional           :: tx,ty,tz,tq,tr,ts
    real(wp),dimension(:),intent(in),optional           :: f1,bcoef1
    real(wp),dimension(:,:),intent(in),optional         :: f2,bcoef2
    real(wp),dimension(:,:,:),intent(in),optional       :: f3,bcoef3
    real(wp),dimension(:,:,:,:),intent(in),optional     :: f4,bcoef4
    real(wp),dimension(:,:,:,:,:),intent(in),optional   :: f5,bcoef5
    real(wp),dimension(:,:,:,:,:,:),intent(in),optional :: f6,bcoef6
    logical,intent(out)                                 :: status_ok

    logical :: error

    status_ok = .false.

    if ((iknot < 0_ip) .or. (iknot > 1_ip)) then

        iflag = 2_ip ! iknot is out of range

    else

        call check('x',nx,kx,x,tx,[3_ip,  4_ip, 5_ip, 6_ip,706_ip,712_ip],iflag,error); if (error) return
        call check('y',ny,ky,y,ty,[7_ip,  8_ip, 9_ip,10_ip,707_ip,713_ip],iflag,error); if (error) return
        call check('z',nz,kz,z,tz,[11_ip,12_ip,13_ip,14_ip,708_ip,714_ip],iflag,error); if (error) return
        call check('q',nq,kq,q,tq,[15_ip,16_ip,17_ip,18_ip,709_ip,715_ip],iflag,error); if (error) return
        call check('r',nr,kr,r,tr,[19_ip,20_ip,21_ip,22_ip,710_ip,716_ip],iflag,error); if (error) return
        call check('s',ns,ks,s,ts,[23_ip,24_ip,25_ip,26_ip,711_ip,717_ip],iflag,error); if (error) return

        if (present(x) .and. present(f1) .and. present(bcoef1)) then
            if (size(x,kind=ip)/=size(f1,1_ip,kind=ip))     then; iflag = 700_ip; return; end if
            if (size(x,kind=ip)/=size(bcoef1,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(f2) .and. present(bcoef2)) then
            if (size(x,kind=ip)/=size(f2,1_ip,kind=ip))     then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f2,2_ip,kind=ip))     then; iflag = 701_ip; return; end if
            if (size(x,kind=ip)/=size(bcoef2,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)/=size(bcoef2,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(f3) .and. &
            present(bcoef3)) then
            if (size(x,kind=ip)/=size(f3,1_ip,kind=ip))     then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f3,2_ip,kind=ip))     then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f3,3_ip,kind=ip))     then; iflag = 702_ip; return; end if
            if (size(x,kind=ip)/=size(bcoef3,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)/=size(bcoef3,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)/=size(bcoef3,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(f4) .and. present(bcoef4)) then
            if (size(x,kind=ip)/=size(f4,1_ip,kind=ip))     then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f4,2_ip,kind=ip))     then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f4,3_ip,kind=ip))     then; iflag = 702_ip; return; end if
            if (size(q,kind=ip)/=size(f4,4_ip,kind=ip))     then; iflag = 703_ip; return; end if
            if (size(x,kind=ip)/=size(bcoef4,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)/=size(bcoef4,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)/=size(bcoef4,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
            if (size(q,kind=ip)/=size(bcoef4,4_ip,kind=ip)) then; iflag = 803_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(f5) .and. present(bcoef5)) then
            if (size(x,kind=ip)/=size(f5,1_ip,kind=ip))     then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f5,2_ip,kind=ip))     then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f5,3_ip,kind=ip))     then; iflag = 702_ip; return; end if
            if (size(q,kind=ip)/=size(f5,4_ip,kind=ip))     then; iflag = 703_ip; return; end if
            if (size(r,kind=ip)/=size(f5,5_ip,kind=ip))     then; iflag = 704_ip; return; end if
            if (size(x,kind=ip)/=size(bcoef5,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)/=size(bcoef5,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)/=size(bcoef5,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
            if (size(q,kind=ip)/=size(bcoef5,4_ip,kind=ip)) then; iflag = 803_ip; return; end if
            if (size(r,kind=ip)/=size(bcoef5,5_ip,kind=ip)) then; iflag = 804_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(s) .and. present(f6) .and. present(bcoef6)) then
            if (size(x,kind=ip)/=size(f6,1_ip,kind=ip))     then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f6,2_ip,kind=ip))     then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f6,3_ip,kind=ip))     then; iflag = 702_ip; return; end if
            if (size(q,kind=ip)/=size(f6,4_ip,kind=ip))     then; iflag = 703_ip; return; end if
            if (size(r,kind=ip)/=size(f6,5_ip,kind=ip))     then; iflag = 704_ip; return; end if
            if (size(s,kind=ip)/=size(f6,6_ip,kind=ip))     then; iflag = 705_ip; return; end if
            if (size(x,kind=ip)/=size(bcoef6,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)/=size(bcoef6,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)/=size(bcoef6,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
            if (size(q,kind=ip)/=size(bcoef6,4_ip,kind=ip)) then; iflag = 803_ip; return; end if
            if (size(r,kind=ip)/=size(bcoef6,5_ip,kind=ip)) then; iflag = 804_ip; return; end if
            if (size(s,kind=ip)/=size(bcoef6,6_ip,kind=ip)) then; iflag = 805_ip; return; end if
        end if

        status_ok = .true.
        iflag = 0_ip

    end if

    contains

        pure subroutine check(s,n,k,x,t,ierrs,iflag,error)  !! check `t`,`x`,`n`,`k` for validity

        implicit none

        character(len=1),intent(in)               :: s     !! coordinate string: 'x','y','z','q','r','s'
        integer(ip),intent(in),optional           :: n     !! size of `x`
        integer(ip),intent(in),optional           :: k     !! order
        real(wp),dimension(:),intent(in),optional :: x     !! abcissae vector
        real(wp),dimension(:),intent(in),optional :: t     !! knot vector `size(n+k)`
        integer(ip),dimension(:),intent(in)       :: ierrs !! int error codes for `n`,`k`,`x`,`t`,
                                                           !! `size(x)`,`size(t)` checks
        integer(ip),intent(out)                   :: iflag !! status return code
        logical,intent(out)                       :: error !! true if there was an error

        integer(ip),dimension(2) :: itmp !! temp integer array

        if (present(n) .and. present(k) .and. present(x) .and. present(t)) then
            itmp = [ierrs(1_ip),ierrs(5)]
            call check_n('n'//s,n,x,itmp,iflag,error);     if (error) return
            call check_k('k'//s,k,n,ierrs(2),iflag,error); if (error) return
            call check_x(s,n,x,ierrs(3),iflag,error);      if (error) return
            if (iknot /= 0_ip) then
                itmp = [ierrs(4),ierrs(6)]
                call check_t('t'//s,n,k,t,itmp,iflag,error); if (error) return
            end if
        end if

        end subroutine check

        pure subroutine check_n(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)         :: s
        integer(ip),intent(in)              :: n
        real(wp),dimension(:),intent(in)    :: x     !! abcissae vector
        integer(ip),dimension(2),intent(in) :: ierr  !! [n<3 check, size(x)==n check]
        integer(ip),intent(out)             :: iflag !! status return code
        logical,intent(out)                 :: error

        if (n < 3_ip) then
            iflag = ierr(1_ip)
            error = .true.
        else
            if (size(x)/=n) then
                iflag = ierr(2)
                error = .true.
            else
                error = .false.
            end if
        end if

        end subroutine check_n

        pure subroutine check_k(s,k,n,ierr,iflag,error)

        implicit none

        character(len=*),intent(in) :: s
        integer(ip),intent(in)      :: k
        integer(ip),intent(in)      :: n
        integer(ip),intent(in)      :: ierr
        integer(ip),intent(out)     :: iflag !! status return code
        logical,intent(out)         :: error

        if ((k < 2_ip) .or. (k >= n)) then
            iflag = ierr
            error = .true.
        else
            error = .false.
        end if

        end subroutine check_k

        pure subroutine check_x(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer(ip),intent(in)            :: n
        real(wp),dimension(:),intent(in)  :: x
        integer(ip),intent(in)            :: ierr
        integer(ip),intent(out)           :: iflag !! status return code
        logical,intent(out)               :: error

        integer(ip) :: i

        error = .true.
        do i=2_ip,n
            if (x(i) <= x(i-1_ip)) then
                iflag = ierr
                return
            end if
        end do
        error = .false.

        end subroutine check_x

        pure subroutine check_t(s,n,k,t,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)         :: s
        integer(ip),intent(in)              :: n
        integer(ip),intent(in)              :: k
        real(wp),dimension(:),intent(in)    :: t
        integer(ip),dimension(2),intent(in) :: ierr  !! [non-decreasing check, size check]
        integer(ip),intent(out)             :: iflag !! status return code
        logical,intent(out)                 :: error

        integer(ip) :: i

        error = .true.

        if (size(t)/=(n+k)) then
            iflag = ierr(2)
            return
        end if

        do i=2_ip,n + k
            if (t(i) < t(i-1_ip))  then
                iflag = ierr(1_ip)
                return
            end if
        end do
        error = .false.

        end subroutine check_t

    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbknot chooses a knot sequence for interpolation of order k at the
!  data points x(i), i=1,..,n.  the n+k knots are placed in the array
!  t.  k knots are placed at each endpoint and not-a-knot end
!  conditions are used.  the remaining knots are placed at data points
!  if n is even and between data points if n is odd.  the rightmost
!  knot is shifted slightly to the right to insure proper interpolation
!  at x(n) (see page 350 of the reference).
!
!### History
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbknot(x,n,k,t)

    implicit none

    integer(ip),intent(in)             :: n  !! dimension of `x`
    integer(ip),intent(in)             :: k
    real(wp),dimension(:),intent(in)   :: x
    real(wp),dimension(:),intent(out)  :: t

    integer(ip) :: i, j, ipj, npj, ip1, jstrt
    real(wp) :: rnot

    !put k knots at each endpoint
    !(shift right endpoints slightly -- see pg 350 of reference)
    rnot = x(n) + 0.1_wp*( x(n)-x(n-1_ip) )
    do j=1_ip,k
        t(j)   = x(1_ip)
        npj    = n + j
        t(npj) = rnot
    end do

    !distribute remaining knots

    if (mod(k,2_ip) == 1_ip)  then

        !case of odd k --  knots between data points

        i = (k-1_ip)/2_ip - k
        ip1 = i + 1_ip
        jstrt = k + 1_ip
        do j=jstrt,n
            ipj = i + j
            t(j) = 0.5_wp*( x(ipj) + x(ipj+1_ip) )
        end do

    else

        !case of even k --  knots at data points

        i = (k/2_ip) - k
        jstrt = k+1_ip
        do j=jstrt,n
            ipj = i + j
            t(j) = x(ipj)
        end do

    end if

    end subroutine dbknot
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbtpcf computes b-spline interpolation coefficients for nf sets
!  of data stored in the columns of the array fcn. the b-spline
!  coefficients are stored in the rows of bcoef however.
!  each interpolation is based on the n abcissa stored in the
!  array x, and the n+k knots stored in the array t. the order
!  of each interpolation is k.
!
!### History
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag)

    integer(ip),intent(in)                :: n  !! dimension of `x`
    integer(ip),intent(in)                :: nf
    integer(ip),intent(in)                :: ldf
    integer(ip),intent(in)                :: k
    real(wp),dimension(:),intent(in)      :: x
    real(wp),dimension(ldf,nf),intent(in) :: fcn
    real(wp),dimension(:),intent(in)      :: t
    real(wp),dimension(nf,n),intent(out)  :: bcoef
    real(wp),dimension(*),intent(out)     :: work   !! work array of size >= `2*k*(n+1)`
    integer(ip),intent(out)               :: iflag  !!   0: no errors
                                                    !! 301: n should be >0

    integer(ip) :: i, j, m1, m2, iq, iw

    ! check for null input

    if (nf > 0_ip)  then

        ! partition work array
        m1 = k - 1_ip
        m2 = m1 + k
        iq = 1_ip + n
        iw = iq + m2*n+1_ip

        ! compute b-spline coefficients

        ! first data set

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw),iflag)
        if (iflag == 0_ip) then
            do i=1_ip,n
                bcoef(1_ip,i) = work(i)
            end do

            !  all remaining data sets by back-substitution

            if (nf == 1_ip)  return
            do j=2_ip,nf
                do i=1_ip,n
                    work(i) = fcn(i,j)
                end do
                call dbnslv(work(iq),m2,n,m1,m1,work)
                do i=1_ip,n
                    bcoef(j,i) = work(i)
                end do
            end do
        end if

    else
        !write(error_unit,'(A)') 'dbtpcf - n should be >0'
        iflag = 301_ip
    end if

    end subroutine dbtpcf
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbintk produces the b-spline coefficients, bcoef, of the
!  b-spline of order k with knots t(i), i=1,...,n+k, which
!  takes on the value y(i) at x(i), i=1,...,n.  the spline or
!  any of its derivatives can be evaluated by calls to [[dbvalu]].
!
!  the i-th equation of the linear system a*bcoef = b for the
!  coefficients of the interpolant enforces interpolation at
!  x(i), i=1,...,n.  hence, b(i) = y(i), for all i, and a is
!  a band matrix with 2k-1 bands if a is invertible.  the matrix
!  a is generated row by row and stored, diagonal by diagonal,
!  in the rows of q, with the main diagonal going into row k.
!  the banded system is then solved by a call to dbnfac (which
!  constructs the triangular factorization for a and stores it
!  again in q), followed by a call to dbnslv (which then
!  obtains the solution bcoef by substitution).  dbnfac does no
!  pivoting, since the total positivity of the matrix a makes
!  this unnecessary.  the linear system to be solved is
!  (theoretically) invertible if and only if
!          t(i) < x(i) < t(i+k),        for all i.
!  equality is permitted on the left for i=1 and on the right
!  for i=n when k knots are used at x(1) or x(n).  otherwise,
!  violation of this condition is certain to lead to an error.
!
!# Error conditions
!
!  * improper input
!  * singular system of equations
!
!### History
!  * splint written by carl de boor [5]
!  * dbintk author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations. (jec)
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbintk(x,y,t,n,k,bcoef,q,work,iflag)

    implicit none

    integer(ip),intent(in)            :: n      !! number of data points, n >= k
    real(wp),dimension(n),intent(in)  :: x      !! vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real(wp),dimension(n),intent(in)  :: y      !! corresponding vector of length n containing data
                                                !! point ordinates.
    real(wp),dimension(*),intent(in)  :: t      !! knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer(ip),intent(in)            :: k      !! order of the spline, k >= 1
    real(wp),dimension(n),intent(out) :: bcoef  !! a vector of length n containing the b-spline coefficients
    real(wp),dimension(*),intent(out) :: q      !! a work vector of length (2*k-1)*n, containing
                                                !! the triangular factorization of the coefficient
                                                !! matrix of the linear system being solved.  the
                                                !! coefficients for the interpolant of an
                                                !! additional data set (x(i),yy(i)), i=1,...,n
                                                !! with the same abscissa can be obtained by loading
                                                !! yy into bcoef and then executing
                                                !! call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
    real(wp),dimension(*),intent(out) :: work   !! work vector of length 2*k
    integer(ip),intent(out)           :: iflag  !! *   0: no errors.
                                                !! * 100: k does not satisfy k>=1.
                                                !! * 101: n does not satisfy n>=k.
                                                !! * 102: x(i) does not satisfy x(i)<x(i+1) for some i.
                                                !! * 103: some abscissa was not in the support of the.
                                                !! corresponding basis function and the system is singular.
                                                !! * 104: the system of solver detects a singular system.
                                                !! although the theoretical conditions for a solution were satisfied.

    integer(ip) :: iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real(wp) :: xi
    logical :: found

    if (k<1_ip) then
        !write(error_unit,'(A)') 'dbintk - k does not satisfy k>=1'
        iflag = 100_ip
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbintk - n does not satisfy n>=k'
        iflag = 101_ip
        return
    end if

    jj = n - 1_ip
    if (jj/=0_ip) then
        do i=1_ip,jj
            if (x(i)>=x(i+1_ip)) then
                !write(error_unit,'(A)') 'dbintk - x(i) does not satisfy x(i)<x(i+1) for some i'
                iflag = 102_ip
                return
            end if
        end do
    end if

    np1 = n + 1_ip
    km1 = k - 1_ip
    kpkm2 = 2_ip*km1
    left = k
    ! zero out all entries of q
    lenq = n*(k+km1)
    do i=1_ip,lenq
        q(i) = 0.0_wp
    end do

    ! loop over i to construct the n interpolation equations
    do i=1_ip,n

        xi = x(i)
        ilp1mx = min(i+k,np1)
        ! find left in the closed interval (i,i+k-1_ip) such that
        !         t(left) <= x(i) < t(left+1_ip)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
            !             ' corresponding basis function and the system is singular'
            iflag = 103_ip
            return
        end if
        found = .false.
        do
            found = (xi<t(left+1_ip))
            if (found) exit
            left = left + 1_ip
            if (left>=ilp1mx) exit
        end do
        if (.not. found) then
            left = left - 1_ip
            if (xi>t(left+1_ip)) then
                !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
                !             ' corresponding basis function and the system is singular'
                iflag = 103_ip
                return
            end if
        end if
        ! the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
        call dbspvn(t, k, k, 1_ip, xi, left, bcoef, work, iwork, iflag)
        if (iflag/=0_ip) return

        ! we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
        ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
        ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
        ! as a two-dim. array , with  2*k-1  rows (see comments in
        ! dbnfac). in the present program, we treat  q  as an equivalent
        ! one-dimensional array (because of fortran restrictions on
        ! dimension statements) . we therefore want  bcoef(j) to go into
        ! entry
        !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
        !            = i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
        ! of q.
        jj = i - left + 1_ip + (left-k)*(k+km1)
        do j=1_ip,k
            jj = jj + kpkm2
            q(jj) = bcoef(j)
        end do

    end do

    ! obtain factorization of a, stored again in q.
    call dbnfac(q, k+km1, n, km1, km1, iflag)

    if (iflag==1) then !success
        ! solve  a*bcoef = y  by backsubstitution
        do i=1_ip,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
        iflag = 0_ip
    else  !failure
        !write(error_unit,'(A)') 'dbintk - the system of solver detects a singular system'//&
        !             ' although the theoretical conditions for a solution were satisfied'
        iflag = 104_ip
    end if

    end subroutine dbintk
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns in w the LU-factorization (without pivoting) of the banded
!  matrix a of order nrow with (nbandl + 1 + nbandu) bands or diagonals
!  in the work array w .
!
!  gauss elimination without pivoting is used. the routine is
!  intended for use with matrices a which do not require row inter-
!  changes during factorization, especially for the totally
!  positive matrices which occur in spline calculations.
!  the routine should not be used for an arbitrary banded matrix.
!
!### Work array
!
! **Input**
!
!        w array of size (nroww,nrow) contains the interesting
!        part of a banded matrix  a , with the diagonals or bands of  a
!        stored in the rows of  w , while columns of  a  correspond to
!        columns of  w . this is the storage mode used in  linpack  and
!        results in efficient innermost loops.
!           explicitly,  a  has  nbandl  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   nbandu  bands above the diagonal
!        and thus, with    middle = nbandu + 1,
!          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
!                                              j=1,...,nrow .
!        for example, the interesting entries of a (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  w
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        all other entries of  w  not identified in this way with an en-
!        try of  a  are never referenced .
!
! **Output**
!
!  * if  iflag = 1, then
!        w contains the lu-factorization of  a  into a unit lower triangu-
!        lar matrix  l  and an upper triangular matrix  u (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  a . this makes it possible to solve any particular linear
!        system  a*x = b  for  x  by a
!              call dbnslv ( w, nroww, nrow, nbandl, nbandu, b )
!        with the solution x  contained in  b  on return .
!  * if  iflag = 2, then
!        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  a  does not have an lu-factorization. this implies that
!        a  is singular in case it is totally positive .
!
!### History
!  * banfac written by carl de boor [5]
!  * dbnfac from CMLIB [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)

    integer(ip),intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
    integer(ip),intent(in) :: nrow    !! matrix order
    integer(ip),intent(in) :: nbandl  !! number of bands of a below the main diagonal
    integer(ip),intent(in) :: nbandu  !! number of bands of a above the main diagonal
    integer(ip),intent(out) :: iflag  !! indicating success(=1) or failure (=2)
    real(wp),dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

    integer(ip) :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
    real(wp) :: factor, pivot

    iflag = 1_ip
    middle = nbandu + 1_ip   ! w(middle,.) contains the main diagonal of a.
    nrowm1 = nrow - 1_ip

    if (nrowm1 < 0_ip) then
        iflag = 2_ip
        return
    else if (nrowm1 == 0_ip) then
        if (w(middle,nrow)==0.0_wp) iflag = 2_ip
        return
    end if

    if (nbandl<=0_ip) then
        ! a is upper triangular. check that diagonal is nonzero .
        do i=1_ip,nrowm1
            if (w(middle,i)==0.0_wp) then
                iflag = 2_ip
                return
            end if
        end do
        if (w(middle,nrow)==0.0_wp) iflag = 2_ip
        return
    end if

    if (nbandu<=0_ip) then
        ! a is lower triangular. check that diagonal is nonzero and
        ! divide each column by its diagonal.
        do i=1_ip,nrowm1
            pivot = w(middle,i)
            if (pivot==0.0_wp) then
                iflag = 2_ip
                return
            end if
            jmax = min(nbandl,nrow-i)
            do j=1_ip,jmax
                w(middle+j,i) = w(middle+j,i)/pivot
            end do
        end do
        return
    end if

    ! a is not just a triangular matrix. construct lu factorization
    do i=1_ip,nrowm1
        ! w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0_wp) then
            iflag = 2_ip
            return
        end if
        ! jmax is the number of (nonzero) entries in column i
        ! below the diagonal.
        jmax = min(nbandl,nrow-i)
        ! divide each entry in column i below diagonal by pivot.
        do j=1_ip,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
        end do
        ! kmax is the number of (nonzero) entries in row i to
        ! the right of the diagonal.
        kmax = min(nbandu,nrow-i)
        ! subtract a(i,i+k)*(i-th column) from (i+k)-th column
        ! (below row i).
        do k=1_ip,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1_ip,jmax
                w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
        end do
    end do

    ! check the last diagonal entry.
    if (w(middle,nrow)==0.0_wp) iflag = 2_ip

    end subroutine dbnfac
!*****************************************************************************************

!*****************************************************************************************
!>
!  Companion routine to [[dbnfac]]. it returns the solution x of the
!  linear system a*x = b in place of b, given the lu-factorization
!  for a in the work array w from dbnfac.
!
!  (with \( a = l*u \), as stored in w), the unit lower triangular system
!  \( l(u*x) = b \) is solved for \( y = u*x \), and y stored in b. then the
!  upper triangular system \(u*x = y \) is solved for x. the calculations
!  are so arranged that the innermost loops stay within columns.
!
!### History
!  * banslv written by carl de boor [5]
!  * dbnslv from SLATEC library [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)

    integer(ip),intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer(ip),intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer(ip),intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer(ip),intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    real(wp),dimension(nroww,nrow),intent(in) :: w !! describes the lu-factorization of a banded matrix a of
                                                   !! order `nrow` as constructed in [[dbnfac]].
    real(wp),dimension(nrow),intent(inout) :: b  !! * **in**: right side of the system to be solved
                                                 !! * **out**: the solution x, of order nrow

    integer(ip) :: i, j, jmax, middle, nrowm1

    middle = nbandu + 1_ip
    if (nrow/=1_ip) then

        nrowm1 = nrow - 1_ip
        if (nbandl/=0_ip) then

            ! forward pass
            ! for i=1,2,...,nrow-1, subtract right side(i)*(i-th column of l)
            !                       from right side (below i-th row).
            do i=1_ip,nrowm1
                jmax = min(nbandl,nrow-i)
                do j=1_ip,jmax
                    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
                end do
            end do

        end if

        ! backward pass
        ! for i=nrow,nrow-1,...,1, divide right side(i) by i-th diagonal
        !                          entry of u, then subtract right side(i)*(i-th column
        !                          of u) from right side (above i-th row).
        if (nbandu<=0_ip) then
            ! a is lower triangular.
            do i=1_ip,nrow
                b(i) = b(i)/w(1_ip,i)
            end do
            return
        end if

        i = nrow
        do
            b(i) = b(i)/w(middle,i)
            jmax = min(nbandu,i-1_ip)
            do j=1_ip,jmax
                b(i-j) = b(i-j) - b(i)*w(middle-j,i)
            end do
            i = i - 1_ip
            if (i<=1_ip) exit
        end do

    end if

    b(1_ip) = b(1_ip)/w(middle,1_ip)

    end subroutine dbnslv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculates the value of all (possibly) nonzero basis
!  functions at x of order max(jhigh,(j+1)*(index-1)), where t(k)
!  <= x <= t(n+1) and j=iwork is set inside the routine on
!  the first call when index=1.  ileft is such that t(ileft) <=
!  x < t(ileft+1).  a call to dintrv(t,n+1,x,ilo,ileft,mflag)
!  produces the proper ileft.  dbspvn calculates using the basic
!  algorithm needed in dbspvd.  if only basis functions are
!  desired, setting jhigh=k and index=1 can be faster than
!  calling dbspvd, but extra coding is required for derivatives
!  (index=2) and dbspvd is set up for this purpose.
!
!  left limiting values are set up as described in dbspvd.
!
!### Error Conditions
!
!  * improper input
!
!### History
!  * bsplvn written by carl de boor [5]
!  * dbspvn author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork,iflag)

    implicit none

    real(wp),dimension(*),intent(in)  :: t        !! knot vector of length n+k, where
                                                  !! n = number of b-spline basis functions
                                                  !! n = sum of knot multiplicities-k
                                                  !! dimension t(ileft+jhigh)
    integer(ip),intent(in)            :: jhigh    !! order of b-spline, 1 <= jhigh <= k
    integer(ip),intent(in)            :: k        !! highest possible order
    integer(ip),intent(in)            :: index    !! index = 1 gives basis functions of order jhigh
                                                  !!       = 2 denotes previous entry with work, iwork
                                                  !!         values saved for subsequent calls to
                                                  !!         dbspvn.
    real(wp),intent(in)               :: x        !! argument of basis functions, t(k) <= x <= t(n+1)
    integer(ip),intent(in)            :: ileft    !! largest integer such that t(ileft) <= x < t(ileft+1)
    real(wp),dimension(k),intent(out) :: vnikx    !! vector of length k for spline values.
    real(wp),dimension(*),intent(out) :: work     !! a work vector of length 2*k
    integer(ip),intent(out)           :: iwork    !! a work parameter.  both work and iwork contain
                                                  !! information necessary to continue for index = 2.
                                                  !! when index = 1 exclusively, these are scratch
                                                  !! variables and can be used for other purposes.
    integer(ip),intent(out)           :: iflag    !!   0: no errors
                                                  !! 201: k does not satisfy k>=1
                                                  !! 202: jhigh does not satisfy 1<=jhigh<=k
                                                  !! 203: index is not 1 or 2
                                                  !! 204: x does not satisfy t(ileft)<=x<=t(ileft+1)

    integer(ip) :: imjp1, ipj, jp1, jp1ml, l
    real(wp) :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i),
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1_ip) then
        !write(error_unit,'(A)') 'dbspvn - k does not satisfy k>=1'
        iflag = 201_ip
        return
    end if
    if (jhigh>k .or. jhigh<1_ip) then
        !write(error_unit,'(A)') 'dbspvn - jhigh does not satisfy 1<=jhigh<=k'
        iflag = 202_ip
        return
    end if
    if (index<1_ip .or. index>2_ip) then
        !write(error_unit,'(A)') 'dbspvn - index is not 1 or 2'
        iflag = 203_ip
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1_ip)) then
        !write(error_unit,'(A)') 'dbspvn - x does not satisfy t(ileft)<=x<=t(ileft+1)'
        iflag = 204_ip
        return
    end if

    iflag = 0_ip

    if (index==1_ip) then
        iwork = 1_ip
        vnikx(1_ip) = 1.0_wp
        if (iwork>=jhigh) return
    end if

    do
        ipj = ileft + iwork
        work(iwork) = t(ipj) - x
        imjp1 = ileft - iwork + 1_ip
        work(k+iwork) = x - t(imjp1)
        vmprev = 0.0_wp
        jp1 = iwork + 1_ip
        do l=1_ip,iwork
            jp1ml = jp1 - l
            vm = vnikx(l)/(work(l)+work(k+jp1ml))
            vnikx(l) = vm*work(l) + vmprev
            vmprev = vm*work(k+jp1ml)
        end do
        vnikx(jp1) = vmprev
        iwork = jp1
        if (iwork>=jhigh) exit
    end do

    end subroutine dbspvn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the b-representation (`t`,`a`,`n`,`k`) of a b-spline
!  at `x` for the function value on `ideriv=0` or any of its
!  derivatives on `ideriv=1,2,...,k-1`.  right limiting values
!  (right derivatives) are returned except at the right end
!  point `x=t(n+1)` where left limiting values are computed.  the
!  spline is defined on `t(k)` \( \le \) `x` \( \le \) `t(n+1)`.
!  dbvalu returns a fatal error message when `x` is outside of this
!  interval.
!
!  To compute left derivatives or left limiting values at a
!  knot `t(i)`, replace `n` by `i-1` and set `x=t(i), i=k+1,n+1`.
!
!### Error Conditions
!
!  * improper input
!
!### History
!  * bvalue written by carl de boor [5]
!  * dbvalu author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbvalu(t,a,n,k,ideriv,x,inbv,work,iflag,val,extrap)

    implicit none

    real(wp),intent(out)             :: val     !! the interpolated value
    integer(ip),intent(in)           :: n       !! number of b-spline coefficients.
                                                !! (sum of knot multiplicities-`k`)
    real(wp),dimension(:),intent(in) :: t       !! knot vector of length `n+k`
    real(wp),dimension(n),intent(in) :: a       !! b-spline coefficient vector of length `n`
    integer(ip),intent(in)           :: k       !! order of the b-spline, `k >= 1`
    integer(ip),intent(in)           :: ideriv  !! order of the derivative, `0 <= ideriv <= k-1`.
                                                !! `ideriv = 0` returns the b-spline value
    real(wp),intent(in)              :: x       !! argument, `t(k) <= x <= t(n+1)`
    integer(ip),intent(inout)        :: inbv    !! an initialization parameter which must be set
                                                !! to 1 the first time [[dbvalu]] is called.
                                                !! `inbv` contains information for efficient processing
                                                !! after the initial call and `inbv` must not
                                                !! be changed by the user.  distinct splines require
                                                !! distinct `inbv` parameters.
    real(wp),dimension(:),intent(inout) :: work !! work vector of length at least `3*k`
    integer(ip),intent(out)          :: iflag   !! status flag:
                                                !!
                                                !! * 0: no errors
                                                !! * 401: `k` does not satisfy `k` \( \ge \) 1
                                                !! * 402: `n` does not satisfy `n` \( \ge \) `k`
                                                !! * 403: `ideriv` does not satisfy 0 \( \le \) `ideriv` \(<\) `k`
                                                !! * 404: `x` is not greater than or equal to `t(k)`
                                                !! * 405: `x` is not less than or equal to `t(n+1)`
                                                !! * 406: a left limiting value cannot be obtained at `t(k)`
    logical,intent(in),optional :: extrap   !! if extrapolation is allowed
                                            !! (if not present, default is False)

    integer(ip) :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real(wp) :: fkmj
    real(wp) :: xt
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    val = 0.0_wp

    if (k<1_ip) then
        iflag = 401_ip  ! dbvalu - k does not satisfy k>=1
        return
    end if

    if (n<k) then
        iflag = 402_ip  ! dbvalu - n does not satisfy n>=k
        return
    end if

    if (ideriv<0_ip .or. ideriv>=k) then
        iflag = 403_ip  ! dbvalu - ideriv does not satisfy 0<=ideriv<k
        return
    end if

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    ! make a temp copy of x (for computing the
    ! interval) in case extrapolation is allowed
    if (extrapolation_allowed) then
        if (x<t(1_ip)) then
            xt = t(1_ip)
        else if(x>t(n+k)) then
            xt = t(n+k)
        else
            xt = x
        end if
    else
        xt = x
    end if

    kmider = k - ideriv

    ! find *i* in (k,n) such that t(i) <= x < t(i+1)
    ! (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1_ip
    call dintrv(t, n+1, xt, inbv, i, mflag)
    if (xt<t(k)) then
        iflag = 404_ip  ! dbvalu - x is not greater than or equal to t(k)
        return
    end if

    if (mflag/=0_ip) then

        if (xt>t(i)) then
            iflag = 405_ip  ! dbvalu - x is not less than or equal to t(n+1)
            return
        end if

        do
            if (i==k) then
                iflag = 406_ip  ! dbvalu - a left limiting value cannot be obtained at t(k)
                return
            end if
            i = i - 1_ip
            if (xt/=t(i)) exit
        end do

    end if

    ! difference the coefficients *ideriv* times
    ! work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

    imk = i - k
    do j=1_ip,k
        imkpj = imk + j
        work(j) = a(imkpj)
    end do

    if (ideriv/=0_ip) then
        do j=1_ip,ideriv
            kmj = k - j
            fkmj = real(kmj,wp)
            do jj=1_ip,kmj
                ihi = i + jj
                ihmkmj = ihi - kmj
                work(jj) = (work(jj+1_ip)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
            end do
        end do
    end if

    ! compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    ! given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

    if (ideriv/=km1) then
        ip1 = i + 1_ip
        kpk = k + k
        j1 = k + 1_ip
        j2 = kpk + 1_ip
        do j=1_ip,kmider
            ipj = i + j
            work(j1) = t(ipj) - x
            ip1mj = ip1 - j
            work(j2) = x - t(ip1mj)
            j1 = j1 + 1_ip
            j2 = j2 + 1_ip
        end do
        iderp1 = ideriv + 1_ip
        do j=iderp1,km1
            kmj = k - j
            ilo = kmj
            do jj=1_ip,kmj
                work(jj) = (work(jj+1_ip)*work(kpk+ilo)+work(jj)*&
                            work(k+jj))/(work(kpk+ilo)+work(k+jj))
                ilo = ilo - 1
            end do
        end do
    end if

    iflag = 0_ip
    val = work(1_ip)

    end subroutine dbvalu
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the largest integer `ileft` in 1 \( \le \) `ileft` \( \le \) `lxt`
!  such that `xt(ileft)` \( \le \) `x` where `xt(*)` is a subdivision of
!  the `x` interval.
!  precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   mflag=0
!         if xt(lxt) <= x           then ileft=lxt, mflag=-2
!```
!
!  that is, when multiplicities are present in the break point
!  to the left of `x`, the largest index is taken for `ileft`.
!
!### History
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 3/4/2017 : added extrapolation option.

    pure subroutine dintrv(xt,lxt,xx,ilo,ileft,mflag,extrap)

    implicit none

    integer(ip),intent(in)             :: lxt    !! length of the `xt` vector
    real(wp),dimension(:),intent(in)   :: xt     !! a knot or break point vector of length `lxt`
    real(wp),intent(in)                :: xx     !! argument
    integer(ip),intent(inout)          :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  distinct splines
                                                 !! require distinct `ilo` parameters.
    integer(ip),intent(out)            :: ileft  !! largest integer satisfying `xt(ileft)` \( \le \) `x`
    integer(ip),intent(out)            :: mflag  !! signals when `x` lies out of bounds
    logical,intent(in),optional        :: extrap !! if extrapolation is allowed
                                                 !! (if not present, default is False)

    integer(ip) :: ihi, istep, middle
    real(wp) :: x

    x = get_temp_x_for_extrap(xx,xt(1_ip),xt(lxt),extrap)

    ihi = ilo + 1_ip
    if ( ihi>=lxt ) then
        if ( x>=xt(lxt) ) then
            mflag = -2_ip
            ileft = lxt
            return
        end if
        if ( lxt<=1 ) then
            mflag = -1_ip
            ileft = 1_ip
            return
        end if
        ilo = lxt - 1_ip
        ihi = lxt
    end if

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1_ip
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=lxt ) then
                if ( x>=xt(lxt) ) then
                    mflag = -2_ip
                    ileft = lxt
                    return
                end if
                ihi = lxt
            else if ( x>=xt(ihi) ) then
                istep = istep*2_ip
                cycle
            end if
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0_ip
            ileft = ilo
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1_ip
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1_ip ) then
                ilo = 1_ip
                if ( x<xt(1_ip) ) then
                    mflag = -1_ip
                    ileft = 1_ip
                    return
                end if
            else if ( x<xt(ilo) ) then
                istep = istep*2_ip
                cycle
            end if
            exit
        end do

    end if

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        middle = (ilo+ihi)/2_ip
        if ( middle==ilo ) then
            mflag = 0_ip
            ileft = ilo
            return
        end if
        ! note. it is assumed that middle = ilo in case ihi = ilo+1
        if ( x<xt(middle) ) then
            ihi = middle
        else
            ilo = middle
        end if
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBSQAD computes the integral on `(x1,x2)` of a `k`-th order
!  b-spline using the b-representation `(t,bcoef,n,k)`.  orders
!  `k` as high as 20 are permitted by applying a 2, 6, or 10
!  point gauss formula on subintervals of `(x1,x2)` which are
!  formed by included (distinct) knots.
!
!  If orders `k` greater than 20 are needed, use [[dbfqad]] with
!  `f(x) = 1`.
!
!### Note
!  * The maximum number of significant digits obtainable in
!    DBSQAD is the smaller of ~300 and the number of digits
!    carried in `real(wp)` arithmetic.
!
!### References
!  * D. E. Amos, "Quadrature subroutines for splines and
!    B-splines", Report SAND79-1825, Sandia Laboratories,
!    December 1979.
!
!### History
!  * Author: Amos, D. E., (SNLA)
!  * 800901  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jacob Williams, 9/6/2017 : refactored to modern Fortran.
!    Added higher precision coefficients.
!
!@note Extrapolation is not enabled for this routine.

    pure subroutine dbsqad(t,bcoef,n,k,x1,x2,bquad,work,iflag)

    implicit none

    real(wp),dimension(:),intent(in)    :: t       !! knot array of length `n+k`
    real(wp),dimension(:),intent(in)    :: bcoef   !! b-spline coefficient array of length `n`
    integer(ip),intent(in)              :: n       !! length of coefficient array
    integer(ip),intent(in)              :: k       !! order of b-spline, `1 <= k <= 20`
    real(wp),intent(in)                 :: x1      !! end point of quadrature interval
                                                   !! in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: x2      !! end point of quadrature interval
                                                   !! in `t(k) <= x <= t(n+1)`
    real(wp),intent(out)                :: bquad   !! integral of the b-spline over (`x1`,`x2`)
    real(wp),dimension(:),intent(inout) :: work    !! work vector of length `3*k`
    integer(ip),intent(out)             :: iflag   !! status flag:
                                                   !!
                                                   !! * 0: no errors
                                                   !! * 901: `k` does not satisfy `1<=k<=20`
                                                   !! * 902: `n` does not satisfy `n>=k`
                                                   !! * 903: `x1` or `x2` or both do
                                                   !!   not satisfy `t(k)<=x<=t(n+1)`

    integer(ip) :: i,il1,il2,ilo,inbv,jf,left,m,mf,mflag,npk,np1
    real(wp) :: a,aa,b,bb,bma,bpa,c1,gx,q,ta,tb,y1,y2
    real(wp),dimension(5) :: s  !! sum

    real(wp),dimension(9),parameter :: gpts = [ &
        &0.577350269189625764509148780501957455647601751270126876018602326483977&
        &67230293334569371539558574952522520871380513556767665664836499965082627&
        &05518373647912161760310773007685273559916067003615583077550051041144223&
        &01107628883557418222973945990409015710553455953862673016662179126619796&
        &4892168_wp,&
        &0.238619186083196908630501721680711935418610630140021350181395164574274&
        &93427563984224922442725734913160907222309701068720295545303507720513526&
        &28872175189982985139866216812636229030578298770859440976999298617585739&
        &46921613621659222233462641640013936777894532787145324672151888999339900&
        &0945406150514997832_wp,&
        &0.661209386466264513661399595019905347006448564395170070814526705852183&
        &49660714310094428640374646145642988837163927514667955734677222538043817&
        &23198010093367423918538864300079016299442625145884902455718821970386303&
        &22362011735232135702218793618906974301231555871064213101639896769013566&
        &1651261150514997832_wp,&
        &0.932469514203152027812301554493994609134765737712289824872549616526613&
        &50084420019627628873992192598504786367972657283410658797137951163840419&
        &21786180750210169211578452038930846310372961174632524612619760497437974&
        &07422632089671621172178385230505104744277222209386367655366917903888025&
        &2326771150514997832_wp,&
        &0.148874338981631210884826001129719984617564859420691695707989253515903&
        &61735566852137117762979946369123003116080525533882610289018186437654023&
        &16761969968090913050737827720371059070942475859422743249837177174247346&
        &21691485290294292900319346665908243383809435507599683357023000500383728&
        &0634351_wp,&
        &0.433395394129247190799265943165784162200071837656246496502701513143766&
        &98907770350122510275795011772122368293504099893794727422475772324920512&
        &67741032822086200952319270933462032011328320387691584063411149801129823&
        &14148878744320432476641442157678880770848387945248811854979703928792696&
        &4254222_wp,&
        &0.679409568299024406234327365114873575769294711834809467664817188952558&
        &57539507492461507857357048037949983390204739931506083674084257663009076&
        &82741718202923543197852846977409718369143712013552962837733153108679126&
        &93254495485472934132472721168027426848661712101171203022718105101071880&
        &4444161_wp,&
        &0.865063366688984510732096688423493048527543014965330452521959731845374&
        &75513805556135679072894604577069440463108641176516867830016149345356373&
        &92729396890950011571349689893051612072435760480900979725923317923795535&
        &73929059587977695683242770223694276591148364371481692378170157259728913&
        &9322313_wp,&
        &0.973906528517171720077964012084452053428269946692382119231212066696595&
        &20323463615962572356495626855625823304251877421121502216860143447777992&
        &05409587259942436704413695764881258799146633143510758737119877875210567&
        &06745243536871368303386090938831164665358170712568697066873725922944928&
        &4383797_wp]

    real(wp),dimension(9),parameter :: gwts = [ &
        &1.0_wp,&
        &0.467913934572691047389870343989550994811655605769210535311625319963914&
        &20162039812703111009258479198230476626878975479710092836255417350295459&
        &35635592733866593364825926382559018030281273563502536241704619318259000&
        &99756987095900533474080074634376824431808173206369174103416261765346292&
        &7888917150514997832_wp,&
        &0.360761573048138607569833513837716111661521892746745482289739240237140&
        &03783726171832096220198881934794311720914037079858987989027836432107077&
        &67872114085818922114502722525757771126000732368828591631602895111800517&
        &40813685547074482472486101183259931449817216402425586777526768199930950&
        &3106873150514997832_wp,&
        &0.171324492379170345040296142172732893526822501484043982398635439798945&
        &76054234015464792770542638866975211652206987440430919174716746217597462&
        &96492293180314484520671351091683210843717994067668872126692485569940481&
        &59429327357024984053433824182363244118374610391205239119044219703570297&
        &7497812150514997832_wp,&
        &0.295524224714752870173892994651338329421046717026853601354308029755995&
        &93821715232927035659579375421672271716440125255838681849078955200582600&
        &19363424941869666095627186488841680432313050615358674090830512706638652&
        &87483901746874726597515954450775158914556548308329986393605934912382356&
        &670244_wp,&
        &0.269266719309996355091226921569469352859759938460883795800563276242153&
        &43231917927676422663670925276075559581145036869830869292346938114524155&
        &64658846634423711656014432259960141729044528030344411297902977067142537&
        &53480628460839927657500691168674984281408628886853320804215041950888191&
        &6391898_wp,&
        &0.219086362515982043995534934228163192458771870522677089880956543635199&
        &91065295128124268399317720219278659121687281288763476662690806694756883&
        &09211843316656677105269915322077536772652826671027878246851010208832173&
        &32006427348325475625066841588534942071161341022729156547776892831330068&
        &8702802_wp,&
        &0.149451349150580593145776339657697332402556639669427367835477268753238&
        &65472663001094594726463473195191400575256104543633823445170674549760147&
        &13716011937109528798134828865118770953566439639333773939909201690204649&
        &08381561877915752257830034342778536175692764212879241228297015017259084&
        &2897331_wp,&
        &0.066671344308688137593568809893331792857864834320158145128694881613412&
        &06408408710177678550968505887782109005471452041933148750712625440376213&
        &93049873169940416344953637064001870112423155043935262424506298327181987&
        &18647480566044117862086478449236378557180717569208295026105115288152794&
        &421677_wp]

    iflag = 0_ip
    bquad = 0.0_wp

    if ( k<1_ip .or. k>20_ip ) then

        iflag = 901_ip ! error return

    else if ( n<k ) then

        iflag = 902_ip ! error return

    else

        aa = min(x1,x2)
        bb = max(x1,x2)
        if ( aa>=t(k) ) then
            np1 = n + 1_ip
            if ( bb<=t(np1) ) then
            if ( aa==bb ) return
            npk = n + k
            ! selection of 2, 6, or 10 point gauss formula
            jf = 0_ip
            mf = 1_ip
            if ( k>4_ip ) then
                jf = 1_ip
                mf = 3_ip
                if ( k>12_ip ) then
                    jf = 4_ip
                    mf = 5_ip
                end if
            end if
            do i = 1_ip , mf
                s(i) = 0.0_wp
            end do
            ilo = 1_ip
            inbv = 1_ip
            call dintrv(t,npk,aa,ilo,il1,mflag)
            call dintrv(t,npk,bb,ilo,il2,mflag)
            if ( il2>=np1 ) il2 = n
            do left = il1 , il2
                ta = t(left)
                tb = t(left+1_ip)
                if ( ta/=tb ) then
                    a = max(aa,ta)
                    b = min(bb,tb)
                    bma = 0.5_wp*(b-a)
                    bpa = 0.5_wp*(b+a)
                    do m = 1_ip , mf
                        c1 = bma*gpts(jf+m)
                        gx = -c1 + bpa
                        call dbvalu(t,bcoef,n,k,0_ip,gx,inbv,work,iflag,y2)
                        if (iflag/=0_ip) return
                        gx = c1 + bpa
                        call dbvalu(t,bcoef,n,k,0_ip,gx,inbv,work,iflag,y1)
                        if (iflag/=0_ip) return
                        s(m) = s(m) + (y1+y2)*bma
                    end do
                end if
            end do
            q = 0.0_wp
            do m = 1_ip , mf
                q = q + gwts(jf+m)*s(m)
            end do
            if ( x1>x2 ) q = -q
                bquad = q
                return
            end if
        end if

        iflag = 903_ip ! error return

    end if

    end subroutine dbsqad
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the value of `x` to use for computing the interval
!  in `t`, depending on if extrapolation is allowed or not.
!
!  If extrapolation is allowed and x is < tmin or > tmax, then either
!  `tmin` or `tmax - 2.0_wp*spacing(tmax)` is returned.
!  Otherwise, `x` is returned.

    pure function get_temp_x_for_extrap(x,tmin,tmax,extrap) result(xt)

    implicit none

    real(wp),intent(in) :: x    !! variable value
    real(wp),intent(in) :: tmin !! first knot vector element for b-splines
    real(wp),intent(in) :: tmax !! last knot vector element for b-splines
    real(wp)            :: xt   !! The value returned (it will either
                                !! be `tmin`, `x`, or `tmax`)
    logical,intent(in),optional :: extrap  !! if extrapolation is allowed
                                           !! (if not present, default is False)

    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    if (extrapolation_allowed) then
        if (x<tmin) then
            xt = tmin
        else if (x>tmax) then
            ! Put it just inside the upper bound.
            ! This is sort of a hack to get
            ! extrapolation to work.
            xt = tmax - 2.0_wp*spacing(tmax)
        else
            xt = x
        end if
    else
        xt = x
    end if

    end function get_temp_x_for_extrap
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a message string associated with the status code.

    pure function get_status_message(iflag) result(msg)

    implicit none

    integer(ip),intent(in)       :: iflag  !! return code from one of the routines
    character(len=:),allocatable :: msg    !! status message associated with the flag

    character(len=10) :: istr   !! for integer to string conversion
    integer(ip)       :: istat  !! for write statement

    select case (iflag)

    case(  0_ip); msg='Successful execution'

    case( -1_ip); msg='Error in dintrv: x < xt(1_ip)'
    case( -2_ip); msg='Error in dintrv: x >= xt(lxt)'

    case(  1_ip); msg='Error in evaluate_*d: class is not initialized'

    case(  2_ip); msg='Error in db*ink: iknot out of range'
    case(  3_ip); msg='Error in db*ink: nx out of range'
    case(  4_ip); msg='Error in db*ink: kx out of range'
    case(  5_ip); msg='Error in db*ink: x not strictly increasing'
    case(  6_ip); msg='Error in db*ink: tx not non-decreasing'
    case(  7_ip); msg='Error in db*ink: ny out of range'
    case(  8_ip); msg='Error in db*ink: ky out of range'
    case(  9_ip); msg='Error in db*ink: y not strictly increasing'
    case( 10_ip); msg='Error in db*ink: ty not non-decreasing'
    case( 11_ip); msg='Error in db*ink: nz out of range'
    case( 12_ip); msg='Error in db*ink: kz out of range'
    case( 13_ip); msg='Error in db*ink: z not strictly increasing'
    case( 14_ip); msg='Error in db*ink: tz not non-decreasing'
    case( 15_ip); msg='Error in db*ink: nq out of range'
    case( 16_ip); msg='Error in db*ink: kq out of range'
    case( 17_ip); msg='Error in db*ink: q not strictly increasing'
    case( 18_ip); msg='Error in db*ink: tq not non-decreasing'
    case( 19_ip); msg='Error in db*ink: nr out of range'
    case( 20_ip); msg='Error in db*ink: kr out of range'
    case( 21_ip); msg='Error in db*ink: r not strictly increasing'
    case( 22_ip); msg='Error in db*ink: tr not non-decreasing'
    case( 23_ip); msg='Error in db*ink: ns out of range'
    case( 24_ip); msg='Error in db*ink: ks out of range'
    case( 25_ip); msg='Error in db*ink: s not strictly increasing'
    case( 26_ip); msg='Error in db*ink: ts not non-decreasing'
    case(700_ip); msg='Error in db*ink: size(x) /= size(fcn,1)'
    case(701_ip); msg='Error in db*ink: size(y) /= size(fcn,2)'
    case(702_ip); msg='Error in db*ink: size(z) /= size(fcn,3)'
    case(703_ip); msg='Error in db*ink: size(q) /= size(fcn,4)'
    case(704_ip); msg='Error in db*ink: size(r) /= size(fcn,5)'
    case(705_ip); msg='Error in db*ink: size(s) /= size(fcn,6)'
    case(706_ip); msg='Error in db*ink: size(x) /= nx'
    case(707_ip); msg='Error in db*ink: size(y) /= ny'
    case(708_ip); msg='Error in db*ink: size(z) /= nz'
    case(709_ip); msg='Error in db*ink: size(q) /= nq'
    case(710_ip); msg='Error in db*ink: size(r) /= nr'
    case(711_ip); msg='Error in db*ink: size(s) /= ns'
    case(712_ip); msg='Error in db*ink: size(tx) /= nx+kx'
    case(713_ip); msg='Error in db*ink: size(ty) /= ny+ky'
    case(714_ip); msg='Error in db*ink: size(tz) /= nz+kz'
    case(715_ip); msg='Error in db*ink: size(tq) /= nq+kq'
    case(716_ip); msg='Error in db*ink: size(tr) /= nr+kr'
    case(717_ip); msg='Error in db*ink: size(ts) /= ns+ks'
    case(800_ip); msg='Error in db*ink: size(x) /= size(bcoef,1)'
    case(801_ip); msg='Error in db*ink: size(y) /= size(bcoef,2)'
    case(802_ip); msg='Error in db*ink: size(z) /= size(bcoef,3)'
    case(803_ip); msg='Error in db*ink: size(q) /= size(bcoef,4)'
    case(804_ip); msg='Error in db*ink: size(r) /= size(bcoef,5)'
    case(805_ip); msg='Error in db*ink: size(s) /= size(bcoef,6)'

    case(100_ip); msg='Error in dbintk: k does not satisfy k>=1'
    case(101_ip); msg='Error in dbintk: n does not satisfy n>=k'
    case(102_ip); msg='Error in dbintk: x(i) does not satisfy x(i)<x(i+1) for some i'
    case(103_ip); msg='Error in dbintk: some abscissa was not in the support of the '//&
                      'corresponding basis function and the system is singular'
    case(104_ip); msg='Error in dbintk: the system of solver detects a singular system '//&
                      'although the theoretical conditions for a solution were satisfied'

    case(201_ip); msg='Error in dbspvn: k does not satisfy k>=1'
    case(202_ip); msg='Error in dbspvn: jhigh does not satisfy 1<=jhigh<=k'
    case(203_ip); msg='Error in dbspvn: index is not 1 or 2'
    case(204_ip); msg='Error in dbspvn: x does not satisfy t(ileft)<=x<=t(ileft+1)'

    case(301_ip); msg='Error in dbtpcf: n should be > 0'

    case(401_ip); msg='Error in dbvalu: k does not satisfy k>=1'
    case(402_ip); msg='Error in dbvalu: n does not satisfy n>=k'
    case(403_ip); msg='Error in dbvalu: ideriv does not satisfy 0<=ideriv<k'
    case(404_ip); msg='Error in dbvalu: x is not greater than or equal to t(k)'
    case(405_ip); msg='Error in dbvalu: x is not less than or equal to t(n+1)'
    case(406_ip); msg='Error in dbvalu: a left limiting value cannot be obtained at t(k)'

    case(501_ip); msg='Error in initialize_*d_specify_knots: tx is not the correct size (kx+nx)'
    case(502_ip); msg='Error in initialize_*d_specify_knots: ty is not the correct size (ky+ny)'
    case(503_ip); msg='Error in initialize_*d_specify_knots: tz is not the correct size (kz+nz)'
    case(504_ip); msg='Error in initialize_*d_specify_knots: tq is not the correct size (kq+nq)'
    case(505_ip); msg='Error in initialize_*d_specify_knots: tr is not the correct size (kr+nr)'
    case(506_ip); msg='Error in initialize_*d_specify_knots: ts is not the correct size (ks+ns)'

    case(601_ip); msg='Error in db*val: x value out of bounds'
    case(602_ip); msg='Error in db*val: y value out of bounds'
    case(603_ip); msg='Error in db*val: z value out of bounds'
    case(604_ip); msg='Error in db*val: q value out of bounds'
    case(605_ip); msg='Error in db*val: r value out of bounds'
    case(606_ip); msg='Error in db*val: s value out of bounds'

    case(901_ip); msg='Error in dbsqad: k does not satisfy 1<=k<=20'
    case(902_ip); msg='Error in dbsqad: n does not satisfy n>=k'
    case(903_ip); msg='Error in dbsqad: x1 or x2 or both do not satisfy t(k)<=x<=t(n+1)'

    case(1001_ip); msg='Error in dbfqad: k does not satisfy k>=1'
    case(1002_ip); msg='Error in dbfqad: n does not satisfy n>=k'
    case(1003_ip); msg='Error in dbfqad: d does not satisfy 0<=id<k'
    case(1004_ip); msg='Error in dbfqad: x1 or x2 or both do not satisfy t(k)<=x<=t(n+1)'
    case(1005_ip); msg='Error in dbfqad: tol is less than dtol or greater than 0.1'

    case(1101_ip); msg='Warning in dbsgq8: a and b are too nearly equal to allow normal integration.'
    case(1102_ip); msg='Error in dbsgq8: ans is probably insufficiently accurate.'

    case default
        write(istr,fmt='(I10)',iostat=istat) iflag
        msg = 'Unknown status flag: '//trim(adjustl(istr))
    end select

    end function get_status_message
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_sub_module
!*****************************************************************************************
