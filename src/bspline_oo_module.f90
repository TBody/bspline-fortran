!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!  date: 12/6/2015
!
!  Object-oriented style wrappers to [[bspline_sub_module]].
!  This module provides classes ([[bspline_1d(type)]], [[bspline_2d(type)]],
!  [[bspline_3d(type)]], [[bspline_4d(type)]], [[bspline_5d(type)]], and [[bspline_6d(type)]])
!  which can be used instead of the main subroutine interface.

    module bspline_oo_module

    use bspline_kinds_module, only: wp, ip
    use,intrinsic :: iso_fortran_env, only: error_unit
    use bspline_sub_module

    implicit none

    private

    integer(ip),parameter :: int_size     = storage_size(1_ip,kind=ip)   !! size of a default integer [bits]
    integer(ip),parameter :: logical_size = storage_size(.true.,kind=ip) !! size of a default logical [bits]
    integer(ip),parameter :: real_size    = storage_size(1.0_wp,kind=ip) !! size of a `real(wp)` [bits]

    type,public,abstract :: bspline_class
        !! Base class for the b-spline types
        private
        integer(ip) :: inbvx = 1_ip  !! internal variable used by [[dbvalu]] for efficient processing
        integer(ip) :: iflag = 1_ip  !! saved `iflag` from the list routine call.
        logical :: initialized = .false. !! true if the class is initialized and ready to use
        logical :: extrap = .false. !! if true, then extrapolation is allowed during evaluation
    contains
        private
        procedure,non_overridable :: destroy_base  !! destructor for the abstract type
        procedure,non_overridable :: set_extrap_flag !! internal routine to set the `extrap` flag
        procedure(destroy_func),deferred,public :: destroy  !! destructor
        procedure(size_func),deferred,public :: size_of !! size of the structure in bits
        procedure,public,non_overridable :: status_ok  !! returns true if the last `iflag` status code was `=0`.
        procedure,public,non_overridable :: status_message => get_bspline_status_message  !! retrieve the last
                                                                                          !! status message
        procedure,public,non_overridable :: clear_flag => clear_bspline_flag  !! to reset the `iflag` saved in the class.
    end type bspline_class

    abstract interface

        pure subroutine destroy_func(me)
        !! interface for bspline destructor routines
        import :: bspline_class
        implicit none
        class(bspline_class),intent(inout) :: me
        end subroutine destroy_func

        pure function size_func(me) result(s)
        !! interface for size routines
        import :: bspline_class,ip
        implicit none
        class(bspline_class),intent(in) :: me
        integer(ip) :: s !! size of the structure in bits
        end function size_func

    end interface

    type,extends(bspline_class),public :: bspline_2d
        !! Class for 2d b-spline interpolation.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: ny = 0_ip  !! Number of \(y\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        integer(ip) :: ky = 0_ip  !! The order of spline pieces in \(y\)
        real(wp),dimension(:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        integer(ip) :: inbvy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloy = 1_ip  !! internal variable used for efficient processing
        real(wp),dimension(:),allocatable :: work_val_1  !! [[db2val] work array of dimension `ky`
        real(wp),dimension(:),allocatable :: work_val_2  !! [[db2val] work array of dimension `3_ip*max(kx,ky)`
        contains
        private
        generic,public :: initialize => initialize_2d_auto_knots,initialize_2d_specify_knots
        procedure :: initialize_2d_auto_knots
        procedure :: initialize_2d_specify_knots
        procedure,public :: evaluate => evaluate_2d
        procedure,public :: destroy => destroy_2d
        procedure,public :: size_of => size_2d
        final :: finalize_2d
    end type bspline_2d

    interface bspline_2d
        !! Constructor for [[bspline_2d(type)]]
        procedure :: bspline_2d_constructor_empty,&
                     bspline_2d_constructor_auto_knots,&
                     bspline_2d_constructor_specify_knots
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  This routines returns true if the `iflag` code from the last
!  routine called was `=0`. Maybe of the routines have output `iflag`
!  variables, so they can be checked explicitly, or this routine
!  can be used.
!
!  If the class is initialized using a function constructor, then
!  this is the only way to know if it was properly initialized,
!  since those are pure functions with not output `iflag` arguments.
!
!  If `status_ok=.false.`, then the error message can be
!  obtained from the [[get_bspline_status_message]] routine.
!
!  Note: after an error condition, the [[clear_bspline_flag]] routine
!  can be called to reset the `iflag` to 0.

    elemental function status_ok(me) result(ok)

    implicit none

    class(bspline_class),intent(in) :: me
    logical                         :: ok

    ok = ( me%iflag == 0_ip )

    end function status_ok
!*****************************************************************************************

!*****************************************************************************************
!>
!  This sets the `iflag` variable in the class to `0`
!  (which indicates that everything is OK). It can be used
!  after an error is encountered.

    elemental subroutine clear_bspline_flag(me)

    implicit none

    class(bspline_class),intent(inout) :: me

    me%iflag = 0_ip

    end subroutine clear_bspline_flag
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get the status message from a [[bspline_class]] routine call.
!
!  If `iflag` is not included, then the one in the class is used (which
!  corresponds to the last routine called.)
!  Otherwise, it will convert the
!  input `iflag` argument into the appropriate message.
!
!  This is a wrapper for [[get_status_message]].

    pure function get_bspline_status_message(me,iflag) result(msg)

    implicit none

    class(bspline_class),intent(in) :: me
    character(len=:),allocatable    :: msg    !! status message associated with the flag
    integer(ip),intent(in),optional :: iflag  !! the corresponding status code

    if (present(iflag)) then
        msg = get_status_message(iflag)
    else
        msg = get_status_message(me%iflag)
    end if

    end function get_bspline_status_message
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_2d]] structure in bits.

    pure function size_2d(me) result(s)

    implicit none

    class(bspline_2d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 6_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,1_ip,kind=ip)*&
                                                    size(me%bcoef,2_ip,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%ty))         s = s + real_size*size(me%ty,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,kind=ip)
    if (allocated(me%work_val_2)) s = s + real_size*size(me%work_val_2,kind=ip)

    end function size_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for contents of the base [[bspline_class]] class.
!  (this routine is called by the extended classes).

    pure subroutine destroy_base(me)

    implicit none

    class(bspline_class),intent(inout) :: me

    me%inbvx = 1_ip
    me%iflag = 1_ip
    me%initialized = .false.
    me%extrap = .false.

    end subroutine destroy_base
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_2d]] class.

    pure subroutine destroy_2d(me)

    implicit none

    class(bspline_2d),intent(inout) :: me

    call me%destroy_base()

    me%nx    = 0_ip
    me%ny    = 0_ip
    me%kx    = 0_ip
    me%ky    = 0_ip
    me%inbvy = 1_ip
    me%iloy  = 1_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%ty))         deallocate(me%ty)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)
    if (allocated(me%work_val_2)) deallocate(me%work_val_2)

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for [[bspline_2d]] class. Just a wrapper for [[destroy_2d]].
    pure elemental subroutine finalize_2d(me)
        type(bspline_2d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sets the `extrap` flag in the class.

    pure subroutine set_extrap_flag(me,extrap)

    implicit none

    class(bspline_class),intent(inout) :: me
    logical,intent(in),optional :: extrap  !! if not present, then False is used

    if (present(extrap)) then
        me%extrap = extrap
    else
        me%extrap = .false.
    end if

    end subroutine set_extrap_flag
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_2d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    elemental function bspline_2d_constructor_empty() result(me)

    implicit none

    type(bspline_2d) :: me

    end function bspline_2d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_2d]] type (auto knots).
!  This is a wrapper for [[initialize_2d_auto_knots]].

    pure function bspline_2d_constructor_auto_knots(x,y,fcn,kx,ky,extrap) result(me)

    implicit none

    type(bspline_2d)                   :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    logical,intent(in),optional      :: extrap  !! if true, then extrapolation is allowed
                                                !! (default is false)

    call initialize_2d_auto_knots(me,x,y,fcn,kx,ky,me%iflag,extrap)

    end function bspline_2d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_2d]] type (user-specified knots).
!  This is a wrapper for [[initialize_2d_specify_knots]].

    pure function bspline_2d_constructor_specify_knots(x,y,fcn,kx,ky,tx,ty,extrap) result(me)

    implicit none

    type(bspline_2d)                   :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)   :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)   :: ty    !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    logical,intent(in),optional      :: extrap  !! if true, then extrapolation is allowed
                                                !! (default is false)

    call initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,me%iflag,extrap)

    end function bspline_2d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_2d]] type (with automatically-computed knots).
!  This is a wrapper for [[db2ink]].

    pure subroutine initialize_2d_auto_knots(me,x,y,fcn,kx,ky,iflag,extrap)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(out)            :: iflag !! status flag (see [[db2ink]])
    logical,intent(in),optional        :: extrap !! if true, then extrapolation is allowed
                                                 !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx,ny

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)

    me%nx = nx
    me%ny = ny

    me%kx = kx
    me%ky = ky

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%bcoef(nx,ny))
    allocate(me%work_val_1(ky))
    allocate(me%work_val_2(3_ip*max(kx,ky)))

    iknot = 0_ip         !knot sequence chosen by db2ink

    call db2ink(x,nx,y,ny,fcn,kx,ky,iknot,me%tx,me%ty,me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_2d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_2d]] type (with user-specified knots).
!  This is a wrapper for [[db2ink]].

    pure subroutine initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,iflag,extrap)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)   :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)   :: ty    !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    integer(ip),intent(out)            :: iflag !! status flag (see [[db2ink]])
    logical,intent(in),optional      :: extrap  !! if true, then extrapolation is allowed
                                                !! (default is false)

    integer(ip) :: nx,ny

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,&
                                  ny=ny,ky=ky,ty=ty,&
                                  iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%ny = ny

        me%kx = kx
        me%ky = ky

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%bcoef(nx,ny))
        allocate(me%work_val_1(ky))
        allocate(me%work_val_2(3_ip*max(kx,ky)))

        me%tx = tx
        me%ty = ty

        call db2ink(x,nx,y,ny,fcn,kx,ky,1_ip,me%tx,me%ty,me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_2d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_2d]] interpolate.  This is a wrapper for [[db2val]].

    pure subroutine evaluate_2d(me,xval,yval,idx,idy,f,iflag)

    implicit none

    class(bspline_2d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)             :: yval  !! \(y\) coordinate of evaluation point.
    integer(ip),intent(in)           :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)           :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db2val]])

    if (me%initialized) then
        call db2val(xval,yval,&
                    idx,idy,&
                    me%tx,me%ty,&
                    me%nx,me%ny,&
                    me%kx,me%ky,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%iloy,&
                    me%work_val_1,me%work_val_2,&
                    extrap=me%extrap)
    else
        iflag = 1_ip
    end if

    me%iflag = iflag

    end subroutine evaluate_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Error checks for the user-specified knot vector sizes.
!
!@note If more than one is the wrong size, then the `iflag` error code will
!      correspond to the one with the highest rank.

    pure subroutine check_knot_vectors_sizes(nx,ny,nz,nq,nr,ns,&
                                             kx,ky,kz,kq,kr,ks,&
                                             tx,ty,tz,tq,tr,ts,iflag)

    implicit none

    integer(ip),intent(in),optional           :: nx
    integer(ip),intent(in),optional           :: ny
    integer(ip),intent(in),optional           :: nz
    integer(ip),intent(in),optional           :: nq
    integer(ip),intent(in),optional           :: nr
    integer(ip),intent(in),optional           :: ns
    integer(ip),intent(in),optional           :: kx
    integer(ip),intent(in),optional           :: ky
    integer(ip),intent(in),optional           :: kz
    integer(ip),intent(in),optional           :: kq
    integer(ip),intent(in),optional           :: kr
    integer(ip),intent(in),optional           :: ks
    real(wp),dimension(:),intent(in),optional :: tx
    real(wp),dimension(:),intent(in),optional :: ty
    real(wp),dimension(:),intent(in),optional :: tz
    real(wp),dimension(:),intent(in),optional :: tq
    real(wp),dimension(:),intent(in),optional :: tr
    real(wp),dimension(:),intent(in),optional :: ts
    integer(ip),intent(out)                   :: iflag  !! 0 if everything is OK

    iflag = 0_ip

    if (present(nx) .and. present(kx) .and. present(tx)) then
        if (size(tx,kind=ip)/=(nx+kx)) then
            iflag = 501_ip  ! tx is not the correct size (nx+kx)
        end if
    end if

    if (present(ny) .and. present(ky) .and. present(ty)) then
        if (size(ty,kind=ip)/=(ny+ky)) then
            iflag = 502_ip  ! ty is not the correct size (ny+ky)
        end if
    end if

    if (present(nz) .and. present(kz) .and. present(tz)) then
        if (size(tz,kind=ip)/=(nz+kz)) then
            iflag = 503_ip  ! tz is not the correct size (nz+kz)
        end if
    end if

    if (present(nq) .and. present(kq) .and. present(tq)) then
        if (size(tq,kind=ip)/=(nq+kq)) then
            iflag = 504_ip  ! tq is not the correct size (nq+kq)
        end if
    end if

    if (present(nr) .and. present(kr) .and. present(tr)) then
        if (size(tr,kind=ip)/=(nr+kr)) then
            iflag = 505_ip  ! tr is not the correct size (nr+kr)
        end if
    end if

    if (present(ns) .and. present(ks) .and. present(ts)) then
        if (size(ts,kind=ip)/=(ns+ks)) then
            iflag = 506_ip  ! ts is not the correct size (ns+ks)
        end if
    end if

    end subroutine check_knot_vectors_sizes
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_oo_module
!*****************************************************************************************
