! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read EQUILIBIUM namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (QC, NU, PC, MU, EPSA,&
     CFLAG,&
     HSFT, VSFT,&
     RO, ZO, DRO, WO,&
     R1, Z1, DR1, W1,&
     R2, Z2, DR2, W2,&
     R3, Z3, DR3, W3,&
     R4, Z4, DR4, W4,&
     R5, Z5, DR5, W5,&
     R6, Z6, DR6, W6,&
     PPED, WPED, RPED,&
     EPS, NS, NR, NF, NW,&
     ACC, H0, HMIN, HMAX)&
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: QC
  real    (kind = c_double), intent (inout) :: NU
  real    (kind = c_double), intent (inout) :: PC
  real    (kind = c_double), intent (inout) :: MU
  real    (kind = c_double), intent (inout) :: EPSA

  integer (kind = c_int),    intent (inout) :: CFLAG

  real    (kind = c_double), intent (inout) :: HSFT
  real    (kind = c_double), intent (inout) :: VSFT

  real    (kind = c_double), intent (inout) :: RO
  real    (kind = c_double), intent (inout) :: ZO
  real    (kind = c_double), intent (inout) :: DRO
  real    (kind = c_double), intent (inout) :: WO

  real    (kind = c_double), intent (inout) :: R1
  real    (kind = c_double), intent (inout) :: Z1
  real    (kind = c_double), intent (inout) :: DR1
  real    (kind = c_double), intent (inout) :: W1

  real    (kind = c_double), intent (inout) :: R2
  real    (kind = c_double), intent (inout) :: Z2
  real    (kind = c_double), intent (inout) :: DR2
  real    (kind = c_double), intent (inout) :: W2
  
  real    (kind = c_double), intent (inout) :: R3
  real    (kind = c_double), intent (inout) :: Z3
  real    (kind = c_double), intent (inout) :: DR3
  real    (kind = c_double), intent (inout) :: W3

  real    (kind = c_double), intent (inout) :: R4
  real    (kind = c_double), intent (inout) :: Z4
  real    (kind = c_double), intent (inout) :: DR4
  real    (kind = c_double), intent (inout) :: W4

  real    (kind = c_double), intent (inout) :: R5
  real    (kind = c_double), intent (inout) :: Z5
  real    (kind = c_double), intent (inout) :: DR5
  real    (kind = c_double), intent (inout) :: W5

  real    (kind = c_double), intent (inout) :: R6
  real    (kind = c_double), intent (inout) :: Z6
  real    (kind = c_double), intent (inout) :: DR6
  real    (kind = c_double), intent (inout) :: W6

  real    (kind = c_double), intent (inout) :: PPED
  real    (kind = c_double), intent (inout) :: WPED
  real    (kind = c_double), intent (inout) :: RPED

  real    (kind = c_double), intent (inout) :: EPS
  integer (kind = c_int),    intent (inout) :: NS
  integer (kind = c_int),    intent (inout) :: NR
  integer (kind = c_int),    intent (inout) :: NF
  integer (kind = c_int),    intent (inout) :: NW

  real    (kind = c_double), intent (inout) :: ACC
  real    (kind = c_double), intent (inout) :: H0
  real    (kind = c_double), intent (inout) :: HMIN
  real    (kind = c_double), intent (inout) :: HMAX
 
  namelist /EQUILIBRIUM_CONTROL/ QC, NU, PC, MU, EPSA,&
       CFLAG,&
       HSFT, VSFT,&
       RO, ZO, DRO, WO,&
       R1, Z1, DR1, W1,&
       R2, Z2, DR2, W2,&
       R3, Z3, DR3, W3,&
       R4, Z4, DR4, W4,&
       R5, Z5, DR5, W5,&
       R6, Z6, DR6, W6,&
       PPED, WPED, RPED,&
       EPS, NS, NR, NF, NW,&
       ACC, H0, HMIN, HMAX
       
  open  (unit = 100, file = 'Inputs/Equilibrium.nml', status = 'old')
  read  (unit = 100, nml  = EQUILIBRIUM_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
