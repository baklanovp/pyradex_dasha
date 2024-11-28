!      main.f

!  This file is part of the RADEX software package
!  to calculate molecular excitation and radiative
!  transfer in a homogeneous medium.

!  Documentation for the program is posted at
!  https://sron.rug.nl/~vdtak/radex/index.shtml

!  Although this program has been thoroughly tested, the
!  authors do not claim that it is free of errors and
!  gives correct results in all situations.

!  Publications using this program should make a reference
!  to our paper: A&A 468, 627 (2007).

!      ---------------------------------------------------------

PROGRAM radex
   use mRadexInc
   use cla,   only: cla_init, cla_register, cla_get, cla_help, cla_key_present &
      , cla_flag, cla_int, cla_char, STRLEN;
   implicit none
   !      Main program: controls program flow and drives subroutines

   ! include 'radex.inc'

   integer niter   ! iteration counter
   integer imore   ! are we running again?
   logical conv    ! are we converged?

   !! !   ***************  UNIT    ************************
   call cla_init();
   call cla_register('-info',   'Print extanded information',  cla_flag, 'f');
   call cla_register('-h',     'Print this help',  cla_flag, 'f');
   !!!   *************** Read command line arguments    ************************
   !!! !Check if any arguments are found
!    narg = command_argument_count();
   !!!Loop over the arguments
   if( cla_key_present('-h') ) then;
      call cla_help();
      stop;
   endif;
   
   is_debug = cla_key_present('-info');


   !      Begin executable statements
   print*
   print*,'   Welcome to Radex, version of '//version
   print*

   !      Get input parameters
   if (is_debug) print*,'calling getinputs'
21 call getinputs

   !      Read data file
   if (is_debug) print*,'calling readdata'
   call readdata

   !      Calculate background radiation field
   if (is_debug) print*,'calling backrad'
   call backrad

   niter = 0
   conv  = .false.

   !      Set up rate matrix, splitting it in radiative and collisional parts
   !      Invert rate matrix to get `thin' starting condition
   if (is_debug) print*,'calling matrix'
   call matrix(niter,conv)

   !      Start iterating
   do 22 niter=1,maxiter

      !      Invert rate matrix using escape probability for line photons
      call matrix(niter,conv)
      if (conv) then
         print*,'Finished in ',niter,' iterations.'
         go to 23
      endif
22 continue

   print*,'   Warning: Calculation did not converge in ',maxiter,' iterations.'

   !      Write output
   if (is_debug) print*,'calling output'
23 call output(niter)

   !      See if user wants more, else call it a day
51 format(A,$)
   write(*,51) '  Another calculation [0/1] ? '
   read(*,*) imore
   write(13,52) imore
52 format(i2)
   if (imore.eq.1) go to 21
   write(*,*) '   Have a nice day.'
   !      Done! Now close log file ...
   close(13)
   !      ...and output file.
   close(8)
   stop
end program radex
