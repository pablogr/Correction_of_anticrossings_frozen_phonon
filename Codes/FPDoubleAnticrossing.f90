
Program ConcavityWithdoubleAnticrossing

! The frozen-phonon method for renormalization of electronic eigenvalues [Capaz et al. PRL 2005] 
! presents one important drawback: even for small nuclear displacements (h), there exist vibrational
! modes where there is a severe inaccuracy due to anticrossing of electronic states. This phenomenon
! completely distorts the eigenvalue_vs_h curve, resulting in severe inaccuracies in the frozen-phonon
! calculation. For example, the typical contribution of a mode to the renormalization of HOMO and LUMO
! in diamondoids is a few meV. However, for modes with anticrossing, this can be above 50 meV. Hence,
! a few modes with anticrossing are responsible for severe dips in accuracy of the total renormalizations.
! This program intends to solve this problem. To this end, we solve the equation of double anticrossing.  
! Proceeding this way, we find the BARE_eigenvalues_vs_h, discounting the anticrossing effect of mixing of
! states. We do it by solving the following equation:
!  
!    | E1-E	  g   |
!	 | g     E2-E | = 0

!
! where E1,E2 are the bare eigenvalues (unknowns for us), E are the dressed eigenvalues
! (E:=E+,E-; they are known from DFT), and g is the coupling (which are calculated
! as half the minimum distance between E2 and E1). IMPORTANT: Note that in our
! notation, the state whose concavity is sought for the calculation of frozen-phonon renormalization is E2
! (e.g. E2=HOMO), and E2 has an anticrossing with E1 (e.g. E1=HOMO-1).
!
! The program below has 3 parts: 
! 1) Initialization: Reading the file with input parameters (FPanticrossing.in) and the files containing the eigenvalues
!    as a function of h (eigenvalue1.dat, eigenvalue2.dat )
! 2) Calculation of couplings: g is calculated as half the distance between the curves of the dressed 
!    eigenvalues; to find it, we build such curves by polynomial interpolation using the data from DFT (from the
!    files which contain the eigenvalues).
! 3) Solution of the equation above. Since the equation is quadratic, it can be analytically solved
!
!
! An example of the input file (FPanticrossing.in) to run the code is the following:
! 
! # Input file for the program for the calculation of the
! # concavity of the bare eigenvalue of an eigenvalue subject
! # to double anticrossing
! 
! Nh = 15
! NhToSolveAnticrossingEquation = 9
! 
! PolynomialDegreeg = 6
! PolynomialDegreeConcavity = 4
! 
! Nintervals = 4097
!
!  ---------------------------------------------------------------------------------------
!
! Being eigenvalue1.dat:
! -16     -6.272565915
! -12     -6.232954997
! -8      -6.198708612
! -6      -6.183874772
! -4      -6.171929038
! -2      -6.164014667
! -1      -6.161926171
! 0       -6.161218334
! 1       -6.161923138
! 2       -6.164008985
! 4       -6.171918496
! 6       -6.18386064
! 8       -6.198692282
! 12      -6.232943587
! 16      -6.27254571
! 
! eigenvalue2.dat
! -16     -5.870199204
! -12     -5.909417627
! -8      -5.944644141
! -6      -5.959556964
! -4      -5.971558114
! -2      -5.979505951
! -1      -5.981603715
! 0       -5.982315736
! 1       -5.981608561
! 2       -5.979515958
! 4       -5.971576763
! 6       -5.959583067
! 8       -5.944676912
! 12      -5.909460719
! 16      -5.870252685
! 
! 



implicit none

integer :: N, Nintervals, NhToSolveAnticrossingEquation, polynomial_degree_g, polynomial_degree_concavity 
real(8) :: Rsq, g
real(8), allocatable :: X(:), Y(:), Yfit(:), h(:), eigv1(:), eigv2(:) 
real(8), allocatable :: fitting_coefs_eigv1(:), fitting_coefs_eigv2(:) 




call read_input_file ( N, polynomial_degree_g, polynomial_degree_concavity, NhToSolveAnticrossingEquation, Nintervals, g )

allocate( Yfit(N), h(N), eigv1(N), eigv2(N) )
allocate( fitting_coefs_eigv1(polynomial_degree_g+1), fitting_coefs_eigv2(polynomial_degree_g+1)  )


call read_eigenvalue_files ( N, h, eigv1, eigv2 )


write(*,*) '  Now calculating the coupling coefficients as a function of the distances between eigenvalue curves.'
call least_squares_Lapack ( N, h, eigv1, polynomial_degree_g, fitting_coefs_eigv1, Yfit, Rsq )
call least_squares_Lapack ( N, h, eigv2, polynomial_degree_g, fitting_coefs_eigv2, Yfit, Rsq )

if (g .lt. -0.1d0) call find_g (N,polynomial_degree_g, Nintervals,h,fitting_coefs_eigv1,fitting_coefs_eigv2, g)


! SOLUTION OF THE double ANTICROSSING EQUATION (to find the concavity of the bare E2)

call solve_double_anticrossing_equation ( N, NhToSolveAnticrossingEquation, polynomial_degree_concavity, & 
                                        & g, h, eigv1, eigv2)



deallocate(Yfit, h, eigv1, eigv2, fitting_coefs_eigv1, fitting_coefs_eigv2  )


end

! =============================================================


subroutine read_input_file (N,polynomial_degree_g,polynomial_degree_concavity, NhToSolveAnticrossingEquation, Nintervals,g  )
                         
  ! System parameters
  integer, Intent(Out) :: N              				 ! The number of points of the eigenvalue files (number of h's)
  integer, Intent(Out) :: polynomial_degree_g 	 		 ! Chosen degree of the interpolation polynomials in the calculation of g
  integer, Intent(Out) :: polynomial_degree_concavity 	 ! Chosen degree of the interpolation polynomial at the final fit for the concavity
  integer, Intent(Out) :: NhToSolveAnticrossingEquation  ! The number of intervals for the division of h (ordinates) in the calculation of g
  integer, Intent(Out) :: Nintervals            		 ! The number of intervals for the division of h used for solving the double anticrossing equation and to find the concavity of the bare eigenvalue
  real(8), Intent(Out) :: g       						 ! Coupling coefficients for the double anticrossing. 
  

 
  ! Inner parameters
  integer :: i,j,k, stat, iostat
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character :: letter 
  logical :: exist 
  
   write(*,*) 
   write(*,*) '  *******************************************************************************' 
   write(*,*) '   Now calculating the concavity of the bare eivenvalue in a 2-state anticrossing '
   write(*,*) '  *******************************************************************************' 
  
   inquire(file='FPanticrossing.in', exist=exist)      
   if (.not. exist) then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************' 
	   write(*,*) '    ERROR: The file "FPanticrossing.in" does not exist. Please, provide it.'
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1) 
   end if
   
  NhToSolveAnticrossingEquation=0 
  polynomial_degree_g=6
  polynomial_degree_concavity=4
  g=-1.d0 

  open(345,file='FPanticrossing.in',status='unknown')

   do k=1,50

      read(345,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 1117

      auxstr1=''
      auxstr2=''
      auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '=') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 1235
        end if
      end do ! i
 
      1235 continue
  
      do j=i+1,100
        letter = inputline(j:j)
        auxstr2 = trim(auxstr2)//letter 
      end do ! j


     if ((auxstr1 .eq. 'Nh').or.(auxstr1 .eq. 'Npoints').or.(auxstr1 .eq. 'N'))   then
          Read(auxstr2, '(I5)' ) N
     end if  
     
     if ((auxstr1 .eq. 'Nintervals').or.(auxstr1 .eq. 'nintervals').or.(auxstr1 .eq. 'NIntervals'))   then
          Read(auxstr2, '(I9)' ) Nintervals
     end if  

     if ((auxstr1 .eq. 'NhToSolveAnticrossingEquation').or.(auxstr1 .eq. 'Nhtosolveanticrossingequation') )   then
          Read(auxstr2, '(I4)' ) NhToSolveAnticrossingEquation
     end if
     
     if ((auxstr1 .eq. 'PolynomialDegreeg').or.(auxstr1 .eq. 'polynomialdegreeg').or.(auxstr1 .eq. 'Polynomialdegreeg'))   then
          Read(auxstr2, '(I3)' ) polynomial_degree_g
     end if  
 
     if ((auxstr1 .eq. 'PolynomialDegreeConcavity').or.(auxstr1 .eq. 'polynomialdegreeconcavity') )   then
          Read(auxstr2, '(I3)' ) polynomial_degree_concavity
     end if  

     if ((auxstr1 .eq. 'g').or.(auxstr1 .eq. 'G') )   then
          Read(auxstr2, '(G18.8)' ) g
     end if  


     
  end do ! k
   

1117 continue
   
   if ( (NhToSolveAnticrossingEquation) .eq. 0 ) NhToSolveAnticrossingEquation=N
   if ( mod(NhToSolveAnticrossingEquation,2) .eq. 0 ) NhToSolveAnticrossingEquation=NhToSolveAnticrossingEquation+1
   
   write(*,*)       
   write(*,*)  '  The number of h values in the input files is ',N
   write(*,*)  '  The number of h values used to solve the threefold anticrossing equation is ', NhToSolveAnticrossingEquation  
   write(*,*)  '  The degree of the interpolating polynomials to get g is ',polynomial_degree_g
   write(*,*)  '  The degree of the interpolating polynomial for the final calculation of the concavity is ', &
                        & polynomial_degree_concavity  
   if ( g .gt. -0.1d0 ) write(*,*) '  The coupling g is ',g
   write(*,*)  
   
!     if ((auxstr1 .eq. 'Omega_units').or.(auxstr1 .eq. 'omega_units').or.(auxstr1 .eq. 'Omega_Units') .or. &
!         & (auxstr1 .eq. 'phonon_frequency_units').or.(auxstr1 .eq. 'Phonon_frequency_units').or. &
!         &   (auxstr1 .eq. 'phonon_Frequency_Units')  )   then
!          units_omegaorig = trim(auxstr2)
!          call StripSpaces (units_omegaorig)
!      end if         


end subroutine read_input_file 


!-----------------------------------------------------------------------------------------                         

subroutine read_eigenvalue_files ( N, h, eigv1, eigv2  )

! This subroutine reads the files with the eigenvalues. Each has N rows with 2 numbers in each: the displacement (h) and the eigenvalue

integer, Intent(In)  :: N    				 ! Size of the mathematical function to fit (number of points)
real(8), Intent(Out) :: h(N)  		   		 ! Vector with the displacements used for the calculation of the eigenvalues
real(8), Intent(Out) :: eigv1(N)  		   	 ! Vector with the first eigenvalue (e.g. HOMO-1a)
real(8), Intent(Out) :: eigv2(N)  		   	 ! Vector with the second eigenvalue (e.g. HOMO)
 

  ! Inner parameters
  integer :: i, j
  real(8) :: auxr0, auxr1, auxr2

  h(:) = 0.d0; eigv1(:)=0.d0; eigv2(:)=0.d0 

  

  open(100,file='eigenvalue1.dat',status='unknown')
  !read(100,*)  
  do i=1,N
    read(100,*) auxr0, auxr1
    !write(*,*) auxr0, auxr1 xxx
    j = nint( auxr0 )
    h(i) = dble(j)
    eigv1(i) = auxr1
  end do
  close(100)

  open(101,file='eigenvalue2.dat',status='unknown')
  !read(101,*)  
  do i=1,N
    read(101,*) auxr0, auxr1
    eigv2(i) = auxr1
  end do
  close(101)
  


end subroutine read_eigenvalue_files 

!-----------------------------------------------------------------------------------------                         


subroutine least_squares_Lapack ( N, X, Y, polynomial_degree, fitting_coefs, Yfit, Rsq )

! This subroutine finds the least squares fitting to a given function Y(X). It provides results
! amazingly similar to those provided by Microsoft Excel.

integer, Intent(In) :: N    				 ! Size of the mathematical function to fit (number of points)
real(8), Intent(In) :: X(N)  		   		 ! Abscise of the mathematical function to fit
real(8), Intent(In) :: Y(N)  				 ! Ordinate of the mathematical function to fit
integer, Intent(In) :: polynomial_degree     ! Degree of the polynomial of the fit
real(8), Intent(Out) :: fitting_coefs(polynomial_degree+1) 	 ! Coefficients of the fitting
real(8), Intent(Out) :: Yfit(N)  			 ! Ordinate of the fitted mathematical function
real(8), Intent(Out) :: Rsq					 ! R-squared coefficient of the fit


! Inner variables
integer :: INFO, i, j
real(8) :: auxr1, auxr2
real(8), allocatable :: A(:,:), B(:,:), WORK(:) 

allocate(A(N,polynomial_degree+1), B(N,1), WORK(4*N*(polynomial_degree+1)) )

! We build matrix A. For every value of X (i.e. for every index of X), the polynomial_degree columns
! of A are X^0, X^1, X^2, ..., X^polynomial_degree

do i=1,N
  do j=1,polynomial_degree+1
    A(i,j) = X(i)**(j-1)
  end do !j
end do !i

    
! We build matrix B, which contains the Y's
do i=1,N
    B(i,1) = Y(i)
  !  write(*,*) i,B(i,1)
end do !i


! Lapack least squares solver (for real arguments, double precision)

call dgels (  'N',       			 &  !
               N ,       			 &  ! Number of rows of the matrix (A) so that Ax=b is solved by minimum squares
               polynomial_degree+1,  &  ! Number of columns of the matrix (A)
               1,					 &  ! Number of rhs'
           	   A,				     &  ! On entry, the M-by-N matrix A. On exit, if M >= N, A is overwritten by details of its QR factorization as returned by DGEQRF;
               N,					 &  ! The leading dimension of the array A.  LDA >= max(1,M).
			   B,					 &  ! On entry, the matrix B of right hand side vectors, stored columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS if TRANS = 'N' and m >= n, rows 1 to n of B contain the least squares solution vectors; the residual sum of squares for the solution in each column is given by the sum of squares of elements N+1 to M in that column;
               N,					 &  ! Dimension of the first index of B			
               WORK,				 &  ! On exit, if INFO = 0, WORK(1) returns the optimal LWORK.			
			   4*N*(polynomial_degree+1), &  ! The dimension of the array WORK.
 			   INFO  ) 				! INFO is INTEGER: = 0:  successful exit
									!	 < 0:  if INFO = -i, the i-th argument had an illegal value
          							!  	 > 0:  if INFO =  i, the i-th diagonal element of the
              						!    triangular factor of A is zero, so that A does not have
           							!     full rank; the least squares solution could not be computed.


Do i=1,polynomial_degree+1
  fitting_coefs(i) = B(i,1)
end do


do i=1,N
  auxr1=0.d0
  do j=1,polynomial_degree+1
    auxr1 = auxr1 + B(j,1) * (X(i)**(j-1))
  end do !j
  Yfit(i)=auxr1
  ! write(*,*) X(i), Y(i), auxr1, (Y(i)-auxr1)
end do !i

call calculate_Rsquared (N, Y, Yfit, Rsq)

deallocate(A,B,WORK)

end subroutine least_squares_Lapack 

! ----------------------------------------------------------------------------------------


subroutine calculate_Rsquared (N, Ytrue, Yfit, Rsq)

integer, Intent(In) :: N    		 ! Size of the fitted mathematical function 
real(8), Intent(In) :: Ytrue(N) 	 ! Ordinate of the fitted mathematical function 
real(8), Intent(In) :: Yfit(N) 		 ! Ordinate of the fitted mathematical function 
real(8), Intent(Out) :: Rsq   		 ! Degree of the polynomial of the fit


! Internal variables
integer :: i, j
real(8) :: auxr1, auxr2, SStot, SSres, Yavg

Yavg=0.d0
do i=1,N
  Yavg=Yavg+Ytrue(i)
end do
Yavg=Yavg/(dble(N))

SStot=0.d0
do i=1,N
  SStot=SStot+ ( Ytrue(i) - Yavg )**2
end do

SSres=0.d0
do i=1,N
  SSres=SSres + ( Ytrue(i) - Yfit(i) )**2
end do


Rsq = 1.d0 - (SSres/SStot)
write(*,*) '  R-squared = ',Rsq


end subroutine calculate_Rsquared



!-------------------------------------------------------------------------------------------

subroutine least_squares_kutre (N, X, Y, m, Yfit)

integer, Intent(In) :: N     ! Size of the mathematical function to fit
real(8), Intent(In) :: X(N)  ! Abscise of the mathematical function to fit
real(8), Intent(In) :: Y(N)  ! Ordinate of the mathematical function to fit
integer, Intent(In) :: m     ! Degree of the polynomial of the fit
real(8), Intent(Out) :: Yfit(N)  ! Ordinate of the fitted mathematical function


! PGR: This subroutine is much less accurate than the Lapack-based one; for example, for the data below,
! the Lapack-based solver has R^2=0.99999330, but this subroutine has R^2=0.99904084.

! X(1) = 0.0
! X(2) = 1.0
! X(3) = 2.0
! X(4) = 3.0
! X(5) = 4.0
! X(6) = 5.0
! 
! Y(1) = 0.0
! Y(2) = 3.0
! Y(3) = 6.0
! Y(4) = 11.0
! Y(5) = 18.0
! Y(6) = 25.9

!*********************************************************
!*    Approximation of a discrete real function F(x) by  *
!*    least squares                                      *
!* ----------------------------------------------------- *
!* Ref.: "Méthodes de calcul numérique, Tome 2 By Claude *
!*        Nowakowski, PSI Edition, 1984" [BIBLI 04].     *
!* ----------------------------------------------------- *
!* SAMPLE RUN:                                           *
!*                                                       *
!* Number of points    : 11                              *
!*                                                       *
!* Degree of polynomial: 3                               *
!*                                                       *
!* Function to approximate:                              *
!*  X(1), Y(1) = 0 0                                     *
!*  X(2), Y(2) = 0.1 0.2                                 *
!*  X(3), Y(3) = 0.2 0.4                                 *
!*  X(4), Y(4) = 0.3 0.6                                 *
!*  X(5), Y(5) = 0.4 0.8                                 *
!*  X(6), Y(6) = 0.5 1                                   *
!*  X(7), Y(7) = 0.6 0.8                                 *
!*  X(8), Y(8) = 0.7 0.6                                 *
!*  X(9), Y(9) = 0.8 0.4                                 *
!*  X(10), Y(10) = 0.9 0.2                               *
!*  X(11), Y(11) = 1 0                                   *
!*                                                       *
!* Polynomial approximation of degree 3 (11 points)      *
!* Coefficients of polynomial:                           *
!*  A(0) =    -0.069930070                               *
!*  A(1) =     3.496503497                               *
!*  A(2) =    -3.496503497                               *                      
!*  A(3) =     0.000000000                               *
!*                                                       *
!* Approximated function:                                *
!*        X           Y                                  *
!*    0.000000   -0.069930                               *
!*    0.100000    0.244755                               *
!*    0.200000    0.489510                               *
!*    0.300000    0.664336                               *
!*    0.400000    0.769231                               *
!*    0.500000    0.804196                               *
!*    0.600000    0.769231                               *
!*    0.700000    0.664336                               *
!*    0.800000    0.489510                               *
!*    0.900000    0.244755                               *
!*    1.000000   -0.069930                               *
!*                                                       *
!*                    F90 Version By J-P Moreau, Paris.  *
!*                           (www.jpmoreau.fr)           *
!*********************************************************



! Inner variables

integer i,ij,j,k,n1,m1,m2
real*8  C(N,N)
real*8  A(N),B(N),Xc(N),Yx(N)
real*8  p,xx,s,yc


 n1=n
 m1=m+1; m2=m+2

  do k=1, m2
    Xc(k)=0.d0
    do i=1, n1
	  Xc(k) = Xc(k) + X(i)**k
    end do
  end do
  yc=0.d0
  do i=1, n1  
    yc = yc + Y(i)
  end do
  do k=1, m
    Yx(k)=0.d0
    do i=1, n1
	  Yx(k) =  Yx(k) + Y(i)*X(i)**k
    end do 
  end do
  do i=1, m1
	do j=1, m1
      ij=i+j-2
      if (i==1.and.j==1)  then
	    C(1,1) = n1
      else 
	    C(i,j)=Xc(ij)
      end if 
    end do
  end do
  B(1)=yc; 
  do i=2,m1
    B(i)=Yx(i-1)
  end do 
  do k=1, m
    do i=k+1, m1
      B(i) = B(i) - C(i,k)/C(k,k)*B(k)
      do j=k+1, m1
        C(i,j) = C(i,j) - C(i,k)/C(k,k)*C(k,j)
      end do
    end do
  end do
  A(m1)=B(m1)/C(m1,m1)
  do i=m, 1, -1
    s=0.d0
    do k=i+1, m1  
	  s = s + C(i,k)*A(k)
    end do 
    A(i) = (B(i)-s)/C(i,i)
  end do
  print *,' '
  write(*,40)  m, n
  print *,' Coefficients of polynomial:'
  do i=1, m1
    write(*,50)  i-1, A(i)
  end do
  print *,' ' 
  print *,' Approximated function:'
  print *,'       X           Y  '
  do i=1, n1
    xx=X(i); p=0.d0
    do k=1, m1
	  p = p*xx + A(m1+1-k)
    end do 
    write(*,60)  xx, p
    Yfit(i)=p
  end do
  print *,' '
  print *,' '


30 format('   X(',I1,'), Y(',I1,') = ')
31 format('   X(',I2,'), Y(',I2,') = ')

40 format('  Polynomial approximation of degree ',I2,' (',I2,' points)'/)
50 format('    A(',I4,') = ',F15.9)
60 format(2F12.6)

end subroutine least_squares_kutre


!-------------------------------------------------------------------------------------------


 subroutine StripSpaces(string)
  
  ! This subroutine removes blank spaces (' ') from the name of the files
  
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

  end subroutine

!-------------------------------------------------------------------------------------------


subroutine find_g ( N, polynomial_degree, Nintervals, h, fitting_coefs1, fitting_coefs2, g )

!     This program calculates the electron-phonon coupling (g) as half the distance between two curves: 
!     a1 + b1*x + c1*x^2 + d1*x^3 ... and a2 + b2*x + c2*x^2 + d2*x^3 ... (from "fitting_coefs1" and "fitting_coefs2").

integer, Intent(In) :: N    				 ! Size of the mathematical function to fit (number of points)
integer, Intent(In) :: polynomial_degree     ! Degree of the polynomial of the fit
integer, Intent(InOut) :: Nintervals         ! Number of intervals used in the calculation of g
real(8), Intent(In) :: h(N)  		   		 ! Displacements with respect to the equilibrium position used for the calculation of eigenvalues
real(8), Intent(In) :: fitting_coefs1(polynomial_degree+1) 	 ! Coefficients of the fitting
real(8), Intent(In) :: fitting_coefs2(polynomial_degree+1) 	 ! Coefficients of the fitting
real(8), Intent(Out) :: g					 ! coupling coefficient

! Local variables
  logical :: exist
  integer :: i, j, k 
  real(8) :: x1, x2, y1, y2, minx1, maxx1, minx2, maxx2, sqdistsofar, auxr1, auxr2, x1mindist, x2mindist
  real(8) :: a1,b1,c1,d1,a2,b2,c2,d2
!  real(8), allocatable :: X0vec(:,:,:), Uvec(:,:,:) ! These are the vectors with the relaxed positions (X0) and with the U-displacements for every vibrational mode

  
! Parameters  
  if ( ( Nintervals .lt. 100 ) .or. (Nintervals .gt. 99999999) ) Nintervals=4097
  minx1=h(1)
  minx2=minx1
  maxx1=h(N)
  maxx2=maxx1
  
 
   sqdistsofar=999999999999.d0
 
   do i=1,Nintervals
    
     x1=minx1+(maxx1-minx1)*dble(i-1)/(dble(Nintervals-1))
     
     y1=0.d0
     do k=1,polynomial_degree+1
       y1=y1+fitting_coefs1(k)*(x1**(k-1))
     end do

   
     do j=1,Nintervals
     
       x2=minx2+(maxx2-minx2)*dble(j-1)/(dble(Nintervals-1))
       y2=0.d0
       do k=1,polynomial_degree+1
         y2=y2+fitting_coefs2(k)*(x2**(k-1))
       end do

       
       auxr1=(x1-x2)**2 + (y1-y2)**2
       
       if (auxr1 .lt. sqdistsofar)  then
          sqdistsofar=auxr1
          x1mindist=x1
          x2mindist=x2 
       end if   
       
      ! write(*,*) i,j,sqrt(auxr1) 
     
     end do
   end do   
   
   g=sqrt(sqdistsofar)/2.d0
   

   write(*,*) '  The minimum distance is ',sqrt(sqdistsofar), ' at x1=',x1mindist, ', x2=',x2mindist
   write(*,*) '  Hence    g = ',g


   
   if   ( ( abs(( x1mindist-h(1))/h(1)) .lt. 0.0001 ) .or. ( abs((x1mindist-h(N))/h(N)) .lt. 0.0001 ) ) then
     write(*,*) ' WARNING: The h for the calculation of g is at the limit of the range given for h. '
     write(*,*) ' Include larger ranges of h in the eigenvalue files (eigenvalue1.dat, ...) for . '
     write(*,*) ' better accuracy.'; write(*,*) 
   end if

   if   ( ( abs(( x2mindist-h(1))/h(1)) .lt. 0.0001 ) .or. ( abs((x2mindist-h(N))/h(N)) .lt. 0.0001 ) ) then
     write(*,*) ' WARNING: The h for the calculation of g is at the limit of the range given for h. '
     write(*,*) ' Include larger ranges of h in the eigenvalue files (eigenvalue1.dat, ...) for . '
     write(*,*) ' better accuracy.'; write(*,*) 
   end if   


end subroutine find_g

!-----------------------------------------------------------------------------------------                         
! 
! 
! subroutine find_initial_guess ( N, h, eigv1, eigv2,initialguesseigv1, initialguesseigv2 )
! 
! ! This subroutine calculates the initial values for the iterative solution of the double anticrossing equation.
! ! E2 and E1 are assumed to be straight lines (with nonzero slope).
! 
! integer, Intent(In)  :: N    				 ! Size of the mathematical function to fit (number of points)
! real(8), Intent(In) :: h(N)  		   		 ! Vector with the displacements used for the calculation of the eigenvalues
! real(8), Intent(In) :: eigv1(N)  		   	 ! Vector with the first eigenvalue (e.g. HOMO-1a)
! real(8), Intent(In) :: eigv2(N)  		   	 ! Vector with the second eigenvalue (e.g. HOMO)
! real(8), Intent(Out) :: initialguesseigv1(N) ! Initial guess for the first eigenvalue (e.g. HOMO-1a)
! real(8), Intent(Out) :: initialguesseigv2(N) ! Initial guess for the second eigenvalue (e.g. HOMO)
! 
! 
!  ! Local variables
!  integer :: i,j
!  real(8) :: auxr1 
!  
!  i=(N+1)/2 
!  
!  ! E1
!  auxr1 =  ( eigv2(N) - eigv1(1) )/( h(N) - h(1) )
!  do j=1,N
!    initialguesseigv1(j) = eigv1(1) + auxr1 * ( h(j) - h(1) )
!  end do
!  
!  write(*,*)
!  
!  ! E2
!  auxr1 =  ( eigv1(N) - eigv2(1) )/( h(N) - h(1) )
!  do j=1,N
!    initialguesseigv2(j) = eigv2(1) + auxr1 * ( h(j) - h(1) )
!  end do
!  
! 
! end subroutine find_initial_guess 
!
!
!-----------------------------------------------------------------------------------------                         


subroutine solve_double_anticrossing_equation ( N, NhToSolveAnticrossingEquation, polynomial_degree_concavity, &
           & g, h, eigv1, eigv2)

! This subroutine calculates the initial values for the iterative solution of the double anticrossing equation.
! E1 and E2 are assumed to be straight lines (with nonzero slope).
!  
!    | E1-E	  g	  |
!	 | g    E2-E  | = 0
!
! where E1,E2 are the bare eigenvalues (unknowns for us), E are the dressed eigenvalues
! (E:=E+,E-; they are known from DFT), and g is the coupling (which is calculated
! as half the minimum distance between E2 and E1). IMPORTANT: Note that in our
! notation, the state whose concavity is sought for the calculation of frozen-phonon renormalization is E2
! (e.g. E2=HOMO), and E2 has an anticrossing with E1 (e.g. E1=HOMO-1).
!
! The equations are solved with the Newton-Raphson iterative method: $\bf{x}_{n+1} = \bf{x}_n - ( J^{-1}(\bf{x}_n) )·( \bf{f}(\bf{x}_n))$
! being J the Jacobian matrix: J_{11}(\bf{x}_n) := \partial f_1 (\bf{x}_n) / \partial x_1 ;  J_{12}(\bf{x}_n) := \partial f_1 (\bf{x}_n) / \partial x_2 ; 
!                              J_{NN}(\bf{x}_n) := \partial f_N (\bf{x}_n) / \partial x_N .
! The method of successive over-relaxation is much simpler to implement, but it presents serious convergence problems. 


integer, Intent(In)  :: N    				 ! Size of the mathematical function to fit (number of points)
integer, Intent(In) :: NhToSolveAnticrossingEquation  ! 
integer, Intent(In) :: polynomial_degree_concavity     ! Degree of the polynomial of the final fit to find concavity
real(8), Intent(In) :: g  			     	 ! Coupling parameter
real(8), Intent(In) :: h(N)  		   		 ! Vector with the displacements used for the calculation of the eigenvalues
real(8), Intent(In) :: eigv1(N)  		   	 ! Vector with the first eigenvalue (e.g. HOMO-1a); E+, dressed, from DFT
real(8), Intent(In) :: eigv2(N)  		   	 ! Vector with the second eigenvalue (e.g. HOMO); E0, dressed, from DFT
!real(8), Intent(Out) :: concavityeigv2       ! The concavity of the bare eigenvalue 2 (E2), which will be used in the frozen-phonon calculation.


 ! Local variables
 integer :: i,j,k,l,m,minhindex,maxhindex,Nmaxiter,poldeg, info
 real(8) :: auxr1, auxr2, E1, E2, E1old, E2old, Eplus, Eminus, tolerance
 real(8), allocatable :: E1vec(:,:),E2vec(:,:),X(:),Y(:) !,XX(:),YY(:),YYfit(:),coefs(:)



 minhindex=1+(N-NhToSolveAnticrossingEquation)/2
 maxhindex=minhindex+NhToSolveAnticrossingEquation-1
 
 write(*,*); write(*,'(A,G11.4,A,G11.4,A)') '   Now solving the double anticrossing equation ( between h = ', & 
         & h(minhindex), ' and h = ',h(maxhindex),')' ; write(*,*) 
 
 allocate( E1vec(maxhindex-minhindex+1,2), E2vec(maxhindex-minhindex+1,2)  )        

! call find_initial_guess ( N, h, eigv1, eigv2,  initialguesseigv1, initialguesseigv2 )

 Nmaxiter=500
 tolerance = 0.0000001d0
 E1vec(:,:)=9999999.d0;  E2vec(:,:)=9999999.d0


 do j=minhindex,maxhindex   
     
    Eplus = eigv1(j) ;         Eminus = eigv2(j) 
 
    i = 1+j-minhindex
    
    E1vec(i,1) =  h(j)
    auxr1 = (Eplus - Eminus)**2 - 4*g*g 
    if (auxr1 .ge. 0.d0) E1vec(i,2) =  ( Eplus + Eminus )/2.d0 + (sqrt( auxr1 ))/2.d0
    !write(*,*) 'j=',j,'i=',i, 'h=',E1vec(i,1), 'E1=',E1vec(i,2)  
    E2vec(i,1) =  h(j)
    auxr1=(Eplus - Eminus)**2 - 4*g*g
    if (auxr1 .ge. 0.d0) E2vec(i,2)  = ( Eplus + Eminus )/2.d0 - (sqrt( (Eplus - Eminus)**2 - 4*g*g )/2.d0)
 
 end do ! j


 
!  do i=1,maxhindex-minhindex+1
!    write(*,*) E2vec(i,1), E2vec(i,2)
!  end do
!  write(*,*)
 
 
 ! We now know the bare eigenvalue, stored in E2vec. We make a polynomial fit to find the concavity.
!  
!  ! First, we count the number of points
!  j=0
!  do i=1,maxhindex-minhindex+1
!    if ( E2vec(i,2) .lt. 9999998.0 ) j=j+1
!  end do
! 
!  ! Second, we build the functions to interpolate
!  allocate(X(j),Y(j))

 open(200,file='eigenvalue1-bare.dat',status='unknown')
 do i=1,maxhindex-minhindex+1
   if ( E1vec(i,2) .lt. 9999998.0 ) then
     write(200,*) E1vec(i,1),E1vec(i,2)   
   end if  
 end do
 close(200)

 open(201,file='eigenvalue2-bare.dat',status='unknown')
 do i=1,maxhindex-minhindex+1
   if ( E2vec(i,2) .lt. 9999998.0 ) then
     write(201,*) E2vec(i,1),E2vec(i,2)
   end if  
 end do
 close(201)
 
 write(*,*) '  The bare eigenvalues were calculated. Find the files to extract the concavity  &
    at eigenvalue1-bare.dat and eigenvalue2-bare.dat. '  
 write(*,*)
 write(*,*) " ** The calculations finished satisfactorily. **"
 write(*,*)
 write(*,*)
 
!  l=(j)/2
!  
! 
!  allocate(XX(l),YY(l),YYfit(l),coefs(l)) 
!  
!  j=l+1
!  do i=maxhindex-minhindex+1,maxhindex-minhindex+1-l+1,-1
!    j=j-1
!    XX(j)=X(i)
!    YY(j)=Y(i)
!  end do
!  
!  do i=1,l
!   write(*,*) XX(i),YY(i)
!  end do
! 
!  poldeg = min(l-1,polynomial_degree_concavity)
! 
!  call least_squares_Lapack ( l, XX, YY, poldeg, coefs, YYfit, auxr1 )
!  
!  write(*,*) ; write(*,'(A,G16.8)') '   RESULT: The concavity of the bare eigenvalue (E2) is: ',coefs(3) ; write(*,*) 
!  write(*,*) '  Please, check the eigenvalue2-bare.dat file, select the appropriate range of points'
!  write(*,*) '  and find the final result from the corresponding fit (coefficient of h-squared).' ; write(*,*) 

 ! deallocate(XX,YY,YYfit,coefs )
 


 deallocate(E1vec,E2vec) !X,Y
 
 end subroutine solve_double_anticrossing_equation 
 
!-------------------------------------------------------------------------------------------


