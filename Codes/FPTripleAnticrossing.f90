
Program ConcavityWithTripleAnticrossing

! This is the program to perform the calculations for correcting the distortions of electron-vibrational
! renormalization due to anticrossings of 3 states, as explained in Risueno, Han, Kumar, Bester PRB 2021.
! Author: Pablo Garcia-Risueño (garcia.risuenoATgmailDOTcom), 2018-2021.
! 
!
! The frozen-phonon method for renormalization of electronic eigenvalues [Capaz et al. PRL 2005] 
! presents one important drawback: even for small nuclear displacements (h), there exist vibrational
! modes where there can be severe inaccuracies due to anticrossing of electronic states. This phenomenon
! completely distorts the eigenvalue_vs_h curve. For example, the typical contribution of a mode to the 
! renormalization of HOMO and LUMO in diamondoids is a few meV. 
! However, for modes with anticrossing, this can be above 50 meV. Hence,
! a few modes with anticrossing are responsible for severe dips in accuracy of the total renormalizations.
! This program intends to solve this problem. To this end, we solve the equation of triple anticrossing.  
! Proceeding this way, we find the BARE_eigenvalues_vs_h, discounting the anticrossing effect of mixing of
! states. We do it by solving the following equation (eq. (9) of Risueno, Han, Kumar, Bester PRB 2021):
!  
!    | Delta + a·h -E	  g	          g3       |
!	 |  g                -E	          g        | = 0
!	 |  g3                g      Delta-a·h -E  |
!
!
! The program below has 3 parts: 
! 1) Initialization: Reading the file with input parameters (FPanticrossing.in) and the files containing the eigenvalues
!    as a function of h (eigenvalue1.dat, eigenvalue2.dat, eigenvalue3.dat)
! 2) Calculation of couplings: g1, g2 are calculated as half the distance between the curves of the dressed 
!    eigenvalues; to find them, we build such curves by polynomial interpolation using the data from DFT (from the
!    files which contain the eigenvalues). g:=(g1+g2)/2. g3 is found from the condition that the gap between the dressed eigenvalues at h=0 is zero.
! 3) Solution of the equation above. Since the unknowns are E1,E2,E3 we have an inverse eigenvalue problem.

! An example of the input file to run the code is the following:
! 
! # Input file for the program for the calculation of the
! # concavity of the bare eigenvalue of an eigenvalue subject
! # to double anticrossing
! 
!  quantity_to_correct = HOMO
!  Nh = 33
!  PolynomialDegreeg = 6
!  Nintervals = 4097
!
!
! Being eigenvalue1.dat (e.g. HOMO-1):   
! -16     -6.430498025
! -15     -6.4231426
! -14     -6.415793552
! -13     -6.40845435
! -12     -6.401128748
! -11     -6.393822139
! -10     -6.386539463
! -9      -6.379286659
! -8      -6.372070707
! -7      -6.364899931
! -6      -6.357784497
! -5      -6.350734742
! -4      -6.343764929
! -3      -6.336890261
! -2      -6.330129825
! -1      -6.323505329
! 0       -6.316975203
! 1       -6.310767603
! 2       -6.304717467
! 3       -6.298925489
! 4       -6.293428443
! 5       -6.288261937
! 6       -6.283456079
! 7       -6.27903268
! 8       -6.275002831
! 9       -6.271362776
! 10      -6.2680949
! 11      -6.265171832
! 12      -6.262556968
! 13      -6.260211297
! 14      -6.258094162
! 15      -6.256168043
! 16      -6.254399585
! 
! Being eigenvalue2.dat:  (IMPORTANT: eigenvalue2.dat MUST contain the eigenvalue to correct, e.g. the HOMO)
! -16     -6.153912529
! -15     -6.160406631
! -14     -6.166678797
! -13     -6.172690795
! -12     -6.178401177
! -11     -6.183766085
! -10     -6.188739865
! -9      -6.193281394
! -8      -6.197354539
! -7      -6.200934406
! -6      -6.204007754
! -5      -6.206573541
! -4      -6.208640367
! -3      -6.210222402
! -2      -6.21133788
! -1      -6.212000524
! 0       -6.212202682
! 1       -6.212000524
! 2       -6.21133788
! 3       -6.210222402
! 4       -6.208640367
! 5       -6.206573541
! 6       -6.204007765
! 7       -6.200934406
! 8       -6.197354539
! 9       -6.193281394
! 10      -6.188739865
! 11      -6.183766085
! 12      -6.178401177
! 13      -6.172690795
! 14      -6.166678797
! 15      -6.160406631
! 16      -6.153912529
! 
! and being eigenvalue3.dat (e.g. HOMO-2)
! -16     -6.254399585
! -15     -6.256168043
! -14     -6.258094162
! -13     -6.260211297
! -12     -6.262556968
! -11     -6.265171832
! -10     -6.2680949
! -9      -6.271362776
! -8      -6.275002831
! -7      -6.27903268
! -6      -6.283456058
! -5      -6.288261937
! -4      -6.293428443
! -3      -6.298925489
! -2      -6.304717467
! -1      -6.310767603
! 0       -6.317127656
! 1       -6.323505329
! 2       -6.330129825
! 3       -6.336890261
! 4       -6.343764929
! 5       -6.350734742
! 6       -6.357784505
! 7       -6.364899931
! 8       -6.372070707
! 9       -6.379286659
! 10      -6.386539463
! 11      -6.393822139
! 12      -6.401128748
! 13      -6.40845435
! 14      -6.415793552
! 15      -6.4231426
! 16      -6.430498025
! 

! For further details in the equations see the articles:
! * "Frozen-phonon method for state anticrossing situations and its application to zero-point motion effects in diamondoids",
! by Risueno, Han, Kumar, Bester, PRB 2021.
! * "Generation of entangled states and error 
! protection from adiabatic avoided level crossings"", by Nicole F. Bell, R. F. Sawyer, Raymond R. Volkas and Yvonne Y. Y. Wong.
! Note that in the latter article the coupling between HOMO-1 and HOMO-2 (i.e. g3, the entry [1,3] of the matrix) is neglected.


implicit none

integer :: N, Nintervals, NhToSolveAnticrossingEquation, polynomial_degree_g, polynomial_degree_concavity, i, j, k 
real(8) :: auxr1, auxr2, auxr3, Rsq, g1, g2, g3, concavityeigv2, check_values, hmindist1, hmindist2, slope_a, Delta
real(8), allocatable :: X(:), Y(:), Yfit(:), h(:), eigv1(:), eigv2(:), eigv3(:)
real(8), allocatable :: eigval_dressed_just_from_coupling(:), eigval_less_dressed_just_from_coupling(:)
real(8), allocatable :: fitting_coefs_eigv1(:), fitting_coefs_eigv2(:), fitting_coefs_eigv3(:)
character(4) :: quantity_to_correct 



call read_input_file ( N, polynomial_degree_g, polynomial_degree_concavity, & 
    NhToSolveAnticrossingEquation, Nintervals, g1, g2, quantity_to_correct )

allocate( X(N), Y(N), Yfit(N), h(N), eigv1(N), eigv2(N), eigv3(N), fitting_coefs_eigv1(polynomial_degree_g+1) )
allocate( eigval_dressed_just_from_coupling(N), eigval_less_dressed_just_from_coupling(N)   )
allocate( fitting_coefs_eigv2(polynomial_degree_g+1), fitting_coefs_eigv3(polynomial_degree_g+1) )


call read_eigenvalue_files ( N, h, eigv1, eigv2, eigv3 )



! CALCULATION OF THE COUPLING COEFFICIENTS (g1, g2)

write(*,*) '  **   Now finding the value of g, coupling between the eigenvalue to analyse **'
write(*,*) '  ** and the two other eigenvalues which participate in an anticrossing with it **'
write(*,*)
write(*,*) "    Fitting of a polynomial to eigenvalue_1:"
call least_squares_Lapack ( N, h, eigv1, polynomial_degree_g, fitting_coefs_eigv1, Yfit, Rsq )
write(*,*) ; write(*,*) "    Fitting of a polynomial to eigenvalue_2:"
call least_squares_Lapack ( N, h, eigv2, polynomial_degree_g, fitting_coefs_eigv2, Yfit, Rsq )
write(*,*) ; write(*,*) "    Fitting of a polynomial to eigenvalue_3:"
call least_squares_Lapack ( N, h, eigv3, polynomial_degree_g, fitting_coefs_eigv3, Yfit, Rsq )
write(*,*)

write(*,*)  "    The values of the coupling of the three states are:"; write(*,*)
write(*,*)  "      g_{21}:"
if (g1 .lt. -0.1d0) call find_g (N,  polynomial_degree_g, Nintervals, h, fitting_coefs_eigv1, fitting_coefs_eigv2, g1, hmindist1)
write(*,*); write(*,*)  "      g_{23}:"
if (g2 .lt. -0.1d0) call find_g ( N, polynomial_degree_g, Nintervals, h, fitting_coefs_eigv3, fitting_coefs_eigv2, g2, hmindist2 )
write(*,*)
check_values = abs(g1/g2)
if ( (check_values .lt. 0.09) .or. (check_values .gt. 1.1)  ) then
    write(*,*) '  WARNING: The values of g1 and g2 differ more than expected for an anticrossing of 3 states.'
    write(*,*) '  Check the shape of the curves to make sure that this is not an anticrossing of 2 states or'
    write(*,*) '  no anticrossing at all. Check also that the curves in eigenvalue1.dat and eigenvalue3.dat are smooth;'
    write(*,*) '  If they are not, exchange rows among them (e.g. all the rows with positive h) until they are smooth.'
    write(*,*)
end if



! CALCULATION OF THE COUPLING COEFFICIENT g3

call findg3(quantity_to_correct, N, h, eigv1, eigv2, eigv3, hmindist1, hmindist2, g1, g2, g3, slope_a, Delta)



! CALCULATION OF EIGENVALUES WITHOUT THE DISTORTION DUE TO ANTICROSSING EFFECT OF THREE STATES

call eigenvaluesalabell( quantity_to_correct, N, h, g1, g2, g3, slope_a, Delta, eigval_dressed_just_from_coupling)



! FINAL OUTPUT TO SCREEN

write(*,*); write(*,*) "  --> The ",quantity_to_correct," as a function of the displacement parameter (h) &
 WITHOUT the effect of anticrossings is what follows:"; write(*,*)

write(*,*) "             h,   ",quantity_to_correct," without distorting effect of anticrossings"
do k=1,N
   eigval_less_dressed_just_from_coupling(k) = eigv2(k) - eigval_dressed_just_from_coupling(k) 
   write(*,1001) "           ",h(k), ",",eigval_less_dressed_just_from_coupling(k)
   1001 format (A,f6.1,A,f12.8)
end do

write(*,*); write(*,*) "      Fit its central region to a parabola using e.g. Excel or OpenOffice. To this, &
 take a few values of h; make sure that you take a number of h's so that it looks still like a parabola. "
write(*,*) "      The concavity of this parabola (this is its coefficient of h^2) is what you have to use in eq. &
 (5) of Risueno et al 2021 to calculate the contribution of the "
write(*,*) "      analysed vibrational mode to the electron-phonon renormalization &
(discard the contribution of that mode which appears at the output of frozen-phonon.x, and use the new one instead). "; write(*,*)

! Don't do what follows; it includes many values of h; you must do it manually selecting so many h's that it is still a parabola.
!write(*,*) "    Eigenvalue_1:"
!call least_squares_Lapack ( N, h, eigval_less_dressed_just_from_coupling, polynomial_degree_g, fitting_coefs_eigv1, Yfit, Rsq )
!write(*,*) "Coefs",fitting_coefs_eigv1
!write(*,*) "Yfit",Yfit


write(*,*) ; write(*,*) '  ** The calculations finished satisfactorily. **'; write(*,*)

deallocate(X,Y, Yfit, h, eigv1, eigv2, eigv3, fitting_coefs_eigv1, fitting_coefs_eigv2, fitting_coefs_eigv3 )


end

! =============================================================


subroutine read_input_file (N,polynomial_degree_g,polynomial_degree_concavity, &
   NhToSolveAnticrossingEquation, Nintervals,g1,g2,quantity_to_correct )
                         
  ! System parameters
  integer, Intent(Out) :: N              				 ! The number of points of the eigenvalue files (number of h's)
  integer, Intent(Out) :: polynomial_degree_g 	 		 ! Chosen degree of the interpolation polynomials in the calculation of g1,g2
  integer, Intent(Out) :: polynomial_degree_concavity 	 ! Chosen degree of the interpolation polynomial at the final fit for the concavity
  integer, Intent(Out) :: NhToSolveAnticrossingEquation  ! The number of intervals for the division of h (ordinates) in the calculation of g
  integer, Intent(Out) :: Nintervals            		 ! The number of intervals for the division of h used for solving the triple anticrossing equation and to find the concavity of the bare eigenvalue
  real(8), Intent(Out) :: g1, g2  						 ! Coupling coefficients for the triple anticrossing. If there is indeed a double anticrossing (i.e. one of the eigenvalues is approximately a straight line), then
   														 ! you must make one of them be 0 (to know if g1=0 or g2=0, first run a calculation without making them 0, and so you will know which one corresponds to the 
   														 ! positive values of h and which one corresponds to the negative values of h).
  character(4), Intent(Out) ::quantity_to_correct      ! Either HOMO or LUMO; it informs which of both must be corrected to remove the effect of anticrossings

 
  ! Inner parameters
  integer :: i,j,k, stat, iostat
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character :: letter 
  logical :: exist 
  
   write(*,*) 
   write(*,*) '  ********************************************************************************' 
   write(*,*) '   Now calculating the eigenvalues without the distorting effects of anticrossings '
   write(*,*) '  ********************************************************************************' 
  
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
  g1=-1.d0; g2=-1.d0

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

     if ((auxstr1 .eq. 'g1').or.(auxstr1 .eq. 'G1') )   then
          Read(auxstr2, '(G18.8)' ) g1
     end if  

     if ((auxstr1 .eq. 'g2').or.(auxstr1 .eq. 'G2') )   then
          Read(auxstr2, '(G18.8)' ) g2
     end if  

     if ((auxstr1 .eq. 'quantity_to_correct').or.(auxstr1 .eq. 'Quantity_to_correct'))   then
         quantity_to_correct = trim(auxstr2)
         call StripSpaces (quantity_to_correct)
         if ( (quantity_to_correct .ne. "HOMO") .and. (quantity_to_correct .ne. "LUMO")) then
            write(*,*) "   ERROR: Invalid value for 'quantity_to_correct' in the input FPanticrossing.in file",&
                quantity_to_correct,"; Please, set it to 'HOMO' or 'LUMO'."
            call exit(1)    
         end if    
         
     end if

     
  end do ! k
   

1117 continue
   
   if ( (NhToSolveAnticrossingEquation) .eq. 0 ) NhToSolveAnticrossingEquation=N
   if ( mod(NhToSolveAnticrossingEquation,2) .eq. 0 ) NhToSolveAnticrossingEquation=NhToSolveAnticrossingEquation+1
   
   write(*,*)       
   write(*,*)  '  The number of h values in the input files is ',N
   !write(*,*)  '  The number of h values used to solve the threefold anticrossing equation is ', NhToSolveAnticrossingEquation  
   write(*,*)  '  The degree of the interpolating polynomials to get g is ',polynomial_degree_g
   !write(*,*)  '  The degree of the interpolating polynomial for the final calculation of the concavity is ', &
   !                     & polynomial_degree_concavity  
   if ( g1 .gt. -0.1d0 ) write(*,*) '  The first coupling g1 is ',g1 
   if ( g2 .gt. -0.1d0 ) write(*,*) '  The second coupling g2 is ',g2
   write(*,*)  
   
!     if ((auxstr1 .eq. 'Omega_units').or.(auxstr1 .eq. 'omega_units').or.(auxstr1 .eq. 'Omega_Units') .or. &
!         & (auxstr1 .eq. 'phonon_frequency_units').or.(auxstr1 .eq. 'Phonon_frequency_units').or. &
!         &   (auxstr1 .eq. 'phonon_Frequency_Units')  )   then
!          units_omegaorig = trim(auxstr2)
!          call StripSpaces (units_omegaorig)
!      end if         


end subroutine read_input_file 


!-----------------------------------------------------------------------------------------                         

subroutine read_eigenvalue_files ( N, h, eigv1, eigv2, eigv3 )

! This subroutine reads the files with the eigenvalues. Each has N rows with 2 numbers in each: the displacement (h) and the eigenvalue

integer, Intent(In)  :: N    				 ! Size of the mathematical function to fit (number of points)
real(8), Intent(Out) :: h(N)  		   		 ! Vector with the displacements used for the calculation of the eigenvalues
real(8), Intent(Out) :: eigv1(N)  		   	 ! Vector with the first eigenvalue (e.g. HOMO-1a)
real(8), Intent(Out) :: eigv2(N)  		   	 ! Vector with the second eigenvalue (e.g. HOMO)
real(8), Intent(Out) :: eigv3(N)  		   	 ! Vector with the third eigenvalue (e.g. HOMO-1b)


  ! Inner parameters
  integer :: i, j
  real(8) :: auxr1, auxr2

  h(:) = 0.00000000; eigv1(:)=0.00000000; eigv2(:)=0.00000000; eigv3(:)=0.00000000





  open(100,file='eigenvalue1.dat',status='unknown')
  do i=1,N
    read(100,*) j, auxr1
    h(i) = dble(j)
    eigv1(i) = auxr1
  end do
  close(100)

  open(101,file='eigenvalue2.dat',status='unknown')
  do i=1,N
    read(101,*) j, auxr1
    eigv2(i) = auxr1
  end do
  close(101)
  
  open(102,file='eigenvalue3.dat',status='unknown')
  do i=1,N
    read(102,*) j, auxr1
    eigv3(i) = auxr1
  end do
  close(102)    


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
  auxr1=0.00000000
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

Yavg=0.0000000000
do i=1,N
  Yavg=Yavg+Ytrue(i)
end do
Yavg=Yavg/(dble(N))

SStot=0.00000000000
do i=1,N
  SStot=SStot+ ( Ytrue(i) - Yavg )**2
end do

SSres=0.00000000000
do i=1,N
  SSres=SSres + ( Ytrue(i) - Yfit(i) )**2
end do


Rsq = 1.0000000000 - (SSres/SStot)
write(*,*) '    R-squared = ',Rsq

if ( Rsq .lt. 0.95 ) then
    write(*,*) "    WARNING: The value of R-squared is too low. Make sure that the data of the .dat file is appropriate"
    write(*,*) "    (e.g. that the functions are a smooth function of h)."
end if

end subroutine calculate_Rsquared


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


subroutine find_g ( N, polynomial_degree, Nintervals, h, fitting_coefs1, fitting_coefs2, g, hmindist   )

!     This program calculates the electron-phonon coupling (g) as half the distance between two curves: 
!     a1 + b1*x + c1*x^2 + d1*x^3 ... and a2 + b2*x + c2*x^2 + d2*x^3 ... (from "fitting_coefs1" and "fitting_coefs2").

integer, Intent(In) :: N    				 ! Size of the mathematical function to fit (number of points)
integer, Intent(In) :: polynomial_degree     ! Degree of the polynomial of the fit
integer, Intent(InOut) :: Nintervals         ! Number of intervals used in the calculation of g
real(8), Intent(In) :: h(N)  		   		 ! Displacements with respect to the equilibrium position used for the calculation of eigenvalues
real(8), Intent(In) :: fitting_coefs1(polynomial_degree+1) 	 ! Coefficients of the fitting
real(8), Intent(In) :: fitting_coefs2(polynomial_degree+1) 	 ! Coefficients of the fitting
real(8), Intent(Out) :: g					 ! coupling coefficient
real(8), Intent(Out) :: hmindist             ! h which gives the minimum distance between the coupled curves

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
   
   g=sqrt(sqdistsofar)/2.0000000000
   

   write(*,*) '      The minimum distance is ',sqrt(sqdistsofar), ' at x1=',x1mindist, ', x2=',x2mindist
   write(*,*) '      Hence    g = ',g
   
   hmindist = ( abs(x1mindist)+abs(x2mindist))/2.000000


   
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


subroutine findg3( quantity_to_correct , N, h_vec, eigv1, eigv2, eigv3, hmindist1,hmindist2, g1, g2, g3, slope_a, Delta)

! This subroutine finds the value of g3, i.e. the coupling between the two {eigenvalues whose
! eigenvectors do couple to the state_to_analyse} between themselves.

character(4), Intent(In) :: quantity_to_correct  ! Either HOMO or LUMO; it informs which quantity must be corrected
integer, Intent(In)  :: N    				 ! Size of the mathematical function to fit (number of points)
real(8), Intent(In) :: h_vec(N)  		     ! Vector with the displacements used for the calculation of the eigenvalues
real(8), Intent(In) :: eigv1(N)  		   	 ! Vector with the first eigenvalue (e.g. HOMO-1a)
real(8), Intent(In) :: eigv2(N)  		   	 ! Vector with the second eigenvalue (e.g. HOMO)
real(8), Intent(In) :: eigv3(N)  		   	 ! Vector with the third eigenvalue (e.g. HOMO-1b)
real(8), Intent(In)  :: hmindist1, hmindist2 ! Minimal distance between the coupled curves (e.g. HOMO and HOMO-1, HOMO and HOMO-2)
real(8), Intent(In)  :: g1                   ! Coupling parameter between HOMO and HOMO-1
real(8), Intent(In)  :: g2                   ! Coupling parameter between HOMO and HOMO-2
real(8), Intent(Out) :: g3                   ! Coupling parameter between HOMO-1 and HOMO-2
real(8), Intent(Out) :: slope_a
real(8), Intent(Out) :: Delta


integer :: i, j, k, maxh, info, iwork(1), sweeper, Nsweep 
real(8) :: auxr1, auxr2, f, g, g3_essaied, lowestgap = 99999
real(8) :: matrixdiag(3), matrixsubdiag(2), matrixin(3,3), matrixout(3,3), ab(3,3), work(6), eigenvalues(3)

g3=9999999 

 if ( eigv1(1) < eigv1(N) ) then
    a = eigv1(1) ! Value of HOMO-1 for lowest h
    b = eigv2(N) ! Value of HOMO for highest h; a & b give the slope (positive) of the diagonal asymptote of the HOMO-1
 else
    a = eigv1(N) ! Value of HOMO-1 for lowest h
    b = eigv2(1) ! Value of HOMO for highest h; a & b give the slope (negative) of the diagonal asymptote of the HOMO-1 
 end if   


 
 !write(*,*) "slope", abs(a-b)/(abs(h_vec(1))+abs(h_vec(N))) 


 ! slope_a is the quotient between the f (coupling energy) of eq. (7) of Bell's paper and our h (displacement size in problems of frozen-phonon); this is the slope of the straight lines which are HOMO-1 and HOMO-2 for high |h|; the factor 2 comes from the fact that in Bell's paper the slope is 2
 slope_a = (abs(a-b)/(abs(h_vec(1))+abs(h_vec(N))) )   !old slope_a = 0.014829 / 2.00000000000
 
 write(*,*) '      The slope ("a" in eq. (9) of Risueno et al. PRB) of the oblique asymptotes is ', slope_a," eV/(a.u.)"
 
 ! Delta is the value of f where the anticoupling takes place (value of the horizontal axis which gives minimal distance between curves of coupled eigenvalues)
 if (quantity_to_correct .eq. "HOMO") then  
    Delta = -((hmindist1+hmindist2)/2.000000000)*slope_a    !!old Delta = ((hmindist1+hmindist2)/2.d0)*slope_a     old Delta = ((4.9389648437+4.935546875)/2.d0)*slope_a
 else
    Delta =  ((hmindist1+hmindist2)/2.000000000)*slope_a  ! "Delta" is like the "Delta" of eq. (9) of Risueno et al. PRB.
 end if
 

 
 write(*,*)'      The gap between asymptotes ("Delta" in eq. (9) of Risueno et al. PRB) of the oblique asymptotes is ',Delta," eV"
 write(*,*)
 
 ! sqrt(2)·omega from Bell's paper: We assume it to be equal to our "g", this is, half the distance between coupled curves
 g =  ( ( g1 + g2 ) /2.0000000000 )  


 h=0.00000000000000
 Nsweep = 100000

 do sweeper=1,Nsweep
 
     g3_essaied = (10.0000000*g/Nsweep)*dble(-Nsweep/2.000 + sweeper)
 
     ah = real(h)*slope_a
     
	 ! We define the matrix:
	 matrixdiag(1) =  Delta + ah
	 matrixdiag(2) =   0.0000000000
	 matrixdiag(3) =  Delta -ah 
 
 
	 matrixsubdiag(1) = g 
	 matrixsubdiag(2) = g 

     write(103,*) h, matrixdiag(1), matrixdiag(2), matrixdiag(3)
     write(104,*) ah, matrixdiag(1), matrixdiag(2), matrixdiag(3)
  
! 	 call dsteqr(  'N',       		    &  ! To compute just eigenvalues (not eigenvectors)
! 					3,				    &  ! Order of the matrix
! 					matrixdiag,			&  ! IN: Diagonal entries of the matrix; OUT: Solution (eigenvalues, ordered)
! 					matrixsubdiag,		&  ! Subiagonal entries of the matrix				
! 					matrixout,			&  ! Output matrix (useless if eigenvectors are not computed)
! 					3,				    &  ! Leading dimension of matrixout
! 					work,			    &  ! Auxiliar array (useless if eigenvectors are not computed)
!					info )				   ! Output which informs if the calculations were properly performed			

    matrixin(:,:)=0.00000000000000; 
    matrixin(1,1)=matrixdiag(1) ; matrixin(2,2)=matrixdiag(2) ; matrixin(3,3)=matrixdiag(3) 
    matrixin(1,2)=matrixsubdiag(1) ;  matrixin(2,3)=matrixsubdiag(2) ; 
    matrixin(2,1)=matrixsubdiag(1) ;  matrixin(3,2)=matrixsubdiag(2) ; 
    matrixin(1,3)= g3_essaied 
    matrixout(:,:)=matrixin(:,:)
    
    
    !ab( 3+i-j,j) = a(i,j)
    ab(3,1)=matrixin(1,1)
    ab(2,2)=matrixin(1,2)
    ab(1,3)=matrixin(1,3)
    ab(3,2)=matrixin(2,2)
    ab(2,3)=matrixin(2,3)
    ab(3,3)=matrixin(3,3)


    call dsbevd	(	'N',       		    &  ! To compute just eigenvalues (not eigenvectors)
					'U',       		    &  ! The input matrix is stored in U (upper) form
					 3,				    &  ! Order of the matrix
					 2,				    &  ! Number of superdiagonals
					 ab,		        &  ! Matrix input cuyos eigvals se calculan; OJO: en notacion rara: if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
					 3,				    &  ! Leading dimension of the matrix
					 eigenvalues,       &  ! Output eigenvalues (in order, from lowest to highest)
					 matrixout,			&  ! Output matrix (useless if eigenvectors are not computed)
					 3,				    &  ! Leading dimension of matrixout
					 work,			    &  ! Auxiliar array  
					 6,				    &  ! Dimension of vector work (2N=6)
					 iwork,			    &  ! Dimension of vector work (2N=6)
					 1,				    &  ! Dimension of iwork
				 	 info )				   ! Output which informs if the calculations were properly performed		
	
	 if (abs(eigenvalues(1) - eigenvalues(2)) .lt. lowestgap) then
	     lowestgap = abs(eigenvalues(1) - eigenvalues(2))
	     g3 = g3_essaied
     end if


 end do  
 
 if ( g3 .eq. 9999999) then
     write(*,*) "   ERROR: Unable to find an appropriate value for g3. Please, check your data."; write(*,*)
     call exit(1)
 end if
 
 write(*,*) "      g_{13}:"
 write(*,*)  "               g =    ",g3
 write(*,*)  "      This value of g_{13} makes the gap between the two eigenavalues which couple to the eigenv. &
  to analyse (at h=0) be",lowestgap," eV."; write(*,*)
 
 

end subroutine findg3

! =============================================================



subroutine eigenvaluesalabell( quantity_to_correct, N, h_vec, g1, g2, g3, slope_a, Delta, eigval_dressed_just_from_coupling)

! This subroutine finds the value of g3, i.e. the coupling between the two {eigenvalues whose
! eigenvectors do couple to the state_to_analyse} between themselves.

character(4), Intent(In) :: quantity_to_correct  ! Either HOMO or LUMO; it informs which quantity must be corrected
integer, Intent(In) :: N    				 ! Size of the mathematical function to fit (number of points)
real(8), Intent(In) :: h_vec(N)  		     ! Vector with the displacements used for the calculation of the eigenvalues
real(8), Intent(In) :: g1                   ! Coupling parameter between HOMO and HOMO-1
real(8), Intent(In) :: g2                   ! Coupling parameter between HOMO and HOMO-2
real(8), Intent(In) :: g3                   ! Coupling parameter between HOMO-1 and HOMO-2
real(8), Intent(In) :: slope_a
real(8), Intent(In) :: Delta
real(8), Intent(Out):: eigval_dressed_just_from_coupling(N) ! This is the eigenvalue_dressed_just_by_anticrossing_effects as a function of h
                                                            ! We will subtract it to the eigenvalue from DFT to find the actual concavity to use in the electron-phonon calculation.

! This program reproduces the scheme of bare and dressed eigenvalues shown by eq. (7) of 
! Bell et al., "Generation of entangled states and error protection from adiabatic avoided level crossings" (PRA),
! which we consider the standard method of triple anticrossing equations.
!
! The parameters must be written at the beginning of the program. 

 

integer :: h, i, j, k, maxh, info, iwork(1)
real(8) :: auxr1, auxr2, f,  g 
real(8) :: matrixdiag(3), matrixsubdiag(2), matrixin(3,3), matrixout(3,3), ab(3,3), work(6), eigenvalues(3)



 open(101,file='eigenvalue-dressed_vs_h.dat',status='unknown') ! These eigenvalues contain the dressing due to the couplings (not to other effects). We will later discount this dressing from the observed dressing, to find the eigenvalues_without_anticrossing_effects -vs- h.
 open(102,file='eigenvalue-dressed_vs_f.dat',status='unknown') !
 open(103,file='eigenvalue-bare_vs_h.dat',status='unknown')    ! These are fully straight lines, asymptotes; they correspond to the states <<bare with respect to the anticrossing effect>>
 open(104,file='eigenvalue-bare_vs_f.dat',status='unknown')    ! These are fully straight lines, asymptotes; they correspond to the states <<bare with respect to the anticrossing effect>>

 do i=1,N
 
     
     h = h_vec(i)
     ah = real(h)*slope_a
     
	 ! We define the matrix:
	 matrixdiag(1) =  Delta + ah
	 matrixdiag(2) =      0.00000000
	 matrixdiag(3) =  Delta - ah
 
     g =  ( ( g1 + g2) /2.d0 )  ! sqrt(2)·omega from Bell's paper: We assume it to be equal to our "g", this is, half the distance between coupled curves
	 matrixsubdiag(1) = g 
	 matrixsubdiag(2) = g 

     write(103,*) h, matrixdiag(1), matrixdiag(2), matrixdiag(3)
     write(104,*) ah, matrixdiag(1), matrixdiag(2), matrixdiag(3)

    matrixin(:,:)=0.000000000; 
    matrixin(1,1)=matrixdiag(1) ; matrixin(2,2)=matrixdiag(2) ; matrixin(3,3)=matrixdiag(3) 
    matrixin(1,2)=matrixsubdiag(1) ;  matrixin(2,3)=matrixsubdiag(2) ; 
    matrixin(2,1)=matrixsubdiag(1) ;  matrixin(3,2)=matrixsubdiag(2) ; 
    matrixin(1,3)=  g3 !0.0036066 ! 2.d0*0.038050477448703823*0.038048050981807080 / (-6.276387301+6.384030242 )  ! g3: condic de que E+=E1: g3=2*g1*g2/(E2-E-)
    matrixout(:,:)=matrixin(:,:)
    
    
    !ab( 3+i-j,j) = a(i,j)
    ab(3,1)=matrixin(1,1)
    ab(2,2)=matrixin(1,2)
    ab(1,3)=matrixin(1,3)
    ab(3,2)=matrixin(2,2)
    ab(2,3)=matrixin(2,3)
    ab(3,3)=matrixin(3,3)


    call dsbevd	(	'N',       		    &  ! To compute just eigenvalues (not eigenvectors)
					'U',       		    &  ! The input matrix is stored in U (upper) form
					 3,				    &  ! Order of the matrix
					 2,				    &  ! Number of superdiagonals
					 ab,		        &  ! Matrix input cuyos eigvals se calculan; OJO: en notacion rara: if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
					 3,				    &  ! Leading dimension of the matrix
					 eigenvalues,       &  ! Output eigenvalues (in order)
					 matrixout,			&  ! Output matrix (useless if eigenvectors are not computed)
					 3,				    &  ! Leading dimension of matrixout
					 work,			    &  ! Auxiliar array  
					 6,				    &  ! Dimension of vector work (2N=6)
					 iwork,			    &  ! Dimension of vector work (2N=6)
					 1,				    &  ! Dimension of iwork
				 	 info )				   ! Output which informs if the calculations were properly performed		
	

     !write(*,*) h,")", eigenvalues(1), eigenvalues(2), eigenvalues(3)
     
     if (quantity_to_correct .eq. "HOMO") then ! HOMO; the dressed eigenvalue is the highest of all three
        eigval_dressed_just_from_coupling(i) =  eigenvalues(3)
     else                                     ! LUMO; the dressed eigenvalue is the lowest of all three
        eigval_dressed_just_from_coupling(i) =  eigenvalues(1)
     end if   
        
     write(101,*) h, eigenvalues(1), eigenvalues(2), eigenvalues(3)
     write(102,*) f, eigenvalues(1), eigenvalues(2), eigenvalues(3)


 end do ! h=-maxh,maxh
 
 close(101)
 close(102)

end subroutine eigenvaluesalabell

! =============================================================







!-----------------------------------------------------------------------------------------                         
