program GenerateDataToPlotCrossovers

! This program reads data from the outputs for frozen-phonon calculations and generates text files
! (as many ones as vibration modes) with the eigenvalues as a function of the displacement parameter.

 integer :: i, j, k, Nmodes, Ncols,  sizeepsilon,number_of_atoms, number_of_species, dimnu, band_occupation
 integer :: NumberofDPs, Numberofmodestoanalyse, readfromxmlfile 
 real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, auxr10
 real(8) :: auxr11, auxr12, auxr13, auxr14, auxr15, auxr16, auxr17, auxr18, auxr19, auxr20
 real(8) :: auxr21, auxr22, auxr23, auxr24, auxr25, auxr26, auxr27, auxr28, auxr29, auxr30
 real(8), allocatable :: matrix(:,:), matrixpm(:,:),epsilonorig(:),epsilon(:),mass(:,:)
 real(8), allocatable :: DPs(:),displacedeigvals(:,:)
 integer, allocatable :: modes_to_analyse(:) 
 character(100) :: auxstrA
 character(150) :: elec_eigval_file_name, xml_eigval_file_name!, mass_file_name
 character(120) :: auxstrC
 character(40) :: auxstrB, units_epsilonorig
 logical :: exist


  call count_Number_of_DPs (NumberofDPs )
  call count_Number_of_modes ( Numberofmodestoanalyse )

    
  allocate(DPs(NumberofDPs), modes_to_analyse (Numberofmodestoanalyse) ) 
  call read_input_file (NumberofDPs, Numberofmodestoanalyse, number_of_atoms, number_of_species, &
        & dimnu, band_occupation,  elec_eigval_file_name, DPs, modes_to_analyse, readfromxmlfile)   
 
        
  call read_QE_occ_unocc_numbers (elec_eigval_file_name, band_occupation, Ne, sizeepsilon)  
 
  allocate(epsilonorig(sizeepsilon),epsilon(sizeepsilon),mass(number_of_atoms,2))
  allocate(displacedeigvals(-NumberofDPs:NumberofDPs,sizeepsilon))
  displacedeigvals(:,:)=0.d0
  
  inquire(file='ToPlot', exist=exist)      
  if (.not. exist)  call system('mkdir ToPlot/')

  ! We read the file with unperturbed (not displaced) eigenvalues 
  xml_eigval_file_name='eigenval.xml'
  call read_QE_output (readfromxmlfile,sizeepsilon, Ne, number_of_atoms, number_of_species, band_occupation, &
                      & elec_eigval_file_name, xml_eigval_file_name, epsilon, mass) 
  do i=1,sizeepsilon
    displacedeigvals(0,i) = epsilon(i)
  end do
  

  do nu=7,dimnu
  
    flag_do_calculation = 0
    do m=1,Numberofmodestoanalyse
        if (nu == modes_to_analyse(m))    flag_do_calculation = flag_do_calculation + 1
    end do
    
  
    if (flag_do_calculation > 0) then
       

        do i=1,NumberofDPs
     
         ! NEGATIVE DISPLACEMENT
           epsilon(:)=0.d0
           ! We write the name of the file to read:
           if ( DPs(i) .ge. 10 ) then
              write(elec_eigval_file_name,*) 'DP',nint(DPs(i)),'/displaced-/out_scf-mode',nu,'.out'
              write(xml_eigval_file_name,*) 'DP',nint(DPs(i)),'/displaced-/eigenval-mode',nu,'.xml'
           else   
              write(elec_eigval_file_name,*) 'DP0',nint(DPs(i)),'/displaced-/out_scf-mode',nu,'.out'
              write(xml_eigval_file_name,*) 'DP0',nint(DPs(i)),'/displaced-/eigenval-mode',nu,'.xml'
           end if   
           call StripSpaces(elec_eigval_file_name);call StripSpaces(xml_eigval_file_name)

    	   ! We read the file with perturbed (displaced) eigenvalues 
    	   call read_QE_output (readfromxmlfile, sizeepsilon, Ne, number_of_atoms, number_of_species, &
    						   & band_occupation, elec_eigval_file_name, xml_eigval_file_name, epsilon, mass) 
    	   do j=1,sizeepsilon
    		 displacedeigvals(-i,j) = epsilon(j)
	    	 !write(*,*) displacedeigvals(-i,j)*27.211385056d0
	       end do       


         ! POSITIVE DISPLACEMENT
           epsilon(:)=0.d0
           ! We write the name of the file to read:
           if ( DPs(i) .ge. 10 ) then
              write(elec_eigval_file_name,*) 'DP',nint(DPs(i)),'/displaced+/out_scf-mode',nu,'.out'
              write(xml_eigval_file_name,*) 'DP',nint(DPs(i)),'/displaced+/eigenval-mode',nu,'.xml'
           else   
              write(elec_eigval_file_name,*) 'DP0',nint(DPs(i)),'/displaced+/out_scf-mode',nu,'.out'
              write(xml_eigval_file_name,*) 'DP0',nint(DPs(i)),'/displaced+/eigenval-mode',nu,'.xml'
           end if   
           call StripSpaces(elec_eigval_file_name);call StripSpaces(xml_eigval_file_name)

    	   ! We read the file with perturbed (displaced) eigenvalues 
    	   call read_QE_output (readfromxmlfile, sizeepsilon, Ne, number_of_atoms, number_of_species, band_occupation, &
    						   & elec_eigval_file_name, xml_eigval_file_name, epsilon, mass) 
    	   do j=1,sizeepsilon
    		 displacedeigvals(i,j) = epsilon(j)
    	   end do  
 
        end do !  do i=1,NumberofDPs
        
 
        write(auxstrA,*) 'ToPlot/mode',nu,'/'; call StripSpaces(auxstrA) 
        inquire(file=auxstrA, exist=exist)       
    	if (.not. exist) then
    	  write(auxstrC,*) 'mkdir  ',auxstrA 
    	  call system(auxstrC)
    	end if  
    
        do j=1,sizeepsilon
    	  write(auxstrB,*) 'state',j,'.out'; call StripSpaces(auxstrB)
    	  write(elec_eigval_file_name,*) auxstrA//auxstrB  ; call StripSpaces(elec_eigval_file_name)
    	  !write(*,*) trim(elec_eigval_file_name)
    	  open(111,file=elec_eigval_file_name,status='unknown')
    	  write(111,*) '# Displacement Param.        eigenvalue (eV)'
    	  do i=-NumberofDPs,NumberofDPs
    		if (i .lt. 0) then
    		   write(111,*) -DPs(abs(i)), '  ', displacedeigvals(i,j)*27.211385056
    		else if (i .gt. 0) then
    		   write(111,*) DPs(i), '  ', displacedeigvals(i,j) *27.211385056
    		else if (i .eq. 0) then
    		   write(111,*) 0.0000d0, '  ', displacedeigvals(i,j) *27.211385056
    		end if 
    	  end do ! j=1,sizeepsilon
    	  close(111)
        end do 	 
        
   end if !(flag_do_calculation == true)       
    
 end do ! nu
 
 
 write(*,*) " ** The calculations finished satisfactorily. Please, find the resulta at the folder called << ToPlot/ >> ** "
 write(*,*)
 write(*,*)

end program GenerateDataToPlotCrossovers


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

                        
 subroutine read_QE_output (readfromxmlfile, dimband, Ne, number_of_atoms,  number_of_species, band_occupation, &
           &  elec_eigval_file_name, xml_eigval_file_name, epsilon, mass)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives electronic eigenvalues and atomic masses.
 
 ! Comment: Ponce2014 considers hbar = m_e = e = 1 (atomic units). This means that, since we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 
 

   ! System parameters
  integer, Intent(In) :: readfromxmlfile                ! If nonzero, read from eigenval.xml files    
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  integer, Intent(In) :: Ne                             ! The number of electrons of the system (i.e. the number of occupied states)
  integer, Intent(In) :: number_of_atoms                ! Total number of atoms of the system
  integer, Intent(In) :: number_of_species              ! Total number of different atomic species of the system
  integer, Intent(In) :: band_occupation                ! Number of electrons (1 or 2) per band in the file output of QE
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  character(100), Intent(In) :: xml_eigval_file_name    ! Name of the  xml file where the electronic eigenvalues are stored
  real(8), Intent(Out) :: epsilon(dimband)				! Electronic eigenvalues
  real(8), Intent(Out) :: mass(number_of_atoms,2)       ! Masses of atoms
  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character(40) :: units_epsilonorig 
  character(12), allocatable :: species_name(:)
  real(8), allocatable :: species_mass(:), epsilonorig(:)
  character(12) :: auxstr4, auxstr5
  character(18) :: auxstr6, auxstr7
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, HOMO, LUMO, fermi_level  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																									!
!                               1st BLOCK: READING EIGENVALUES										!
!																									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(epsilonorig(dimband))
epsilonorig(:)=0.00000000000000000
  
if ( readfromxmlfile .eq. 0) then


  open(333,file=elec_eigval_file_name,status='unknown')
  !write(*,*) '   Reading ' , trim(elec_eigval_file_name)   

  do k=1,20000

      read(333,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 3456

      auxstr1=''; auxstr2='' 
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'Endofself-consistentcalculation') then  
          go to 2235
        end if
      end do ! i
      
   end do !k   
   
 2235 continue    
 
  read(333,'(A)',IOSTAT=stat) inputline; read(333,'(A)',IOSTAT=stat) inputline 
  auxstr1=''
  do i=1,60
     letter = inputline(i:i)
     auxstr1 = trim(auxstr1)//letter 
     auxstr2=''
     j=LNBLNK(auxstr1)   ! Last index of the chain that is not a blank space
     if ( j.ge.5) then
       do k=1,5
          letter=inputline(i-6+k:i-6+k)
          auxstr2 = trim(auxstr2)//letter
       end do
     end if
     
     if ( auxstr2 .eq. 'bands') then 
          units_epsilonorig=''
          do j=i+1,i+7
              letter=inputline(j:j)
              if ( (letter .ne. '(') .and.(letter .ne. ')') .and.(letter .ne. ',') .and. &
                & (letter .ne. '.') .and.(letter .ne. ';') .and.(letter .ne. ':')  ) then
                 units_epsilonorig = trim(units_epsilonorig)//letter
              end if
          end do
          if ( ( trim(units_epsilonorig) .eq. 'CM**-1' ) .or. ( trim(units_epsilonorig) .eq. 'Cm**-1' ) ) units_epsilonorig='cm**-1' 
          if ( ( trim(units_epsilonorig) .eq. 'mev' )    .or. ( trim(units_epsilonorig) .eq. 'MeV' ) )    units_epsilonorig='meV'
          if ( ( trim(units_epsilonorig) .eq. 'ev' )    .or. ( trim(units_epsilonorig) .eq. 'EV' ) )    units_epsilonorig='eV'
          if ( ( trim(units_epsilonorig) .eq. 'A.u.' ) .or. ( trim(units_epsilonorig) .eq. 'A.u.' ) .or. &
           &  ( trim(units_epsilonorig) .eq. 'Au.' )  .or. ( trim(units_epsilonorig) .eq. 'AU' ) .or.   &
            & ( trim(units_epsilonorig) .eq. 'hartree' ) .or. ( trim(units_epsilonorig) .eq. 'HARTREE' )  &
            & .or.  ( trim(units_epsilonorig) .eq. 'Hartree' ) )                                          units_epsilonorig='a.u.'

     end if
  end do ! i        



  read(333,'(A)',IOSTAT=stat) inputline
        
  auxstr2='' 
  i=1             ! i is the band counter
  do k=1,20000 
    read(333,'(A)',IOSTAT=stat) inputline  

    do j=1,100
     
     if ((stat .ne. 0)  .or. (i.gt.dimband)) go to 3456
     letter = inputline(j:j)
     if (letter .eq. 'h') go to 3456
     auxstr2 = trim(auxstr2)//letter
     if ( (auxstr2 .ne. '') .and. (letter .eq. ' ') ) then
        ! write(*,*) ' auxstr2 = ', auxstr2
        read(auxstr2,*) auxr1
        
!        if (band_occupation .eq. 2) then 
           !epsilonorig(i) = auxr1; epsilonorig(i+1) = auxr1;  i = i+2   
!        else 
           epsilonorig(i) = auxr1
           i = i+1         
!        end if  
        
        auxstr2=''
     end if
   end do ! j  
  end do   ! k
     
        
3456 continue

  close(333)
  
  
else ! we read the eigenvalues from eigenval.xml files


      !write(*,*) 'now reading',trim(xml_eigval_file_name)


      units_epsilonorig='a.u.'
      open(333,file=xml_eigval_file_name,status='unknown')
      

	  do k=1,20000

		  read(333,'(A)',IOSTAT=stat) inputline

		  if (stat .ne. 0) go to 6456

		  auxstr1=''; auxstr2='' 
  
		  do i=1,90
			letter = inputline(i:i)
			auxstr1 = trim(auxstr1)//letter 
			if ( auxstr1 .eq. '<EIGENVALUES') then  
			  go to 6235
			end if
		  end do ! i
	  
		end do !k   
   
		 6235 continue    
 

         i=1
         do j=1,dimband
          !  read(333,'(A)',IOSTAT=stat) inputline
          ! write(*,*) inputline
		  !  read(333,'(G30.23)') epsilonorig(j) !epsilon_output(j)
		   read(333,*) auxr1
		   !write(*,*) auxr1 !xxxxx
           epsilonorig(i) = auxr1
           i = i+1         
         end do
	 
	    6456 continue
	    
	    close(333)  

  
end if !if ( readfromxmlfiles .eq. 0) then  



  ! Now we make the 0 of energy in the Fermi level; we establish the Fermi level in the medium point of the band gap
  if ( mod(Ne,2) .eq. 1) then
     HOMO = epsilonorig(Ne);  LUMO = HOMO
  else
     HOMO = epsilonorig(Ne/2);  LUMO = epsilonorig(Ne/2+1)
  end if
  fermi_level = (HOMO+LUMO)/2.000000000000000000
  do i=1,dimband
    epsilonorig(i) = epsilonorig(i) !UUUUU- fermi_level
    !write(*,*) "bbb",i,epsilonorig(i) !xxxxxx
  end do 
  
  
  
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																									!
!                               2nd BLOCK: READING ATOMIC MASSES									!
!																									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(species_name(number_of_species),species_mass(number_of_species))
  
  mass(:,2) = -1.000000000000000
  
 open(333,file=elec_eigval_file_name,status='unknown')

  do k=1,20000

      read(333,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 7456

      auxstr1=''; auxstr2='' 
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'atomicspeciesvalencemass') then  
          go to 2435
        end if
      end do ! i
      
   end do !k   
   
 2435 continue  
  
  
  do i=1,number_of_species
    read(333,*) species_name(i), auxr1, species_mass(i), auxstr1
    call StripSpaces(species_name(i))
  end do

         
   do k=1,20000

      read(333,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 7456

      auxstr1=''; auxstr2='' 
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'siten.atompositions') then  
          go to 1435
        end if
      end do ! i
      
   end do !k   
   
 1435 continue   

  do i=1,number_of_atoms
    read(333,*) j, auxstr4
    call StripSpaces(auxstr4)
    
    do k=1,number_of_species
       if ( auxstr4 .eq. species_name(k)) then
         mass(i,1) = dble(i)
         mass(i,2) = species_mass(k) * 1822.888486192 
       end if
    end do
    
  end do
  
  
  do i=1,number_of_atoms
     if (mass(i,2) .lt. 0.00000000) then
        write(*,*) '  ERROR: Masses are incorrect or they were not properly read.'
        stop
     end if
  end do
  
  
  close(333)
  
 
 7456 continue
 
 
   call epsilon_to_atomic_units (units_epsilonorig,dimband,epsilonorig,epsilon)
  
  deallocate(species_name,epsilonorig)


 end subroutine read_QE_output
 
!-------------------------------------------------------------------------------------------------------

subroutine  count_Number_of_DPs ( NumberofDPs  )

  ! System parameters
  integer, Intent(Out) :: NumberofDPs              ! The number of displacement parameters to be used
 
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxstr4, auxstr5
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l,iostat,unitt !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2

  NumberofDPs = 0
 

  ! Reading parameters from file
  
    
   open(345,file='GenerateDataToPlotCrossovers.in',status='unknown')


   do k=1,50

      read(345,'(A)',IOSTAT=stat) inputline
      !write(*,*) inputline
      if (stat .ne. 0) go to 3317

      auxstr1=''; auxstr2=''; auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '#') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 1235
        end if
      end do ! i
 
      1235 continue
  
          
   
    if ((trim(auxstr1) .eq. 'Displacement_parameters').or.(trim(auxstr1) .eq. 'displacement_parameters').or. &
             & (trim(auxstr1) .eq. 'Displacement_Parameters') .or. ( trim(auxstr1) .eq. 'DisplacementParameters' ))  then
       
      do l=1,100000 
       
		  read(345,'(A)',IOSTAT=stat) inputline
		  auxstr1=''; auxstr2=''; auxstr3=''
   
		  do i=1,20
		   letter = inputline(i:i)
		   if (letter .ne. '#') then
			 auxstr1 = trim(auxstr1)//letter 
		   else 
		     close(345)
			 return
		   end if
		  end do ! i
 
          NumberofDPs = NumberofDPs + 1

      
      end do ! l
  
    end if  
    
  end do ! k
   
    

3317 continue

   close(345)
   


   write(*,*) 'ERROR: Please, specify the displacement parameters used with a block similar to the one '
   write(*,*) ' below in GenerateDataToPlotCrossovers.in'

   write(*,*) ' Displacement_parameters######'
   write(*,*) ' 01'
   write(*,*) ' 04'
   write(*,*) ' 06'
   write(*,*) ' 08'
   write(*,*) ' 12'
   write(*,*) ' 16'
   write(*,*) ' ############################'




end subroutine  count_number_of_DPs 

!-------------------------------------------------------------------------------------------------------


subroutine  count_Number_of_modes ( Numberofmodestoanalyse )

  ! System parameters
  integer, Intent(Out) :: Numberofmodestoanalyse   ! The number of vibrational modes (normal modes) to analyse               
  
    
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxstr4, auxstr5
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l,iostat,unitt !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2

  Numberofmodestoanalyse = 0

  ! Reading parameters from file
  
    
   open(346,file='GenerateDataToPlotCrossovers.in',status='unknown')


   do k=1,50

      read(346,'(A)',IOSTAT=stat) inputline
      !write(*,*) inputline
      if (stat .ne. 0) go to 5517

      auxstr1=''; auxstr2=''; auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '#') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 5235
        end if
      end do ! i
 
      5235 continue
  
          
   
      if ((trim(auxstr1) .eq. 'Modes_to_analyse').or.(trim(auxstr1) .eq. 'Modes_to_Analyse').or. &
             & (trim(auxstr1) .eq. 'modes_to_analyse') )  then
       
      do l=1,100000 
       
		  read(346,'(A)',IOSTAT=stat) inputline
		  auxstr1=''; auxstr2=''; auxstr3=''
   
		  do i=1,20
		   letter = inputline(i:i)
		   if (letter .ne. '#') then
			 auxstr1 = trim(auxstr1)//letter 
		   else 
		     close(346)
			 return
		   end if
		  end do ! i
 
          Numberofmodestoanalyse = Numberofmodestoanalyse + 1

      
      end do ! l
  
    end if  
    
  end do ! k
   
    

5517 continue

   close(346)
   


   write(*,*) 'ERROR: Please, specify the displacement parameters used with a block similar to the one '
   write(*,*) ' below in GenerateDataToPlotCrossovers.in'

   write(*,*) ' Modes_to_analyse######'
   write(*,*) ' 44'
   write(*,*) ' 45'
   write(*,*) ' 55'
   write(*,*) ' 56'
   write(*,*) ' 66'
   write(*,*) ' 67'
   write(*,*) ' ############################'




end subroutine  count_number_of_modes

!-------------------------------------------------------------------------------------------------------


subroutine read_input_file (NumberofDPs,Numberofmodestoanalyse,number_of_atoms, number_of_species, dimnu, band_occupation, &
          &  elec_eigval_file_name, DPs, modes_to_analyse, readfromxmlfile)

 ! This subroutine reads the GenerateDataToPlotCrossovers.in file, which contains the names of the files where one must read the data
 ! as well as other relevant data. 

  ! System parameters
  integer, Intent(In)  :: NumberofDPs 
  integer, Intent(In)  :: Numberofmodestoanalyse
  integer, Intent(Out) :: number_of_atoms                ! The number of atoms of the system
  integer, Intent(Out) :: number_of_species                ! The number of different atomic species of the system
  integer, Intent(Out) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  integer, Intent(Out) :: band_occupation                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  character(100), Intent(Out) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  real(8), Intent(Out) :: DPs (NumberofDPs)               ! The values of the different displacement parameters considered 
  integer, Intent(Out) :: modes_to_analyse (Numberofmodestoanalyse)     ! The indices of modes to analyse (between 7 and 3N_atoms)
  integer, Intent(Out) :: readfromxmlfile                 ! Different from 0 if the eigenvalues are read from eigenval.xml from Quantum Espresso
 
   
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxstr4, auxstr5
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,iostat,unitt !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2


  number_of_atoms = 0; number_of_species = 0
  band_occupation=2; readfromxmlfile = 0
  


  ! Reading parameters from file
  
   write(*,*); write(*,*); 
   write(*,*) '  ************************************************************************************' 
   write(*,*) '  ********************************  NOW RUNNING   ************************************'   
   write(*,*) '  *********************** GenerateDataToPlotCrossovers.x *****************************' 
   write(*,*) '  ************************************************************************************'   
   write(*,*) 
   write(*,*); write(*,*) '  ****** FROM THE INPUT FILE OF GenerateDataToPlotCrossovers.x '
    write(*,*) ' (GenerateDataToPlotCrossovers.in): ***********************'   
  
  
  write(auxstr4,*)  NumberofDPs; call StripSpaces(auxstr4)
  write(*,*);  write(*,*) '         The number of considered displacement paramters is ',trim(auxstr4) 

  
     
   open(349,file='GenerateDataToPlotCrossovers.in',status='unknown')




   do k=1,250

      read(349,'(A)',IOSTAT=stat) inputline
      !write(*,*) inputline
      if (stat .ne. 0) go to 43317

      auxstr1=''; auxstr2=''; auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '#') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 11123
        end if
      end do ! i
 
      11123 continue
  
          
   
    if ((trim(auxstr1) .eq. 'Modes_to_analyse').or.(trim(auxstr1) .eq. 'modes_to_analyse').or. &
             & (trim(auxstr1) .eq. 'Modes_to_Analyse') )  then
      do l=1,Numberofmodestoanalyse
          read(349,*) modes_to_analyse(l)
      end do ! l
      go to 43317
    end if  
    
  end do ! k

43317 continue

   close(349) 
    
  
  
    open(345,file='GenerateDataToPlotCrossovers.in',status='unknown')

 

   do k=1,250

      read(345,'(A)',IOSTAT=stat) inputline
      !write(*,*) inputline
      if (stat .ne. 0) go to 3317

      auxstr1=''; auxstr2=''; auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '#') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 1238
        end if
      end do ! i
 
      1238 continue
  
          
   
    if ((trim(auxstr1) .eq. 'Displacement_parameters').or.(trim(auxstr1) .eq. 'displacement_parameters').or. &
             & (trim(auxstr1) .eq. 'Displacement_Parameters') .or. ( trim(auxstr1) .eq. 'DisplacementParameters' ))  then
      do l=1,NumberofDPs
          read(345,*) DPs(l)
      end do ! l
      go to 3317
    end if  
    
  end do ! k
   

3317 continue

   close(345) 
  
  
   
   open(345,file='GenerateDataToPlotCrossovers.in',status='unknown')


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



     if ((auxstr1 .eq. 'elec_eigval_file_name').or.(auxstr1 .eq. 'Elec_eigval_file_name'))   then
         elec_eigval_file_name = trim(auxstr2)
         call StripSpaces (elec_eigval_file_name)
     end if


   
    if ((auxstr1 .eq. 'Number_of_atoms').or.(auxstr1 .eq. 'number_of_atoms').or. &
             & (auxstr1 .eq. 'Number_atoms'))  Read(auxstr2, '(I3)' ) number_of_atoms 

    if ((auxstr1 .eq. 'band_occupation').or.(auxstr1 .eq. 'Band_Occupation').or. &
             & (auxstr1 .eq. 'Band_occupation'))  Read(auxstr2, '(I3)' ) band_occupation
             
    if ((auxstr1 .eq. 'Number_of_species').or.(auxstr1 .eq. 'number_of_species').or. &
             & (auxstr1 .eq. 'Number_species'))  Read(auxstr2, '(I3)' ) number_of_species           
             
    if ((auxstr1 .eq. 'min_band_index').or.(auxstr1 .eq. 'Min_band_index'))  Read(auxstr2, '(I3)' ) min_band_index 
    if ((auxstr1 .eq. 'max_band_index').or.(auxstr1 .eq. 'Max_band_index'))  Read(auxstr2, '(I3)' ) max_band_index
  
    if ((auxstr1 .eq. 'Read_from_xml_file').or.(auxstr1 .eq. 'read_from_xml_file'))  Read(auxstr2, '(I3)' ) readfromxmlfile
  
    
  end do ! k
  
   

1117 continue




   if ((band_occupation .ne. 1) .and. (band_occupation .ne. 2)) then
    band_occupation = 2
    write(*,*)
    write(*,*) '  <<<<<<<<<<<< WARNING: The number of electrons per band in ',trim(elec_eigval_file_name)
    write(*,*) '  was assumed to be 2. >>>>>>>>>>>'
    write(*,*)   
  end if
 
   dimnu = number_of_atoms*3 ! Strictly speaking it is number_of_atoms*3-6, but we make it bigger to avoid problems of lacks of space (some freqs are 0 and they are anyway read)
 
  

  if (number_of_atoms .le. 0) then
     write(*,*) ' ERROR: Please, specify a valid number of atoms in GenerateDataToPlotCrossovers.in; '
     write(*,*) '        E.g.: Number_of_atoms =  2 '
     stop
  end if

  if (number_of_species .le. 0) then
     write(*,*) ' ERROR: Please, specify a valid number of different atomic species in GenerateDataToPlotCrossovers.in; '
     write(*,*) '        E.g.: Number_of_species =  1 '
     stop
  end if
  
  write(*,*) 
   write(*,*)   '   Name of the file that stores electronic eigenvalues and masses:  ',trim(elec_eigval_file_name),', '
  
  if (band_occupation .eq. 1) then
    write(*,*) '                                                                      with up to 1 electron per band.'
  else  if(band_occupation .eq. 2) then
    write(*,*) '                                                                      with up to 2 electrons per band.'
  else  
    write(*,*) '  ERROR: Unphysical number of electrons per band (',band_occupation,').'; write(*,*)
    call exit(1)     
  end if     
  
  write(*,'(A,I5)') '                               Number of different atomic species = ', number_of_species
  write(*,'(A,I5)') '                                                  Number of atoms = ', number_of_atoms
  write(auxstr4,*) dimnu-6; call StripSpaces(auxstr4)
  write(auxstr5,*) dimnu; call StripSpaces(auxstr5) 
  write(*,*)        '                                       Number of phonon branches =    ', &
           & trim(auxstr5), ' (',trim(auxstr4),' valid)'
  write(*,*)
  

  inquire(file=elec_eigval_file_name, exist=exist)       
  if (.not. exist) then
      write(*,*) 
      write(*,*) '**********************************************************************************' 
      write(*,*) '    ERROR: The file ',trim(elec_eigval_file_name),' does not exist. Please, provide it.'
      write(*,*) '**********************************************************************************'
      write(*,*)
      call exit(1) 
  end if
  
  


end subroutine read_input_file


!-------------------------------------------------------------------------------------------
 
!       
 subroutine read_QE_occ_unocc_numbers (elec_eigval_file_name, band_occupation, Ne, sizeepsilon)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives the numbers of occ and unocc states.

   ! System parameters
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  integer, Intent(In)  :: band_occupation               ! 1 or 2, number of electrons per band in the QE output file
  integer, Intent(Out) :: Ne                            ! The number of electrons of the system (i.e. the number of occupied states)
  integer, Intent(Out) :: sizeepsilon                   ! size of the vector storing the electronic eigenvalues
  
     
    ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  integer :: Nks     ! Number of Kohn-Sham states
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9  
  
    
  write(*,*)      
  write(*,*); write(*,*) '  ****** FROM THE QUANTUM ESPRESSO OUTPUT (', trim(elec_eigval_file_name),'): *******************'
  write(*,*)   



  open(433,file=elec_eigval_file_name,status='unknown')

  do m=1,20000
  
   do k=1,2000

      read(433,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 2117

      auxstr1=''; auxstr2=''; auxstr3=''
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( (auxstr1 .eq. 'numberofelectrons=') .or.  ( auxstr1 .eq. 'Numberofelectrons=')  ) then  
          go to 2235
        end if
      end do ! i
      
   end do !k   
   
 2235 continue      

   do j=i+1,i+18
     letter = inputline(j:j)
     auxstr2 = trim(auxstr2)//letter
   end do
   read(auxstr2,*) auxr1
   Ne=nint(auxr1)
   write(*,'(A,I5)') '                                       Number of electrons = ', Ne

  close(433)
  
end do ! m
  
 2117 continue 
  
  
  close(433)
  
  open(451,file=elec_eigval_file_name,status='unknown')

  do m=1,20000
  
   do k=1,2000

      read(451,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 3117

      auxstr1=''; auxstr2=''; auxstr3=''
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( (auxstr1 .eq. 'numberofKohn-Shamstates=') .or.  ( auxstr1 .eq. 'NumberofKohn-Shamstates=')  ) then  
          go to 2237
        end if
      end do ! i
      
   end do !k   
   
 2237 continue      
      
   do j=i+1,i+18
     letter = inputline(j:j)
     auxstr2 = trim(auxstr2)//letter
   end do
   
   read(auxstr2,*) Nks
  write(*,'(A,I5)') '                                Number of Kohn-Sham states = ', Nks
  write(*,*)
      
  close(451)
  
end do ! m
  
  3117 continue
  
  if (band_occupation .eq. 2) then
     sizeepsilon = Nks   !2*Nks
  else if (band_occupation .eq. 1) then 
     sizeepsilon = Nks
  else   
      write(*,*) '  ERROR: Unphysical number of electrons per band'
      call exit(1)
  end if

 end subroutine read_QE_occ_unocc_numbers
  


!-------------------------------------------------------------------------------------------




subroutine epsilon_to_atomic_units (units_epsilonorig,dimband,epsilonorig,epsilon)
     
  ! We use atomic units as done in Ponce et al. PRB 90, 214304 (2014); this is because our derivation of gDW
  ! depends on the definitions given in that paper.   
  ! This subroutine writes the phonon freqs. in the units of the electronic eigenvectors, multiplying the former
  ! values of the omegas by the 'factor'.   
       
  ! System parameters
  character(40), Intent(In) :: units_epsilonorig
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  real(8), Intent(In) :: epsilonorig(dimband)			! Electronic eigenvalues in their original units
  real(8), Intent(Out):: epsilon(dimband)				! Electronic eigenvalues in atomic units
  
 
  ! Local variable
  integer :: i
  real(8) :: factor  


  
 
  factor = 0.000000000
 
 
 !  write(*,*) '  *** Now converting the units of the the electronic eigenvalues (',trim(units_epsilonorig)
 !  write(*,*) '                                  to atomic units). ***'

  if      (units_epsilonorig.eq.'cm**-1') then
    factor=1.000000000/219474.631370515d0 ! 0.00000455633d0
  else if (units_epsilonorig.eq.'meV') then  
    factor=1.0000000000/27211.385056d0  ! 0.0000367502d0 
  else if (units_epsilonorig.eq.'eV') then  
    factor=1.000000000/27.211385056d0 ! 0.0367502d0
  else if (units_epsilonorig.eq.'a.u.') then  
    factor=1.0000000000           
  else
    write(*,*) '  ERROR: Conversion of units not established; please, use "eV", "meV", "cm**-1", or "a.u.", '
    write(*,*) '  or rewrite the <<epsilon_to_atomic_units>> subroutine'
    call exit(1)    
  end if

   do i=1,dimband 
     epsilon(i) = (epsilonorig(i))*factor
   end do
  
! We write the phonon frequencies with a format that Maat can read:
  
! We write the electronic eigenvalues with a format that Maat can read:
! epsilonorig0.dat: 1st row: number of bands, number of k-points, number of spins; then all three indices and the corresponding unrenormalized electronic eigenvalue
  open(633,file='epsilon0.dat',status='unknown')
  write(633,*) dimband, '  1    1  ' 
  do i=1,dimband
    write(633,*) i, '  1  ', '  1  ', epsilon(i)
  end do
  close(633)
  
        
       
end subroutine epsilon_to_atomic_units   


!-------------------------------------------------------------------------------------------
