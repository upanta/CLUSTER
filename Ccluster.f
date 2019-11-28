      
      module myatoms
      
      TYPE :: atom
      CHARACTER(15) :: name     ! name of object
      REAL, Dimension(3) :: r   ! coordinates
      INTEGER ::  status = 0    ! 0 if inside, -1 if outside
      END TYPE atom
      
      end module myatoms


      Program cluster
      use myatoms
      implicit none
      
      
      TYPE(atom) :: at(8)       ! 8 atoms in the unit cell
      Integer :: i, xi, yi, zi, ki
      Real :: LatPar            ! lattice parameter  
      Real :: Radius            ! radius of the cluster
      Integer :: Ncell, Natom   ! estimate to reserve arays
      TYPE(atom), dimension(:), allocatable :: Atoms 
      Integer :: Total_counter = 0
      Integer :: Atoms_belong = 0
      Integer :: counter_quadrant = 0
      Real :: distance
      Integer :: N_arg
      character (15) :: STR_arg
      
!
! reading arguments
!

      N_arg = IARGC()           ! reads the number of arguments

      if (N_arg .lt. 1) then
         write(*,*) 'this program  needs  one argument'
         write(*,*) N_arg, 'arguments were entered'
         write(*,*) 'enter the desired radius '
         write(*,*) '------------------------------------'
         write(*,*)
         write(*,*) 'Usage'
         write(*,*) 'cluster.exe 4.5'
         write(*,*) 'cluster.exe creates a nanoparticle of rad 4.5 Ang.'
         write(*,*) 'the coordinates in the XYZ format'
         write(*,*)
         write(*,*) 'exiting ...'
         stop
      end if


        CALL GETARG(1, STR_arg) ! returns the argument 1
        read(STR_arg,*) Radius ! read filename first


      


      
      Ncell = 4                 !NINT(Radius/LatPar)
      
!     write (*,*) "Ncell = ", Ncell 
!     if (Ncell .lt. Radius/LatPar) Ncell = Ncell+1
      
      Natom= 8*8*(Ncell**3)
      allocate(Atoms(Natom))


! Defining the crystal structure
      
      LatPar = 3.57             ! in agstroms, but  x 0.5291772 to get a.u. 
      
      at(1)%name = 'C'
      at(2)%name = 'C'
      at(3)%name = 'C'
      at(4)%name = 'C'
      
      at(5)%name = 'C'
      at(6)%name = 'C'
      at(7)%name = 'C'
      at(8)%name = 'C'


      at(1)%r = (/ 0.0, 0.0, 0.0 /) 
      at(2)%r = (/ 0.5, 0.5, 0.0 /)
      at(3)%r = (/ 0.0, 0.5, 0.5 /)
      at(4)%r = (/ 0.5, 0.0, 0.5 /) 
      
      at(5)%r = (/ 0.25, 0.25, 0.25 /)
      at(6)%r = (/ 0.75, 0.75, 0.25 /)
      at(7)%r = (/ 0.75, 0.25, 0.75 /)
      at(8)%r = (/ 0.25, 0.75, 0.75 /)
      
      
      
!     do i=1, 16, 1 
!     write(*,10) at(i)%name, at(i)%r(1:3), at(i)%status
!     end do
!      write(*,*)
      
      
!     (1)
!     generating atoms in the first quadrant
      counter_quadrant = 0
      do zi=1, Ncell, 1
         do yi=1, Ncell, 1
            do xi=1, Ncell, 1
               do i=1, 8, 1
                  counter_quadrant = counter_quadrant + 1 
                  CALL Copy_atom(Atoms( counter_quadrant ), at(i))
                  Atoms( counter_quadrant )%r(1) = at(i)%r(1) + (xi-1)
                  Atoms( counter_quadrant )%r(2) = at(i)%r(2) + (yi-1)  
                  Atoms( counter_quadrant )%r(3) = at(i)%r(3) + (zi-1)
               end do           ! i
            end do              ! xi 
         end do                 ! yi
      end do                    ! zi
      Total_counter = counter_quadrant
      
      
      
!     (2)
!     replicating atoms from the 1-st quadrant to the second one
!     by simpli shifting it along x, Ncell tymes of the lattice parameter   
      do ki=1, counter_quadrant, 1
         Total_counter = Total_counter + 1
         CALL Copy_atom(Atoms( Total_counter ), Atoms(ki))
         Atoms( Total_counter )%r(1) = Atoms(ki)%r(1) - Ncell
      end do
  
      
!     (3)
!     replicating atoms from the 1-st quadrant to the (3) one
!     by simpli shifting it along y, Ncell tymes of the lattice parameter   
      do ki=1, counter_quadrant, 1
         Total_counter = Total_counter + 1
         CALL Copy_atom(Atoms( Total_counter ), Atoms(ki))
         Atoms( Total_counter )%r(2) = Atoms(ki)%r(2) - Ncell
      end do
      
      
      
!     (4)
!     replicating atoms from the 1-st quadrant to the (4) one
!     by simpli shifting it along x,y by Ncell tymes of the lattice parameter   
      do ki=1, counter_quadrant, 1
         Total_counter = Total_counter + 1
         CALL Copy_atom(Atoms( Total_counter ), Atoms(ki))
         Atoms( Total_counter )%r(1) = Atoms(ki)%r(1) - Ncell
         Atoms( Total_counter )%r(2) = Atoms(ki)%r(2) - Ncell
      end do
      
      
      
      
!     (5-8)
!     replicating atoms from the 1-4 quadrants to the 5-8 ones
!     by simpli shifting it along z by Ncell tymes of the lattice parameter   
      do ki=1, counter_quadrant*4, 1
         Total_counter = Total_counter + 1
         CALL Copy_atom(Atoms( Total_counter ), Atoms(ki))
         Atoms( Total_counter )%r(3) = Atoms(ki)%r(3) - Ncell
      end do
      
      
      
      
      do ki=1, Total_counter, 1
         Atoms(ki)%r(1)=Atoms(ki)%r(1)*LatPar
         Atoms(ki)%r(2)=Atoms(ki)%r(2)*LatPar
         Atoms(ki)%r(3)=Atoms(ki)%r(3)*LatPar
         distance=sqrt(Atoms(ki)%r(1)*Atoms(ki)%r(1)+
     1     Atoms(ki)%r(2)*Atoms(ki)%r(2)+
     2     Atoms(ki)%r(3)*Atoms(ki)%r(3))
         if (distance .gt. Radius) then
!        write(*,*) distance 
              Atoms(ki)%status = - 1
         endif
      end do
      
!     select atoms which are inside the sphere 
      


      
      do i=1, Total_counter, 1
         if (Atoms(i)%status .eq. 0) then
            Atoms_belong = Atoms_belong + 1
         endif   
      end do
      write(*,12) Atoms_belong
      write(*,*)
      
      do i=1, Total_counter, 1
         if (Atoms(i)%status .eq. 0) then
            write(*,10)  Atoms(i)%name, Atoms(i)%r(1:3)
         endif   
      end do
      
      
      
      
      deallocate(Atoms)
      
      
      
!..........................................
      
 10   format(a2, 3f12.6, i3)
 11   format(i4, '  ', a2, 3f12.6, i3)
 12   format(i4)
      
      
      stop
      end program cluster
      

      
!-----------------------------------------
      
      
      subroutine Copy_atom(receiver,source)
      use myatoms
      
      TYPE(atom), intent(out) :: receiver
      TYPE(atom), intent(in) :: source
      
      receiver%name = source%name
      receiver%r(1) = source%r(1)
      receiver%r(2) = source%r(2)
      receiver%r(3) = source%r(3)
      receiver%status  = source%status
      
      end subroutine Copy_atom
 
