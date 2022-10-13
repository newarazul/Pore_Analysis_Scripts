!******************************************************************!

! Program to calculate ssvvaf and OH strecth SFG spectra for Pores !
! Adaption by Frederik Zysk
!******************************************************************!

program pore_ssvvaf
        
        implicit none

!****************************************************************!

! Note: Variables (Nstep, Ncfstep, Nifreq, Nifreq, dt and rcut) should
! !

! be initialzed manually before exceuting the program            !

!****************************************************************!
integer                 :: Nstep, idum1 
integer                 :: Natom
integer                 :: istep,iatom
integer                 :: h2o_atoms
integer                 :: Noh
integer                 :: places(1:500)
integer,parameter       :: Ncfstep = 2000
integer,parameter       :: Nifreq =1, Nffreq = 5000
real(kind=8),parameter       :: pi = 3.14159265
integer                 ::counterx,countery

real(8)                 ::grid(0:100,0:100)
real(8)                 ::density2(0:100,0:100)
integer                 :: iradius
Real(kind=8), parameter :: dt = 0.5d0
real(kind=8), parameter :: rcut = 2
real(8), allocatable :: vz(:,:)
real(8), allocatable :: Bxx(:,:), Byy(:,:)
real(8), allocatable :: hcoord(:,:,:)
real(8), allocatable :: ocoord(:,:,:)
real(8), allocatable :: beta(:,:)

real(8)                 ::density(0:600)
real(8)              ::         statistic(0:600)
Real(8), dimension (:,:,:), allocatable :: h2o_coord
Real(8), dimension (:,:,:), allocatable :: coord
character(len=3), dimension (:), allocatable :: atom_name
Real(8), dimension (:,:), allocatable :: com_xycoord
Real(8), dimension (:,:,:), allocatable :: vel
real(8), dimension (:,:), allocatable :: alpha

real(8) vvaf_tot(0:Ncfstep)
real(8) vvaf_ave(0:Ncfstep)
real(8) svvaf_ave(0:Ncfstep)
real(8) vdos(Nifreq:Nffreq), svdos(Nifreq:Nffreq)
real(8) srgvaf(Nifreq:Nffreq), sigvaf(Nifreq:Nffreq)
real(8) rgvaf(Nifreq:Nffreq), igvaf(Nifreq:Nffreq)
real(8)         ::multi,sizing
!determine number of atoms and steps with routine

call stepsandatoms(Natom, Nstep)

!allocate fitting memory according to number of atoms and steps

allocate(coord(1:3,1:Natom,1:Nstep))
allocate(atom_name(1:Natom))
allocate(com_xycoord(1:2,1:Nstep))

!Reading the trajectory

call read_standard_trajectory(Natom,Nstep,atom_name,coord)

!Calculating the center of mass for H20 only

!call findh2o(Natom, Nstep, coord, atom_name, h2o_atoms, places)

!allocate(vel(1:3,1:h2o_atoms,1:Nstep))
!allocate(h2o_coord(1:3,1:h2o_atoms,1:Nstep))


!call fillh2oarray(Nstep, Natom,h2o_atoms, places, coord, h2o_coord)




!call centerofmass(h2o_atoms,Nstep,h2o_coord,com_xycoord)


!calculate actual density



!do istep=1, Nstep
! do iatom=1, h2o_atoms, 3
!   do iradius=1,600
!     if(sqrt((h2o_coord(1,iatom,istep)*h2o_coord(1,iatom,istep))+&
!             &(h2o_coord(2,iatom,istep)*h2o_coord(2,iatom,istep))) &
!             &.lt. (iradius*0.01) .and. sqrt((h2o_coord(1,iatom,istep)*h2o_coord(1,iatom,istep))+&
!             &(h2o_coord(2,iatom,istep)*h2o_coord(2,iatom,istep)))&
!             &.ge. ((iradius-1)*0.01)) then 
!                statistic(iradius)=statistic(iradius)+1
!     end if
!end do
!end do   
!end do




!open(59,file="statistics.dat")
open(58,file="18oh_water.xyz")
do istep=1, Nstep
 write(58,*) 144
 write(58,*)
 do iatom=Natom-143, Natom
        write(58,*) atom_name(iatom),coord(1,iatom,istep),coord(3,iatom,istep),coord(2,iatom,istep)
 end do
end do

close(58)
        
!do iradius=1, 600
!density(iradius)=((18.0152*statistic(iradius)/6.022E+023)/((((iradius*0.01)*&
!&(iradius*0.01))-(((iradius-1)*0.01)*((iradius-1)*0.01)))*3.14*14.32*1.0E-24))
!        write(58,*) 0, 1
!        write(58,*) (iradius-1)*0.01, density(iradius)/Nstep
!        write(59,*) (iradius-1)*0.01, statistic(iradius)/Nstep
!        end do

!close(58)
!close(59)
!sizing=15
!print*,sizing
!open(59,file="flat.dat")
!do istep=1, Nstep
!  print*, "We are at", istep
!        do iatom=1, h2o_atoms, 3
!        print*, iatom
!        do counterx=0,sizing-1
!             if(h2o_coord(1,iatom,istep) .ge. (-14.32+(counterx*28.64/sizing)) .and. &
!                     &h2o_coord(1,iatom,istep) .lt. (-14.32+((counterx+1)*28.64/sizing))) then
!               do countery=0,sizing-1
!               if(h2o_coord(2,iatom,istep) .ge. (-14.32+(countery*28.64/sizing)) .and. &
!                       &h2o_coord(2,iatom,istep) .lt. (-14.32+((countery+1)*28.64/sizing))) then
!                   grid(counterx,countery)=grid(counterx,countery)+1
!                end if

!                end do

!        end if
!        end do 
        
!        end do
!        end do
!multi=(1/sizing)
!print*, multi

!do counterx=0,sizing-1
! do countery=0,sizing-1
!density2(counterx,countery)=(((grid(counterx,countery)/(6.022E+023*Nstep))*18.0152)/(28.64*multi*28.64*multi*14.32*1.0E-024))
!   print*, grid(counterx,countery)
!   write(59,*) (counterx*28.64/sizing)-14.32, (countery*28.64/sizing)-14.32, density2(counterx,countery)
!   end do
!   end do







end program




!*************************************************!
!Routine to read trajectory                       !
!*************************************************!
        subroutine read_standard_trajectory(Natom,Nstep,atom_name,coord)
        
        implicit none
        
        integer, intent(in)     :: Natom, Nstep
        Real(8)                 :: coord(1:3,1:Natom,1:Nstep)
        character(LEN=3)        :: atom_name(1:Natom)

        integer :: iatom, istep, i1
        character(len=10) :: cdum1

        open(51, file="Pore.xyz",status="old")
        do istep=1, Nstep
                read(51,*)  
                read(51,*) cdum1, cdum1, cdum1
                do iatom =1, Natom
                        read(51,*) atom_name(iatom),(coord(i1,iatom,istep),i1=1,3)
                end do
        end do
        close(51)
        end subroutine


!*********************************************!
!Routine to find the number of atoms and steps!
!*********************************************!
        subroutine stepsandatoms(Natom, Nstep)

        implicit none
        
        integer, intent(out) :: Natom, Nstep                               
        integer :: reads                                                       
        integer :: lines
        
        lines = 1
        open(50, file="Pore.xyz")
        read(50,*) Natom
        do  
                read(50,*,IOSTAT=reads)
                if (reads /= 0) exit 
                lines = lines + 1
        end do
        print*, lines
        Nstep= lines/(2+Natom)
        print*, Natom, Nstep    
        close(50)
        Nstep=200000    
        end subroutine



!*************************************************!
!Routine to check for H20 and generate a new array!
!*************************************************!
        subroutine findh2o(Natom, Nstep, coord, atom_name, h2o_atoms, places)
        
        implicit none
        
        integer                 :: Natom, Nstep, iatom, istep, h2o_atoms
        
        Real(8)                 :: coord(1:3,1:Natom,1:Nstep)
        character(LEN=3)        :: atom_name(1:Natom)
        integer                 :: places(1:500) 


        h2o_atoms = 1
        istep = 1
        print*, Nstep
        do iatom = 1, Natom -2
            if ((atom_name(iatom) .eq. "O" .and. atom_name(iatom+1) .eq. "H" &
                & .and. atom_name(iatom+2) .eq. "H")) then
                   if(((coord(1,iatom,1)- coord(1,iatom+1,1)) <= 2.0) .and. &
                      &  ((coord(1,iatom,1)- coord(1,iatom+2,1)) <=2.0)) then
                          places(h2o_atoms) = iatom
                          h2o_atoms = h2o_atoms +1 
                          places(h2o_atoms) = iatom+1
                          h2o_atoms = h2o_atoms +1
                          places(h2o_atoms) = iatom+2
                          h2o_atoms = h2o_atoms +1
                          print*, h2o_atoms
                   end if
            end if
        end do
        h2o_atoms = h2o_atoms -1
        end subroutine

        subroutine fillh2oarray(Nstep, Natom, h2o_atoms, places, coord, h2o_coord)

        implicit none
        
        integer :: istep, iatom, h2o_atoms, Nstep, Natom
        integer :: places(1:500)
        real(8) :: h2o_coord(1:3,1:h2o_atoms,1:Nstep)
        real(8) :: coord(1:3,1:Natom,1:Nstep)

        do istep =1, Nstep !-1
          do iatom =1, h2o_atoms
            h2o_coord(1:3,iatom,istep) = coord(1:3,places(iatom),istep)
          end do
        end do    
        end subroutine
        
        !Subroutine for the center of mass!
        !*********************************!
        
        subroutine centerofmass(h2o_atoms,Nstep,h2o_coord,com_xycoord)
                
        implicit none
        integer                  :: istep, iatom, Nstep, h2o_atoms
        Real(8), dimension (:,:) :: com_xycoord(1:2,1:Nstep)
        Real(8)                  :: h2o_coord(1:3,1:h2o_atoms,1:Nstep)

        com_xycoord = 0.d0
         do istep = 1, Nstep !-1 
           do iatom = 1 , h2o_atoms, 1
             com_xycoord(1:2,istep) = com_xycoord(1:2,istep) + h2o_coord(1:2,iatom,istep)
           end do     
         end do
        com_xycoord(1:2,:) = com_xycoord(1:2,:) / (h2o_atoms)
        do istep = 1, Nstep !-1
          do iatom = 1, h2o_atoms
            h2o_coord(1:2,iatom,istep) = h2o_coord(1:2,iatom,istep)-com_xycoord(1:2,istep)
          end do
        end do
        

        end subroutine

        !subroutine calculate the nearest distance!
        subroutine calc_distan(distan)

        implicit none
        
        real(8) distan(1:5)

        distan(4) = distan(1) * distan(1) + distan(2) * distan(2) + distan(3) *distan(3)
        distan(5) = dsqrt(distan(4))
        
        return
        end subroutine
