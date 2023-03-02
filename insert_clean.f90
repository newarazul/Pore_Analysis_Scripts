program watervolume

implicit none

integer         :: Natom, Nstep,hatom

character(LEN=3), allocatable          :: atom_name(:), atom_name2(:)

real(8), allocatable            :: coord(:,:,:), coord_center(:,:,:), coord2(:,:,:), coord_center2(:,:,:)
real(8), allocatable            :: com_xycoord(:,:)
real(8), allocatable            :: newpore(:,:,:) 
real(8)                         :: positionx, positiony, positionx2, positiony2, positionz,positionz2
real(8)                         :: positionx3, positiony3, positionz3
real(8), dimension(3)           :: shift1, shift2, shift3
real(8)                         :: cellsizex,cellsizey,cellsizez



!determine number of atoms and steps to allocate appropriate memory 
call stepsandatoms(Natom,Nstep)

!allocating arrays
allocate(newpore(1:3,1:1563,1:Nstep))
allocate(atom_name2(1:1563))
allocate(atom_name(1:Natom))
allocate(coord(1:3,1:Natom,1:Nstep))
allocate(com_xycoord(1:3,1:Nstep))
allocate(coord_center(1:3,1:Natom,1:Nstep))
allocate(coord2(1:3,1:1563,1:Nstep))
allocate(coord_center2(1:3,1:1536,1:Nstep))


!read in trajectory. In this case bulk water trajectory
!keep in mind water trajectory is hard coded as "water.xyz"
call read_standard_trajectory(Natom,Nstep,atom_name,coord)

!use center of mass routine
call centerofmass(Nstep,Natom,coord,coord_center)


!cell dimension in Angstrom
cellsizex=25.919
cellsizey=25.919
cellsizez=25.78


!shift for first channel
positionz=cellsizez/2
positionx=4.511
positiony=(sqrt(3.0)/3)*cellsizey

!shift for second channel
positionx2=0.5*cellsizex
positiony2=(sqrt(3.0)/6)*cellsizey
positionz2=cellsizez/2

!shift for third channel
positionx3=0.0
positiony3=0.0
positionz3=cellsizez/2



!create shift vector
shift1=(/positionx, positiony, positionz/)
shift2=(/positionx2, positiony2, positionz2/)
shift3=(/positionx3, positiony3, positionz3/)

hatom=1




!first channel filled
call chosewater(Natom,Nstep,coord_center,cellsizez,shift1,hatom,newpore,atom_name2)
!second channel filled
call chosewater(Natom,Nstep,coord_center,cellsizez,shift2,hatom,newpore,atom_name2)
!third channel filled
call chosewater(Natom,Nstep,coord_center,cellsizez,shift3,hatom,newpore,atom_name2)

!writes water of desired geometry into file "new.xyz"
call writetrajectory(Natom,Nstep,hatom,newpore,atom_name2)

!check density from number of water atoms and provided pore volume
call checkdensity(hatom,cellsizez)








end program

!**************************************************!
!Routine to calculate the density to check for realism                  !
!**************************************************!
        subroutine checkdensity(hatom,cellsizez)
        
        implicit none
        
        integer                 :: hatom
        real(8)                 :: cellsizez,density

        density= ((((hatom/3)/6.022E+23)*18.01528)/((4.5)*(4.5)*3.14*cellsizez*1.00E-24))*0.5
        print*,"Calculated density is",density 

        end subroutine

!**************************************************!
!Routine to write into file                     !
!**************************************************!
        subroutine writetrajectory(Natom,Nstep,hatom,newpore,atom_name2)
        
        implicit none
        
        integer                 :: hatom,Natom,Nstep,iatom
        character(LEN=3)        :: atom_name2(1:Natom)
        real(8)                 :: newpore(1:3,1:Natom,1:Nstep)

        open(54, file="new.xyz")
        do iatom=1, hatom-1
                write(54,*) atom_name2(iatom), newpore(1:3,iatom,1)
        end do
        close(54)        


        end subroutine



!**************************************************!
!Routine to take water with the correct radius                        !
!**************************************************!
        subroutine chosewater(Natom,Nstep,coord_center,cellsizez,shift,hatom,newpore,atom_name2)
        
        implicit none
        
        integer, intent(in)     :: Natom, Nstep
        Real(8)                 :: coord_center(1:3,1:Natom,1:Nstep),newpore(1:3,1:Natom,1:Nstep)
        character(LEN=3)        :: atom_name2(1:Natom)
        Real(8)                 :: cellsizez
        integer :: iatom, hatom
        real(8), dimension(3)           ::shift
        
        
        do iatom=1, Natom-2, 3
        if(sqrt((coord_center(1,iatom,Nstep)*coord_center(1,iatom,Nstep))+(coord_center(2,iatom,Nstep)*&
                &coord_center(2,iatom,Nstep))) .lt. (4.65) .and. abs(coord_center(3,iatom,Nstep)) &
                &.lt. (cellsizez/2)) then
           newpore(1:3,hatom,Nstep)  =coord_center(1:3,iatom,Nstep)+shift(1:3)
           newpore(1:3,hatom+1,Nstep)=coord_center(1:3,iatom+1,Nstep)+shift(1:3)
           newpore(1:3,hatom+2,Nstep)=coord_center(1:3,iatom+2,Nstep)+shift(1:3)
           atom_name2(hatom)="O"
           atom_name2(hatom+1)="H"
           atom_name2(hatom+2)="H"
           hatom=hatom+3
        end if
        end do

        end subroutine

!**************************************************!
!Rouutine to read trajectory                       !
!**************************************************!
        subroutine read_standard_trajectory(Natom,Nstep,atom_name,coord)
        
        implicit none
        
        integer, intent(in)     :: Natom, Nstep
        Real(8)                 :: coord(1:3,1:Natom,1:Nstep)
        character(LEN=3)        :: atom_name(1:Natom)

        integer :: iatom, istep, i1
!        character(len=10) :: cdum1

        open(51, file="water.xyz",status="old")
        do istep=1, Nstep
                read(51,*)  
                do iatom =1, Natom
                        read(51,*) atom_name(iatom),(coord(i1,iatom,istep),i1=1,3)
!                        print*, atom_name(iatom),(coord(i1,iatom,istep),i1=1,3)
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
        open(50, file="water.xyz")
        read(50,*) Natom
        do  
                read(50,*,IOSTAT=reads)
                if (reads /= 0) exit 
                lines = lines + 1
        end do
        print*, lines
        Nstep =1
        !        Nstep= lines/(1+Natom)
        print*, Natom, Nstep    
        close(50)
!        Nstep=1    
        end subroutine

!*********************************************!
!Subroutine for the center of mass            !
!*********************************************!
        
       subroutine centerofmass(Nstep,Natom,coord,coord_center)
                
        implicit none
        integer                  :: iatom, Natom, Nstep
        Real(8), dimension (:,:) :: coord(1:3,1:Natom,1:Nstep)
        Real(8)                  :: coord_center(1:3,1:Natom,1:Nstep)
        Real(8)                  :: center(1:3)

        !print*, Nh2o
        !calculate the center of x and y 
        center(1:3)= 0.d0 
           do iatom = 1 , Natom
             center(1:3) = center(1:3) + coord(1:3,iatom,1)
           end do     
        center(1:3) = center(1:3)/Natom
!        print*, center(1:3)
          do iatom =1, Natom        
            coord_center(1:3,iatom,1)=coord(1:3,iatom,1)-center(1:3)
!            coord_center(3,iatom,1)=coord(3,iatom,1)-center(3)
!        print*, coord_center(1:3,iatom,1)
          end do
        end subroutine
