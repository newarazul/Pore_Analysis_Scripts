!******************************************************************!

! Program to calculate ssvvaf and OH strecth SFG spectra for Pores !
! Adaption by Frederik Zysk
!******************************************************************!

function betrag3(i) result(betrag1) 
        real(8) :: i(1:3)
        real(8) :: betrag1
        betrag1=sqrt((i(1)**2) +(i(2)**2) +(i(3)**2))
end function 

function ortskoordinate(x) result(radius)
        real(8) :: x(1:3)
        real(8) :: radius        
        radius=NINT(sqrt((x(1)*x(1))+(x(2)*x(2)))*10)
END function 

function distance_atomic(x,y) result(distance)
        real(8) :: x(1:3)
        real(8) :: y(1:3)
        real(8) :: distance
        distance=sqrt(((x(1)-y(1))*(x(1)-y(1)))+((x(3)-y(3))*(x(3)-y(3)))+((x(2)-y(2))*(x(2)-y(2))))
END function

!function distance_periodic(x,y) result(distance_p)
!        real(8) :: x(1:3)
!        real(8) :: y(1:3)
!        real(8) :: distance_p
!        real(8) :: distance_atomic
!        real(8) :: box_size=14.32



!function c_vector(x,y) result(vec1)
!        real(8) :: x(1:3)
!        real(8) :: y(1:3)
!        real(8) :: vec1(1:3)
!        vec1(1:3)=y(1:3)-x(1:3)
!end function

function hydrogen_angle(x,y,z) result(angle_h)
        !give O O H
        !anlge between three points
        real(8) :: betrag3
        real(8) :: y(1:3)
        real(8) :: x(1:3)
        real(8) :: z(1:3)
        real(8) :: vec1(1:3), vec2(1:3)
        real(8) :: top,down
        real(8) :: angle_h, mid1
        vec1(1:3)=z-x
        vec2(1:3)=y-x
        top=dot_product(vec1,vec2)
        down=betrag3(vec1)*betrag3(vec2)
        if(down .ne. 0) then
        mid1=top/down
        end if
        if(mid1 .gt. 1.0) then
                mid1=1.0
        end if
        angle_h=dacos(mid1)*57.53
end function        







!angle=0
!do istep=1,Nstep-1
!do iatom=1, Noh
!det=(ocoord(1,iatom,istep)*(hcoord(3,iatom,istep)-ocoord(3,iatom,istep)))-&
!        &((hcoord(1,iatom,istep)-ocoord(1,iatom,istep))*ocoord(3,iatom,istep))
!dot=ocoord(1,iatom,istep)*(hcoord(1,iatom,istep)-ocoord(1,iatom,istep))+&
!        &(ocoord(3,iatom,istep)*(hcoord(3,iatom,istep)-ocoord(3,iatom,istep)))
!angle2(iatom,istep)=atan2(det,dot)
!angle2(iatom,istep)=angle2(iatom,istep)*57.53
!print*, angle2(iatom,istep)
!end do
!end do




program pore_ssvvaf
!use     ::mymath        
        implicit none

!****************************************************************!

! Note: Variables (Nstep, Ncfstep, Nifreq, Nifreq, dt and rcut) should
! !

! be initialzed manually before exceuting the program            !

!****************************************************************!
integer                 :: Nstep, idum1, idum2 
integer                 :: Natom, me,you, irange
integer                 :: istep,iatom
integer                 :: h2o_atoms, know
integer                 :: Noh, ioh, batom, catom
integer                 :: places(1:500)
integer,parameter       :: Ncfstep = 2000
integer,parameter       :: Nifreq =1, Nffreq = 5000
real(kind=8),parameter       :: pi = 3.14159265

Real(kind=8), parameter :: dt = 0.5d0
real(kind=8), parameter :: rcut = 2
real(8), allocatable :: vz(:,:)
real(8), allocatable :: Bxx(:,:), Byy(:,:)
real(8), allocatable :: hcoord(:,:,:)
real(8), allocatable :: ocoord(:,:,:)
real(8), allocatable :: beta(:,:)
real(8), allocatable :: angle2(:,:)
real(8)                 :: det, dot
real(8)                 ::betrag3, hydrogen_angle
!real(8)                 ::betrag1,betrag2
real(8), allocatable    ::bwinkel(:,:)


!real(8)                 ::c_vector
Real(8), dimension (:,:,:), allocatable :: h2o_coord
Real(8), dimension (:,:,:), allocatable :: coord, center_coord, coord_hydroxil_o
character(len=3), dimension (:), allocatable :: atom_name
Real(8), dimension (:,:), allocatable :: com_xycoord
Real(8), dimension (:,:,:), allocatable :: vel
real(8), dimension (:,:), allocatable :: alpha
!real(8)                :: angle(1:360,1:360)
real(8), allocatable   :: vektoroh1(:,:,:)
real(8), allocatable   :: vektoroh2(:,:,:)
real(8),allocatable    ::halfvektor(:,:,:)
real(8)                 :: distan1(1:3), numberofbonds(1:100)
real(8)                 :: distan2(1:3)
real(8)                 :: thecount, speed, angle_total
real(8)                 :: dist
real(8)                 ::a(1:3), b(1:3)
real(8)                 :: distance_atomic, ortskoordinate
real(8)                 :: countbonds(0:500), countwater(0:500), countbondsx(0:500), countbondsy(0:500), countbonds2(0:500)

real(8) vvaf_tot(0:Ncfstep)
real(8) vvaf_ave(0:Ncfstep)
real(8) svvaf_ave(0:Ncfstep)
real(8) vdos(Nifreq:Nffreq), svdos(Nifreq:Nffreq)
real(8) srgvaf(Nifreq:Nffreq), sigvaf(Nifreq:Nffreq)
real(8) rgvaf(Nifreq:Nffreq), igvaf(Nifreq:Nffreq)

!determine number of atoms and steps with routine

call stepsandatoms(Natom, Nstep)

!allocate fitting memory according to number of atoms and steps

allocate(coord(1:3,1:Natom,1:Nstep))
allocate(atom_name(1:Natom))
allocate(com_xycoord(1:3,1:Nstep))

!Reading the trajectory

call read_standard_trajectory(Natom,Nstep,atom_name,coord)


!calculating the distance of OH stretches
Noh=(Natom/3)*2
allocate(Bxx(1:Natom,1:Nstep))
allocate(Byy(1:Natom,1:Nstep))
allocate(vz(1:Natom,1:Nstep))
allocate(hcoord(1:3,1:Natom,1:Nstep))
allocate(ocoord(1:3,1:Natom,1:Nstep))
!allocate(alpha(1:h2o_atoms,1:Nstep))
allocate(beta(1:Natom,1:Nstep))
!allocate(angle2(1:Natom,1:Nstep))
allocate(bwinkel(1:Natom,1:Nstep))
allocate(vektoroh1(1:3,1:Natom,1:Nstep))
allocate(vektoroh2(1:3,1:Natom,1:Nstep))
allocate(halfvektor(1:3,1:Natom,1:Nstep))
allocate(center_coord(1:3,1:Natom,1:Nstep))
allocate(coord_hydroxil_o(1:3,1:Natom,1:Nstep))





!wrap trajectory in z direction (cp2k does not wrap and produces the trajectory for this analysis)

do istep=1, Nstep
do iatom=1, Natom
   if(coord(3,iatom,istep) .gt. 7.16) then
           coord(3,iatom,istep)=coord(3,iatom,istep)-14.32
   end if
end do
end do







com_xycoord = 0.d0
 do istep = 1, Nstep!-1 
    do iatom = 1 , Natom
      com_xycoord(1:3,istep) = com_xycoord(1:3,istep) + coord(1:3,iatom,istep)
           end do     
         end do
        com_xycoord(1:3,:) = com_xycoord(1:3,:)/(Natom)
        do istep = 1, Nstep !-1
          do iatom = 1, Natom
            center_coord(1:3,iatom,istep) = center_coord(1:3,iatom,istep)-com_xycoord(1:3,istep)
!            print*, com_xycoord(1:3,istep)
    end do
 end do

 
!copy all water into hcoord and ocoord array



do istep =1, Nstep-1
  idum1 = 0
  do iatom =607, 753, 3
idum1=idum1+1
hcoord(:,idum1,istep) = coord(:,iatom+1,istep)
ocoord(:,idum1,istep) = coord(:,iatom,istep)
!print*, ocoord(:,idum1,istep)
idum1=idum1+1
hcoord(:,idum1,istep) = coord(:,iatom+2,istep)
ocoord(:,idum1,istep) = coord(:,iatom,istep)
end do
end do



!check for silanols that are close to water and calculate their coordination number for angle < 120 and distance < 3.55 A


do istep=1, Nstep
  idum1=0
  do iatom=551, 607
    idum2=0
    if(atom_name(iatom) .eq. "H") then
      idum1=idum1+1      
      coord_hydroxil_o(1:3,idum1,istep)=coord(1:3,iatom,istep)
      do batom=1, 702
      if(atom_name(batom) .eq. "O") then
        if(distance_atomic(coord_hydroxil_o(1:3,idum1,istep),coord(1:3,batom,istep)) .lt. 1.0) then
                idum2=idum2+1
!                print*, "paare Hydrogen:", iatom, "hydrogen", batom, "vorkommen_gleiche",idum2, "Hs", idum1
                do catom=1,98,2 
                   if(distance_atomic(ocoord(1:3,catom,istep),coord(1:3,batom,istep)) .lt. 3.55) then        
                        print*, "yes", catom, distance_atomic(ocoord(1:3,catom,istep),coord(1:3,batom,istep))        
                        if(hydrogen_angle(ocoord(:,catom,istep),ocoord(:,batom,istep),hcoord(:,iatom,istep)) .lt. 120) then
                        countbonds2(ortskoordinate(ocoord(1:3,catom,istep)))=countbonds2(ortskoordinate(ocoord(1:3,catom,istep)))+1
                        end if
                end if     
        end do
        end if
      end if
      end do
      end if 
  end do
end do 


!do istep=1, Nstep
!  do iatom=1, Natom
!    if(atom_name(iatom) .eq. "O") then
!      coord_hydroxil_o(1:3,iatom,istep)=coord(1:3,iatom,istep)
!    end if 
!  end do
!end do 


!do iatom=1, h2o_atoms
!do istep = 1, Nstep-1
!        print*, h2o_coord(1:3,iatom,Nstep-1), h2o_coord(1:3,iatom,Nstep)
!        vel(:,iatom,istep)=((h2o_coord(:,iatom,istep+1)-h2o_coord(:,iatom,istep))/dt)
        
        
        
        !        if(vel(3,iatom,istep) .ge. 0.1) then
!                print*, vel(1:3,iatom,istep)
!        end if
!enddo
!enddo

!angle between O and OH vektor. 

!angle=0
!do istep=1,Nstep-1
!do iatom=1, Noh
!det=(ocoord(1,iatom,istep)*(hcoord(3,iatom,istep)-ocoord(3,iatom,istep)))-&
!        &((hcoord(1,iatom,istep)-ocoord(1,iatom,istep))*ocoord(3,iatom,istep))
!dot=ocoord(1,iatom,istep)*(hcoord(1,iatom,istep)-ocoord(1,iatom,istep))+&
!        &(ocoord(3,iatom,istep)*(hcoord(3,iatom,istep)-ocoord(3,iatom,istep)))
!angle2(iatom,istep)=atan2(det,dot)
!angle2(iatom,istep)=angle2(iatom,istep)*57.53
!print*, angle2(iatom,istep)
!end do
!end do

!calculate OH hydroxil











do istep=1, Nstep-1
! do irange=1,99
  do iatom=1, 98, 2
    countwater(ortskoordinate(ocoord(:,iatom,istep)))=countwater(ortskoordinate(ocoord(:,iatom,istep)))+1
!    print*, ortskoordinate(ocoord(:,iatom,istep))
    do batom=1, 98, 2
        if(distance_atomic(ocoord(:,iatom,istep),ocoord(:,batom,istep)) .lt. 3.55) then                
                if(iatom .ne. batom) then
                        if(hydrogen_angle(ocoord(:,batom,istep),ocoord(:,batom,istep),hcoord(:,iatom,istep)) .lt. 120) then
                        countbonds(ortskoordinate(ocoord(:,iatom,istep)))=countbonds(ortskoordinate(ocoord(:,iatom,istep)))+1
                        end if
                end if
        end if
        end do
        
        end do
        end do    
        

do irange=0, 200
!     print*, countwater(irange)
     if(countwater(irange) .ne. 0) then
     countbondsx(irange)=countbonds(irange)/countwater(irange)
     countbondsy(irange)=countbonds2(irange)/countwater(irange)
     end if 
end do


!calculate with OH hydroxil groups

!get O and H for the hydroxil groups

!ch33. get C-H-O

















!            &sqrt(ocoord(1,iatom,istep)*ocoord(1,iatom,istep)) .lt. (irange*0.1)) then
!       do batom=1, Noh, 2
!          if(ocoord(1,iatom,istep)-ocoord(1,batom,istep) .lt. 1.3) then
!                numberofbonds(irange)=numberofbonds(irange)+1
!          end if
!       end do
!       end if 
!   end do
!   end do
!   end do

   
!#print*, numberofbonds


!do idum1=1,360
!do idum2=1,360
!        if(angle_total .gt. 0) then
!                angle(idum1,idum2)=angle(idum1,idum2)-(angle_total/(360*360))
!        end if
!end do
!end do
 
 

open(61, file="conformation.dat")
do idum1=1,400
      write(61,*) idum1, countbondsx(idum1), countbonds(idum1),countwater(idum1),countbondsy(idum1), countbonds2(idum1)
      end do
close(61)








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
                read(51,*) cdum1 
                read(51,*) 
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
!        print*, lines
!        Nstep= lines/(2+Natom)
!        print*, Natom, Nstep    
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
        


!Subroutine for length of vector!
subroutine calculatelengthofvector(alpha,Noh, vel,h2o_atoms,Nstep,h2o_coord, hcoord,&
               & ocoord,vz, Bxx)


integer         :: istep, iatom, Nstep, h2o_atoms, Noh, number1
real(8)         :: distan1(1:5), distan2(1:5)
real(8)         :: nvec1(1:3), nvec2(1:3)
real(8)         :: veloh1(1:3), veloh2(1:3), veloh11(1:3), veloh22(1:3)
real(8)         :: vel(1:3,1:h2o_atoms,1:Nstep)
real(8)         :: hcoord(1:3,1:Noh,1:Nstep), ocoord(1:3,1:Noh,1:Nstep)
real(8)         :: vz(1:Noh,1:Nstep)
real(8)         :: velproj1, velproj2
real(8)         :: h2o_coord(1:3,1:h2o_atoms,1:Nstep)
real(8)         :: Bxx(1:Noh,1:Nstep)
real(8)         ::beta(1:Noh,1:Nstep)
real(8)         ::alpha(1:h2o_atoms,1:Nstep)

do istep =1, Nstep-1
  idum1 = 0
  number1 = 0
  do iatom =1 , h2o_atoms-2, 3
    distan1(1:3) = h2o_coord(1:3,iatom+1,istep)- h2o_coord(1:3,iatom,istep)
    distan2(1:3) = h2o_coord(1:3,iatom+2,istep)- h2o_coord(1:3,iatom,istep)    
    idum1=idum1+1
    call calc_distan(distan1)
    nvec1(1:3) = distan1(1:3) / distan1(5)
    call calc_distan(distan2)
    nvec2(1:3) = distan2(1:3) / distan2(5)
    veloh1(1:3) = vel(1:3,iatom+1,istep) - vel(1:3,iatom,istep)
    veloh2(1:3) = vel(1:3,iatom+2,istep) - vel(1:3,iatom,istep)
    velproj1 = dot_product(veloh1,nvec1)
    velproj2 = dot_product(veloh2,nvec2)
    vz(idum1,istep)=veloh1(2)
    Bxx(idum1,istep) = velproj1
    idum1=idum1+1
    vz(idum1,istep)=veloh2(2)
    Bxx(idum1,istep) = velproj2
   end do
end do


end subroutine



!****************************************************************************!
!Subroutine that calculates the angle between the O vector and the OH- vector!
!****************************************************************************!

subroutine calculate_angle(Noh,hcoord,ocoord,Nstep,h2o_atoms,alpha)


implicit none
        
integer         :: istep, idum1, Nstep, h2o_atoms, ioh, Noh,iatom
real(8)         :: hcoord(1:3,1:Noh,1:Nstep)
real(8)         :: ocoord(1:3,1:Noh,1:Nstep)
real(8)         :: alpha(1:h2o_atoms,1:Nstep)
real(8)         :: dot
real(8)         :: det


!calculate angle 
        
do istep =1, Nstep -1
idum1 =1
 do ioh = 1, Noh, 2
 dot=ocoord(1,ioh,istep)
 det=-1*ocoord(2,ioh,istep)
! print*,dot,det
if(det .ne. 0) then
 alpha(idum1,istep)=atan2(det,dot)
 idum1=idum1+1
 alpha(idum1,istep)=atan2(det,dot)
 idum1=idum1+1
 alpha(idum1,istep)=atan2(det,dot)
 idum1=idum1+1
else 
 if(dot .ge. 0) then
 alpha(idum1,istep)=0
 idum1=idum1+1
 alpha(idum1,istep)=0
 idum1=idum1+1
 alpha(idum1,istep)=0
 idum1=idum1+1
 else if(dot .lt. 0) then
 alpha(idum1,istep)=3.14159
 idum1=idum1+1
 alpha(idum1,istep)=3.14159
 idum1=idum1+1
 alpha(idum1,istep)=3.14159
 idum1=idum1+1
 end if
 end if
 end do
end do
print*, atan2(0.0,-0.5)

end subroutine

!******************************************************************************!
!subroutine to define the surface and the correlation. possible to do auto,
!intra and inter. comment the fitting chapter and adapt rcut and the top of the
!program. rcut=2 is intra and anything bigger inter.
!******************************************************************************!

subroutine surface_correlation(rcut,Noh,vvaf_ave,vvaf_tot,Nstep,Ncfstep, hcoord, ocoord,Bxx,vz)

implicit none

integer         :: idum1, istep, Nstep, Ncfstep, Noh, i1, iatom
integer         :: iatom2, H_bond, iiatom

real(8)         ::beta(1:Noh,1:Nstep)
real(kind=8)    :: rcut
real(kind=8)    :: ri(1:5), rj(1:5), rij(1:5), distan_OH(1:5)
real(8)         :: ddistan(1:5)
real(8)         :: ocoord(1:3,1:Noh,1:Nstep),hcoord(1:3,1:Noh,1:Nstep)
real(8)         :: Bxx(1:Noh,1:Nstep)
real(8)         :: vz(1:Noh,1:Nstep)
real(8)         :: vvaf_tot(0:Ncfstep), vvaf_ave(0:Ncfstep)
real(8)         :: alphaoh1(1:Noh,1:Nstep)

print*, "4: Calculating the ssvvaf and ssvvcf"
vvaf_ave=0.d0; vvaf_tot=0.d0; idum1 = 0

do istep = 1, Nstep-Ncfstep-1
!istep=1
idum1 = 0
 do iatom = 1, Noh
  ddistan(1:3) = hcoord(1:3,iatom,istep) - ocoord(1:3,iatom,istep)
  call calc_distan(ddistan)
!  H_bond = 0

! do iiatom = 1, Noh
!    if(ocoord(2,iatom,istep) .ge. 0.0d0) then
!    distan_OH(1:3) = hcoord(1:3,iatom,istep) - ocoord(1:3,iiatom,istep)
!    call calc_distan(distan_OH)
!    else 
!    distan_OH(1:3) = hcoord(1:3,iatom,istep) + ocoord(1:3,iiatom,istep)
!    call calc_distan(distan_OH)
!    end if
!    if(distan_OH(5) > 1.25 .and. distan_OH(5) < 2.45)then
!      H_bond = H_bond + 1
!    end if 
!  end do


 
!  ri(1:3) = (hcoord(1:3,iatom,istep) + ocoord(1:3,iatom,istep))/2
  
  !Define surface region. In this case divided into top and lower part and in
  !rings aroung the zylinder. Adjust the xy value as desired. 
  !***************************************************************************! 
  if(sqrt((ocoord(1,iatom,istep)*ocoord(1,iatom,istep)) +&
     & (ocoord(2,iatom,istep)*ocoord(2,iatom,istep))) .ge. 4.0d0 .and. sqrt((ocoord(1,iatom,istep)&
     &*ocoord(1,iatom,istep)) +&
     & (ocoord(2,iatom,istep)*ocoord(2,iatom,istep))) .le. 7.0d0 .and. ddistan(5) .lt. 1.2) then
         idum1 = idum1 + 1
         do i1 = 0, Ncfstep
!         if(vz(iatom,istep)  .ge. 0) then
!         idum1=idum1+1
         !Autocorellation calculated with projection of vz(here vxy) on Ovector
         vvaf_tot(i1) = vvaf_tot(i1) + vz(iatom,istep)*Bxx(iatom,istep+i1)
!         if(vz(iatom,istep) .ge. 0.1 .or. vz(iatom,istep) .lt. -0.1) then
!                 print*, vz(iatom,istep)
!         end if 
!         end if
!         print*, vz(iatom,Nstep-1)
         !*alphaoh1(iatom,istep+i1)
         !*COS(alphaoh1(iatom,istep))
        
         !Intracorrelation if you dont want to include inter!

!         if(mod(iatom,2) == 1) then
!           vvaf_tot(i1) = vvaf_tot(i1) + vz(iatom+1,istep) * Bxx(iatom,istep+i1)
!         elseif(mod(iatom,2) == 0) then
!           vvaf_tot(i1) = vvaf_tot(i1) + vz(iatom-1,istep) * Bxx(iatom,istep+i1)
!         endif
        
         ! Cross bond coorelations between two different bonds in a same water
         ! molecule and different water molecules. calculating intra and inter
         ! correlation

        
!         do iatom2 = 1,Noh

!          if(iatom == iatom2) then
!          elseif(iatom .ne. iatom2) then
!                rj(1:3)=(hcoord(1:3,iatom2,istep)+ocoord(1:3,iatom2,istep))/2
!                rij(1:3) = ri(1:3) - rj(1:3)
!                call calc_distan(rij)
!                if(rij(5) < rcut) then
!                  vvaf_tot(i1) = vvaf_tot(i1) + vz(iatom,istep) * Bxx(iatom2,istep+i1)*COS(alphaoh1(iatom,istep))
!                endif

!         endif

         end do

 
!         end do
! else if(ocoord(2,iatom,istep) .lt. 0.0d0 .and. sqrt((ocoord(1,iatom,istep)*ocoord(1,iatom,istep)) +&
!     & (ocoord(2,iatom,istep)*ocoord(2,iatom,istep))) .ge. 2.0d0 .and. sqrt((ocoord(1,iatom,istep)&
!     &*ocoord(1,iatom,istep)) +&
!     & (ocoord(2,iatom,istep)*ocoord(2,iatom,istep))) .le. 7.0d0 .and. ddistan(5) .lt. 1.2) then
!             idum1 = idum1 + 1
!             do i1 = 0, Ncfstep
        !Autocorellation calculated with projection of vz(here vxy) on Ovector
!         vvaf_tot(i1) = vvaf_tot(i1) + vz(iatom,istep)*Bxx(iatom,istep+i1)
         !*COS(alphaoh1(iatom,istep))
         
         !lower part of the surface division
!         elseif(ocoord(2,iatom,istep) .lt. 0.0d0 .and. sqrt((ocoord(1,iatom,istep)*ocoord(1,iatom,istep)) +&
!                & (ocoord(2,iatom,istep)*ocoord(2,iatom,istep))) .ge. 7.0d0 .and. sqrt((ocoord(1,iatom,istep)&
!                &*ocoord(1,iatom,istep)) +&
!                & (ocoord(2,iatom,istep)*ocoord(2,iatom,istep))) .le. 13.0d0 .and. (ddistan(5) .le. 1.2d0) .and. H_bond .eq. 0) then
!          idum1 = idum1 + 1
!          do i1 = 0, Ncfstep
         !Autocorellation. vz as a negative contribution
!             vvaf_tot(i1)=vvaf_tot(i1)-vz(iatom,istep)*Bxx(iatom,istep+i1)*COS(alphaoh1(iatom,istep))
         
!          do iatom2 = 1,Noh

!          if(iatom == iatom2) then
!          elseif(iatom .ne. iatom2) then
!                rj(1:3)=(hcoord(1:3,iatom2,istep)+ocoord(1:3,iatom2,istep))/2
!                rij(1:3) = ri(1:3) - rj(1:3)
!                call calc_distan(rij)
!                if(rij(5) < rcut) then
!                  vvaf_tot(i1) = vvaf_tot(i1) - vz(iatom,istep) * Bxx(iatom2,istep+i1)*COS(alphaoh1(iatom,istep))
!                endif

!          endif

!end do
        end if
        end do
end do 
!print*, idum1
vvaf_ave = vvaf_tot/dble(idum1)


!vvaf_ave=vvaf_ave/(MAXVAL(vvaf_ave))



end subroutine

!***********************************************************************
!subroutine to smooth and write ssvvaf.dat file
!***********************************************************************
subroutine smoothing(dt, pi, Ncfstep, vvaf_ave, svvaf_ave)

implicit none

real(kind=8)    :: dt 
real(kind=8)    :: pi
integer         :: i1, Ncfstep
real(8)         :: svvaf_ave(0:Ncfstep), vvaf_ave(0:Ncfstep)

open(56, file="new.dat", STATUS="replace",action="write")
do i1 = 0, Ncfstep
 svvaf_ave(i1) = vvaf_ave(i1) *dcos(pi*dble(i1)*dt/(dble(Ncfstep)*dt*2.d0)) *dcos(pi*dble(i1)*dt/(dble(Ncfstep)*dt*2.d0)) 
 write(56,*) i1*dt, vvaf_ave(i1), svvaf_ave(i1)
enddo
close(56) 
end subroutine

!**************************************************************************
!FFT subroutine
!**************************************************************************

subroutine fft (vaf,vdos,rgvaf,igvaf,dt,Nstep,Nifreq,Nffreq)


implicit none



integer Nstep, Nifreq, Nffreq



real(8) vaf(0:Nstep), dt
real(8) rgvaf(Nifreq:Nffreq), igvaf(Nifreq:Nffreq)
real(8) agvaf(Nifreq:Nffreq), vdos(Nifreq:Nffreq)

real(8) light_vel !cm/fs
parameter(light_vel = 299792458.d-13)
real(8) pi
parameter(pi = 3.14159265)

integer ifreq, istep
rgvaf = 0.d0; igvaf = 0.d0; agvaf = 0.d0

do ifreq = Nifreq, Nffreq
  do istep = 0, Nstep-1
  if((vaf(istep) *sin(2.d0*pi*light_vel*istep*(dt)*ifreq) * dt) .ge. 0) then
    rgvaf(ifreq) = rgvaf(ifreq) + vaf(istep) *cos(2.d0*pi*light_vel*istep*(dt)*ifreq) * dt
    igvaf(ifreq) = igvaf(ifreq) + vaf(istep) *sin(2.d0*pi*light_vel*istep*(dt)*ifreq) * dt
    endif 
  enddo
enddo
print*, cos(2.d0*pi*light_vel*istep*(dt)*4000)
do ifreq = Nifreq, Nffreq
   vdos(ifreq) = rgvaf(ifreq)*rgvaf(ifreq) + igvaf(ifreq)*igvaf(ifreq)
   vdos(ifreq) = sqrt(vdos(ifreq))
enddo

end subroutine

















