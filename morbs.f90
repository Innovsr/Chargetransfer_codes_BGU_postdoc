!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***************************************************************!
!***************************************************************!
!*********This code is written by DR. SOURAV ROY****************!
!*****************Post-Doctoral Fellow *************************!
!**********In the department of Chemistry***********************!
!******Under the Supervision of Prof. Sebastian Kozuch**********!
!*********Ben-Gurion University of the Negev, Israle ***********!
!*******************JUNE,2016 **********************************!
!***************************************************************!
!***************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!****************************************************************************************************
!!! ************************************Format of the
!inputfile**************************************
!****************************************************************************************************
!code words||| acc_don_thr: acceptor, doner, third ||| dimer_ch_sp:
!charge and spin of the dimer    *
!||| mono1_ch_sp: charge and spin of momoner-1 ||| mono2_ch_sp: charge
!and spin of monomer-2        *
!||| inp: name of inputfile  ||| opt_inp: if u want to optimize the
!input .com file then write "yes"*
!|||mono_option: preparation of two monomers inputfile : put 'distance':
!to short out the atoms for *
!two monomerfiles with the help of distance among them atoms or put
!'radius': use radius of atoms   *
!for shorting or 'user':for manual insert the number of the atoms  |||
!mono1_atom: insert the serial*
!number of atoms for monomer1 where first number should be total number
!of atoms in monomer1        *
!||| mono2_atom: similar for monomer2 ||| block starting and ending by
!'@'specifying the first four *
!line of dimer and monomer input file ||| run_gauss: is the command to
!run the gaussian             *
!||| run_cube: command to generate .cube files
!*
!                                                                                                   *
!acc_don_thr  1 5 3
!*
!dimer_ch_sp  0 1
!*
!mono1_ch_sp  -1 1
!*
!mono2_ch_sp  1 1
!*
!inp          ethilin.com
!*
!opt_inp      no
!*
!mono_option    distance
!*
!mono1_atom   4 1 2 3 4
!*
!mono2_atom   4 5 6 7 8
!*
!@@@@@
!*
!%chk=ethilin_dim.chk
!*
!%nproc=8
!*
!%mem=50gb
!*
!#p m062x 6-311g(d,p) opt freq
!*
!@@@@@@@
!*
!run_gauss    g09
!*
!run_cube     cubegen
!*
!****************************************************************************************************
!****************************************************************************************************
!****************************************************************************************************


module commondat
implicit none

real,public::dongroup,accgroup,at_rad
real::dummy,min_grd,grd_sp,ocdi,ocm1,ocm2
integer::mono1,mono2,at_covrad,mno_at,nacc,ndon,third,fg9,fg8,fg7,fg6,fg5,fg3,fg2,fg1,mnopt,mon1,mon2&
,mnm1,mnm2,steps,MNXY,MNZ,MN,fg4,noorbdi,noorbm1,noorbm2,NMODI,NMOM1,NMOM2,orbdi(350),orbm1(350),orbm2(350)&
,fg10,fg11,fg12,fg13,fg14,fg15,fg16,fg17,fg18,flg20,flg21,flg22,nmon1,nmon2,&
sdimo(1000),pdimo(1000),ddimo(1000),sdi,pdi,ddi,sm1,sm2,pm1,pm2,dm1,dm2,sm1mo(1000),pm1mo(1000),dm1mo(1000),&
sm2mo(1000),pm2mo(1000),dm2mo(1000),ddi1,ddimo1(1000),dm11,dm1mo1(1000),dm21,dm2mo1(1000),noinp,lopt(5)
character(len=18)::wordstr(40)
character(15):: basis,gauss
character(75):: cubcom,fchk_com,cubcommo
character(10):: funal
character(len=2):: stt,atmdi(15),atmm1(15),atmm2(15)
character(len=2),public :: accst, donst,at_list(88)
character(len=35),public:: infname, inputfilename,lname
character(len=70)::line1,line2,line3
character(len=90)::line13,line27,line28,line29
character(len=110)::line35,line36,line37
character(len=110)::inp_log(5)
character(len=35)::molname,cubname,logname,mononame1,mononame2,dimname,chknamedi,chknamem1,chknamem2&
,fchknamedi,fchknamem1,fchknamem2,cubnamedi,cubnamem1,cubnamem2,lognamem2,lognamem1,lognamedi&
,substract,integral,cubnamemodi,cubnamemom1,cubnamemom2,substract_oct,dimer_oct,mono1_oct,mono2_oct
character(len=2),public:: ch1, ch2, sp1, sp2, disc, ch0,sp0
dimension at_covrad(88),at_rad(88)

common /coord/dummy(300,3),ocdi,ocm1,ocm2
common /atom/accst(100),donst(100),stt(300),min_grd(3),steps(3),grd_sp
common /options/basis,funal,ch0,sp0,ch1,sp1,ch2,sp2,infname,gauss,cubcom
common /allname/molname,cubname,logname,mononame1,mononame2,dimname,chknamedi,chknamem1,chknamem2&
,fchknamedi,fchknamem1,fchknamem2,cubnamedi,cubnamem1,cubnamem2,lognamem2,lognamem1,lognamedi&
,substract,integral,line13,cubnamemodi,cubnamemom1,cubnamemom2,cubcommo,substract_oct,dimer_oct,&
mono1_oct,mono2_oct,inp_log
common /input/line1,line2,line3,wordstr,fchk_com,line27,line28,line29,line35,line36,line37
common /mono/mono1,mono2,MNXY,MNZ,MN,NMODI,NMOM1,NMOM2
common /group/dongroup(300,3),accgroup(300,3),mon1(200),mon2(200),orbdi,orbm1,orbm2
common /dat/mno_at,nacc,ndon,third,fg9,fg8,fg7,fg6,fg5,fg4,fg3,fg2,fg1,mnopt,mnm1,mnm2,noorbdi,noorbm1,noorbm2&
,fg10,fg11,fg12,fg13,fg14,fg15,fg16,fg17,fg18,flg20,flg21,flg22,&
sdimo,pdimo,ddimo,sdi,pdi,ddi,sm1,sm2,pm1,pm2,dm1,dm2,sm1mo,pm1mo,dm1mo,&
sm2mo,pm2mo,dm2mo,ddi1,ddimo1,dm11,dm1mo1,dm21,dm2mo1,noinp,lopt,nmon1,nmon2


DATA at_covrad /37, 32, 134, 90, 82, 77, 75, 73, 71, 69, 154, 130, 118, 111,&
106, 102, 99,97, 196, 174, 144, 136, 125, 127, 139, 125, 126, 121, 138, 131,&
126, 122, 119, 116, 114, 110, 211, 192, 162, 148, 137, 145, 156, 126, 135, 131, 153,&
148, 144, 141, 138, 135, 133, 130, 225, 198, 169, 204, 203, 201, 199, 198, 198,&
196, 194, 192, 192, 189, 190, 187, 160, 150, 138, 146, 159, 128, 137,&
128, 144, 149, 148, 147, 146, 140, 150, 145, 260, 221/

DATA at_list /'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',&
'Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',&
'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr',&
'Nd','Pm','SM','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',&
'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra'/


DATA at_rad /3.5, 3.5, 11.1, 8.7, 7.2, 6.5, 5.5, 5.0, 4.6, 4.3, 13.5, 12.0, 10.9,&
9.9, 9.3, 8.6, 8.1, 7.7, 17.4, 14.9, 12.3, 11.3, 10.5, 10.1, 10.0, 10.0, 9.9, 9.9,&
10.0, 10.7, 10.8, 10.5, 10.3, 9.9, 9.8, 9.6, 18.5, 16.4, 13.9, 12.4, 11.5, 11.1,&
10.9, 10.7, 10.7, 11.0, 11.5, 12.7, 12.3, 12.1, 12, 11.7, 11.4, 11.2, 20.1, 17.0,&
14.5, 14.1, 14.1, 14.1, 14.0, 13.9, 15.9, 13.8, 13.6, 13.6, 13.5, 13.5, 13.4,&
13.4, 13.4, 12.3, 11.5, 11.1, 11.0, 10.8, 10.9, 11.1, 11.5, 12.8, 12.7, 12.6,&
12.5, 12.5, 12.4, 12.4, 26.0, 22.1/

save
end module commondat
!!***************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main

use commondat
implicit none

real::finish,start,time1
integer:: argnum,k,l,ndate,j,i,l3
logical :: fileexists
character(len=18) :: orbin
character(len=35)::iname,name,name2,name1
!character(len=35)::mononame1, mononame2, dimname
character(len=55)::namo
character(len=2) :: atomin,natom
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
character(len=50)  :: inplogname
integer,dimension(8) :: values



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif
argnum=iargc()
if(argnum.eq.0)then
PRINT*,'input file not exist please put the inputfile name'
stop
endif

if(argnum.gt.1)then
print*,'please put only one argument here : inputfile name'
stop
endif

call getarg(1,inputfilename)
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=21,file=TRIM(inputfilename),status='old')
ELSE
PRINT*,'SORRY This input file does not exist or you may not provide the &
filename at all'
stop
ENDIF

call readinp


iname=inp_log(1)
infname=iname
!print*,infname
call prepname
write(*,104)trim(iname)
104 Format(// a,' is taken as your initial gaussian input file'//)
l=len(TRIM(iname))
lname=iname(l-3:l)

if(lname.ne.'.log') then
print*,'There is no gaussian log file present in the input file'
endif
!if(lname.eq.'.xyz') call readcom(iname)

call cubegrid_read
if(noorbdi.ne.0.or.noorbm1.ne.0.or.noorbm2.ne.0)then
call molorb
endif
do i=1,noinp
namo=inp_log(i)
if(lname.eq.'.log') call readlog(namo)
call readfchk(namo)
if(lopt(i).eq.1) call readorb_lin(namo)
if(lopt(i).eq.2) call readorb_nlin(namo)
enddo

!400 print*,'sourav1'
!if (fg4.eq.1)then
!call molorb
!endif

!print*,'sourav2'
!400 stop
stop
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readlog(iname)
!reading of .log file for molecular coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
integer ::  i, j, k, nargs ,num_args,p,ndate, jj, kk,&
lenth, l,  io, MDP, ii,iii,cn,an,at,mno_at1
character(len=56)::lines
character(len=35)::iname,logn
character(len=70)::lines2(300)
dimension an(300)
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
integer,dimension(8) :: values



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif

mnm1=0
mnm2=0
do i=1,300
do j=1,3
dummy(i,j)=0.0
enddo
enddo
open(unit=23,file=iname,status='old')
MDP=0
do
read(23,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(23)

i=0
do iii=1,MDP
read(23,'(a)')lines
if(lines(26:46).eq.'Standard orientation:'.or.lines(27:44).eq.'Input orientation:')then
if(iii.gt.i)i=iii
endif
enddo
rewind(23)

do iii=1,i
read(23,*)
enddo
do iii=i+1,MDP
read(23,'(a)')lines
if(lines(2:28).eq.'Rotational constants (GHZ):'.or.lines(2:41).eq.'Symmetry turned off by external request.'&
.or.lines(21:48).eq.'Distance matrix (angstroms):')goto 500
enddo
500 k=iii-i
rewind(23)

j=0
!print*,i,i+4,i+5,i+k-2,i+k-1,MDP
do ii=1,i+4
read(23,*)
enddo
do ii=i+5,i+k-2
j=j+1
!read(23,'(a,a,a,a,a,a)')(ann(j,k),k=1,6)
read(23,'(a)')lines2(j)
enddo
do iii=i+k-1,MDP
read(23,*)
enddo
open(unit=25,file='temp_2',status='unknown')
do i=1,j
write(25,'(a)')lines2(i)
enddo
rewind(25)
!print*,'j',j
mno_at=j
do j=1,mno_at
read(25,*)cn,an(j),at,(dummy(j,k),k=1,3)
stt(j)=at_list(an(j))
enddo

return
end subroutine readlog


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubegrid_read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
integer::i,j,div,extra,ndate
real::angtobohr,greatest,dummy_max,dummy_min,max_grd
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
integer,dimension(8) :: values
dimension dummy_max(300,3),dummy_min(300,3),max_grd(3)



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif



! change the unit from angstrom to Bohr

open(unit=3,file=cubname,status='unknown')
read(3,*)
read(3,*)
read(3,*)(min_grd(i),i=1,3)
read(3,*)steps(1),grd_sp
read(3,*)steps(2)
read(3,*)steps(3)
rewind(3)

extra=0
div=steps(3)/6
if(mod(steps(3),6).ne.0)extra=1
MN=div+extra
MNXY=steps(1)*steps(2)*MN
MNZ=steps(3)
return
end subroutine cubegrid_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readfchk(inplogname)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,i,i1,j,dialphae,m1alphae,m2alphae,alphaedi,betaedi,alphaem1,betaem1,alphaem2,betaem2,&
di1,p1,l1,l2,k,kk,c(100),d(100),f1,cc(50,100),dd(50,100),l3,j1,j2,nlast,mnnum,ndate,loop,j10,fg35,fg36
real::momat1(6,5),momat2(6,5),start
character(len=50)::lines,XYZ1(5),XYZ2(5),XYZ3(5),a1,inplogname,fchkname
character(len=3)::a,b,e,na,nd,matorb(5)
character(len=8)::fmt1,fmt2,f2,aa,bb,ans
character(len=75)::logn
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
integer,dimension(8) :: values

call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif


l=len(trim(inplogname))
fchkname=trim(inplogname(1:l-3))//trim('fchk')
lines='Sourav_Roy'
open(unit=43,file=trim(fchkname),status='old')

INQUIRE(file=trim(fchkname),EXIST=file1)
IF (file1) THEN


open(unit=43,file=trim(fchkname),status='old')
MDP=0
do
read(43,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(1:26).eq.'Number of alpha electrons')dialphae=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(43)

else
print*,trim(fchkname),'file does not exist'
stop
endif
read(43,*)
read(43,*)
read(43,*)(XYZ3(j),j=1,4),mno_at
do i=1,dialphae-3
read(43,*)
enddo
read(43,*)(XYZ3(j),j=1,5),alphaedi
read(43,*)(XYZ3(j),j=1,5),betaedi
read(43,*)(XYZ3(j),j=1,5),NMODI
rewind(43)

return
end subroutine readfchk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readorb_lin(inplogname)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,i,i1,j,dialphae,m1alphae,m2alphae,alphaedi,betaedi,alphaem1,betaem1,alphaem2,betaem2,&
di1,p1,l1,l2,k,kk,c(100),d(100),f1,cc(50,100),dd(50,100),l3,j1,j2,nlast,mnnum,ndate,loop,j10,fg35,fg36
real::momat1(6,5),momat2(6,5),start
character(len=50)::lines,XYZ1(5),XYZ2(5),XYZ3(5),a1,inplogname,lognm
character(len=3)::a,b,e,na,nd,matorb(5)
character(len=8)::fmt1,fmt2,f2,aa,bb,ans
character(len=75)::logn
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
integer,dimension(8) :: values

sdi=0
pdi=0
ddi=0
sm2=0
pm2=0
dm2=0
do i=1,1000
sdimo(i)=0
pdimo(i)=0
ddimo(i)=0
sm2mo(i)=0
pm2mo(i)=0
dm2mo(i)=0
enddo
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif
open(unit=9,file='Molecular_Orbitals_lin.dat',status='unknown')

l=len(trim(inplogname))
lognm=inplogname(l-8:l-5)
logn=inplogname(l-8:l-4)
if(lognm.eq.'mono')goto 800
fg35=0
fg36=0
INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:36)
if(lines(6:36).eq.'Molecular Orbital Coefficients:')fg35=1
if(lines(6:44).eq.'Alpha Molecular Orbital Coefficients:')fg36=1
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif

if (fg35.eq.1)loop=1
if (fg36.eq.1)loop=2
do j10=1,loop
INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:36)
if (fg35.eq.1)then
if(lines(6:36).eq.'Molecular Orbital Coefficients:')di1=MDP
endif
if (fg36.eq.1.and.j10.eq.1)then
if(lines(6:43).eq.'Alpha Molecular Orbital Coefficients:')di1=MDP
endif
if (fg36.eq.1.and.j10.eq.2)then
if(lines(6:42).eq.'Beta Molecular Orbital Coefficients:')di1=MDP
endif
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif

write(9,*)'Molecular orbital file of',molname
write(9,*)
if (fg35.eq.1) write(9,*)'Molecular orbitals of Dimer'
if (fg36.eq.1.and.j10.eq.1) write(9,*)'Alpha Molecular orbitals of Dimer'
if (fg36.eq.1.and.j10.eq.2) write(9,*)'Beta Molecular orbitals of Dimer'
write(9,*)
l=mod(NMODI,5)
nlast=l
k=(NMODI-l)/5
if(l.eq.0)kk=k
if(l.ne.0)kk=k+1
l1=0
l2=0
do i=1,di1+1
read(45,*)
enddo
do j=1,kk
do i=1,3
read(45,*)
enddo
do i=1,NMODI
!read(45,'(a,a,a)')a,b,e
read(45,*)a,b,e
!print*,a,b,e

if(nacc.lt.10)fmt1='(I1)'
if(nacc.gt.9.and.nacc.lt.100)fmt1='(I2)'
if(nacc.gt.99.and.nacc.lt.1000)fmt1='(I3)'
if(ndon.lt.10)fmt2='(I1)'
if(ndon.ge.10.and.ndon.lt.100)fmt2='(I2)'
if(ndon.ge.100.and.ndon.lt.1000)fmt2='(I3)'
write(na,fmt1)nacc
write(nd,fmt2)ndon
!print*,'na,nd',trim(na),nd,b
!read(e,'(I10)')f1
if(b.eq.na)then
l1=l1+1
c(l1)=i
!print*,'c',c(l1),nacc,f1,e
endif
if(b.eq.nd)then
l2=l2+1
d(l2)=i
!print*,'d',d(l2),f1,ndon,e
endif
enddo
enddo
rewind(45)

do i=1,di1+1
read(45,*)
enddo
do j=1,kk
l1=0
l2=0
l3=0
do i=1,3
read(45,*)
enddo
do i=1,c(j)
l3=l3+1
read(45,*)
enddo
do i=1,d(j)-1
l3=l3+1
!print*,i
read(45,*)aa,bb
if(bb.eq.'2PX')then
l1=l1+1
cc(j,l1)=l3
!print*,'cc',cc(j,l1),bb,l1
endif
enddo
do i=1,NMODI-(d(j)-1)-c(j)
l3=l3+1
!print*,i
read(45,*)aa,bb
if(bb.eq.'2PX')then
l2=l2+1
dd(j,l2)=l3
!print*,'dd',dd(j,l2),bb,l2
endif
enddo
!print*,'**************************',NMODI
enddo
rewind(45)

l1=0
do i=1,di1+1
read(45,*)
enddo
do j=1,k
j1=0
j2=0
do i=1,6
do l=1,5
momat1(i,l)=0.0
momat2(i,l)=0.0
enddo
enddo
read(45,*)
read(45,*)(matorb(i1),i1=1,5)
!print*,(matorb(i1),i1=1,5)
read(45,*)
do i=1,cc(j,1)-1
read(45,*)
enddo
do i=1,4
j1=j1+1
read(45,*)aa,bb,(momat1(j1,l),l=1,5)
!print*,(momat1(j1,l),l=1,5)
enddo
do i=1,dd(j,1)-(cc(j,1)+4)
read(45,*)
enddo
do i=1,3
j2=j2+1
read(45,*)aa,bb,(momat2(j2,l),l=1,5)
!print*,(momat2(j2,l),l=1,5),j2
enddo
do i=1,NMODI-(dd(j,1)+2)
read(45,*)
enddo
do i=1,5
l1=l1+1
write(9,*)'****************************MO=',l1,matorb(i)

if(momat1(1,i).eq.0.0.and.momat2(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat2(2,i).eq.0.0)then
if(momat1(3,i)*momat2(3,i).gt.0.0.and.momat1(3,i).ne.0.0)then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is sigma Anti-bonding orbital'
endif
if(momat1(3,i)*momat2(3,i).lt.0.0.and.momat1(3,i).ne.0.0)then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is sigma Bonding orbital'
endif
if(momat1(3,i)*momat2(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
ddi=ddi+1
ddimo(ddi)=l1
write(9,*)'Dimer MO=',l1,'is Delta orbital'
endif
if(momat1(3,i)*momat2(3,i).eq.0.0.and.momat1(4,i).ne.0.0)then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is Sigma orbital'
endif
else
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif

if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Anti-Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Anti-Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
endif
enddo
enddo

do i=1,6
do j=1,5
momat1(i,j)=0.0
momat2(i,j)=0.0
enddo
enddo

j1=0
j2=0
if(nlast.ne.0)then
read(45,*)
read(45,*)(matorb(i1),i1=1,nlast)
read(45,*)
!print*,i,cc(j,1)-1
do i=1,cc(1,1)-1
read(45,*)
enddo
do i=1,4
j1=j1+1
read(45,*)aa,bb,(momat1(j1,l),l=1,nlast)
!print*,aa,bb,(momat1(j1,l),l=1,2)
enddo
do i=1,dd(1,1)-(cc(1,1)+4)
read(45,*)
enddo
do i=1,3
j2=j2+1
read(45,*)aa,bb,(momat2(j2,l),l=1,nlast)
!print*,(momat2(j2,l),l=1,nlast)
enddo
do i=1,NMODI-(dd(1,1)+2)
read(45,*)
enddo
do i=1,nlast
l1=l1+1
write(9,*)'****************************MO=',l1,matorb(i)
if(momat1(1,i).eq.0.0.and.momat2(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat2(2,i).eq.0.0)then
if(momat1(3,i)*momat2(3,i).gt.0.0.and.momat1(3,i).ne.0.0)then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is sigma Anti-bonding orbital'
endif
if(momat1(3,i)*momat2(3,i).lt.0.0.and.momat1(3,i).ne.0.0)then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is sigma Bonding orbital'
endif
if(momat1(3,i)*momat2(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
ddi=ddi+1
ddimo(ddi)=l1
write(9,*)'Dimer MO=',l1,'is Delta orbital'
endif
if(momat1(3,i)*momat2(3,i).eq.0.0.and.momat1(4,i).ne.0.0)then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is Sigma orbital'
endif
else
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).gt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).gt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif

if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Anti-Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)).or.abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)).or.abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).eq.0.0.and.momat1(2,i)*momat2(2,i).lt.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is Pi Anti-Bonding orbital'
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Anti-Bonding orbital'
endif
endif
if(momat1(1,i)*momat2(1,i).lt.0.0.and.momat1(2,i)*momat2(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sdi=sdi+1
sdimo(sdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Sigma Anti-Bonding orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pdi=pdi+1
pdimo(pdi)=l1
write(9,*)'Dimer MO=',l1,'is deformed Pi Bonding orbital'
endif
endif
endif
enddo

endif
enddo
if(lognm.eq.'dime')goto 801
!!!! Orbitals of MONOMER 1

800 do mnnum=1,1
fg36=0
p1=0
!if (mnnum.eq.1)NMODI=NMOM1*2
!if (mnnum.eq.2)NMODI=NMOM2*2
if (logn.eq.'mono1')then
!logn=lognamem1
write(9,*)
write(9,*)'Molecular Orbitals for monomer 1'
write(9,*)
endif
if (logn.eq.'mono2')then
!logn=lognamem2
write(9,*)
write(9,*)'Molecular Orbitals for monomer 2'
write(9,*)
endif
INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=46,file=trim(inplogname),status='old')
MDP=0
do
read(46,'(a)',iostat=io)lines
!print*,lines(6:42)
if(lines(6:44).eq.'Alpha Molecular Orbital Coefficients:')di1=MDP
if(lines(6:36).eq.'Molecular Orbital Coefficients:')fg36=1
if(lines(6:36).eq.'Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(46)

else
print*,trim(inplogname),'file does not exist'
stop
endif

sm1=0
pm1=0
dm1=0
do i=1,1000
sm1mo(i)=0
pm1mo(i)=0
dm1mo(i)=0
enddo
!print*,di1
if(fg36.eq.0) then
write(9,*)'Alpha Molecular Orbitals'
write(9,*)
endif
if(fg36.eq.1) then
write(9,*)
endif
l=mod(NMODI,5)
nlast=l
k=(NMODI-l)/5
if(l.eq.0)kk=k
if(l.ne.0)kk=k+1
l1=0
l2=0
do i=1,di1+1
read(46,*)
enddo
do j=1,kk
do i=1,3
read(46,*)
enddo
do i=1,NMODI
read(46,'(a,a,a)')a,b,e
read(e,'(I10)')f1
if(f1.eq.nmon2)then
l1=l1+1
c(l1)=i
endif
enddo
enddo
rewind(46)

do i=1,di1+1
read(46,*)
enddo
do j=1,kk
l1=0
l2=0
l3=0
do i=1,3
read(46,*)
enddo
do i=1,c(j)
l3=l3+1
read(46,*)
enddo
do i=1,NMODI-c(j)
l3=l3+1
read(46,*)aa,bb
!print*,aa,bb
if(bb.eq.'2PX')then
l2=l2+1
cc(j,l2)=l3
endif
enddo
enddo
rewind(46)

do i=1,di1+1
read(46,*)
enddo
do j=1,k
j1=0
do i=1,6
do l=1,5
momat1(i,l)=0.0
enddo
enddo
read(46,*)
read(46,*)(matorb(i1),i1=1,5)
read(46,*)
do i=1,cc(j,1)-1
read(46,*)
enddo
do i=1,4
j1=j1+1
read(46,*)aa,bb,(momat1(j1,l),l=1,5)
!print*,(momat1(j1,l),l=1,5)
enddo
do i=1,(NMODI)-(cc(j,1)+3)
read(46,*)
enddo
do i=1,5
p1=p1+1
write(9,*)'************************ Monomer ',mnnum,',MO=',p1,matorb(i)
if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
if(momat1(4,i).eq.0.0)then
dm1=dm1+1
dm1mo(dm1)=p1
write(9,*)'Monomer MO=',p1,'is Delta orbital'
else
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0)then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is a Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is a Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif

enddo
enddo

do i=1,6
do j=1,5
momat1(i,j)=0.0
enddo
enddo

j1=0
j2=0
if(nlast.ne.0)then
read(46,*)
read(46,*)(matorb(i1),i1=1,nlast)
read(46,*)
do i=1,cc(1,1)-1
read(46,*)
enddo
do i=1,4
j1=j1+1
read(46,*)aa,bb,(momat1(j1,l),l=1,nlast)
!print*,(momat1(j1,l),l=1,nlast)
enddo
do i=1,(NMODI)-(cc(1,1)+3)
read(46,*)
enddo
do i=1,nlast
p1=p1+1
write(9,*)'************************ Monomer ',mnnum,'MO=',p1,matorb(i)

if(momat1(3,i).eq.0.0.and.momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0)then
if(momat1(4,i).eq.0.0)then
dm1=dm1+1
dm1mo(dm1)=p1
write(9,*)'Monomer MO=',p1,'is Delta orbital'
else
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0)then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
!if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
!sm1=sm1+1
!sm1mo(sm1)=p1
!write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
!endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif

enddo

endif

rewind(46)

if(fg36.eq.1)goto 801

MDP=0
do
read(46,'(a)',iostat=io)lines
!print*,lines(6:43)
if(lines(6:43).eq.'Beta Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(46)


!print*,di1
write(9,*)'Beta Molecular Orbitals'
write(9,*)
l=mod(NMODI,5)
nlast=l
k=(NMODI-l)/5
if(l.eq.0)kk=k
if(l.ne.0)kk=k+1
l1=0
l2=0
do i=1,di1+1
read(46,*)
enddo
do j=1,kk
do i=1,3
read(46,*)
enddo
do i=1,NMODI
read(46,'(a,a,a)')a,b,e
read(e,'(I10)')f1
if(f1.eq.nmon2)then
l1=l1+1
c(l1)=i
endif
enddo
enddo
rewind(46)

do i=1,di1+1
read(46,*)
enddo
do j=1,kk
l1=0
l2=0
l3=0
do i=1,3
read(46,*)
enddo
do i=1,c(j)
l3=l3+1
read(46,*)
enddo
do i=1,(NMODI)-c(j)
l3=l3+1
read(46,*)aa,bb
if(bb.eq.'2PX')then
l2=l2+1
cc(j,l2)=l3
endif
enddo
enddo
rewind(46)

l1=0
do i=1,di1+1
read(46,*)
enddo
do j=1,k
j1=0
j2=0
do i=1,3
do l=1,5
momat1(i,l)=0.0
enddo
enddo
read(46,*)
read(46,*)(matorb(i1),i1=1,5)
!print*,(matorb(i1),i1=1,5)
read(46,*)
do i=1,cc(j,1)-1
read(46,*)
enddo
do i=1,4
j1=j1+1
read(46,*)aa,bb,(momat1(j1,l),l=1,5)
!print*,(momat1(j1,l),l=1,5)
enddo
do i=1,NMODI-(cc(j,1)+3)
read(46,*)
enddo
do i=1,5
p1=p1+1
write(9,*)'************************ Monomer ',mnnum,'MO=',p1,matorb(i)
if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).eq.0.0)then
if(momat1(4,i).eq.0.0)then
dm1=dm1+1
dm1mo(dm1)=p1
write(9,*)'Monomer  MO=',p1,'is Delta orbital'
else
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer  MO=',p1,'is Sigma orbital'
endif
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0)then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
!if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
!sm1=sm1+1
!sm1mo(sm1)=p1
!write(9,*)'Monomer MO=',p1,'is Sigma orbital'
!endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
enddo
enddo
do i=1,3
do j=1,5
momat1(i,j)=0.0
enddo
enddo

j1=0
j2=0
if(nlast.ne.0)then
read(46,*)
read(46,*)(matorb(i1),i1=1,nlast)
!print*,(matorb(i1),i1=1,nlast)
read(46,*)
do i=1,cc(1,1)-1
read(46,*)
enddo
do i=1,4
j1=j1+1
read(46,*)aa,bb,(momat1(j1,l),l=1,nlast)
!print*,(momat1(j1,l),l=1,nlast)
enddo
do i=1,(NMODI)-(cc(1,1)+3)
read(46,*)
enddo
do i=1,nlast
p1=p1+1
write(9,*)'************************ Monomer ',mnnum,'MO=',p1,matorb(i)

!write(9,*)p1
!write(9,*)momat1(1,i),momat1(2,i),momat1(3,i)
if(momat1(3,i).eq.0.0.and.momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0)then
if(momat1(4,i).eq.0.0)then
write(9,*)'Monomer MO=',p1,'is Delta orbital'
else
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0)then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is Sigma orbital'
endif
!if(momat1(1,i).eq.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
!sm1=sm1+1
!sm1mo(sm1)=p1
!write(9,*)'Monomer MO=',p1,'is Sigma orbital'
!endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).eq.0.0.and.momat1(4,i).eq.0.0)then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is Pi orbital'
endif
if(momat1(1,i).eq.0.0.and.momat1(2,i).ne.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(2,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(2,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).eq.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif
if(momat1(1,i).ne.0.0.and.momat1(2,i).eq.0.0.and.momat1(3,i).ne.0.0.and.momat1(4,i).ne.0.0)then
if(abs(momat1(3,i)).gt.abs(momat1(1,i)))then
sm1=sm1+1
sm1mo(sm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Sigma orbital'
endif
if(abs(momat1(3,i)).lt.abs(momat1(1,i)))then
pm1=pm1+1
pm1mo(pm1)=p1
write(9,*)'Monomer MO=',p1,'is deformed Pi orbital'
endif
endif

enddo
endif
rewind(46)
if(mnnum.eq.2)then
sm2=sm1
pm2=pm1
dm2=dm1
do i=1,sm2
sm2mo(i)=sm1mo(i)
enddo
do i=1,pm2
pm2mo(i)=pm1mo(i)
enddo
do i=1,dm2
dm2mo(i)=dm1mo(i)
enddo

endif
enddo


801 return
end subroutine readorb_lin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readorb_nlin(inplogname)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,i,j,dialphae,m1alphae,m2alphae,alphaedi,betaedi,alphaem1,betaem1,alphaem2,betaem2,&
di1,l1,l2,k,kk,c(500),d(100),f1,cc(50,100),dd(50,100),l3,j1,j2,nlast,mnnum,i1,i2,i3,i4,i6,last,mateigen(6),ll&
,ro(3),j3,j5,lst,j6,j7,j8,ndate,loop,j10,fg35,fg36
real::momat1(5000,5),start,st,momat4(6)
character(len=50)::lines,XYZ1(5),XYZ2(5),XYZ3(5),inplogname,lognm
character(len=3)::a,b,e,na,matorb(6),th,nm
character(len=3)::cch1,cch2
character(len=5)::momat2(5000)
character(len=8)::fmt1,fmt2,fmt,f2,aa,bb,ans,nd
character(len=75)::logn
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
integer,dimension(8) :: values



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif

!print*,mno_at
lognamedi=inplogname

open(unit=9,file='Molecular_Orbitals_nlin.dat',status='unknown')

l=len(trim(inplogname))
lognm=inplogname(l-8:l-5)
logn=inplogname(l-8:l-4)
if(lognm.eq.'mono')goto 800

fg35=0
fg36=0

INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:36)
if(lines(6:36).eq.'Molecular Orbital Coefficients:')fg35=1
if(lines(6:44).eq.'Alpha Molecular Orbital Coefficients:')fg36=1
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif
if (fg35.eq.1)loop=1
if (fg36.eq.1)loop=2
do j10=1,loop
INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:36)
if (fg35.eq.1)then
if(lines(6:36).eq.'Molecular Orbital Coefficients:')di1=MDP
endif
if (fg36.eq.1.and.j10.eq.1)then
if(lines(6:43).eq.'Alpha Molecular Orbital Coefficients:')di1=MDP
endif
if (fg36.eq.1.and.j10.eq.2)then
if(lines(6:42).eq.'Beta Molecular Orbital Coefficients:')di1=MDP
endif
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif

write(9,*)'Molecular orbital file of',molname
write(9,*)
if (fg35.eq.1) write(9,*)'Molecular orbitals of Dimer'
if (fg36.eq.1.and.j10.eq.1) write(9,*)'Alpha Molecular orbitals of Dimer'
if (fg36.eq.1.and.j10.eq.2) write(9,*)'Beta Molecular orbitals of Dimer'

write(9,*)

j3=0
if(NMODI.ge.5)then
l=mod(NMODI,5)
nlast=l
k=(NMODI-l)/5
if(l.eq.0)kk=k
if(l.ne.0)kk=k+1
lst=5
endif
if(NMODI.lt.5)then
nlast=0
kk=1
k=1
lst=NMODI
endif
l1=0
l2=0
!print*,mno_at
do j1=1,mno_at
!print*,j1
do i=1,di1+1
read(45,*)
enddo
do j=1,kk
do i=1,3
read(45,*)
enddo
do i=1,NMODI
!read(45,'(a,a,a)')a,b,e
read(45,*)a,b,e
!print*,a,b,e

if(j1.lt.10)fmt1='(I1)'
if(j1.gt.9.and.j1.lt.100)fmt1='(I2)'
if(j1.gt.99.and.j1.lt.1000)fmt1='(I3)'
write(na,fmt1)j1
!print*,'na,nd',trim(na),nd,b
!read(e,'(I10)')f1
if(b.eq.na)then
l1=l1+1
c(l1)=i
!print*,'cccc',c(l1),i,l1,j
goto 430
endif
enddo
enddo
430 rewind(45)
enddo
ll=0
i6=0
do i=1,di1+1
read(45,*)
enddo
do j=1,k
l1=0
!print*,'sourav'
j1=0
j2=0
do i=1,5000
do l=1,5
momat1(i,l)=0.0
enddo
enddo
read(45,*)(mateigen(i1),i1=1,lst)
read(45,*)(matorb(i1),i1=1,lst)
!read(45,*)st,nd,(mateigen(i1),i1=1,5)
read(45,*)
!print*,(matorb(i1),i1=1,lst),(mateigen(i1),i1=1,lst)
do i1=1,mno_at
!print*,i1,mno_at
l1=l1+1
read(45,*)a,b,e,momat2(l1),(momat1(l1,i2),i2=1,lst)
!print*,l1,(momat1(l1,i2),i2=1,5)
last=c(i1+1)-1
if(i1.eq.mno_at)last=NMODI
!print*,'last',last
do i3=c(i1)+1,last
l1=l1+1
read(45,*)th,momat2(l1),(momat1(l1,i2),i2=1,lst)
!print*,l1,(momat1(l1,i2),i2=1,lst),momat1(l1,2)
enddo
enddo
do i1=1,lst
j3=j3+1
do i=1,6
momat4(6)=0.0
enddo
do i2=c(nacc)+1,c(nacc+1)-1
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
if(nm.eq.'PX')then
i6=i6+1
!print*,i2,i1,momat1(i2,i1)
momat4(1)=momat1(i2,i1)
endif
if(nm.eq.'PY')then
i6=i6+1
momat4(2)=momat1(i2,i1)
endif
if(nm.eq.'PZ')then
i6=i6+1
momat4(3)=momat1(i2,i1)
endif
enddo
!print*,ndon,c(ndon)
last=c(ndon+1)-1
if(ndon.eq.mno_at)last=NMODI
do i2=c(ndon)+1,last
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
!print*,nm
if(nm.eq.'PX')then
i6=i6+1
momat4(4)=momat1(i2,i1)
endif
if(nm.eq.'PY')then
i6=i6+1
momat4(5)=momat1(i2,i1)
endif
if(nm.eq.'PZ')then
i6=i6+1
momat4(6)=momat1(i2,i1)
endif
enddo
if(i6.lt.3)then
if(i6.eq.0)write(9,*)'All MO are Sigma orbitals'
goto 702
endif
!print*,momat4(1),momat4(2),momat4(3),momat4(4),momat4(5),momat4(6)
do i=1,3
l=0
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
if (ro(1).eq.1.and.ro(2).eq.0.and.ro(3).eq.2)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.1.and.ro(3).eq.2)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.2)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.1.and.ro(3).eq.0)then
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.2.and.ro(3).eq.0)then
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.0.and.ro(3).eq.0)then
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.2.and.ro(3).eq.0)then
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.0)then
write(9,*)'The MO is a DELTA orbital'
endif
!if (ro(1).eq.2.and.ro(2).eq.2.and.ro(3).ne.2)then
!write(9,*)'The MO is a DELTA orbital'
!endif
!if (ro(2).eq.2.and.ro(3).eq.2.and.ro(1).ne.2)then
!write(9,*)'The MO is a DELTA orbital'
!endif
!if (ro(1).eq.2.and.ro(3).eq.2.and.ro(2).ne.2)then
!write(9,*)'The MO is a DELTA orbital'
!endif
goto 759
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).ge.6.and.abs(momat4(1))/abs(momat4(3)).ge.6)then
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(1).eq.3)then
!if(abs(momat4(1))/abs(momat4(2)).lt.9.and.abs(momat4(1))/abs(momat4(3)).lt.9)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a  DELTA orbital'
endif
if(abs(momat4(1))/abs(momat4(2)).gt.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a  DELTA  orbital'
!endif
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).ge.6.and.abs(momat4(2))/abs(momat4(3)).ge.6)then
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(2).eq.3)then
!if(abs(momat4(2))/abs(momat4(1)).lt.9.and.abs(momat4(2))/abs(momat4(3)).lt.9)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(2))/abs(momat4(1)).gt.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif

if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(2)).ge.6.and.abs(momat4(3))/abs(momat4(1)).ge.6)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
endif
if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
endif
if(ro(3).eq.3)then
!if(abs(momat4(3))/abs(momat4(1)).lt.9.and.abs(momat4(3))/abs(momat4(2)).lt.9)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(3))/abs(momat4(1)).gt.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif


759 enddo
enddo
l1=0
do i=1,mno_at
do j=1,6
momat1(i,j)=0.0
enddo
enddo
if (nlast.ne.0)then
read(45,*)
read(45,*)(matorb(i1),i1=1,nlast)
read(45,*)
do i1=1,mno_at
!print*,i1,mno_at
l1=l1+1
read(45,*)a,b,e,momat2(l1),(momat1(i1,i2),i2=1,nlast)
!print*,a,b,e,momat2(l1),(momat1(l1,i2),i2=1,nlast)
last=c(i1+1)-1
if(i1.eq.mno_at)last=NMODI
do i3=c(i1)+1,last
l1=l1+1
read(45,*)th,momat2(l1),(momat1(i1,i2),i2=1,nlast)
!print*,momat2(l1),(momat1(l1,i2),i2=1,nlast)
enddo
enddo

do i1=1,nlast
j3=j3+1
do i=1,6
momat4(6)=0.0
enddo
do i2=c(nacc)+1,c(nacc+1)-1
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
if(nm.eq.'PX')then
!print*,i2,i1,momat1(i2,i1)
momat4(1)=momat1(i2,i1)
endif
if(nm.eq.'PY')then
momat4(2)=momat1(i2,i1)
endif
if(nm.eq.'PZ')then
momat4(3)=momat1(i2,i1)
endif
enddo
last=c(ndon+1)-1
if(ndon.eq.mno_at)last=NMODI
do i2=c(ndon)+1,last
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
if(nm.eq.'PX')then
momat4(4)=momat1(i2,i1)
endif
if(nm.eq.'PY')then
momat4(5)=momat1(i2,i1)
endif
if(nm.eq.'PZ')then
momat4(6)=momat1(i2,i1)
endif
enddo
!print*,momat4(1),momat4(2),momat4(3),momat4(4),momat4(5),momat4(6)
do i=1,3
l=0
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
if (ro(1).eq.1.and.ro(2).eq.0.and.ro(3).eq.2)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.1.and.ro(3).eq.2)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.2)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.1.and.ro(3).eq.0)then
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.2.and.ro(3).eq.0)then
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.0.and.ro(3).eq.0)then
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.2.and.ro(3).eq.0)then
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.0)then
write(9,*)'The MO is a DELTA orbital'
endif
goto 751
if (ro(1).eq.1.and.ro(2).eq.1.and.ro(3).ne.1)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
if (ro(2).eq.1.and.ro(3).eq.1.and.ro(1).ne.1)then
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.1.and.ro(3).eq.1.and.ro(2).ne.1)then
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.1.and.ro(3).eq.1)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.2.and.ro(3).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(2).eq.2.and.ro(3).eq.2.and.ro(1).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(1).eq.2.and.ro(3).eq.2.and.ro(2).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).ge.6.and.abs(momat4(1))/abs(momat4(3)).ge.6)then
if (momat4(1)*momat4(4).gt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(1)*momat4(4).lt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(1).eq.3)then
!if(abs(momat4(1))/abs(momat4(2)).lt.9.and.abs(momat4(1))/abs(momat4(3)).lt.9)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a  DELTA orbital'
endif
if(abs(momat4(1))/abs(momat4(2)).gt.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a  DELTA  orbital'
!endif
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).ge.6.and.abs(momat4(2))/abs(momat4(3)).ge.6)then
if (momat4(2)*momat4(5).gt.0.0)write(9,*)'The MO is a PI Bonding orbital'
if (momat4(2)*momat4(5).lt.0.0)write(9,*)'The MO is a PI Anti-Bonding orbital'
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(2).eq.3)then
!if(abs(momat4(2))/abs(momat4(1)).lt.9.and.abs(momat4(2))/abs(momat4(3)).lt.9)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(2))/abs(momat4(1)).gt.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif

if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(2)).ge.6.and.abs(momat4(3))/abs(momat4(1)).ge.6)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
endif
if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
if (momat4(3)*momat4(6).gt.0.0)write(9,*)'The MO is a SIGMA Anti-Bonding orbital'
if (momat4(3)*momat4(6).lt.0.0)write(9,*)'The MO is a SIGMA Bonding orbital'
endif
endif
if(ro(3).eq.3)then
!if(abs(momat4(3))/abs(momat4(1)).lt.9.and.abs(momat4(3))/abs(momat4(2)).lt.9)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(3))/abs(momat4(1)).gt.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif


751 enddo

endif
702 rewind(45)
enddo


if(lognm.eq.'dime')goto 801

!!! MONOMER ORBITALS

800 INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=46,file=trim(inplogname),status='old')
MDP=0
do
read(46,'(a)',iostat=io)lines
!print*,lines(6:42)
if(lines(6:36).eq.'Molecular Orbital Coefficients:')fg36=1
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(46)

else
print*,trim(inplogname),'file does not exist'
stop
endif

j6=2
if(fg36.eq.1)j6=1
do j5=1,j6
do i=1,6
momat4(i)=0
enddo
do i=1,3
ro(i)=0
enddo
do i=1,500
c(i)=0
d(i)=0
enddo
if(j6.eq.1)then
INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:42)
if(lines(6:36).eq.'Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif
endif


if(j5.eq.1.and.j6.eq.2)then
do i=1,6
momat4(i)=0
enddo
do i=1,3
ro(i)=0
enddo
do i=1,500
c(i)=0
d(i)=0
enddo
INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:42)
if(lines(6:42).eq.'Alpha Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif
endif

if(j5.eq.2.and.j6.eq.2)then
do i=1,6
momat4(i)=0
enddo
do i=1,3
ro(i)=0
enddo
do i=1,500
c(i)=0
d(i)=0
enddo

INQUIRE(file=trim(inplogname),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(inplogname),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:42)
if(lines(6:41).eq.'Beta Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(inplogname),'file does not exist'
stop
endif
endif

write(9,*)
write(9,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(9,*)
if(j5.eq.1.and.j6.eq.2.and.logn.eq.'mono1')write(9,*)'Alpha Molecular orbitals of mnomer1'
if(j5.eq.1.and.j6.eq.2.and.logn.eq.'mono2')write(9,*)'Alpha Molecular orbitals of mnomer2'
if(j5.eq.2.and.j6.eq.2.and.logn.eq.'mono1')write(9,*)'Beta Molecular orbitals of mnomer1'
if(j5.eq.2.and.j6.eq.2.and.logn.eq.'mono2')write(9,*)'Beta Molecular orbitals of mnomer2'
if(j6.eq.1.and.logn.eq.'mono1')write(9,*)'Molecular orbitals of mnomer1'
if(j6.eq.1.and.logn.eq.'mono2')write(9,*)'Molecular orbitals of mnomer2'
write(9,*)
if(logn.eq.'mono2')nmon1=nmon2
if(logn.eq.'mono2')mnm1=mnm2

j3=0
if(NMODI.ge.5)then
l=mod(NMODI,5)
lst=5
nlast=l
k=(NMODI-l)/5
if(l.eq.0)kk=k
if(l.ne.0)kk=k+1
endif
if(NMODI.lt.5)then
lst=NMODI
nlast=0
kk=1
k=1
endif
l1=0
l2=0
!print*,nlast,mno_at
do j1=1,mno_at
!print*,j1
do i=1,di1+1
read(45,*)
enddo
do j=1,kk
do i=1,3
read(45,*)
enddo
do i=1,NMODI
!read(45,'(a,a,a)')a,b,e
read(45,*)a,b,e
!print*,'***',b,nmon1,i,NMODI

if(j1.lt.10)fmt1='(I1)'
if(j1.gt.9.and.j1.lt.100)fmt1='(I2)'
if(j1.gt.99.and.j1.lt.1000)fmt1='(I3)'
write(na,fmt1)j1
!print*,'na,nd',trim(na),nd,b
!read(e,'(I10)')f1
if(b.eq.na)then
l1=l1+1
c(l1)=i
!print*,'c',c(l1),i,l1
goto 431
endif
enddo
enddo
431 rewind(45)
enddo
ll=0
do i=1,di1+1
read(45,*)
enddo
do j=1,k
l1=0
i6=0
!print*,'sourav'
j1=0
j2=0
do i=1,5000
do l=1,5
momat1(i,l)=0.0
enddo
enddo
read(45,*)(mateigen(i1),i1=1,lst)
!read(45,*)
read(45,*)(matorb(i1),i1=1,lst)
!read(45,*)
!read(45,*)st,nd,(mateigen(i1),i1=1,5)
!print*,(matorb(i1),i1=1,lst),(mateigen(i1),i1=1,lst)
read(45,*)
!print*,mnm1,mno_at
do i1=1,mno_at
l1=l1+1
read(45,*)a,b,e,momat2(l1),(momat1(l1,i2),i2=1,lst)
!print*,l1,(momat1(l1,i2),i2=1,lst)
last=c(i1+1)-1
if(i1.eq.mno_at)last=NMODI
!print*,'last',last
do i3=c(i1)+1,last
l1=l1+1
read(45,*)th,momat2(l1),(momat1(l1,i2),i2=1,lst)
!print*,l1,(momat1(l1,i2),i2=1,5),momat1(l1,2)
enddo
enddo
!print*,lst
do i1=1,lst
j3=j3+1
do i=1,6
momat4(6)=0.0
enddo
last=c(nmon1+1)-1
!print*,last
if(mno_at.eq.1)last=NMODI
do i2=c(nmon1)+1,last
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
!print*,'*****',nm
if(nm.eq.'PX')then
i6=i6+1
momat4(1)=momat1(i2,i1)
endif
if(nm.eq.'PY')then
i6=i6+1
momat4(2)=momat1(i2,i1)
endif
if(nm.eq.'PZ')then
i6=i6+1
momat4(3)=momat1(i2,i1)
endif
enddo
!write(9,*)momat4(1),momat4(2),momat4(3)
if(i6.lt.3)then
if(i6.eq.0)write(9,*)'*** All MO are Sigma orbitals ****'
goto 701
endif
do i=1,3
l=0
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
if (ro(1).eq.1.and.ro(2).eq.0.and.ro(3).eq.2)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.1.and.ro(3).eq.2)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.2)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.1.and.ro(3).eq.0)then
write(9,*)'The MO is a PI orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.2.and.ro(3).eq.0)then
write(9,*)'The MO is a PI  orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.0.and.ro(3).eq.0)then
write(9,*)'The MO is a PI  orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.2.and.ro(3).eq.0)then
write(9,*)'The MO is a PI  orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.0)then
write(9,*)'The MO is a DELTA orbital'
endif
goto 752
if (ro(1).eq.1.and.ro(2).eq.1.and.ro(3).ne.1)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(2).eq.1.and.ro(3).eq.1.and.ro(1).ne.1)then
write(9,*)'The MO is a PI orbital'
endif
if (ro(1).eq.1.and.ro(3).eq.1.and.ro(2).ne.1)then
write(9,*)'The MO is a PI orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.2.and.ro(3).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(2).eq.2.and.ro(3).eq.2.and.ro(1).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(1).eq.2.and.ro(3).eq.2.and.ro(2).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.1.and.ro(3).eq.1)then
write(9,*)'The MO is a DELTA orbital'
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).ge.6.and.abs(momat4(1))/abs(momat4(3)).ge.6)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(1).eq.3)then
!if(abs(momat4(1))/abs(momat4(2)).lt.9.and.abs(momat4(1))/abs(momat4(3)).lt.9)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a  DELTA orbital'
endif
if(abs(momat4(1))/abs(momat4(2)).gt.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a  DELTA  orbital'
!endif
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).ge.6.and.abs(momat4(2))/abs(momat4(3)).ge.6)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(2).eq.3)then
!if(abs(momat4(2))/abs(momat4(1)).lt.9.and.abs(momat4(2))/abs(momat4(3)).lt.9)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(2))/abs(momat4(1)).gt.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif

if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(2)).ge.6.and.abs(momat4(3))/abs(momat4(1)).ge.6)then
write(9,*)'The MO is a SIGMA orbital'
endif
endif
if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
write(9,*)'The MO is a SIGMA orbital'
endif
endif
if(ro(3).eq.3)then
!if(abs(momat4(3))/abs(momat4(1)).lt.9.and.abs(momat4(3))/abs(momat4(2)).lt.9)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(3))/abs(momat4(1)).gt.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif


752 enddo
enddo
l1=0
do i=1,5000
do j=1,6
momat1(i,j)=0.0
enddo
enddo
if (nlast.ne.0)then
read(45,*)
read(45,*)(matorb(i1),i1=1,nlast)
read(45,*)
do i1=1,mno_at
!print*,i1,mno_at
l1=l1+1
read(45,*)a,b,e,momat2(l1),(momat1(l1,i2),i2=1,nlast)
!print*,a,b,e,momat2(l1),(momat1(l1,i2),i2=1,nlast)
last=c(i1+1)-1
if(i1.eq.mno_at)last=NMODI
do i3=c(i1)+1,last
l1=l1+1
read(45,*)th,momat2(l1),(momat1(l1,i2),i2=1,nlast)
!print*,momat2(l1),(momat1(l1,i2),i2=1,nlast)
enddo
enddo

do i1=1,nlast
j3=j3+1
do i=1,6
momat4(6)=0.0
enddo

last=c(nmon1+1)-1
if(mno_at.eq.1)last=NMODI
do i2=c(nmon1)+1,last
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
if(nm.eq.'PX')then
!print*,i2,i1,momat1(i2,i1)
momat4(1)=momat1(i2,i1)
endif
if(nm.eq.'PY')then
momat4(2)=momat1(i2,i1)
endif
if(nm.eq.'PZ')then
momat4(3)=momat1(i2,i1)
endif
enddo
print*,momat4(1),momat4(2),momat4(3),momat4(4),momat4(5),momat4(6)
do i=1,3
l=0
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
if (ro(1).eq.1.and.ro(2).eq.0.and.ro(3).eq.2)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.1.and.ro(3).eq.2)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.2)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.1.and.ro(3).eq.0)then
write(9,*)'The MO is a PI orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.2.and.ro(3).eq.0)then
write(9,*)'The MO is a PI  orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.0.and.ro(3).eq.0)then
write(9,*)'The MO is a PI  orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.2.and.ro(3).eq.0)then
write(9,*)'The MO is a PI  orbital'
endif
if (ro(1).eq.0.and.ro(2).eq.0.and.ro(3).eq.0)then
write(9,*)'The MO is a DELTA orbital'
endif
goto 753
if (ro(1).eq.1.and.ro(2).eq.1.and.ro(3).ne.1)then
write(9,*)'The MO is a SIGMA orbital'
endif
if (ro(2).eq.1.and.ro(3).eq.1.and.ro(1).ne.1)then
write(9,*)'The MO is a PI orbital'
endif
if (ro(1).eq.1.and.ro(3).eq.1.and.ro(2).ne.1)then
write(9,*)'The MO is a PI orbital'
endif
if (ro(1).eq.2.and.ro(2).eq.2.and.ro(3).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(2).eq.2.and.ro(3).eq.2.and.ro(1).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(1).eq.2.and.ro(3).eq.2.and.ro(2).ne.2)then
write(9,*)'The MO is a DELTA orbital'
endif
if (ro(1).eq.1.and.ro(2).eq.1.and.ro(3).eq.1)then
write(9,*)'The MO is a DELTA orbital'
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).ge.6.and.abs(momat4(1))/abs(momat4(3)).ge.6)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(1).eq.3)then
if(abs(momat4(1))/abs(momat4(2)).le.5.and.abs(momat4(1))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a  DELTA orbital'
endif
if(abs(momat4(1))/abs(momat4(2)).gt.5.and.abs(momat4(1))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a  DELTA  orbital'
!endif
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).ge.6.and.abs(momat4(2))/abs(momat4(3)).ge.6)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(2).eq.3)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a PI orbital'
endif
endif
if(ro(2).eq.3)then
!if(abs(momat4(2))/abs(momat4(1)).lt.9.and.abs(momat4(2))/abs(momat4(3)).lt.9)then
if(abs(momat4(2))/abs(momat4(1)).le.5.and.abs(momat4(2))/abs(momat4(3)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(2))/abs(momat4(1)).gt.5.and.abs(momat4(2))/abs(momat4(3)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif

if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(2)).ge.6.and.abs(momat4(3))/abs(momat4(1)).ge.6)then
write(9,*)'The MO is a SIGMA orbital'
endif
endif
if(ro(3).eq.3)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
write(9,*)'The MO is a SIGMA orbital'
endif
endif
if(ro(3).eq.3)then
!if(abs(momat4(3))/abs(momat4(1)).lt.9.and.abs(momat4(3))/abs(momat4(2)).lt.9)then
if(abs(momat4(3))/abs(momat4(1)).le.5.and.abs(momat4(3))/abs(momat4(2)).gt.5)then
write(9,*)'The MO is a distorted DELTA orbital'
endif
if(abs(momat4(3))/abs(momat4(1)).gt.5.and.abs(momat4(3))/abs(momat4(2)).le.5)then
write(9,*)'The MO is a distorted DELTA  orbital'
!endif
endif
endif


753 enddo

endif
701 rewind(45)
enddo
801 return
end subroutine readorb_nlin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
character(len=70)::infnam,str1,line5,line6,line8,line9,line10,line11,line12,line14,runn,line15
character(len=12)::abc,sttr(100),str2,charst(100),bochtrns,lpe_mono,lin
character(len=10)::string
character(len=18)::lettstr(40)
!character(len=:), allocatable::lettstr(:)
integer::i,MDP,io,l,j,m,k,ll,kk,jj,stn,stn1,a,b,p,q,t,i3,i4,i2,i5,c,d,m5,m6
integer::st_num1(100),st_num2(100),st_num3(100),st_num4(100),number(500),ndate
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
character(len=1)  :: t1, t2, t3
character(len=110)  :: line34
integer,dimension(8) :: values



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif




!open(unit=21,file='input.dat',status='unknown')


MDP=0
i3=0
i4=0
a=0
b=0
fg9=0
fg10=0
fg11=0
fg12=0
fg13=0
fg14=0
fg15=0
fg16=0
fg17=0
fg18=0

noorbdi=0
noorbm1=0
noorbm2=0
do k=1,350
orbdi(k)=0
orbm1(k)=0
orbm2(k)=0
enddo

do
read(21,'(a)',iostat=io)

if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(21)

do i5=1,MDP
read(21,'(a)')charst(i5)
enddo
do i5=1,MDP
ll=len(trim(charst(i5)))
do k=1,ll+1
line8=charst(i5)(k:k)
if(line8.eq.'')then
i4=k
goto 433
endif
enddo
433 a=1
 b=i4
sttr(i5)=charst(i5)(a:b)
!print*,sttr(i5),charst(i5)(a:b)
enddo

rewind(21)

301 do i=1,MDP
kk=0
q=0
do j=1,100
st_num1(j)=0
st_num2(j)=0
st_num3(j)=0
st_num4(j)=0
enddo
do j=1,500
number(j)=0
enddo
str2=sttr(i)
if(str2.eq.'oct_di')then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 439
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 440
439 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
440 enddo
stn=0
do k=1,500
number(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(fg10.eq.1)goto 389
if(line34(a:b).ne.'oct_di')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')number(stn)
!print*,number(stn),stn
endif
endif
endif
389 enddo
m5=0
do k=1,stn
if(number(k).eq.0)then
do m6=number(k-1)+1,number(k+1)-1
m5=m5+1
orbdi(m5)=m6
enddo
endif
if(number(k).ne.0)then
m5=m5+1
orbdi(m5)=number(k)
endif
enddo
noorbdi=m5
!print*,noorbdi,(orbdi(k),k=1,noorbdi)
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'oct_m1')then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line36=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 539
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 540
539 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
540 enddo
stn=0
do k=1,500
number(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(fg10.eq.1)goto 388
if(line34(a:b).ne.'oct_m1')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')number(stn)
!print*,number(stn),stn
endif
endif
endif
388 enddo
m5=0
do k=1,stn
if(number(k).eq.0)then
do m6=number(k-1)+1,number(k+1)-1
m5=m5+1
orbm1(m5)=m6
enddo
endif
if(number(k).ne.0)then
m5=m5+1
orbm1(m5)=number(k)
endif
enddo
noorbm1=m5
!print*,noorbm1,(orbm1(k),k=1,noorbm1)
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'oct_m2')then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line37=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 639
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 640
639 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
640 enddo
stn=0
do k=1,500
number(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(fg10.eq.1)goto 387
if(line34(a:b).ne.'oct_m2')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')number(stn)
!print*,number(stn),stn
endif
endif
endif
387 enddo
m5=0
do k=1,stn
if(number(k).eq.0)then
do m6=number(k-1)+1,number(k+1)-1
m5=m5+1
orbm2(m5)=m6
enddo
endif
if(number(k).ne.0)then
m5=m5+1
orbm2(m5)=number(k)
endif
enddo
noorbm2=m5
!print*,noorbm2,(orbm2(k),k=1,noorbm2)
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'IhaveMOs')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,lpe_mono
if(lpe_mono.eq.'yes')fg9=1
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'cent_at_no')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,nacc,ndon,nmon1,nmon2
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'mol')then
stn=0
stn1=0
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 739
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 740
739 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
740 enddo
stn=0
do k=1,500
number(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).eq.'l:'.or.line34(a:b).eq.'nl:')then
stn1=stn1+1
if(line34(a:b).eq.'l:')fg6=1
if(line34(a:b).eq.'nl:')fg6=2
lopt(stn1)=fg6
endif
if(line34(a:b).ne.'mol')then
if(line34(a:b).ne.'l:')then
if(line34(a:b).ne.'nl:')then
stn=stn+1
inp_log(stn)=line34(a:b)
!print*,number(stn),stn
endif
endif
endif
endif
enddo
noinp=stn
rewind(21)
endif

enddo

return
end subroutine readinp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prepname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!preparation of name of output file
use commondat
implicit none
integer::l,ndate
character(len=35):: name,mname1, mname2,mname3,mname4,mname5,mname6,mname7,mname8,diname,&
m1name,m2name, filename,mname9,mname10
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
integer,dimension(8) :: values



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif

l=len(TRIM(infname))
name=infname(1:l-10)
molname=name
mname1='_mono1.com'
mname2='_mono2.com'
mname3='_dimer.com'
mname4='.cubgrid'
mname5='.chk'
mname6='.fchk'
mname7='.cub'
mname8='.log'
mname9='.out'
mname10='dens_sub.cub'
logname=trim(name)//trim(mname8)
mononame1=trim(name)//trim(mname1)
mononame2=trim(name)//trim(mname2)
dimname=trim(name)//trim(mname3)
cubname=trim(name)//trim(mname4)
integral=trim(name)//trim(mname9)
substract=trim(name)//trim(mname10)
substract_oct=trim(name)//trim('_mo_sub.cub')
dimer_oct=trim(name)//trim('_dimer_mo.cub')
mono1_oct=trim(name)//trim('_mono1_mo.cub')
mono2_oct=trim(name)//trim('_mono2_mo.cub')
l=len(TRIM(mononame1))
m1name=mononame1(1:l-4)
l=len(TRIM(mononame2))
m2name=mononame2(1:l-4)
l=len(TRIM(dimname))
diname=dimname(1:l-4)
chknamedi=trim(diname)//trim(mname5)
chknamem1=trim(m1name)//trim(mname5)
chknamem2=trim(m2name)//trim(mname5)
fchknamedi=trim(diname)//trim(mname6)
fchknamem1=trim(m1name)//trim(mname6)
fchknamem2=trim(m2name)//trim(mname6)
cubnamedi=trim(diname)//trim(mname7)
cubnamem1=trim(m1name)//trim(mname7)
cubnamem2=trim(m2name)//trim(mname7)
lognamedi=trim(diname)//trim(mname8)
lognamem1=trim(m1name)//trim(mname8)
lognamem2=trim(m2name)//trim(mname8)
cubnamemodi=trim(diname)//trim('_mo.cub')
cubnamemom1=trim(m1name)//trim('_mo.cub')
cubnamemom2=trim(m2name)//trim('_mo.cub')

return
end subroutine prepname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine molorb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none


!real::start,finish,time,dival(NMODI,MNXY*MNZ),m1val(NMOM1,MNXY*MNZ),m2val(NMOM2,MNXY*MNZ),grd1(5)&
!m1vall(NMOM1*MNXY*MNZ),m2vall(NMOM2*MNXY*MNZ),
real::start,finish,time,divall(1*MNXY*6),dival(MNXY*MNZ),m1val(MNXY*MNZ),m2val(MNXY*MNZ),grd1(5)&
,grd2(3,4),at_grd(10,300),divalsq(MNXY*MNZ),m1valsq(MNXY*MNZ),char1,char2,char3,char4,&
m2valsq(MNXY*MNZ),divallsq(MNXY*MNZ),di,m1,m2,at_grdd(200,6),noz,int,norm_gr,sumxy(MNZ),ox,oy,oz
logical::file1
integer::io,MDP,n1,n2,l,l1,l2,l3,i,k,x,mdp1(3000),nline,nlast,nlstlin,int1,int2,int3(3),ll,noe,norbmo,norbdi,j,tndp,ndate&
,uplim,i1,galp,i3,i2,li,lin,intd1,noed(200),noem1(100),noem2(100),nx,ny,nz
character(75)::line,name,opt,cubnm,integral_mo,namo,substract_oct_new,substract_oct1
!character(len=3)::name
character(len=8)::fmt
character(len=8)  :: date
character(len=5)  :: zone
character(len=150)  :: ln10
integer,dimension(8) :: values
character(len=200)::cubndi,cubnm1,cubnm2


fmt='(I3.3)'
open(unit=43,file='script1',status='unknown')
do i=1,noorbdi
x=orbdi(i)
write(name,fmt)x
cubndi=trim('MO_dim_')//trim(name)//trim('.cub')
write(43,111)trim(cubcommo),trim(name),trim(fchknamedi),trim(cubndi),trim(cubname)
enddo
open(unit=53,file='script2',status='unknown')
do i=1,noorbm1
write(name,fmt)orbm1(i)
cubnm1=trim('MO_mon1_')//trim(name)//trim('.cub')
write(53,111)trim(cubcommo),trim(name),trim(fchknamem1),trim(cubnm1),trim(cubname)
enddo
open(unit=63,file='script3',status='unknown')
do i=1,noorbm2
write(name,fmt)orbm2(i)
cubnm2=trim('MO_mon2_')//trim(name)//trim('.cub')
write(63,111)trim(cubcommo),trim(name),trim(fchknamem2),trim(cubnm2),trim(cubname)
enddo
111 Format (a,a,1x,a,1x,a,1x,'-1',1x,'h',1x,a)

close(43)
close(53)
close(63)
call systemqq ("chmod +x script1")
call systemqq ("chmod +x script2")
call systemqq ("chmod +x script3")
CALL SYSTEMQQ ("./script1")
CALL SYSTEMQQ ("./script2")
CALL SYSTEMQQ ("./script3")
write(*,101)
101 Format(//'********************************************************************************')
write(*,100)trim('MO files of dimer'),trim('monomer1'),trim('and monomer2')
100 Format( a,1x,a,1x,a,' are preparing')
write(*,*)'********************************************************************************'

!do i=1,noorbdi
i=101
write(name,fmt)orbdi(noorbdi)
cubndi=trim('MO_dim_')//trim(name)//trim('.cub')
448 INQUIRE(file=trim(cubndi),EXIST=file1)
IF (file1) THEN

k=0

open(unit=i,file=trim(cubndi),status='old')
447 MDP=0
k=k+1
time=time+10
do
read(i,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(i)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 548
202 call cpu_time(start)
if(start.lt.time)goto 202

goto 447
else
goto 448
endif
548 close(i)
!enddo

!do i=1,noorbm1
i=103
write(name,fmt)orbm1(noorbm1)
cubnm1=trim('MO_mon1_')//trim(name)//trim('.cub')
348 INQUIRE(file=trim(cubnm1),EXIST=file1)
IF (file1) THEN

k=0

open(unit=i,file=trim(cubnm1),status='old')
347 MDP=0
k=k+1
time=time+10
do
read(i,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(i)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 648
302 call cpu_time(start)
if(start.lt.time)goto 302

goto 347
else
goto 348
endif
648 close(i)
!enddo

!do i=1,noorbm2
i=105
write(name,fmt)orbm2(noorbm2)
cubnm2=trim('MO_mon2_')//trim(name)//trim('.cub')
248 INQUIRE(file=trim(cubnm2),EXIST=file1)
IF (file1) THEN

k=0

open(unit=i,file=trim(cubnm2),status='old')
247 MDP=0
k=k+1
time=time+10
do
read(i,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(i)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 449
402 call cpu_time(start)
if(start.lt.time)goto 402

goto 247
else
goto 248
endif

449 close(i)
!enddo


!print*,'sourav2'

!print*,'sourav2'
CALL SYSTEMQQ ("rm script1")
CALL SYSTEMQQ ("rm script2")
CALL SYSTEMQQ ("rm script3")
CALL SYSTEMQQ ("mkdir dump_MO")
CALL SYSTEMQQ ("mv MO* dump_MO/")

return
end subroutine molorb

