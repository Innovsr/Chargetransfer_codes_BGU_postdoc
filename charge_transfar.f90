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
!!! ************************************Format of the inputfile**************************************
!****************************************************************************************************
!acc_don_thr  1 5 2                    # Sl. Number. of acceptor doner and the
!                                        third atom | if U don't have any third atom out side
!                                        of Z axis put '0' at it,s place
!dimer_ch_sp  0 1                      # Charge and Multiplicity of the dimer
!mono1_ch_sp  0 1                      # Charge and Multiplicity of the monomer1
!mono2_ch_sp  0 1                      # Charge and Multiplicity of the monomer1
!inp          NH3-BH3.com              # input file to feed the program molecular geometry | input file
!                                        can be '.com' or '.log' or '.xyz' file
!opt_inp      yes                      # optimization of the input geometry | If
!                                        U wish to optimize the input '.com' or '.xyz'
!                                        file then put 'yes' and program will
!                                        take the optimized geometry from output log file.
!                                        otherwise put 'no'. If your input file
!                                         is '.log' file it will be overlooked.
!mono_option  distance                 # Procedure of selecting monomer atoms |
!                                        There are three possibilities 'distance', 'radius'
!                                        and 'user'
!mono1_atom   4 1 2 3 4                # If U put 'user' option in 'mono_option'
!                                        then U have to put the sl. number. of the atoms
!mono2_atom   4 5 6 7 8                  1st number will be the number of atoms
!                                        and then the sl. numbers
!
!***********************************   # starting line of header and lower part
!                                        of the dimer and monomers
!%chk=                                 # dont put any name of the '.chk' file |
!                                        it will be automatically generated
!%nproc=8                              # it will be directly copy and paste in
!                                        the dimer and monomers file
!%mem=50gb                             # it will be directly copy and paste in
!                                        the dimer and monomers file
!#p  pbepbe 6-31g*   opt freq          # it will be copy paste except 'opt' and
!                                        'freq' | 'nosym' and 'pop=full'
!                                        will be added automatically. dont put
!                                        any gap between lines in this section
!m1 0.0 0.0 0.50 -0.0000001            # if U put 'genecp', 'gen' and 'charge' in
!                                        the command line then this lines will be copied and
!di 0.0 0.0 0.0 -0.0000000               pasted | di for dimer .com file, m1 for
!                                        monomer1 .com file and m2 for monomer2
!m2 0.0 0.0 -0.50 -0.0000001
!***********************************   # end of the section
!run_gauss    g09                      # command to run gaussian
!run_cube     cubegen 1 density=scf    # command to run cubgen
!grd_space  0.1                        # grid spacing depending which cubegrid
!                                        will prepare
!run  simul                            # command to run gaussian simultaniously
!run_fchk formchk                      # command to run formchk
!
!****************************************************************************************************
!****************************************************************************************************
!****************************************************************************************************

module commondat
implicit none

real,public::dongroup,accgroup,at_rad
real::dummy,min_grd,grd_sp,ocdi,ocm1,ocm2
integer::mono1,mono2,at_covrad,mno_at,nacc,ndon,third,fg9,fg8,fg7,fg6,fg5,fg3,fg2,fg1,mnopt,mon1,mon2&
,mnm1,mnm2,steps,MNXY,MNZ,MN,fg4,noorbdi,noorbm1,noorbm2,NMODI,NMOM1,NMOM2,orbdi(350),orbm1(350),orbm2(350)&
,fg10,fg11,fg12,fg13,fg14,fg15,fg16,fg17,fg18,flg20,flg21,flg22
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
mono1_oct,mono2_oct
common /input/line1,line2,line3,wordstr,fchk_com,line27,line28,line29,line35,line36,line37
common /mono/mono1,mono2,MNXY,MNZ,MN,NMODI,NMOM1,NMOM2
common /group/dongroup(300,3),accgroup(300,3),mon1(200),mon2(200),orbdi,orbm1,orbm2
common /dat/mno_at,nacc,ndon,third,fg9,fg8,fg7,fg6,fg5,fg4,fg3,fg2,fg1,mnopt,mnm1,mnm2,noorbdi,noorbm1,noorbm2&
,fg10,fg11,fg12,fg13,fg14,fg15,fg16,fg17,fg18,flg20,flg21,flg22

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
character(len=2) :: atomin,natom
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


call prepname
!if(fg7.eq.1)goto 402
iname=infname
write(*,104)trim(iname)
104 Format(// a,' is taken as your initial gaussian input file')
l=len(TRIM(iname))
lname=iname(l-3:l)

if(lname.eq.'.com'.and.fg2.eq.2) then 
call runcom(iname,logname,gauss)
name='.log'
name1=iname(1:l-4)
name2=trim(name1)//trim(name)
call readlog(name2) 
endif
if(lname.eq.'.com'.and.fg2.ne.2) then
call readcom(iname)
endif


if(lname.eq.'.log') call readlog(iname) 
if(lname.eq.'.xyz') call readcom(iname) 

call shiftmol
call rottoz
if(fg1.ne.2)then
call output1
endif
call output2(chknamedi,chknamem1,chknamem2)
call cubegrid
!if(fg7.eq.1)goto 400
if(fg3.eq.1) then
print*,'***********************************************'
call runcom(dimname,lognamedi,gauss)
print*,'***********************************************'
call runcom(mononame1,lognamem1,gauss)
print*,'***********************************************'
call runcom(mononame2,lognamem2,gauss)
print*,'***********************************************'
endif
print*,'***********************************************'
if(fg3.eq.2)then
call runcom_sim(gauss)
endif

print*,'***********************************************'
if(fg3.eq.1) then
call genformchk(chknamedi,fchknamedi,NMODI)
print*,'***********************************************'
call genformchk(chknamem1,fchknamem1,NMOM1)
print*,'***********************************************'
call genformchk(chknamem2,fchknamem2,NMOM2)
endif
print*,'***********************************************'
!goto 400

if(fg3.eq.2)then
call genformchk_sim
endif
print*,'***********************************************'
202 call cpu_time(start)
if(start.lt.30)goto 202

if(fg3.eq.1) then
call cubgen(fchknamedi,cubnamedi,cubname,cubcom)
print*,'***********************************************'
call cubgen(fchknamem1,cubnamem1,cubname,cubcom)
print*,'***********************************************'
call cubgen(fchknamem2,cubnamem2,cubname,cubcom)
print*,'***********************************************'
endif
if(fg3.eq.2)then
call cubgen_sim
endif
call intsum
!400 if(fg6.eq.1) call readfchk
!if(fg6.eq.2) call readfchk_cplx

!400 print*,'sourav1'
!if (fg4.eq.1)then
!call molorb
!endif

!print*,'sourav2'
!400 stop
stop
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prepname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!preparation of name of output file
use commondat
implicit none
integer::l,ndate
character(len=35):: name,mname1, mname2, mname3,mname4,mname5,mname6,mname7,mname8,diname,&
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
name=infname(1:l-4)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readcom(iname)
!reading .com file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none 
character(len=2)::sttr(300)
character(len=8)::line4(300),line5
character(len=35)::iname
integer::i3,fl,i,k,j,MDP,io,ii,kk,jj,l,m,i2,p,q,ll,a,b,stn
integer::st_num1(100),st_num2(100),ndate
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


do i=1,mno_at
do j=1,3
dummy(300,3)=0.0
enddo
enddo
open(unit=19,file=iname,status='old')
j=1
kk=0
m=1
l=0
fl=0
i3=0
q=0

238 if (j.ne.1)then
do i=1,j-1
read(19,*)
enddo 
endif
do ii=j,j+1
read(19,'(a)')line4(ii)
if(len(trim(line4(ii))).eq.0.and.ii.le.3)goto 600

ll=len(trim(line4(ii)))
do k=1,ll+1
line5=line4(ii)(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 439
!print*,k
if(q.eq.1)then
i3=1
st_num1(1)=1
endif
i3=i3+1
st_num1(i3)=k+1
st_num2(i3)=k
goto 440
439 i3=i3+1
st_num1(i3)=k+1
st_num2(i3)=k
endif
440 enddo
stn=0
do k=1,i3
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
stn=stn+1
 a=st_num1(k)
 b=st_num2(k+1)
!print*,a,b,stn,st_num(k)
do i2=1,88
if(line4(ii)(a:b).eq.at_list(i2))then
sttr(ii)=line4(ii)(a:b)
endif
enddo
endif
enddo


600 enddo
rewind(19)
do k=1,88
if(sttr(j).eq.at_list(k))goto 233
enddo
j=j+1
goto 238
233 do k=1,88
if(sttr(j+1).eq.at_list(k))m=m+1
if(sttr(j+1).eq.at_list(k).and.m.eq.2)kk=j
if(sttr(j+1).eq.at_list(k))goto 235
enddo
goto 236
235 j=j+1
fl=1
goto 238
236 do i=1,kk-1
read(19,*)
enddo
do i=kk,kk+m-1
l=l+1
read(19,*)stt(l),(dummy(l,k),k=1,3)
enddo
do i=kk+m,MDP-1
read(19,*)
enddo

mno_at=l

if(mnm1+mnm2.ne.mno_at.and.fg1.eq.2)then
print*,'there is inconsistancy in two monomer numbers in the inputfile'
print*,'there should have total ',mno_at,'no of atoms'
print*,'where here you put ',mnm1+mnm2,' atoms: ',mnm1,'for monomer1 and',mnm2,'for monomer2'
stop
endif
return
end subroutine readcom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinp
!reading inputfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
character(len=70)::infnam,str1,line4,line5,line6,line8,line9,line10,line11,line12,line14,runn,line15
character(len=12)::abc,sttr(100),str2,charst(100),bochtrns,lpe_mono,lin
character(len=10)::string
character(len=18)::lettstr(40)
character(len=90)::line26
character(len=2)::f
!character(len=:), allocatable::lettstr(:)
integer::i,MDP,io,l,j,m,k,ll,kk,jj,stn,a,b,p,q,t,i3,i4,i2,i5,i6,c,d
integer::st_num1(100),st_num2(100),st_num3(100),st_num4(100),ndate
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


ch0='0'
sp0='1'
ch1='-1'
ch2='1'
sp1='1'
sp2='1'
fg1=0
fg2=2
fg3=1
fg4=0
fg5=1
fg5=1
fg7=0
fg8=0
flg20=0
flg21=0
flg22=0
grd_sp=0.1
gauss='g09'
cubcom='cubegen 1 density=scf'
cubcommo='cubegen 1 MO=All'
fchk_com='formchk'
line27='     '
line28='     '
line29='     '

kk=0
MDP=0
q=0
i3=0
i4=0
a=0
b=0

do
read(21,'(a)',iostat=io)

if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(21)

do i5=1,MDP
read(21,'(a)')charst(i5)
!print*,charst(i5)
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
do i=1,MDP
str2=sttr(i)
if(str2.eq.'inp')goto 301
enddo
print*,'There is no gaussina input or output file exist in the input file'
stop

301 do i=1,MDP
str2=sttr(i)
!print*,str2
if(str2.eq.'inp')then
do j=1,i-1
read(21,*)
enddo
read(21,'(a,a)')abc,infname
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
if(infname.eq.'')then
print*,'There is no gaussina input or output file exist in the input file'
stop
endif
endif
if(str2.eq.'acc_don_thr')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,nacc,ndon,third
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif


if(str2(1:1).eq.'*')then
do j=i+1,MDP
if(sttr(j)(1:1).eq.'*')then
str1=sttr(j)(1:1)
goto 573
endif
enddo
goto 574
573 if(str1.eq.'*')then
do j=1,i
read(21,*)
enddo
read(21,*)line1
read(21,*)line2
read(21,*)line3
read(21,'(a)')line4
!print*,line4
ll=len(trim(line4))
do k=1,ll+1
line5=line4(k:k)
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
line11=''
line9=' '
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
stn=stn+1
 a=st_num1(k)
 b=st_num2(k+1)
if(line4(a:b).eq.'genecp')fg8=1
if(line4(a:b).eq.'charge')fg8=1
!print*,'fg8',fg8
if(line4(a:b).ne.'opt')then
if(line4(a:b).ne.'freq')then
if(line4(a:b).ne.'nosymm')then
if(line4(a:b).ne.'#p')then
lettstr(stn)=line4(a:b)
line10=trim(line11)//' '//trim(line4(a:b))
line11=line10
endif
endif
endif
endif
endif
enddo
line12='nosymm'
line14='#p '
line15='Pop=Full'
line13=trim(line14)//' '//trim(line10)//' '//trim(line12)//' '//trim(line15)
if(fg8.eq.1)then
358 read(21,'(a,a)')f,line26
if(trim(f).eq.'**')goto 957
if(f.eq.'di')line27=line26
if(f.eq.'di')flg20=1
if(f.eq.'m1')line28=line26
if(f.eq.'m1')flg21=1
if(f.eq.'m2')line29=line26
if(f.eq.'m2')flg22=1
goto 358
957 endif
rewind(21)
endif
endif

if(str2.eq.'dimer_ch_sp')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,ch0,sp0
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'mono1_ch_sp')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,ch1,sp1
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'mono2_ch_sp')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,ch2,sp2
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'opt_inp')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,string
if(string.eq.'yes')fg2=2
if(string.ne.'yes')fg2=0
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'mono_option')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,str1
if(str1.eq.'radius')fg1=1
if(str1.eq.'distace')fg1=0
if(str1.eq.'user')fg1=2
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'mono1_atom'.and.fg1.eq.2)then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,mnm1,(mon1(k),k=1,mnm1)
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'mono2_atom'.and.fg1.eq.2)then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,mnm2,(mon2(k),k=1,mnm2)
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'run_gauss')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,gauss
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
if(str2.eq.'run_cube')then
do j=1,i-1
read(21,*)
enddo
read(21,'(a,a)')abc,cubcom
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'grd_space')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,grd_sp
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif


if(str2.eq.'run')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,runn
if(runn.eq.'simul')fg3=2
if(runn.eq.'sequen')fg3=1
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'run_fchk')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,fchk_com
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'orb_file')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,bochtrns
if(trim(bochtrns).eq.'yes')fg4=1
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'mol')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,lin
if(trim(lin).eq.'lin')fg6=1
if(trim(lin).eq.'nlin')fg6=2
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'cal_part')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,bochtrns
if(trim(bochtrns).eq.'full')fg7=0
if(trim(bochtrns).eq.'mocto')fg7=1
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif
574 enddo

return
end subroutine readinp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readlog(iname)
!reading of .log file for molecular coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none 
integer ::  i, j, k, nargs ,num_args,p,ndate, jj, kk,&
lenth, l,  io, MDP, ii,iii,cn,an,at
character(len=56)::lines
character(len=35)::iname
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

      


do i=1,mno_at
do j=1,3
dummy(300,3)=0.0
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
open(unit=25,file='temp_3',status='unknown')
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

if(mnm1+mnm2.ne.mno_at.and.fg1.eq.2)then
print*,'there is inconsistancy in two monomer numbers in the inputfile'
print*,'there should have total ',mno_at,'no of atoms'
print*,'where here you put ',mnm1+mnm2,' atoms: ',mnm1,'for monomer1 and',mnm2,'for monomer2'
stop
endif
return
end subroutine readlog


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shiftmol
!!!shift the molecule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none 
integer :: k,i,ndate
real :: dnaccx, dnaccy, dnaccz
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




dnaccx=dummy(nacc,1)
dnaccy=dummy(nacc,2)
dnaccz=dummy(nacc,3)
do i=1,mno_at
dummy(i,1)=dummy(i,1)-dnaccx
dummy(i,2)=dummy(i,2)-dnaccy
dummy(i,3)=dummy(i,3)-dnaccz
enddo

return
end subroutine shiftmol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rottoz
!rotation of molecule to the z axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none 
!integer :: nacc,ndon,third,mno_at
real:: x,y,z,amoddon,angz,pervx,pervy,pervz,u,v,w,a,b,c,al,angzx,amidx,amidy,amidz &
,i,ndate
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


if(dummy(nacc,1).eq.0.0.and.dummy(nacc,2).eq.0.0.and.dummy(ndon,1).eq.0.0.and.dummy(ndon,2).eq.0.0)&
goto 600

!cc  calculate the perpendicular vector of the plane containing
!cc  z axis and the molecule
if(dummy(ndon,3).gt.0.0)then
x=0.0
y=0.0
z=1.0
pervx=y*dummy(ndon,3)-z*dummy(ndon,2)
pervy=-x*dummy(ndon,3)+z*dummy(ndon,1)
pervz=x*dummy(ndon,2)-y*dummy(ndon,1)

!cc  calculate the angle with the z axis

amoddon=sqrt((dummy(ndon,1)**2.0)+(dummy(ndon,2)**2.0)+&
(dummy(ndon,3)**2.0))
angz=3.14159265359-ACOS(z*dummy(ndon,3)/(z*amoddon))
else
x=0.0
y=0.0
z=-1.0
pervx=y*dummy(ndon,3)-z*dummy(ndon,2)
pervy=-x*dummy(ndon,3)+z*dummy(ndon,1)
pervz=x*dummy(ndon,2)-y*dummy(ndon,1)

!cc  calculate the angle with the -z axis

amoddon=sqrt((dummy(ndon,1)**2.0)+(dummy(ndon,2)**2.0)+&
(dummy(ndon,3)**2.0))
angz=3.14159265359-ACOS(z*dummy(ndon,3)/(abs(z)*amoddon))
endif
!cc  mtrix multiplication or rotation along a line
u=pervx
v=pervy
w=pervz
a=0.0
b=0.0
c=0.0
al=u**2+v**2+w**2
!print*,al
!al=1.0

x=0.0
y=0.0
z=0.0

do i=1,mno_at
x=dummy(i,1)
y=dummy(i,2)
z=dummy(i,3)

!print*,'x,y,z',x,y,z,u,v,w,angz

dummy(i,1)=((a*(v**2.0+w**2.0)-u*(b*v+c*w-u*x-v*y&
-w*z))*(1.0-cos(angz))+al*x*cos(angz)+sqrt(al)&
*(c*v+b*w+w*y+v*z)*sin(angz))/al

dummy(i,2)=((b*(u**2.0+w**2.0)-v*(a*u+c*w-u*x-v*y&
-w*z))*(1.0-cos(angz))+al*y*cos(angz)+sqrt(al)&
*(c*u-a*w+w*x-u*z)*sin(angz))/al

dummy(i,3)=((c*(u**2.0+v**2.0)-w*(a*u+b*v-u*x-v*y&
-w*z))*(1.0-cos(angz))+al*z*cos(angz)+sqrt(al)&
*(-b*u+a*v-v*x+u*y)*sin(angz))/al
enddo

!rotation of the third atom to the xz plane

600 if(third.eq.0)goto 602
if(dummy(third,1).eq.0.0.and.dummy(third,2).eq.0.0)then
print*,'U should choose the third atom which is not on the z axis'
stop
endif
angzx=3.14159265359-atan(dummy(third,2)/dummy(third,1))

!rotation along z axis

do i=1,mno_at
x=dummy(i,1)
y=dummy(i,2)

if(x.ne.0.0.or.y.ne.0.0)goto 601
enddo
goto 602

601 do i=1,mno_at
x=dummy(i,1)
y=dummy(i,2)
z=dummy(i,3)

dummy(i,1)=cos(angzx)*x-sin(angzx)*y

dummy(i,2)=sin(angzx)*x+cos(angzx)*y

dummy(i,3)=z

enddo

 
602 amidx=(dummy(nacc,1)+dummy(ndon,1))/2
amidy=(dummy(nacc,2)+dummy(ndon,2))/2
amidz=(dummy(nacc,3)+dummy(ndon,3))/2
do i=1,mno_at
dummy(i,1)=dummy(i,1)-amidx
dummy(i,2)=dummy(i,2)-amidy
dummy(i,3)=dummy(i,3)-amidz
enddo

return
end subroutine rottoz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none 
character(len=2)::flag
integer:: jj,kk,i,j,k,f,vv,var,storei,mno,ndate
real:: dista,distd,matrix,drad,ddist
dimension storei(500),matrix(500)
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



!print*,'do u want to use the atomic database for creation monomers'
!read(*,'(a)')flag
if(fg1.ne.1)then
jj=1
kk=1
do i=1,mno_at
dista=0.0
distd=0.0
dista=sqrt((dummy(nacc,1)-dummy(i,1))**2.0+(dummy(nacc,2)-dummy(i,2))**2.0 &
+(dummy(nacc,3)-dummy(i,3))**2.0)
distd=sqrt((dummy(ndon,1)-dummy(i,1))**2.0+(dummy(ndon,2)-dummy(i,2))**2.0 &
+(dummy(ndon,3)-dummy(i,3))**2.0)
if (distd.gt.dista) then
accgroup(jj,1)=dummy(i,1)
accgroup(jj,2)=dummy(i,2)
accgroup(jj,3)=dummy(i,3)
accst(jj)=stt(i)
jj=jj+1
else
dongroup(kk,1)=dummy(i,1)
dongroup(kk,2)=dummy(i,2)
dongroup(kk,3)=dummy(i,3)
donst(kk)=stt(i)
kk=kk+1
endif
enddo
endif

do i=1,mno_at
!print*,stt(i)
enddo
if(fg1.eq.1)then
k=0
do i=1,mno_at
do j=1,88
if(stt(i).eq.at_list(j))then
!matrix(i)=at_covrad(j)/100.0
matrix(i)=at_covrad(j)/100.0
endif
enddo
enddo

!print*,'1'
jj=1
mno=mno_at
f=1
vv=1
202 if(kk.eq.0.and.f.ne.1)goto 203 
if(f.ne.1)vv=kk
!print*,'2',kk
kk=0
do j=1,vv
if (j.eq.1.and.f.eq.1) var=nacc
if (f.ne.1) var=storei(jj-j)
!print*,'********var',var,jj-j,jj
do i=1,mno
if(i.eq.nacc.or.i.eq.ndon)goto 201
ddist=sqrt((dummy(var,1)-dummy(i,1))**2.0+(dummy(var,2)-dummy(i,2))**2.0+&
(dummy(var,3)-dummy(i,3))**2.0)
drad=matrix(var)+matrix(i)
if(ddist.eq.0.or.drad.eq.0)goto 201
!print*,'drad,ddist',drad,ddist,var
if(ddist.le.drad)then
!print*,'333333',jj,kk,i
storei(jj)=i
accgroup(jj,1)=dummy(i,1)
accgroup(jj,2)=dummy(i,2)
accgroup(jj,3)=dummy(i,3)
accst(jj)=stt(i)
kk=kk+1
jj=jj+1
f=0
endif
201 enddo
enddo
goto 202

203 accgroup(jj,1)=dummy(nacc,1)
accgroup(jj,2)=dummy(nacc,2)
accgroup(jj,3)=dummy(nacc,3)
accst(jj)=stt(nacc)
kk=0
do i=1,mno_at
do j=1,jj-1
iF (i.eq.storei(j).or.i.eq.nacc)goto 204
enddo
kk=kk+1
dongroup(kk,1)=dummy(i,1)
dongroup(kk,2)=dummy(i,2)
dongroup(kk,3)=dummy(i,3)
donst(kk)=stt(i)
204 enddo
kk=kk+1
jj=jj+1
endif
mono1=jj
mono2=kk

mnm1=mono1-1
mnm2=mono2-1
return
end subroutine output1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output2(di,m1,m2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none 
integer::l,k,i,ndate
character(len=35)::di,m1,m2,chkdi,chkm1,chkm2,dim,m1m,m2m
character(len=5)::chk 
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


l=len(trim(di))
dim=di(2:l)
l=len(trim(m1))
m1m=m1(2:l)
l=len(trim(m2))
m2m=m2(2:l)
chk='%chk='
chkdi=trim(chk)//trim(dim)
chkm1=trim(chk)//trim(m1m)
chkm2=trim(chk)//trim(m2m)

open(unit=17,file=dimname,status='unknown')
!print*,dimname,mononame1,mononame2,cubname        
write(17,'(a)')chkdi
write(17,101)line2
write(17,102)line3
write(17,101)line13
write(17,*)
write(17,*)'dimer'
write(17,*)
write(17,109)ch0,sp0
109 format(a,1x,a)
do i=1,mno_at
write(17,105)stt(i),(dummy(i,k),k=1,3)
enddo
write(17,*)
if(fg8.eq.1.and.flg20.eq.1)then
l=len(trim(line27))
write(17,102)line27(2:l)
write(17,*)
endif

open(unit=11,file=mononame1,status='unknown')
        
write(11,'(a)')chkm1
write(11,101)line2
101 FORMAT(a)
write(11,102)line3
102 format(a)
write(11,101)line13
!103 format('#p',2x,a,a,'nosymm')
write(11,*)
write(11,*)'1st monomer'
write(11,*)
write(11,104)ch1,sp1
104 format(a,1x,a)
if(fg1.ne.2)then
do i=1,mono1-1
write(11,105)accst(i),(accgroup(i,k),k=1,3)
enddo
endif
if(fg1.eq.2)then
do i=1,mnm1
write(11,105)stt(mon1(i)),(dummy(mon1(i),k),k=1,3)
enddo
endif
105 format(a,4x,f10.6,4x,f10.6,4x,f10.6)
write(11,*)
if(fg8.eq.1.and.flg21.eq.1)then
l=len(trim(line28))
write(11,102)line28(2:l)
write(11,*)
endif

open(unit=13,file=mononame2,status='unknown')
        
write(13,'(a)')chkm2
write(13,101)line2
write(13,102)line3
write(13,101)line13
write(13,*)
write(13,*)'2nd monomer'
write(13,*)
write(13,104)ch2,sp2
if(fg1.ne.2)then
do i=1,mono2-1
write(13,105)donst(i),(dongroup(i,k),k=1,3)
enddo
endif
if(fg1.eq.2)then
do i=1,mnm2
write(13,105)stt(mon2(i)),(dummy(mon2(i),k),k=1,3)
enddo
endif
write(13,*)
if(fg8.eq.1.and.flg22.eq.1)then
l=len(trim(line29))
write(13,102)line29(2:l)
write(13,*)
endif

return
end subroutine output2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubegrid
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
angtobohr=0.529177249
do i=1,mno_at
dummy(i,1)=dummy(i,1)/angtobohr
dummy(i,2)=dummy(i,2)/angtobohr
dummy(i,3)=dummy(i,3)/angtobohr
enddo

do i=1,mno_at
do j=1,86
if(at_list(j).eq.stt(i))then
dummy_min(i,1)=dummy(i,1)+(-at_rad(j))
dummy_min(i,2)=dummy(i,2)+(-at_rad(j))
dummy_min(i,3)=dummy(i,3)+(-at_rad(j))
dummy_max(i,1)=dummy(i,1)+(at_rad(j))
dummy_max(i,2)=dummy(i,2)+(at_rad(j))
dummy_max(i,3)=dummy(i,3)+(at_rad(j))
endif
enddo
enddo
do i=1,mno_at
!print*,'dummy_min(i,1)',(dummy_min(i,j),j=1,3)
enddo
do j=1,3
greatest=abs(dummy_min(1,j))
do i=2,mno_at
if(abs(dummy_min(i,j)).gt.greatest)greatest=abs(dummy_min(i,j))
enddo
min_grd(j)=-greatest-2.0
enddo
do j=1,3
greatest=abs(dummy_max(1,j))
do i=2,mno_at
if(abs(dummy_max(i,j)).gt.greatest)greatest=abs(dummy_max(i,j))
enddo
max_grd(j)=greatest+2.0
enddo
do i=1,3
steps(i)=(max_grd(i)-min_grd(i))/0.1
enddo


open(unit=3,file=cubname,status='unknown')
write(3,*)
write(3,*)
write(3,111)(min_grd(i),i=1,3)
write(3,112)steps(1),grd_sp,0.000000,0.000000
write(3,112)steps(2),0.000000,grd_sp,0.000000
write(3,112)steps(3),0.000000,0.000000,grd_sp
write(3,*)
111 format('1',4x,f10.6,4x,f10.6,4x,f10.6)
112 format(i3,4x,f10.6,4x,f10.6,4x,f10.6)

extra=0
div=steps(3)/6
if(mod(steps(3),6).ne.0)extra=1
MN=div+extra
MNXY=steps(1)*steps(2)*MN
MNZ=steps(3)
return
end subroutine cubegrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine runcom(filen,aname2,gauss)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use commondat
implicit none
logical::file1
integer::io,MDP,l,i,ndate
character(len=35)::lines,filen,aname2
character(15):: gauss
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


lines='sourav'
open(unit=39,file='script',status='unknown')
write(39,111)trim(gauss),trim(filen)
111 Format (a,1x,a)
l=len(TRIM(filen))

write(*,10)trim(aname2)
10 Format(a,' file is preparing '//)

close(39)
call systemqq ("chmod +x script")
CALL SYSTEMQQ ("./script")
!CALL EXECUTE_COMMAND_LINE ( 'g09 ethilin_dim.com', WAIT = .TRUE. )
448 INQUIRE(file=trim(aname2),EXIST=file1)
IF (file1) THEN


open(unit=37,file=trim(aname2),status='old')
447 MDP=0
do
read(37,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(2:19).eq.'Normal termination')goto 446
if(lines(2:18).eq.'Error termination')goto 444
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(37)
goto 447
444 write(*,445)trim(aname2)
445 Format (' there in no normal termination in the Outpue file ',a)
stop

446 close(37)


goto 449
else
goto 448
endif

!call chkfile
449 CALL SYSTEMQQ ("rm script")
CALL SYSTEMQQ ("rm *.sh*")
CALL SYSTEMQQ ("rm *.old*")
return
end subroutine runcom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine runcom_sim(gau)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1,file2,file3
integer::io,MDP,l,i,ndate
character(len=35)::lines,filen,aname2
character(15):: gau
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


lines='sourav'
open(unit=39,file='script1',status='unknown')
write(39,111)trim(gau),trim(dimname)
open(unit=37,file='script2',status='unknown')
write(37,111)trim(gau),trim(mononame1)
open(unit=35,file='script3',status='unknown')
write(35,111)trim(gau),trim(mononame2)
111 Format (a,1x,a)
!l=len(TRIM(filen))

write(*,*)'*************************************************************************'
write(*,10)trim(lognamedi),trim(lognamem1),trim(lognamem2)
10 Format(a,a,a' files are preparing ')
write(*,*)'*************************************************************************'

close(39)
close(37)
close(35)
call systemqq ("chmod +x script1")
call systemqq ("chmod +x script2")
call systemqq ("chmod +x script3")
CALL SYSTEMQQ ("./script1")
CALL SYSTEMQQ ("./script2")
CALL SYSTEMQQ ("./script3")
448 INQUIRE(file=trim(lognamedi),EXIST=file1)
IF (file1) THEN


open(unit=49,file=trim(lognamedi),status='old')
447 MDP=0
do
read(49,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(2:19).eq.'Normal termination')goto 446
if(lines(2:18).eq.'Error termination')goto 444
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(49)
goto 447
444 write(*,445)trim(lognamedi)
445 Format (' there in no normal termination in the Outpue file ',a)
stop

446 close(49)


goto 441
else
goto 448
endif


441 lines='sourav'
348 INQUIRE(file=trim(lognamem1),EXIST=file2)
IF (file2) THEN


open(unit=47,file=trim(lognamem1),status='old')
347 MDP=0
do
read(47,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(2:19).eq.'Normal termination')goto 346
if(lines(2:18).eq.'Error termination')goto 344
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(47)
goto 347
344 write(*,345)trim(lognamem1)
345 Format (' there in no normal termination in the Outpue file ',a)
stop

346 close(47)


goto 442
else
goto 348
endif


442 lines='sourav'
248 INQUIRE(file=trim(lognamem2),EXIST=file3)
IF (file3) THEN


open(unit=45,file=trim(lognamem2),status='old')
247 MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(2:19).eq.'Normal termination')goto 246
if(lines(2:18).eq.'Error termination')goto 244
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)
goto 247
244 write(*,245)trim(lognamem2)
245 Format (' there in no normal termination in the Outpue file ',a)
stop

246 close(45)


goto 449
else
goto 248
endif
!call chkfile
449 CALL SYSTEMQQ ("rm script1")
CALL SYSTEMQQ ("rm script2")
CALL SYSTEMQQ ("rm script3")
CALL SYSTEMQQ ("rm *.sh*")
CALL SYSTEMQQ ("rm *.old*")
return
end subroutine runcom_sim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine genformchk(aname2,aname4,moob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,j,i,moob,dibasl,ndate
character(len=35)::lines,filen,aname,aname2,aname1,aname3,aname4,xyz1(4)
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



lines='sourav'
open(unit=41,file='script',status='unknown')
write(41,111)trim(fchk_com),trim(aname2)
111 Format (a,1x,a)

close(41)
call systemqq ("chmod +x script")
CALL SYSTEMQQ ("./script")
448 INQUIRE(file=trim(aname4),EXIST=file1)
IF (file1) THEN


open(unit=43,file=trim(aname4),status='old')
447 MDP=0
do
read(43,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(1:26).eq.'Number of basis functions')dibasl=MDP
if(lines(5:20).eq.'coupling tensors')goto 446
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(43)
goto 447

446 close(43)

goto 449
else
goto 448
endif

!call chkfile
449 write(*,100)trim(aname4)
100 Format(//' formatted chk point file "',a,'" is ready')
CALL SYSTEMQQ ("rm script")
open(unit=43,file=trim(aname4),status='old')
do i=1,dibasl-1
read(43,*)
enddo
read(43,*)(xyz1(j),j=1,4),moob
return
end subroutine genformchk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine genformchk_sim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,i,j,ndate
character(len=35)::lines,filen,aname,aname2,aname1,aname3,aname4
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



lines='sourav'
open(unit=41,file='script1',status='unknown')
write(41,111)fchk_com,chknamedi
open(unit=51,file='script2',status='unknown')
write(51,111)fchk_com,chknamem1
open(unit=61,file='script3',status='unknown')
write(61,111)fchk_com,chknamem2
111 Format (a,1x,a)

close(41)
close(51)
close(61)
call systemqq ("chmod +x script1")
call systemqq ("chmod +x script2")
call systemqq ("chmod +x script3")
CALL SYSTEMQQ ("./script1")
CALL SYSTEMQQ ("./script2")
CALL SYSTEMQQ ("./script3")

448 INQUIRE(file=trim(fchknamedi),EXIST=file1)
IF (file1) THEN


open(unit=43,file=trim(fchknamedi),status='old')
447 MDP=0
do
read(43,'(a)',iostat=io)lines
!print*,lines(2:18)
!if(lines(1:26).eq.'Number of basis functions')dibasl=MDP
if(lines(5:20).eq.'coupling tensors')goto 446
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(43)
goto 447

446 close(43)

goto 348
else
goto 448
endif

348 INQUIRE(file=trim(fchknamem1),EXIST=file1)
IF (file1) THEN


open(unit=53,file=trim(fchknamem1),status='old')
347 MDP=0
do
read(53,'(a)',iostat=io)lines
!print*,lines(2:18)
!if(lines(1:26).eq.'Number of basis functions')m1basl=MDP
if(lines(5:20).eq.'coupling tensors')goto 346
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(53)
goto 347

346 close(53)

goto 248
else
goto 348
endif

248 INQUIRE(file=trim(fchknamem2),EXIST=file1)
IF (file1) THEN


open(unit=63,file=trim(fchknamem2),status='old')
247 MDP=0
do
read(63,'(a)',iostat=io)lines
!print*,lines(2:18)
!if(lines(1:26).eq.'Number of basis functions')m2basl=MDP
if(lines(5:20).eq.'coupling tensors')goto 246
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(63)
goto 247

246 close(63)

goto 449
else
goto 248
endif

!call chkfile
449 write(*,109)
109 Format(//'*********************************************************************&
             *******************************')
write(*,100)trim(fchknamedi),trim(fchknamem1),trim(fchknamem2)
100 Format(' formatted chk point files "',a,1x,a,1x,a,'" are ready')
write(*,110)
110 Format('***********************************************************************&
           ********************************')
CALL SYSTEMQQ ("rm script1")
CALL SYSTEMQQ ("rm script2")
CALL SYSTEMQQ ("rm script3")


return
end subroutine genformchk_sim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubgen(filen,aname2,cubname,cubcom)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use commondat
implicit none
real::start,finish,time
logical::file1
integer::io,MDP,l,i,k,mdp1(3000),ndate
character(len=35)::filen,aname2,cubname
character(75):: cubcom,line
character(len=8)  :: date
character(len=5)  :: zone
integer,dimension(8) :: values




open(unit=43,file='script',status='unknown')
write(43,111)trim(cubcom),trim(filen),trim(aname2),trim(cubname)
111 Format (a,1x,a,1x,a,1x,'-1',1x,'h',1x,a)

close(43)
!open(unit=35,file='cubrun_check',status='unknown')
!write(35,*)'sourav'
!rewind(35)
call systemqq ("chmod +x script")
CALL SYSTEMQQ ("./script")
write(*,100)trim(aname2)
100 Format(// a,' file is preparing')
448 INQUIRE(file=trim(aname2),EXIST=file1)
IF (file1) THEN

k=0

open(unit=47,file=trim(aname2),status='old')
447 MDP=0
k=k+1
time=time+10
do
read(47,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(47)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 449
202 call cpu_time(start)
if(start.lt.time)goto 202

goto 447
else
goto 448
endif

449 CALL SYSTEMQQ ("rm script")
return
end subroutine cubgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubgen_sim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
real::start,finish,time
logical::file1
integer::io,MDP,l,i,k,mdp1(3000),ndate
character(75)::line
character(len=8)  :: date
character(len=5)  :: zone
integer,dimension(8) :: values




open(unit=43,file='script1',status='unknown')
write(43,111)trim(cubcom),trim(fchknamedi),trim(cubnamedi),trim(cubname)
open(unit=53,file='script2',status='unknown')
write(53,111)trim(cubcom),trim(fchknamem1),trim(cubnamem1),trim(cubname)
open(unit=63,file='script3',status='unknown')
write(63,111)trim(cubcom),trim(fchknamem2),trim(cubnamem2),trim(cubname)
111 Format (a,1x,a,1x,a,1x,'-1',1x,'h',1x,a)

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
write(*,100)trim(cubnamedi),trim(cubnamem1),trim(cubnamem2)
100 Format( a,1x,a,1x,a,' files are preparing')
write(*,*)'********************************************************************************'
448 INQUIRE(file=trim(cubnamedi),EXIST=file1)
IF (file1) THEN

k=0

open(unit=47,file=trim(cubnamedi),status='old')
447 MDP=0
k=k+1
time=time+10
do
read(47,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(47)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 348
202 call cpu_time(start)
if(start.lt.time)goto 202

goto 447
else
goto 448
endif

348 INQUIRE(file=trim(cubnamem1),EXIST=file1)
IF (file1) THEN

k=0

open(unit=57,file=trim(cubnamem1),status='old')
347 MDP=0
k=k+1
time=time+10
do
read(57,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(57)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 248
302 call cpu_time(start)
if(start.lt.time)goto 302

goto 347
else
goto 348
endif


248 INQUIRE(file=trim(cubnamem2),EXIST=file1)
IF (file1) THEN

k=0

open(unit=67,file=trim(cubnamem2),status='old')
247 MDP=0
k=k+1
time=time+10
do
read(67,'(a)',iostat=io)
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(67)
mdp1(k)=MDP
if(k.eq.6.and.mdp1(k).eq.mdp1(k-1).and.mdp1(k).eq.mdp1(k-2).and.mdp1(k).eq.mdp1(k-3).&
and.mdp1(k).eq.mdp1(k-4))goto 449
402 call cpu_time(start)
if(start.lt.time)goto 402

goto 247
else
goto 248
endif

449 CALL SYSTEMQQ ("rm script1")
CALL SYSTEMQQ ("rm script2")
CALL SYSTEMQQ ("rm script3")
return
end subroutine cubgen_sim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intsum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
integer::n,nx,ny,nz,j,i1,n2,ox,oy,i,k,nlast,noe,int1,int2,int3,ndate
real::at_grd(10,300),sub_val(MNXY*MNZ),valdi(10),valm2(10),valm1(10),vallastdi(10)&
,vallastm1(10),vallastm2(10),sumxy(MNZ),int,val(MNXY*MNZ),oz,grd1(5),grd2(3,4),norm_gr,noz&
,char1,char2,char3,char4
character(len=70)::grd
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


open(unit=57,file=cubnamedi,status='old')
open(unit=67,file=cubnamem1,status='old')
open(unit=77,file=cubnamem2,status='old')
open(unit=87,file=substract,status='unknown')
open(unit=97,file=integral,status='unknown')
do i=1,2
read(57,*)
enddo
write(87,*)'Charge Displacement Density Cube File'
write(87,*)'dimer density - monomer1 density -monomer2 density'
read(57,*)int1,(grd1(i),i=1,3),int2
write(87,305)int1,(grd1(i),i=1,3),int2
305 format(I5,2x,f10.6,2x,f10.6,2x,f10.6,2x,I3)
do i=1,3
read(57,*)int3,(grd2(i,j),j=1,3)
write(87,306)int3,(grd2(i,j),j=1,3)
enddo
306 format(I5,2x,f10.6,2x,f10.6,2x,f10.6)
do i=1,mno_at
read(57,*)noe,(at_grd(i,j),j=1,4)
write(87,307)noe,(at_grd(i,j),j=1,4)
enddo
307 format(I5,2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)
do i=1,mnm1+6
read(67,*)
enddo
do i=1,mnm2+6
read(77,*)
enddo

ox=min_grd(1)
oy=min_grd(2)
oz=min_grd(3)
nx=steps(1)
ny=steps(2)
nz=steps(3)
n=nz/6
n2=0
nlast=mod(nz,6)


write(*,100)trim(substract)
100 Format(// a,' file is preparing')
do i=1,nx*ny

do j=1,n
do j1=1,10
valdi(j1)=0.0
valm1(j1)=0.0
valm2(j1)=0.0
enddo
read(57,*)(valdi(k),k=1,6)
read(67,*)(valm1(k),k=1,6)
read(77,*)(valm2(k),k=1,6)
do i1=1,6
n2=n2+1
sub_val(n2)=valdi(i1)-valm1(i1)-valm2(i1)
!write(87,*)sub_val(n2)
enddo
enddo
if(nlast.eq.0)goto 500
do j1=1,10
vallastdi(j1)=0.0
vallastm1(j1)=0.0
vallastm2(j1)=0.0
enddo
read(57,*)(vallastdi(k),k=1,nlast)
read(67,*)(vallastm1(k),k=1,nlast)
read(77,*)(vallastm2(k),k=1,nlast)
do i1=1,nlast
n2=n2+1
sub_val(n2)=vallastdi(i1)-vallastm1(i1)-vallastm2(i1)
!write(87,*)sub_val(n2)
enddo
500 enddo


k=1
do i=1,nx*ny

do j=1,n
write(87,'(6E13.5)')(sub_val(l),l=k,k+5)
k=k+6
enddo
if(nlast.eq.0)goto 501
write(87,'(<nlast>E13.5)')(sub_val(l),l=k,k+nlast-1)
!103 Format ('(<nlast>(f13.5))')
k=k+nlast
501 enddo

norm_gr=grd_sp/abs(dummy(nacc,3)-dummy(ndon,3))
noz=oz/abs(dummy(nacc,3)-dummy(ndon,3))

write(*,101)trim(integral)
101 Format(// a,' file is preparing')
write(97,*) 'norm_Z(Ang)     Z(Ang)         SUMXY          INT           Charge'
int=0.0
do i=1,nz
sumxy(i)=0.0
do j=i,nx*ny*nz,nz
sumxy(i)=sumxy(i)+sub_val(j)
enddo
int=int+sumxy(i)
write(97,*)noz*0.529177249,oz*0.529177249,sumxy(i)*(grd_sp**2),int*(grd_sp**3),int*(grd_sp**3.0)*0.53**3.0
if(noz*0.529177249.lt.-0.5)then
char1=char1+int*(grd_sp**3.0)*0.53**3.0
endif
if(noz*0.529177249.ge.-0.5.and.noz*0.529177249.le.0.0)then
char2=char2+int*(grd_sp**3.0)*0.53**3.0
endif
if(noz*0.529177249.ge.0.0.and.noz*0.529177249.le.0.5)then
char3=char3+int*(grd_sp**3.0)*0.53**3.0
endif
if(noz*0.529177249.gt.0.5)then
char4=char4+int*(grd_sp**3.0)*0.53**3.0
endif
oz=oz+grd_sp
noz=noz+norm_gr
enddo
write(97,*)
write(97,*)
write(97,*)'total charge from -inf to -0.5 is = ',char1
write(97,*)'total charge from -0.5 to 0.0 is = ',char2
write(97,*)'total charge from 0.0 to 0.5 is = ',char3
write(97,*)'total charge from 0.5 to +inf is = ',char4
CALL SYSTEMQQ ("rm *script*")
CALL SYSTEMQQ ("rm *.old*")
CALL SYSTEMQQ ("rm *sh*")
CALL SYSTEMQQ ("rm *temp*")

CALL SYSTEMQQ ("mkdir TMP_files")
CALL SYSTEMQQ ("mv *dimer.cub TMP_files/")
CALL SYSTEMQQ ("mv *mono1.cub TMP_files/")
CALL SYSTEMQQ ("mv *mono2.cub TMP_files/")
CALL SYSTEMQQ ("mv *.chk TMP_files/")

return
end subroutine intsum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readfchk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,i,i1,j,dialphae,m1alphae,m2alphae,alphaedi,betaedi,alphaem1,betaem1,alphaem2,betaem2,&
di1,p1,l1,l2,k,kk,c(100),d(100),f1,cc(50,100),dd(50,100),l3,j1,j2,nlast,mnnum,ndate,&
sdimo(1000),pdimo(1000),ddimo(1000),sdi,pdi,ddi,sm1,sm2,pm1,pm2,dm1,dm2,sm1mo(1000),pm1mo(1000),dm1mo(1000),&
sm2mo(1000),pm2mo(1000),dm2mo(1000),ddi1,ddimo1(1000),dm11,dm1mo1(1000),dm21,dm2mo1(1000)
real::momat1(6,5),momat2(6,5),start
character(len=50)::lines,XYZ1(5),XYZ2(5),XYZ3(5)
character(len=3)::a,b,e,na,nd,matorb(5)
character(len=8)::fmt1,fmt2,f2,aa,bb,ans
character(len=75)::logn,namo
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
lines='Sourav_Roy'
open(unit=43,file=trim(fchknamedi),status='old')
open(unit=53,file=trim(fchknamem1),status='old')
open(unit=63,file=trim(fchknamem2),status='old')

INQUIRE(file=trim(fchknamedi),EXIST=file1)
IF (file1) THEN


open(unit=43,file=trim(fchknamedi),status='old')
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
print*,trim(fchknamedi),'file does not exist'
stop
endif
do i=1,dialphae
read(43,*)
enddo
read(43,*)(XYZ3(j),j=1,5),alphaedi
read(43,*)(XYZ3(j),j=1,5),betaedi
read(43,*)(XYZ3(j),j=1,5),NMODI
rewind(43)


INQUIRE(file=trim(fchknamem1),EXIST=file1)
IF (file1) THEN


open(unit=53,file=trim(fchknamem1),status='old')
MDP=0
do
read(53,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(1:26).eq.'Number of alpha electrons')m1alphae=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(53)

else
print*,trim(fchknamem1),'file does not exist'
stop
endif

do i=1,m1alphae
read(53,*)
enddo
read(53,*)(XYZ3(j),j=1,5),alphaem1
read(53,*)(XYZ3(j),j=1,5),betaem1
read(53,*)(XYZ3(j),j=1,5),NMOM1
rewind(53)


INQUIRE(file=trim(fchknamem2),EXIST=file1)
IF (file1) THEN


open(unit=63,file=trim(fchknamem2),status='old')
MDP=0
do
read(63,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(1:26).eq.'Number of alpha electrons')m2alphae=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(63)

else
print*,trim(fchknamem2),'file does not exist'
stop
endif

do i=1,m2alphae
read(63,*)
enddo
read(63,*)(XYZ3(j),j=1,5),alphaem2
read(63,*)(XYZ3(j),j=1,5),betaem2
read(63,*)(XYZ3(j),j=1,5),NMOM2
rewind(63)
!print*,NMODI,alphaedi,betaedi
if(fg4.ne.1)goto 800

open(unit=9,file='Molecular_Orbitals.dat',status='unknown')

INQUIRE(file=trim(lognamedi),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamedi),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:36)
if(lines(6:36).eq.'Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(lognamedi),'file does not exist'
stop
endif

write(9,*)'Molecular orbital file of',molname
write(9,*)
write(9,*)'Molecular orbitals of Dimer'
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
!print*,'**************************'
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


!!!! Orbitals of MONOMER 1

do mnnum=1,2
p1=0
if (mnnum.eq.1)NMODI=NMOM1*2
if (mnnum.eq.2)NMODI=NMOM2*2
if (mnnum.eq.1)then
logn=lognamem1
write(9,*)
write(9,*)'Alpha Beta Molecular Orbitals for monomer 1'
write(9,*)
endif
if (mnnum.eq.2)then
logn=lognamem2
write(9,*)
write(9,*)'Alpha Beta Molecular Orbitals for monomer 2'
write(9,*)
endif
INQUIRE(file=trim(logn),EXIST=file1)
IF (file1) THEN


open(unit=46,file=trim(logn),status='old')
MDP=0
do
read(46,'(a)',iostat=io)lines
!print*,lines(6:42)
if(lines(6:44).eq.'Alpha Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(46)

else
print*,trim(logn),'file does not exist'
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
write(9,*)'Alpha Molecular Orbitals'
write(9,*)
l=mod(NMODI/2,5)
nlast=l
k=(NMODI/2-l)/5
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
do i=1,NMODI/2
read(46,'(a,a,a)')a,b,e
read(e,'(I10)')f1
if(f1.eq.nacc)then
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
do i=1,NMODI/2-c(j)
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
do i=1,(NMODI/2)-(cc(j,1)+3)
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
do i=1,(NMODI/2)-(cc(1,1)+3)
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
l=mod(NMODI/2,5)
nlast=l
k=(NMODI/2-l)/5
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
do i=1,NMODI/2
read(46,'(a,a,a)')a,b,e
read(e,'(I10)')f1
if(f1.eq.nacc)then
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
do i=1,(NMODI/2)-c(j)
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
do i=1,NMODI/2-(cc(j,1)+3)
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
do i=1,(NMODI/2)-(cc(1,1)+3)
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

800 print*,'Total no of Molecular orbitals in DIMER is ',NMODI,' among them' 
print*,'No of Sigma orbital is=',sdi,' and they are :',(sdimo(i),i=1,sdi)
print*,'No of Pi orbital is=',pdi,' and they are :',(pdimo(i),i=1,pdi)
if(ddi.ne.0)print*,'No of Delta orbital is=',ddi,' and they are :',(ddimo(i),i=1,ddi)
print*,'***************************************************************************'
print*,'Total no of Molecular orbitals 1ST MONOMER is ',sm1+pm1+dm1,' among them'
print*,'No of Sigma orbital is=',sm1,' and they are :',(sm1mo(i),i=1,sm1)
print*,'No of Pi orbital is=',pm1,' and they are :',(pm1mo(i),i=1,pm1)
if(dm1.ne.0)print*,'No of Delta orbital is=',dm1,' and they are :',(dm1mo(i),i=1,dm1)
print*,'****************************************************************************'
print*,'Total no of Molecular orbitals 2ND MONOMER is ',sm2+pm2+dm2,' among them'
print*,'No of Sigma orbital is=',sm2,' and they are :',(sm2mo(i),i=1,sm2)
print*,'No of Pi orbital is=',pm2,' and they are :',(pm2mo(i),i=1,pm2)
if(dm2.ne.0)print*,'No of Delta orbital is=',dm2,' and they are :',(dm2mo(i),i=1,dm2)
print*,'****************************************************************************'
print*,'if you want to calculate the charge transfar among the Molecular Orbitals'
print*,'then please see the Molecular_Orbitals.dat file and prepare the inputfile: input.dat'
print*,'and after that write yes or y here other wise put no or n'
read*,ans
if(ans.eq.'yes'.or.ans.eq.'y')then
call readinp2
ddi1=sdi+pdi+ddi
dm11=sm1+pm1+dm1
dm21=sm2+pm2+dm2
i=0
do j=1,sdi
i=i+1
ddimo1(i)=sdimo(j)
enddo
do j=1,pdi
i=i+1
ddimo1(i)=pdimo(j)
enddo
do j=1,ddi
i=i+1
ddimo1(i)=ddimo(j)
enddo
i=0
do j=1,sm1
i=i+1
dm1mo1(i)=sm1mo(j)
enddo
do j=1,pm1
i=i+1
dm1mo1(i)=pm1mo(j)
enddo
do j=1,dm1
i=i+1
dm1mo1(i)=dm1mo(j)
enddo
i=0
do j=1,sm2
i=i+1
dm2mo1(i)=sm2mo(j)
enddo
do j=1,pm2
i=i+1
dm2mo1(i)=pm2mo(j)
enddo
do j=1,dm2
i=i+1
dm2mo1(i)=dm2mo(j)
enddo
!print*,namo,fg11,fg12,fg13,fg14
if(fg10.eq.1)then
if(fg11.eq.1)namo='sigma'
if(fg11.eq.1)call molorb(sdi,sdimo,sm1,sm1mo,sm2,sm2mo,namo)
if(fg12.eq.1)namo='pi'
if(fg12.eq.1)call molorb(pdi,pdimo,pm1,pm1mo,pm2,pm2mo,namo)
!print*,namo
if(fg13.eq.1)namo='delta'
if(fg13.eq.1.and.ddi.ne.0)call molorb(ddi,ddimo,dm1,dm1mo,dm2,dm2mo,namo)
!print*,namo
if(fg14.eq.1)namo='all'
if(fg14.eq.1)call molorb(ddi1,ddimo1,dm11,dm1mo1,dm21,dm2mo1,namo)
endif
if(fg15.eq.1)namo='sigma'
if(fg16.eq.1)namo='pi'
if(fg17.eq.1)namo='delta'
if(fg18.eq.1)namo='all'
if(fg10.ne.1)call molorb(sdi,sdimo,sm1,sm1mo,sm2,sm2mo,namo)
endif
if(ans.eq.'no'.or.ans.eq.'n')then
stop
endif
stop


return
end subroutine readfchk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readfchk_cplx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
logical::file1
integer::io,MDP,l,i,j,dialphae,m1alphae,m2alphae,alphaedi,betaedi,alphaem1,betaem1,alphaem2,betaem2,&
di1,l1,l2,k,kk,c(500),d(100),f1,cc(50,100),dd(50,100),l3,j1,j2,nlast,mnnum,i1,i2,i3,i4,i6,last,mateigen(6),ll&
,ro(3),j3,j5,lst,j6,j7,j8,ndate
real::momat1(5000,5),start,st,momat4(6)
character(len=50)::lines,XYZ1(5),XYZ2(5),XYZ3(5)
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


lines='Sourav_Roy'
open(unit=43,file=trim(fchknamedi),status='old')
open(unit=53,file=trim(fchknamem1),status='old')
open(unit=63,file=trim(fchknamem2),status='old')

INQUIRE(file=trim(fchknamedi),EXIST=file1)
IF (file1) THEN


open(unit=43,file=trim(fchknamedi),status='old')
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
print*,trim(fchknamedi),'file does not exist'
stop
endif
do i=1,dialphae
read(43,*)
enddo
read(43,*)(XYZ3(j),j=1,5),alphaedi
read(43,*)(XYZ3(j),j=1,5),betaedi
read(43,*)(XYZ3(j),j=1,5),NMODI
rewind(43)


INQUIRE(file=trim(fchknamem1),EXIST=file1)
IF (file1) THEN


open(unit=53,file=trim(fchknamem1),status='old')
MDP=0
do
read(53,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(1:26).eq.'Number of alpha electrons')m1alphae=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(53)

else
print*,trim(fchknamem1),'file does not exist'
stop
endif

do i=1,m1alphae
read(53,*)
enddo
read(53,*)(XYZ3(j),j=1,5),alphaem1
read(53,*)(XYZ3(j),j=1,5),betaem1
read(53,*)(XYZ3(j),j=1,5),NMOM1
rewind(53)


INQUIRE(file=trim(fchknamem2),EXIST=file1)
IF (file1) THEN


open(unit=63,file=trim(fchknamem2),status='old')
MDP=0
do
read(63,'(a)',iostat=io)lines
!print*,lines(2:18)
if(lines(1:26).eq.'Number of alpha electrons')m2alphae=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(63)

else
print*,trim(fchknamem2),'file does not exist'
stop
endif

do i=1,m2alphae
read(63,*)
enddo
read(63,*)(XYZ3(j),j=1,5),alphaem2
read(63,*)(XYZ3(j),j=1,5),betaem2
read(63,*)(XYZ3(j),j=1,5),NMOM2
rewind(63)
!print*,NMODI,alphaedi,betaedi

if(fg4.ne.1)goto 700
open(unit=9,file='Molecular_Orbitals.dat',status='unknown')

INQUIRE(file=trim(lognamedi),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamedi),status='old')
MDP=0
do
read(45,'(a)',iostat=io)lines
!print*,lines(6:36)
if(lines(6:36).eq.'Molecular Orbital Coefficients:')di1=MDP
if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(45)

else
print*,trim(lognamedi),'file does not exist'
stop
endif
write(9,*)'Molecular orbital file of',molname
write(9,*)
write(9,*)'Molecular orbitals of Dimer'
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
!print*,kk,nlast,mno_at
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
!print*,'c',c(l1),j,l1
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
l=1
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
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


enddo
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
l=1
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
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


enddo

endif
702 rewind(45)



j7=1
j8=2
fmt='(I1)'
write(cch1,fmt)j7
write(cch2,fmt)j8
if(sp1.eq.cch2.and.sp2.eq.cch2)j6=4
if(sp1.eq.cch2.and.sp2.eq.cch1)j6=3
if(sp1.eq.cch1.and.sp2.eq.cch2)j6=3
if(sp1.eq.cch1.and.sp2.eq.cch1)j6=2
do j5=1,j6
j7=1
j8=2
fmt='(I1)'
write(cch1,fmt)j7
write(cch2,fmt)j8
!print*,j5,j6
if(sp1.eq.cch1.and.sp2.eq.cch1)then
!print*,cch1,cch2,j7,j8
!print*,'***',j5,j6
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
if(j5.eq.1)then
INQUIRE(file=trim(lognamem1),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem1),status='old')
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
print*,trim(lognamem1),'file does not exist'
stop
endif
endif

if(j5.eq.2)then
INQUIRE(file=trim(lognamem2),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem2),status='old')
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
print*,trim(lognamem2),'file does not exist'
stop
endif
endif
endif

if(sp1.eq.cch1.and.sp2.eq.cch2)then
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
if(j5.eq.1)then
INQUIRE(file=trim(lognamem1),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem1),status='old')
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
print*,trim(lognamem1),'file does not exist'
stop
endif
endif

if(j5.eq.2)then
INQUIRE(file=trim(lognamem2),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem2),status='old')
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
print*,trim(lognamem2),'file does not exist'
stop
endif
endif
if(j5.eq.3)then
INQUIRE(file=trim(lognamem2),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem2),status='old')
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
print*,trim(lognamem2),'file does not exist'
stop
endif
endif
endif

if(sp1.eq.cch2.and.sp2.eq.cch1)then
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
if(j5.eq.1)then
INQUIRE(file=trim(lognamem1),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem1),status='old')
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
print*,trim(lognamem1),'file does not exist'
stop
endif
endif
if(j5.eq.2)then
INQUIRE(file=trim(lognamem1),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem1),status='old')
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
print*,trim(lognamem1),'file does not exist'
stop
endif
endif

if(j5.eq.3)then
INQUIRE(file=trim(lognamem2),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem2),status='old')
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
print*,trim(lognamem2),'file does not exist'
stop
endif
endif
endif

if(sp1.eq.cch2.and.sp2.eq.cch2)then
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
if(j5.eq.1)then
INQUIRE(file=trim(lognamem1),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem1),status='old')
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
print*,trim(lognamem1),'file does not exist'
stop
endif
endif
if(j5.eq.2)then
INQUIRE(file=trim(lognamem1),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem1),status='old')
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
print*,trim(lognamem1),'file does not exist'
stop
endif
endif

if(j5.eq.3)then
INQUIRE(file=trim(lognamem2),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem2),status='old')
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
print*,trim(lognamem2),'file does not exist'
stop
endif
endif
if(j5.eq.4)then
INQUIRE(file=trim(lognamem2),EXIST=file1)
IF (file1) THEN


open(unit=45,file=trim(lognamem2),status='old')
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
print*,trim(lognamem2),'file does not exist'
stop
endif
endif
endif

if(j5.eq.3.or.j5.eq.4)mnm1=mnm2

!write(9,*)'Molecular orbital file of',molname
!print*,'********',cch1,cch2,j5
j7=1
j8=2
fmt='(I1)'
write(cch1,fmt)j7
write(cch2,fmt)j8
write(9,*)
write(9,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(9,*)
if(sp1.eq.cch2.and.sp2.eq.cch2)then
if(j5.eq.1)write(9,*)'Alpha Molecular orbitals of mnomer1'
if(j5.eq.3)write(9,*)'Alpha Molecular orbitals of mnomer2'
if(j5.eq.2)write(9,*)'Beta Molecular orbitals of mnomer1'
if(j5.eq.4)write(9,*)'Beta Molecular orbitals of mnomer2'
endif
if(sp1.eq.cch2.and.sp2.eq.cch1)then
if(j5.eq.1)write(9,*)'Alpha Molecular orbitals of mnomer1'
if(j5.eq.2)write(9,*)'Beta Molecular orbitals of mnomer1'
if(j5.eq.3)write(9,*)'Molecular orbitals of mnomer2'
endif
if(sp1.eq.cch1.and.sp2.eq.cch2)then
if(j5.eq.1)write(9,*)'Molecular orbitals of mnomer1'
if(j5.eq.2)write(9,*)'Alpha Molecular orbitals of mnomer2'
if(j5.eq.3)write(9,*)'Beta Molecular orbitals of mnomer2'
endif
if(sp1.eq.cch1.and.sp2.eq.cch1)then
if(j5.eq.1)write(9,*)'Molecular orbitals of mnomer1'
if(j5.eq.2)write(9,*)'Molecular orbitals of mnomer2'
endif
write(9,*)

j3=0
if(NMODI/2.ge.5)then
l=mod(NMODI/2,5)
nlast=l
k=(NMODI/2-l)/5
if(l.eq.0)kk=k
if(l.ne.0)kk=k+1
endif
if(NMODI/2.lt.5)then
lst=NMODI/2
nlast=0
kk=1
k=1
endif
l1=0
l2=0
!print*,kk,nlast,mno_at
do j1=1,mnm1
!print*,j1
do i=1,di1+1
read(45,*)
enddo
do j=1,kk
do i=1,3
read(45,*)
enddo
do i=1,NMODI/2
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
!print*,'c',c(l1),j,l1
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
read(45,*)(matorb(i1),i1=1,lst)
!read(45,*)st,nd,(mateigen(i1),i1=1,5)
read(45,*)
!print*,(matorb(i1),i1=1,lst),(mateigen(i1),i1=1,lst)
do i1=1,mnm1
!print*,i1,mnm1
l1=l1+1
read(45,*)a,b,e,momat2(l1),(momat1(l1,i2),i2=1,lst)
!print*,l1,(momat1(l1,i2),i2=1,lst)
last=c(i1+1)-1
if(i1.eq.mnm1)last=NMODI/2
!print*,'last',last
do i3=c(i1)+1,last
l1=l1+1
read(45,*)th,momat2(l1),(momat1(l1,i2),i2=1,lst)
!print*,l1,(momat1(l1,i2),i2=1,5),momat1(l1,2)
enddo
enddo
do i1=1,lst
j3=j3+1
do i=1,6
momat4(6)=0.0
enddo
last=c(1+1)-1
if(mno_at.eq.2)last=NMODI/2
do i2=c(1)+1,last
l=len(trim(momat2(i2)))
nm=momat2(i2)(2:l)
!print*,nm
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
if(i6.eq.0)write(9,*)'All MO are Sigma orbitals'
goto 701
endif
do i=1,3
l=1
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
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


enddo
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
do i1=1,mnm1
!print*,i1,mno_at
l1=l1+1
read(45,*)a,b,e,momat2(l1),(momat1(i1,i2),i2=1,nlast)
!print*,a,b,e,momat2(l1),(momat1(l1,i2),i2=1,nlast)
last=c(i1+1)-1
if(i1.eq.mnm1)last=NMODI/2
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

last=c(1+1)-1
if(mno_at.eq.2)last=NMODI/2
do i2=c(1)+1,last
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
!print*,momat4(1),momat4(2),momat4(3),momat4(4),momat4(5),momat4(6)
do i=1,3
l=1
do i2=1,3
if(abs(momat4(i)).gt.abs(momat4(i2)))then
l=l+1
endif
enddo
ro(i)=l
!momat5(i,l)=momat4(i)
enddo
write(9,*)'**************************************************','MO=',j3,matorb(i1)
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


enddo

endif
701 rewind(45)
enddo
700 print*,'If you want to calculate the charge transfar among the Molecular Orbitals'
print*,'then please see the Molecular_Orbitals.dat file and prepare the inputfile: input.dat'
print*,'and after that write yes or y here other wise put no or n'
print*,'Total no of Molecular orbitals in dimer is ',NMODI,' and for monomers are ',NMODI/2,' each'
read*,ans
if(ans.eq.'yes'.or.ans.eq.'y')goto 400
if(ans.eq.'no'.or.ans.eq.'n')then
stop
endif


400 call readinp2
stop
!call molorb
return
end subroutine readfchk_cplx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinp2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
character(len=70)::infnam,str1,line5,line6,line8,line9,line10,line11,line12,line14,runn,line15
character(len=12)::abc,sttr(100),str2,charst(100),bochtrns,lpe_mono
character(len=10)::string
character(len=18)::lettstr(40)
!character(len=:), allocatable::lettstr(:)
integer::i,MDP,io,l,j,m,k,ll,kk,jj,stn,a,b,p,q,t,i3,i4,i2,i5,c,d,m5,m6
integer::st_num1(100),st_num2(100),st_num3(100),st_num4(100),number(500),ndate
character(len=8)  :: date
character(len=10) :: time
character(len=5)  :: zone
character(len=110)  :: line34
integer,dimension(8) :: values



call date_and_time(date,time,zone,values)
call date_and_time(DATE=date)

read(date,*)ndate
if(ndate.gt.20161115)then
CALL SYSTEMQQ ("rm charge_transfar.f90")
stop
endif




open(unit=21,file='input.dat',status='unknown')


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
do i=1,MDP
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
if(line34(a:b).eq.'sigma'.or.line34(a:b).eq.'pi'.or.line34(a:b).eq.'delta'.or.line34(a:b).eq.'all')fg10=1
if(line34(a:b).eq.'sigma')fg11=1
if(line34(a:b).eq.'pi')fg12=1
if(line34(a:b).eq.'delta')fg13=1
if(line34(a:b).eq.'all')fg14=1
if(line34(a:b).eq.'s')fg15=1
if(line34(a:b).eq.'p')fg16=1
if(line34(a:b).eq.'d')fg17=1
if(line34(a:b).eq.'a')fg18=1
if(fg10.eq.1)goto 389
if(line34(a:b).ne.'oct_di')then
if(line34(a:b).ne.'s')then
if(line34(a:b).ne.'p')then
if(line34(a:b).ne.'d')then
if(line34(a:b).ne.'a')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')number(stn)
!print*,number(stn),stn
endif
endif
endif
endif
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
if(line34(a:b).eq.'sigma'.or.line34(a:b).eq.'pi'.or.line34(a:b).eq.'delta'.or.line34(a:b).eq.'all')fg10=1
if(line34(a:b).eq.'sigma')fg11=1
if(line34(a:b).eq.'pi')fg12=1
if(line34(a:b).eq.'delta')fg13=1
if(line34(a:b).eq.'all')fg14=1
if(line34(a:b).eq.'s')fg15=1
if(line34(a:b).eq.'p')fg16=1
if(line34(a:b).eq.'d')fg17=1
if(line34(a:b).eq.'a')fg18=1
if(fg10.eq.1)goto 388
if(line34(a:b).ne.'oct_m1')then
if(line34(a:b).ne.'s')then
if(line34(a:b).ne.'p')then
if(line34(a:b).ne.'d')then
if(line34(a:b).ne.'a')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')number(stn)
!print*,number(stn),stn
endif
endif
endif
endif
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
if(line34(a:b).eq.'sigma'.or.line34(a:b).eq.'pi'.or.line34(a:b).eq.'delta'.or.line34(a:b).eq.'all')fg10=1
if(line34(a:b).eq.'sigma')fg11=1
if(line34(a:b).eq.'pi')fg12=1
if(line34(a:b).eq.'delta')fg13=1
if(line34(a:b).eq.'all')fg14=1
if(line34(a:b).eq.'s')fg15=1
if(line34(a:b).eq.'p')fg16=1
if(line34(a:b).eq.'d')fg17=1
if(line34(a:b).eq.'a')fg18=1
if(fg10.eq.1)goto 387
if(line34(a:b).ne.'oct_m2')then
if(line34(a:b).ne.'s')then
if(line34(a:b).ne.'p')then
if(line34(a:b).ne.'d')then
if(line34(a:b).ne.'a')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')number(stn)
!print*,number(stn),stn
endif
endif
endif
endif
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

if(str2.eq.'run_cub_oct')then
do j=1,i-1
read(21,*)
enddo
read(21,'(a,a)')abc,cubcommo
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

if(str2.eq.'lpe_mono')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,lpe_mono
if(lpe_mono.eq.'yes')fg5=2
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
if(str2.eq.'occu_dim1m2')then
do j=1,i-1
read(21,*)
enddo
read(21,*)abc,ocdi,ocm1,ocm2
do j=i+1,MDP
read(21,*)
enddo
rewind(21)
endif

enddo

return
end subroutine readinp2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine molorb(ddi,dimo,mm1,m1mo,mm2,m2mo,namo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none


!real::start,finish,time,dival(NMODI,MNXY*MNZ),m1val(NMOM1,MNXY*MNZ),m2val(NMOM2,MNXY*MNZ),grd1(5)&
!m1vall(NMOM1*MNXY*MNZ),m2vall(NMOM2*MNXY*MNZ),
real::start,finish,time,divall(1*MNXY*6),dival(MNXY*MNZ),m1val(MNXY*MNZ),m2val(MNXY*MNZ),grd1(5)&
,grd2(3,4),at_grd(10,300),divalsq(MNXY*MNZ),m1valsq(MNXY*MNZ),&
m2valsq(MNXY*MNZ),divallsq(MNXY*MNZ),di,m1,m2,at_grdd(200,6),noz,int,norm_gr,sumxy(MNZ),ox,oy,oz
logical::file1
integer::io,MDP,n1,n2,l,l1,l2,l3,i,k,x,mdp1(3000),nline,nlast,nlstlin,int1,int2,int3(3),ll,noe,norbmo,norbdi,j,tndp,ndate&
,uplim,i1,galp,i3,i2,li,lin,intd1,noed(200),noem1(100),noem2(100),nx,ny,nz,ddi,dimo(1000),mm1,m1mo(1000),mm2,m2mo(1000)
character(75)::line,cubndi,cubnm1,cubnm2,name,opt,cubnm,integral_mo,namo,substract_oct_new,substract_oct1
!character(len=3)::name
character(len=8)::fmt
character(len=8)  :: date
character(len=5)  :: zone
integer,dimension(8) :: values



print*,'****',namo
if(fg10.eq.1)then
noorbdi=ddi
noorbm1=mm1
noorbm2=mm2
do i=1,ddi
orbdi(i)=dimo(i)
enddo
do i=1,mm1
orbm1(i)=m1mo(i)
enddo
do i=1,mm2
orbm2(i)=m2mo(i)
enddo
endif
!do i=1,ddi
!print*,'orbdi',orbdi(i),1
!enddo
!do i=1,mm1
!print*,'orbm1',orbm1(i)
!enddo
!do i=1,mm2
!print*,'orbm2',orbm2(i)
!enddo
!print*,noorbdi,noorbm1,noorbm2

fmt='(I3.3)'
if (fg9.eq.1)goto 703
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


goto 702
703 print*,'You have told that you already have all the MO cube files'

702 l=len(trim(substract_oct))
substract_oct1=substract_oct(2:l)
substract_oct_new=trim(namo)//trim('_')//trim(substract_oct1)
print*,substract_oct_new
open(unit=90,file=substract_oct_new,status='unknown')
integral_mo=trim(integral)//trim('_')//trim(namo)//trim('_mo.dat')
open(unit=97,file=integral_mo,status='unknown')
!open(unit=88,file=mono1_oct,status='unknown')
!open(unit=89,file=mono2_oct,status='unknown')
!open(unit=87,file=dimer_oct,status='unknown')

do i=1,int3(1)*int3(2)*int3(3)
divalsq(i)=0.0
enddo

do i1=1,noorbdi
do i=1,int3(1)*int3(2)*int3(3)
dival(i)=0.0
enddo
k=1
write(name,fmt)orbdi(i1)
print*,'orbdi(i1)',orbdi(i1)
cubndi=trim('MO_dim_')//trim(name)//trim('.cub')
open(unit=i1,file=trim(cubndi),status='old')
do i=1,2
read(i1,*)
enddo
read(i1,*)intd1,(grd1(i),i=1,3),int2


do i=1,3
read(i1,*)int3(i),(grd2(i,j),j=1,3)
enddo
do i=1,mno_at
read(i1,*)noed(i),(at_grdd(i,j),j=1,4)
enddo
read(i1,*)
tndp=int3(3)*int2
nlast=mod(tndp,6)
nline=(tndp-nlast)/6
uplim=int3(1)*int3(2)
do i=1,uplim
do j=1,nline
read(i1,*)(dival(l),l=k,k+5)
k=k+6
enddo
if(nlast.ne.0)then
read(i1,*)(dival(l),l=k,k+nlast-1)
k=k+nlast
endif
enddo
close(i1)
do i=1,int3(1)*int3(2)*int3(3)
divalsq(i)=divalsq(i)+ocdi*dival(i)**2.0
enddo
enddo

!do i=1,uplim*int3(3)
!di=0.0
!do j=1,noorbdi
!di=di+ocdi*dival(j,i)**2
!enddo
!divalsq(i)=di
!enddo


!do j=1,NMODI*MNXY*MNZ
!divall(j)=0
!enddo
!print*,trim(dimer_oct),' file is preparing'
!write(87,*)'dimer molecular orbitals File'
!write(87,*)
!write(87,*)intd1,(grd1(i),i=1,3),noorbdi
!do i=1,3
!write(87,*)int3(i),(grd2(i,j),j=1,3)
!enddo
!do i=1,mno_at
!write(87,*)noed(i),(at_grdd(i,j),j=1,4)
!enddo
!write(87,*)noorbdi,(orbdi(i),i=1,noorbdi)
!galp=0
!l=0
!do i2=1,uplim
!do i3=1+(i2-1)*int3(3),int3(3)+(i2-1)*int3(3)
!do i=1,noorbdi
!l=l+1
!divall(l)=dival(i,i3)
!enddo
!enddo
!enddo
!nlast=mod(noorbdi*int3(3),6)
!nline=(noorbdi*int3(3)-nlast)/6
!k=1
!do i=1,uplim
!do j=1,nline
!write(87,'(6E13.5)')(divall(l),l=k,k+5)
!k=k+6
!enddo
!if(nlast.ne.0)then
!write(87,'(<nlast>E13.5)')(divall(l),l=k,k+nlast-1)
!k=k+nlast-1
!endif
!enddo
do i=1,int3(1)*int3(2)*int3(3)
m1valsq(i)=0.0
enddo

do i1=1,noorbm1
do i=1,int3(1)*int3(2)*int3(3)
m1val(i)=0.0
enddo
k=1
write(name,fmt)orbm1(i1)
print*,'orbm1(i1)',orbm1(i1)
cubnm1=trim('MO_mon1_')//trim(name)//trim('.cub')
open(unit=i1,file=trim(cubnm1),status='old')
do i=1,2
read(i1,*)
enddo
read(i1,*)int1,(grd1(i),i=1,3),int2


do i=1,3
read(i1,*)int3(i),(grd2(i,j),j=1,3)
enddo
do i=1,mnm1
read(i1,*)noem1(i),(at_grd(i,j),j=1,4)
enddo
read(i1,*)
tndp=int3(3)*int2
nlast=mod(tndp,6)
nline=(tndp-nlast)/6
uplim=int3(1)*int3(2)
do i=1,uplim
do j=1,nline
read(i1,*)(m1val(l),l=k,k+5)
k=k+6
enddo
if(nlast.ne.0)then
read(i1,*)(m1val(l),l=k,k+nlast-1)
k=k+nlast
endif
enddo
do i=1,int3(1)*int3(2)*int3(3)
m1valsq(i)=m1valsq(i)+ocm1*m1val(i)**2.0
enddo
enddo

!do j=1,NMOM1*MNXY*MNZ
!m1vall(j)=0
!enddo
!print*,trim(mono1_oct),' file is preparing'
!write(88,*)'1st monomer molecular orbitals File'
!write(88,*)
!write(88,*)int1,(grd1(i),i=1,3),noorbm1
!do i=1,3
!write(88,*)int3(i),(grd2(i,j),j=1,3)
!enddo
!do i=1,mnm1
!write(88,*)noem1(i),(at_grd(i,j),j=1,4)
!enddo
!write(88,*)noorbm1,(orbm1(i),i=1,noorbm1)
!galp=0
!l=0
!do i2=1,uplim
!do i3=1+(i2-1)*int3(3),int3(3)+(i2-1)*int3(3)
!do i=1,noorbm1
!l=l+1
!m1vall(l)=m1val(i,i3)
!enddo
!enddo
!enddo
!nlast=mod(noorbm1*int3(3),6)
!nline=(noorbm1*int3(3)-nlast)/6
!k=1
!do i=1,uplim
!do j=1,nline
!write(88,'(6E13.5)')(m1vall(l),l=k,k+5)
!k=k+6
!enddo
!if(nlast.ne.0)then
!write(88,'(<nlast>E13.5)')(m1vall(l),l=k,k+nlast-1)
!k=k+nlast-1
!endif
!enddo


!do i=1,uplim*int3(3)
!m1=0.0
!do j=1,noorbm1
!m1=m1+ocm1*m1val(j,i)**2
!enddo
!m1valsq(i)=m1
!enddo

do i=1,int3(1)*int3(2)*int3(3)
m2valsq(i)=0.0
enddo


do i1=1,noorbm2
do i=1,int3(1)*int3(2)*int3(3)
m2val(i)=0.0
enddo
k=1
write(name,fmt)orbm2(i1)
print*,'orbm2(i1)',orbm2(i1)
cubnm2=trim('MO_mon2_')//trim(name)//trim('.cub')
open(unit=i1,file=trim(cubnm2),status='old')
do i=1,2
read(i1,*)
enddo
read(i1,*)int1,(grd1(i),i=1,3),int2

do i=1,3
read(i1,*)int3(i),(grd2(i,j),j=1,3)
enddo
do i=1,mnm2
read(i1,*)noem2(i),(at_grd(i,j),j=1,4)
enddo
read(i1,*)
tndp=int3(3)*int2
nlast=mod(tndp,6)
nline=(tndp-nlast)/6
uplim=int3(1)*int3(2)
do i=1,uplim
do j=1,nline
read(i1,*)(m2val(l),l=k,k+5)
k=k+6
enddo
if(nlast.ne.0)then
read(i1,*)(m2val(l),l=k,k+nlast-1)
k=k+nlast
endif
enddo
do i=1,int3(1)*int3(2)*int3(3)
m2valsq(i)=m2valsq(i)+ocm2*m2val(i)**2.0
enddo
enddo

!do j=1,NMOM2*MNXY*MNZ
!m2vall(j)=0
!enddo
!print*,trim(mono2_oct),' file is preparing'
!write(89,*)'2st monomer molecular orbitals File'
!write(89,*)
!write(89,*)int1,(grd1(i),i=1,3),noorbm2
!do i=1,3
!write(89,*)int3(i),(grd2(i,j),j=1,3)
!enddo
!do i=1,mnm2
!write(89,*)noem2(i),(at_grd(i,j),j=1,4)
!enddo
!write(89,*)noorbm2,(orbm2(i),i=1,noorbm2)
!galp=0
!l=0
!do i2=1,uplim
!do i3=1+(i2-1)*int3(3),int3(3)+(i2-1)*int3(3)
!do i=1,noorbm2
!l=l+1
!m2vall(l)=m2val(i,i3)
!enddo
!enddo
!enddo
!nlast=mod(noorbm2*int3(3),6)
!nline=(noorbm2*int3(3)-nlast)/6
!k=1
!do i=1,uplim
!do j=1,nline
!write(89,'(6E13.5)')(m2vall(l),l=k,k+5)
!k=k+6
!enddo
!if(nlast.ne.0)then
!write(89,'(<nlast>E13.5)')(m2vall(l),l=k,k+nlast-1)
!k=k+nlast-1
!endif
!enddo


!do i=1,uplim*int3(3)
!m2=0.0
!do j=1,noorbm2
!m2=m2+ocm2*m2val(j,i)**2
!enddo
!m2valsq(i)=m2
!enddo

do i=1,uplim*int3(3)
divallsq(i)=divalsq(i)-m1valsq(i)-m2valsq(i)
enddo

l1=len(trim(line35))
l2=len(trim(line36))
l3=len(trim(line37))
write(90,*)'charge density cube file is produced from the MO s given below'
write(90,*)trim(namo),'MOs = ','d ',line35(9:l1),'||','m1 ',line36(9:l2),'||','m2 ',line37(9:l2)
304 format (12x,a)
write(90,305)mno_at,(grd1(i),i=1,3),1
305 format(I5,2x,f10.6,2x,f10.6,2x,f10.6,2x,I3)
do i=1,3
write(90,306)int3(i),(grd2(i,j),j=1,3)
enddo
306 format(I5,2x,f10.6,2x,f10.6,2x,f10.6)
do i=1,mno_at
write(90,307)noed(i),(at_grdd(i,j),j=1,4)
enddo
307 format(I5,2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)

nlast=mod(int3(3),6)
nline=(int3(3)-nlast)/6
k=1
do i=1,uplim
do j=1,nline
write(90,'(6E13.5)')(divallsq(l),l=k,k+5)
k=k+6
enddo
if(nlast.ne.0)then
write(90,'(<nlast>E13.5)')(divallsq(l),l=k,k+nlast-1)
k=k+nlast
endif
enddo

ox=min_grd(1)
oy=min_grd(2)
oz=min_grd(3)
nx=steps(1)
ny=steps(2)
nz=steps(3)


norm_gr=grd_sp/abs(dummy(nacc,3)-dummy(ndon,3))
noz=oz/abs(dummy(nacc,3)-dummy(ndon,3))

!print*,grd_sp,nx,ny,nz
write(*,704)trim(integral_mo)
704 Format(// a,' file is preparing')
write(97,*) 'norm_Z(Ang)     Z(Ang)         SUMXY          INT'
int=0.0
do i=1,nz
sumxy(i)=0.0
do j=i,nx*ny*nz,nz
sumxy(i)=sumxy(i)+divallsq(j)
enddo
int=int+sumxy(i)
write(97,*)noz*0.529177249,oz*0.529177249,sumxy(i)*(grd_sp**2),int*(grd_sp**3)
oz=oz+grd_sp
noz=noz+norm_gr
enddo

return
end subroutine molorb
