program main

!&&& part 1

! 1.0 include library
implicit none

real :: beta=0.99
real :: alpha=0.3
real :: delta=0.1
real :: gamma=0.5

! real,parameter :: theta=0.3
integer,parameter :: maxage=65
integer,parameter :: retage=45

real :: kmax=10.0
real :: kmin=0.0
integer,parameter :: kgrid=191

real :: gradkm(9)=(/ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 /)
real :: tol=0.001
! real :: tau=theta/(2.0+theta)	!tau=0, so pen=0 too
real :: L=real(retage-1)/real(maxage)
real :: kinitmax=0.5
real :: kinitmin=0.2
integer,parameter :: kinitgrid=100

integer,parameter :: maxiter1=100
integer,parameter :: maxiter2=200


real kspace(kgrid),kdiff,kinitspace(kinitgrid),kinitstep,kinit
real v(kgrid,maxage),v3(3),d1k(maxage)
real theta,tau

integer d(kgrid,maxage),d1(maxage),iter,initm

real K,w,r,pen,sum,tK,tkv(kinitgrid),kstep,fktr(kgrid),fktw(kgrid),gradk
integer js,jmin,jmax,jl,ju,i,age,kmaxindr(kgrid),kmaxindw(kgrid),gm,vmax

! dynamic
integer iss,it ! age for age index, it for time in transition
integer,parameter :: tr=220 ! total period for transition path

real klong(tr),vlong(kgrid,maxage,tr),d2k(maxage,tr),tklong(tr),kdifflong(tr)
real wlong(tr),rlong(tr),penlong(tr),fktrlong(kgrid,tr),fktwlong(kgrid,tr)
integer kmaxindrlong(kgrid,tr),kmaxindwlong(kgrid,tr)
integer d2(maxage,tr),dlong(kgrid,maxage,tr)


kstep=(kmax-kmin)/real(kgrid-1)
kspace=(/ ( kmin+real(i-1)*kstep, i=1,kgrid ) /)
kinitstep=(kinitmax-kinitmin)/real(kinitgrid-1)
kinitspace=(/ ( kinitmin+real(i-1)*kinitstep, i=1,kinitgrid ) /)

open(unit=7,file='C:\Users\NREM\Desktop\Dropbox\cversion\ckv')

888 format (4(I5,2X))

! kinit loop
! do initm=1,kinitgrid
! 	kinit=kinitspace(initm)
kinit=3.26


! gradk loop
do gm=9,9
	gradk=gradkm(gm)

! steady state loop
do iss=1,2
	if (iss==1) theta=0.0
	if (iss==2) theta=0.3
! 	tau=theta/(2.0+theta)
	tau=(maxage-retage+1)*theta/( retage-1+theta*(maxage-retage+1) )
	
	iter=0
	kdiff=10.0
	K=Kinit

do while( (kdiff>tol).and.(iter<maxiter1) )
	iter=iter+1
	w=(1.0-alpha)*(K**alpha)*(L**(-alpha))
	r=alpha*(K**(alpha-1.0))*(L**(1.0-alpha))-delta
	
	
	pen=theta*(1.0-tau)*w
	
	! upper limit for Kt+1 known after rt, wt and pen are known
	! Ct+Kt+1=(1+rt)Kt+(1-tau)wt 		workers upper limit fktw
	! Ct+Kt+1=(1+rt)Kt+pen			retirees upper limit fktr
	! all workers upper limit for Kt+1 are same
	! all retirees upper limit for Kt+1 are same
	
	do i=1,kgrid
		fktr(i)=(1.0+r)*kspace(i)+pen
		fktw(i)=(1.0+r)*kspace(i)+(1.0-tau)*w
		kmaxindr(i)=count(fktr(i)>kspace)
		kmaxindw(i)=count(fktw(i)>kspace)
		if (kmaxindr(i)==0) kmaxindr(i)=1
		if (kmaxindr(i)==0) print *, 'kmaxindr(i)=0'
		if (kmaxindw(i)==0) kmaxindw(i)=1
		if (kmaxindw(i)==0) print *, 'kmaxindw(i)=0'
	end do
	

	! value fn and decision rule for maxage
	d(:,maxage)=1
	v(:,maxage)=(/( fktr(i)**(1.0-gamma)/(1.0-gamma),i=1,kgrid )/)
	
	! age 64 - age 1
	do age=maxage-1,1,-1
		js=1
		do i=1,kgrid
			jmin=js
			if (age>retage-1)	jmax=kmaxindr(i)
			if (age<retage)	jmax=kmaxindw(i)
			do while ((jmax-jmin)>2)
				jl=floor(real(jmin+jmax)/2.0)
				ju=jl+1
				v3(1)=util1(i,kspace(jl))+beta*v(jl,age+1)
				v3(2)=util1(i,kspace(ju))+beta*v(ju,age+1)
				if (v3(2)>v3(1)) then
					jmin=jl
				else
					jmax=ju
				end if
			end do
			v3(1)=util1(i,kspace(jmin))+beta*v(jmin,age+1)
			v3(3)=util1(i,kspace(jmax))+beta*v(jmax,age+1)
			if (jmax>jmin) then
				v3(2)=util1(i,kspace(jmin+1))+beta*v(jmin+1,age+1)
			else
				v3(2)=v3(1)
			end if
			js=jmin+(maxloc(v3,dim=1))-1
				
			v(i,age)=maxval(v3)
			d(i,age)=js	! d(i,age) is position of optimal kt+1 when kt=kspace(i) at age age
			
			
		end do
	end do
	d1(1)=1	! d1(age) is position of starting kt at each age
	d1k(1)=0.0	! d1k(age) is exact value of starting kt at each age
	sum=0.0
	do age=2,maxage
		d1(age)=d(d1(age-1),age-1)
		d1k(age)=kspace(d1(age))
		sum=sum+d1k(age)
	end do
	tK=sum/real(maxage)
! 	tkv(initm)=tK
	kdiff=abs(K-tK)/K
	
		print *, 'iteration :', iter
		print *, 'K is:', K
		print *, 'tK is:', tK
		print *, 'kdiff is:', kdiff
	
	K=gradk*K+(1.0-gradk)*tK
! 	write (7,*) kinit, tk
end do ! end kdiff loop

! print *, 'gradk= ', gradk, 'kinit= ', kinit

if (kdiff<=tol) then
	print *, ' Steady State Converged, iter= ', iter, 'kdiff= ', kdiff
	print *, 'K= ', K, 'kgrid=', kgrid
	print *, ''
	print *, ''

 	if (any(d1==kgrid)) print *, 'Kt+1 has reached upper bound of state space'


! connecting steady states and transition path
	if (iss==1) then
		klong(1:2)=tk
		d2(:,2)=d1(:)		! initial capital at start of t=2 is same as initial steady state
		! d2(age,it) is beginning-of-period capital position each age at period it
	else if (iss==2) then
		klong(tr)=tk
		vlong(:,:,tr)=v(:,:) ! last stage value in transition is value in final steady state
	end if
	
end if ! end kdiff if





print *, ''
print *, ''
pause

end do ! end iss loop
end do ! end gradk loop
! end do ! end kinit loop





! iteration for transition
! do gm=9,9
    gradk=0.98

iter=0
kdiff=10.0
do it=3,tr-1
    klong(it)=klong(2)+(it-2)*(klong(tr)-klong(2))/(tr-2)
end do

do while( (kdiff>tol).and.(iter<maxiter2) )
	iter=iter+1
	do it=2,tr-1
		wlong(it)=(1.0-alpha)*(Klong(it)**alpha)*(L**(-alpha))
		rlong(it)=alpha*(Klong(it)**(alpha-1.0))*(L**(1.0-alpha))-delta
		penlong(it)=theta*(1.0-tau)*wlong(it)

		! upper limit for Ki+1 known after rt, wt and pen are known
		! here i is position of starting K and t marks the exact time during transition
		! cash at hand depends on both starting K and time specific interest and wage
		! C+Ki+1=(1+rt)Ki+(1-tau)wt 		workers upper limit fktw
		! C+Ki+1=(1+rt)Ki+PENt				retirees upper limit fktr
		do i=1,kgrid
			fktrlong(i,it)=( 1.0+rlong(it) )*kspace(i)+penlong(it)
			fktwlong(i,it)=( 1.0+rlong(it) )*kspace(i)+(1.0-tau)*wlong(it)
			kmaxindrlong(i,it)=count( fktrlong(i,it)>kspace )
			kmaxindwlong(i,it)=count( fktwlong(i,it)>kspace )
			if ( kmaxindrlong(i,it)==0 ) kmaxindrlong(i,it)=1
			if ( kmaxindrlong(i,it)==0 ) print *, 'kmaxindrlong(i,it)==0, i=',i,'it=',it
			if ( kmaxindwlong(i,it)==0 ) kmaxindwlong(i,it)=1
			if ( kmaxindwlong(i,it)==0 ) print *, 'kmaxindwlong(i,it)==0, i=',i,'it=',it
		end do
	end do

	! for the steady states we only need to compute for one cohort and map it to other cohorts
	! for the transition we need to compute separately for each cohort

	! maxage of agents who die during transition from t=2 to t=tr-1
	do it=tr-1,2,-1
		do i=1,kgrid
			dlong(i,maxage,it)=1
			vlong(i,maxage,it)=fktrlong(i,it)**(1.0-gamma)/(1.0-gamma)
		end do
	end do


	! remaining agents
	do it=tr-1,2,-1 ! each loop calculates the behaviors of all agents at that t
		do age=maxage-1,1,-1
			js=1
			do i=1,kgrid
				jmin=js
				if (age>retage-1)	jmax=kmaxindrlong(i,it)
				if (age<retage)		jmax=kmaxindwlong(i,it)
				do while ((jmax-jmin)>2)
					jl=floor(real(jmin+jmax)/2.0)
					ju=jl+1
					v3(1)=util2(i,kspace(jl))+beta*vlong(jl,age+1,it+1)
					v3(2)=util2(i,kspace(ju))+beta*vlong(ju,age+1,it+1)
					if (v3(2)>v3(1)) then
						jmin=jl
					else
						jmax=ju
					end if
				end do
				v3(1)=util2(i,kspace(jmin))+beta*vlong(jmin,age+1,it+1)
				v3(3)=util2(i,kspace(jmax))+beta*vlong(jmax,age+1,it+1)
				if (jmax>jmin) then
					v3(2)=util2(i,kspace(jmin+1))+beta*vlong(jmin+1,age+1,it+1)
				else
					v3(2)=v3(1)
				end if
				js=jmin+(maxloc(v3,dim=1))-1

				vlong(i,age,it)=maxval(v3)
				dlong(i,age,it)=js	! dlong(i,age,it) is position of optimal ki+1 when kt=kspace(i) at age age and time t


			end do ! end i loop
		end do ! end age loop
	end do ! end it loop
	
	d2(1,:)=1	! d2(age,it) is starting position of K at each age and time it
	d2k(1,:)=0.0		! d2k(age,it) is exact value of starting K at each age and time it
	do it=3,tr-1
		sum=0.0
		do age=2,maxage
			d2(age,it)=dlong(d2(age-1,it-1),age-1,it-1)
			d2k(age,it)=kspace(d2(age,it))
			sum=sum+d2k(age,it)
		end do
		tklong(it)=sum/real(maxage)
		kdifflong(it)=( tklong(it)-klong(it) )/klong(it)
	end do
	kdiff=maxval( abs(kdifflong) )
	
    print *, 'iteration :', iter
    print *, 'gradk :', gradk
    print *, 'kdiff is:', kdiff
	print *, ''
	print *, ''
	
    do it=3,tr-1
		klong(it)=gradk*klong(it)+(1.0-gradk)*tklong(it)
	end do

end do ! end kdiff loop

if (kdiff<=tol) then
	print *, 'Transition converged, iter= ', iter, 'kdiff= ', kdiff
    if (any( d2==kgrid )) print *, 'kt+1 has reached upper bound of state space'
end if

! end do ! end gradk loop

do it=1,tr
    write (7,*) klong(it)
end do

contains

real function util1(i1,ktp11)	!i1 is integer, ktp11 is real
implicit none
integer i1
real ct1,ktp11
if (age>retage-1) ct1=fktr(i)-ktp11
if (age<retage)	ct1=fktw(i)-ktp11
if (ct1>0) then
	util1=ct1**(1.0-gamma)/(1.0-gamma)
else	! this occurs when kmaxind=1 & all ct<=0
	util1=-10000.0
end if
end function



real function util2(i2,ktp22)
implicit none
integer i2
real ct2,ktp22
if (age>retage-1) ct2=fktrlong(i,it)-ktp22
if (age<retage) ct2=fktwlong(i,it)-ktp22
if (ct2>0) then
    util2=ct2**(1.0-gamma)/(1.0-gamma)
else
    util2=-10000.0
end if
end function



end program





