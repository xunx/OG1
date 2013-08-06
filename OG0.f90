program main

!&&& part 1

! 1.0 include library
implicit none

real :: beta=0.99
real :: alpha=0.3
real :: delta=0.1
real :: gamma=0.5

real,parameter :: theta=0.3
integer,parameter :: maxage=3
integer,parameter :: retage=3

real :: kmax=5.5
real :: kmin=0.0
integer,parameter :: kgrid=100

real :: gradkm(9)=(/ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 /)
real :: tol=0.001
real :: tau=theta/(2.0+theta)	!tau=0, so pen=0 too
real :: L=real(retage-1)/real(maxage)
real :: kinitmax=0.5
real :: kinitmin=0.2
integer,parameter :: kinitgrid=100

integer,parameter :: maxiter=1000


real kspace(kgrid),kdiff,kinitspace(kinitgrid),kinitstep,kinit
real v(kgrid,maxage),v3(3),kcross(maxage)

integer d(kgrid,maxage),d1(maxage),iter,initm

real K,w,r,pen,sum,tK,tkv(kinitgrid),kstep,fktr(kgrid),fktw(kgrid),gradk
integer js,jmin,jmax,jl,ju,i,t,kmaxindr(kgrid),kmaxindw(kgrid),gm

kstep=(kmax-kmin)/real(kgrid-1)
kspace=(/ ( kmin+real(i-1)*kstep, i=1,kgrid ) /)
kinitstep=(kinitmax-kinitmin)/real(kinitgrid-1)
kinitspace=(/ ( kinitmin+real(i-1)*kinitstep, i=1,kinitgrid ) /)

open(unit=7,file='C:\Users\NREM\Desktop\Dropbox\cversion\ckv')

888 format (4(I5,2X))

! do initm=1,kinitgrid
! 	kinit=kinitspace(initm)
kinit=0.309

do gm=8,8
	gradk=gradkm(gm)

	iter=0
	kdiff=10.0
	K=Kinit

do while( (kdiff>tol).and.(iter<maxiter) )
	sum=0.0
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
	v(:,maxage)=(/( ((1.0+r)*kspace(i)+pen)**(1.0-gamma)/(1.0-gamma),i=1,kgrid )/)
	
	! age 64 - age 1
	do t=maxage-1,1,-1
		js=1
		do i=1,kgrid
			jmin=js
			if (t>retage-1)	jmax=kmaxindr(i)
			if (t<retage)	jmax=kmaxindw(i)
			do while ((jmax-jmin)>2)
				jl=floor(real(jmin+jmax)/2.0)
				ju=jl+1
				v3(1)=util(i,kspace(jl))+beta*v(jl,t+1)
				v3(2)=util(i,kspace(ju))+beta*v(ju,t+1)
				if (v3(2)>v3(1)) then
					jmin=jl
				else
					jmax=ju
				end if
			end do
			v3(1)=util(i,kspace(jmin))+beta*v(jmin,t+1)
			v3(3)=util(i,kspace(jmax))+beta*v(jmax,t+1)
			if (jmax>jmin) then
				v3(2)=util(i,kspace(jmin+1))+beta*v(jmin+1,t+1)
			else
				v3(2)=v3(1)
			end if
			js=jmin+(maxloc(v3,dim=1))-1
				
			v(i,t)=maxval(v3)
			d(i,t)=js	! d(i,t) is position of optimal kt+1 when kt=kspace(i) at age t
			
			if ( (d(i,t)==kgrid).and.(t>1) ) then	! kt+1 reached upperbound for age>1
				print *, 'kt+1 reached upper bound'
				print *, 'i=', i, 't=', t
				pause
			else if ( (d(i,t)==kgrid).and.(t==1).and.(i==1) ) then	! kt+1 reached upperbound for zero initial asset at age 1
				print *, 'kt+1 reached upper bound'
				print *, 'i=', i, 't=', t
				pause
			end if
			
		end do
	end do
	d1(1)=1	! d1(t) is position of starting kt at each age
	kcross(1)=0.0	! kcross(t) is exact value of starting kt at each age
	do t=2,maxage
		d1(t)=d(d1(t-1),t-1)
		kcross(t)=kspace(d1(t))
		sum=sum+kcross(t)
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
end do

! print *, 'gradk= ', gradk, 'kinit= ', kinit

if (kdiff<=tol) then
	print *, 'Converged, iter= ', iter, 'kdiff= ', kdiff
	print *, 'K= ', K, 'kgrid=', kgrid
	print *, ''
	print *, ''

! 	if (any(d==kgrid)) print *, 'Kt+1 has reached upper bound of state space'
	if (any(d(:,1:maxage-1)==1)) print *, 'Kt+1 has reached lower bound of state space'
	! print *, 'd(kgrid,maxage) is:'
	! do i=1,kgrid
	! 	write (7,888) d(i,:)
	! 	print *, d(i,:)
	! 	write (7,*) d(i,:)
	! end do

end if





print *, ''
print *, ''

end do ! end gradk loop
! end do ! end kinit loop


contains

real function util(i1,ktp11)	!i1 is integer, ktp11 is real
implicit none
integer i1
real ct1,ktp11
if (t>retage-1) ct1=fktr(i)-ktp11
if (t<retage)	ct1=fktw(i)-ktp11

if (ct1>0) then
	util=ct1**(1.0-gamma)/(1.0-gamma)
else	! this occurs when kmaxind=1 & all ct<=0
	util=-10000.0
end if

end function



end program





