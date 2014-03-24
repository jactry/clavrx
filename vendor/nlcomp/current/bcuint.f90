	SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
	
	
	IMPLICIT NONE
	REAL, DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
	REAL, INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
	REAL, INTENT(OUT) :: ansy,ansy1,ansy2
	INTEGER :: i
	REAL :: t,u
	REAL, DIMENSION(4,4) :: c
	call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
	if (x1u == x1l .or. x2u == x2l) print*,'error'
	t=(x1-x1l)/(x1u-x1l)
	u=(x2-x2l)/(x2u-x2l)
	ansy=0.0
	ansy2=0.0
	ansy1=0.0
	do i=4,1,-1
		ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
		ansy2=t*ansy2+(3.0*c(i,4)*u+2.0 *c(i,3))*u+c(i,2)
		ansy1=u*ansy1+(3.0*c(4,i)*t+2.0 *c(3,i))*t+c(2,i)
	end do
	ansy1=ansy1/(x1u-x1l)
	ansy2=ansy2/(x2u-x2l)
	END SUBROUTINE bcuint
