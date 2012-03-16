program psim

implicit none
real dt, time
double precision x(10000,3), v(10000,3)
double precision B(3), E(3)
double precision amu, pi, q, mass, radius, ke
integer numofcounts
integer numofions
integer i,j,k,l
double precision xmin, xmax, ymin, ymax, zmin, zmax
double precision vel, dist

integer seed
real randnum
real theta, phi, thetamin, thetamax, thetarange, phimin, phimax, phirange
double precision veldefault

!  Initialize  !
dt = .01
time = 0
numofcounts = 5000
numofions = 10
radius = 252000.0
seed = 0

veldefault = -207954.066
B(1) = 0.0
B(2) = 0.0
B(3) = -325.0d-9
E(1) = -0.01
E(2) = 0.0
E(3) = 0.0
amu = 1.667e-27
pi = 3.14159265358979323846264338
q = 1.6022e-19
mass = 19.0 * amu
xmin = -3*radius
xmax = 3*radius
ymin = -12*radius
ymax = 3*radius
zmin = -3*radius
zmax = 3*radius

thetamin = -pi/6
thetamax = pi/6
thetarange = thetamax - thetamin
phimin = pi/3
phimax = 2*pi/3
phirange = phimax - phimin

do i = 1,numofions
   call random(seed,randnum)
   theta = thetamin + randnum * thetarange
   call random(seed,randnum)
   phi = phimin + randnum * phirange
   v(i,1) = veldefault * cos(theta) * cos(phi)
   v(i,2) = veldefault * cos(theta) * sin(phi)
   v(i,3) = veldefault * sin(theta)
   !  Puts the starting position of the particle randomly on a  !
   !  10000x10000 plane located at y=radius.  !
   call random(seed,randnum)
   x(i,1) = 5000 - randnum * 10000
   x(i,2) = -1 * radius + 1.0
   call random(seed,randnum)
   x(i,3) = 5000 - randnum * 10000
end do

open(10,file="tkexyz")

!  BEGIN MAIN LOOP  !
do i = 1,numofcounts
   time = time + dt
   do j = 1,numofions
      if (x(j,1).eq.0.0) then
         goto 20
      end if
      if ((x(j,1)**2.0 + x(j,2)**2.0 + x(j,3)**2.0).le.radius) then
         if (x(j,1).gt.xmax.or.x(j,1).lt.xmin.or.x(j,2).gt.ymax) then
            if (x(j,2).lt.ymin.or.x(j,3).gt.zmax.or.x(j,3).lt.zmin) then
               x(j,1) = 0.0
            end if
         end if
      end if

      x(j,1) = x(j,1) + v(j,1) * dt + 1/2 * (q/mass) * (E(1) + v(j,2) * B(3) - v(j,3) * B(2)) * dt**2.0
      x(j,2) = x(j,2) + v(j,2) * dt + 1/2 * (q/mass) * (E(2) + v(j,3) * B(1) - v(j,1) * B(3)) * dt**2.0
      x(j,3) = x(j,3) + v(j,3) * dt + 1/2 * (q/mass) * (E(3) + v(j,1) * B(2) - v(j,2) * B(1)) * dt**2.0

      v(j,1) = v(j,1) + (q/mass) * (E(1) + v(j,2) * B(3) - v(j,3) * B(2)) * dt
      v(j,2) = v(j,2) + (q/mass) * (E(2) + v(j,3) * B(1) - v(j,1) * B(3)) * dt
      v(j,3) = v(j,3) + (q/mass) * (E(3) + v(j,1) * B(2) - v(j,2) * B(1)) * dt

      ke = 0.5 * mass * (v(j,1) ** 2.0 + v(j,2) ** 2.0 + v(j,3) ** 2.0)/q

      vel = sqrt(v(j,1) ** 2.0 + v(j,2) ** 2.0 + v(j,3) ** 2.0)
      do k = 1,j-1
         dist = sqrt((x(j,1)-x(k,1)) ** 2.0 + (x(j,2)-x(k,2)) ** 2.0 + (x(j,3)-x(k,3)) ** 2.0)
         if (dist.le.10) then
            print*,"COLLISION:",j,k
            call random(seed,randnum)
            theta = pi/2 - randnum * pi
            call random(seed,randnum)
            phi = randnum * 2 * pi
            v(j,1) = vel * cos(theta) * cos(phi)
            v(j,2) = vel * cos(theta) * sin(phi)
            v(j,3) = vel * sin(theta)
            call random(seed,randnum)
            theta = pi/2 - randnum * pi
            call random(seed,randnum)
            phi = randnum * 2 * pi
            v(k,1) = vel * cos(theta) * cos(phi)
            v(k,2) = vel * cos(theta) * sin(phi)
            v(k,3) = vel * sin(theta)
         end if
      end do

      write (10,*) time, ke, x(j,1), x(j,2), x(j,3), j
      20 continue

   end do
end do

end program psim


subroutine random(seed,randnum)
integer seed
real randnum
randnum = rand(seed)
seed = INT(randnum*12954971)
return
end
