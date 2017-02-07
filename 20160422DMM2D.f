      Program DMM2D
      implicit none
      integer i,j,s,n,te,l,gam,nu,c1,no2,rno2,c,nevents,infin,nmsp1
      integer lower,sp1,sm1,nms,s2,nm1,nm2sp1,nmsp2,upper,left,right
      integer sp(4),ne,c2,c3,no2p1,kl
c                           ,c4
      parameter(s=10,n=s*s,no2=n/2,nmsp1=(n-s)+1,sp1=s+1)
      parameter(nm1=n-1,nm2sp1=nmsp1-s,nmsp2=(n-s)+2,s2=2*s,nms=n-s)
      parameter(sm1=s-1,no2p1=no2+1)
      integer y(n,5),state(n),neighbour(n),neighheap(n,4),cl(n)
      integer cl1(n)
c                   ,cl2(n)
      logical z(n)
c                 ,z1(n)
      real*8 t(n),tau,K,ti1,t2,tto1,tto2,a,b,d(1),Infinity,to
c                                                         ,tr(n)
      real*8 ti2,t3,tto3,k0,k1,k2,k3,k4
      parameter(Infinity=Huge(1.d0))
      logical ch
      rno2=no2-(1-mod(n,2))
      open(unit=10,file='DMM2DOutput')
      do kl=1,100
      K=(kl-1)*2.5d-2
      te=1
      tto1=0.d0
      tto2=0.d0
      tto3=0.d0
      c1=0
      c2=0
      c3=0
c      c4=0
      k0=dexp(-K)
      k1=dexp(-K/2.d0)
      k2=1.d0
      k3=dexp(K/2.d0)
      k4=dexp(K)
      nevents=Huge(c1)
c             Huge(c1)
      do gam=1,n
          call neighbours(gam,upper,left,right,lower)
          y(gam,1)=gam
          y(gam,2)=upper
          y(gam,3)=left
          y(gam,4)=right
          y(gam,5)=lower
          neighheap(gam,1)=upper
          neighheap(gam,2)=left
          neighheap(gam,3)=right
          neighheap(gam,4)=lower
c          tr(gam)=-1.d0
c          cl2(gam)=0
c          z1(gam)=.false.
      end do
c      call cpu_time(ti1)
      do i=1,te

c     Initialize the state of the system.  For the AB model, I choose 1
c     to represent the state A and 0 to represent the state B.

         j=0
c          gam=y(1,1)
          neighbour=0
          do gam=1,n
              call neighbours(gam,upper,left,right,lower)
              if((.5d0+kiss64()/18446744073709551616.d0).ge.0.5d0) then
	          state(gam)=1
                  neighbour(upper)=neighbour(upper)+1
                  neighbour(left)=neighbour(left)+1
                  neighbour(right)=neighbour(right)+1
                  neighbour(lower)=neighbour(lower)+1
              else
                  state(gam)=0
              end if
          end do
          do gam=1,n
              select case(neighbour(gam)+5*state(gam))
                  case(0)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k0
                  case(1)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k1
                  case(2)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k2
                  case(3)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k3
                  case(4)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k4
                  case(5)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k4
                  case(6)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k3
                  case(7)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k2
                  case(8)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k1
                  case(9)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k0
                  case default
                      write(*,*) 'Error in Initialization',j,
     &                    neighbour(gam)
              end select
          end do
          do l=1,no2
              z(l)=.true.
              cl(l)=no2p1-l
          end do
          c=no2
          ch=.true.
c          write(*,*)i,c4
c          write(*,*)cl2,z1
c          do l=1,c4
c              tr(cl2(l))=-1.d0
c              z1(cl2(l))=.false.
c              cl2(l)=0
c          end do
c          c4=0
c          write(*,*)i,j
c          write(*,*)tr,z1
c           do l=1,n
c               tr(l)=-1.d0
c           end do

c      Construct the heap
c      call cpu_time(ti2)
          call heap1
c      call cpu_time(t3)
c      tto2=tto2+t3-ti2
c      c=0
          do while(t(y(1,1))<1000000.d0.and.j.le.nevents)
           j=j+1

c     Generate a new value for the smallest

              gam=y(1,1)
c      call cpu_time(ti1)
c          tau=minval(t)
c      call cpu_time(t2)
              tau=t(gam)
c      tto2=tto2+t3-ti2
c      tto3=tto3+t2-ti1
      to=tau
              upper=y(1,2)
              left=y(1,3)
              right=y(1,4)
              lower=y(1,5)
              if(state(gam)==1) then
                  state(gam)=0
                  neighbour(upper)=neighbour(upper)-1
                  neighbour(left)=neighbour(left)-1
                  neighbour(right)=neighbour(right)-1
                  neighbour(lower)=neighbour(lower)-1
c                  if(neighheap(1,1).le.no2) then
c                  if(z(neighheap(1,1)).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,1)
c                      z(neighheap(1,1))=.true.
c                  end if
c                  end if
c                  if(neighheap(1,2).le.no2) then
c                  if(z(neighheap(1,2)).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,2)
c                      z(neighheap(1,2))=.true.
c                  end if
c                  end if
c                  if(neighheap(1,3).le.no2) then
c                  if(z(neighheap(1,3)).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,3)
c                      z(neighheap(1,3))=.true.
c                  end if
c                  end if
c                  if(neighheap(1,4).le.no2) then
c                  if(z(neighheap(1,4)).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,4)
c                      z(neighheap(1,4))=.true.
c                  end if
c                  end if
                  call update1(upper,1)
                  call update1(left,2)
                  call update1(right,3)
                  call update1(lower,4)
                  select case(neighbour(gam))
                      case(0)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k0+tau
                      case(1)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k1+tau
                      case(2)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k2+tau
                      case(3)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k3+tau
                      case(4)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k4+tau
                      case default
                          write(*,*)'Error in Min Update',neighbour(gam)
                  end select
c                  if(neighbour(gam)>0) then
c                      t(gam)=tau-dlog(.5d0+kiss64()/1844674407370955
c     &                1616.d0)*4.d0/neighbour(gam)
c                  else
c                          if(z1(gam).eqv..false.) then
c                          c4=c4+1
c                          cl2(c4)=gam
c                          z1(gam)=.true.
c                          end if
c                      tr(gam)=-4.d0*dlog(.5d0+kiss64()/1844674407370955
c     &                1616.d0)
c                      t(gam)=Infinity
c                  end if
              else
                  state(gam)=1
                  neighbour(upper)=neighbour(upper)+1
                  neighbour(left)=neighbour(left)+1
                  neighbour(right)=neighbour(right)+1
                  neighbour(lower)=neighbour(lower)+1
c                  if(neighheap(1,1)>1) then
c                  if(z(neighheap(1,1)/2).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,1)/2
c                      z(neighheap(1,1)/2)=.true.
c                  end if
c                  end if
c                  if(neighheap(1,2)>1) then
c                  if(z(neighheap(1,2)/2).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,2)/2
c                      z(neighheap(1,2)/2)=.true.
c                  end if
c                  end if
c                  if(neighheap(1,3)>1) then
c                  if(z(neighheap(1,3)/2).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,3)/2
c                      z(neighheap(1,3)/2)=.true.
c                  end if
c                  end if
c                  if(neighheap(1,4)>1) then
c                  if(z(neighheap(1,4)/2).eqv..false.) then
c                      c=c+1
c                      cl(c)=neighheap(1,4)/2
c                      z(neighheap(1,4)/2)=.true.
c                  end if
c                  end if
                  call update0(upper,1)
                  call update0(left,2)
                  call update0(right,3)
                  call update0(lower,4)
                  select case(neighbour(gam))
                      case(0)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k4+tau
                      case(1)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k3+tau
                      case(2)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k2+tau
                      case(3)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k1+tau
                      case(4)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k0+tau
                      case default
                          write(*,*)'Error in Min Update',neighbour(gam)
                  end select
c                  t(gam)=tau-dlog(.5d0+kiss64()/18446744073709551616
c     &            .d0)
              end if
              if(z(1).eqv..false.) then
                  c=c+1
                  cl(c)=1
                  z(1)=.true.
              end if
              ch=.true.
              
c     Update the heap

c      call cpu_time(ti2)
          call heap1
c      call cpu_time(t3)
c      tto2=tto2+t3-ti2
c      c=0
          end do
      tto1=tto1+to
      c1=c1+j
      end do
c      call cpu_time(t2)
c      write(*,*) c1,tto1,tto3,tto2,t(y(1,1)),minval(t)
      write(10,*)K,sum(state)
      end do
      contains

c     The kiss64 (Keep It Simple Stupid 64-bit) algorithm was downloaded
c     from http://fortranwiki.org/fortran/show/kiss64 on June 26th, 2015
c     .  It is a pseudorandom number generator, as explained in https://
c     groups.google.com/forum/#!topic/comp.lang.fortran/qFv18ql_WlU , it
c     was created by George Marsaglia and consists of a Multiply-With
c     -Carry, a Xorshift, and a Congruential algorithm.  The original
c     code was posted on comp.lang.fortran on February 28th, 2009.

        function kiss64()
          integer*8, save :: x, y, z, c
          integer*8 :: t, k, m, s, kiss64
          data x, y, z, c
     &          / 1234567890987654321_8,
     &            362436362436362436_8,
     &            1066149217761810_8,
     &            123456123456123456_8 /
          m(x,k) = ieor(x, ishft(x,k))  ! statement function
          s(x) = ishft(x, -63)          ! statement function
          t = ishft(x, 58) + c
          if (s(x) .eq. s(t)) then
             c = ishft(x, -6) + s(x)
          else
             c = ishft(x, -6) + 1 - s(x + t)
          endif
          x = t + x
          y = m(m(m(y,13_8),-17_8), 43_8)
          z = 6906969069_8 * z + 1234567
          kiss64 = x + y + z
        end function kiss64
      subroutine neighbours(lq,upper,left,right,lower)
      integer, intent(in) :: lq
      integer, intent(out) :: upper,left,right,lower
      nu=0
      if(lq.le.s) nu=nu+1
      if(mod(lq,s)==1) nu=nu+2
      if(mod(lq,s)==0) nu=nu+4
      if(lq.ge.nmsp1) nu=nu+8
      select case (nu)
          case (0)
              upper=lq-s
              left=lq-1
              right=lq+1
              lower=lq+s
          case (1)
              upper=lq+nms
              left=lq-1
              right=lq+1
              lower=lq+s
          case (2)
              upper=lq-s
              left=lq+sm1
              right=lq+1
              lower=lq+s
          case (3)
              upper=nmsp1
              left=s
              right=2
              lower=sp1
          case (4)
              upper=lq-s
              left=lq-1
              right=lq-sm1
              lower=lq+s
          case (5)
              upper=n
              left=sm1
              right=1
              lower=s2
          case (8)
              upper=lq-s
              left=lq-1
              right=lq+1
              lower=lq-nms
          case (10)
              upper=nmsp1-s
              left=n
              right=nmsp1+1
              lower=1
          case (12)
              upper=nms
              left=n-1
              right=nmsp1
              lower=s
          case default
              write(*,*)"Error in Neighbours"
      end select
      end subroutine neighbours

      subroutine update1(upper,c3)
      integer upper,c3
c                  if(state(upper)==0) then
c                  select case (neighbour(upper))
c                      case(0)
c                          if(z1(upper).eqv..false.) then
c                          c4=c4+1
c                          cl2(c4)=upper
c                          z1(upper)=.true.
c                          end if
c                          tr(upper)=t(upper)-tau
c                          t(upper)=Infinity
c                      case(1)
c                          t(upper)=2.d0*t(upper)-tau
c                      case(2)
c                          t(upper)=1.5d0*t(upper)-.5*tau
c                      case(3)
c                          t(upper)=1.3333333333333333333d0*t(upper)
c     &                            -.33333333333333333333d0*tau
c                      case default
c                          write(*,*)"Error in Update1",neighbour(upper)
c                  end select
c                  end if
      if(state(upper)==0) then
          t(upper)=tau+(t(upper)-tau)*k3
          if(neighheap(1,c3).le.no2) then
              if(z(neighheap(1,c3)).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)
                  z(neighheap(1,c3))=.true.
              end if
          end if
      else
          t(upper)=tau+(t(upper)-tau)*k1
          if(neighheap(1,c3)>1) then
              if(z(neighheap(1,c3)/2).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)/2
                  z(neighheap(1,c3)/2)=.true.
              end if
          end if
      end if
      end subroutine update1

      subroutine update0(upper,c3)
      integer upper,c3
c                  if(state(upper)==0) then
c                  select case (neighbour(upper))
c                      case(1)
c                          if(tr(upper)<0.d0) then
c                              tr(upper)=-4.d0*dlog(.5d0+
c     &                        kiss64()/18446744073709551616.d0)
c                          c4=c4+1
c                          cl2(c4)=upper
c                          z1(upper)=.true.
c                          end if
c                          t(upper)=tr(upper)+tau
c                      case(2)
c                          t(upper)=.5d0*t(upper)+.5d0*tau
c                      case(3)
c                          t(upper)=.66666666666666666666d0*t(upper)
c     &                       +.33333333333333333333d0*tau
c                      case(4)
c                          t(upper)=.75d0*t(upper)+.25d0*tau
c                      case default
c                          write(*,*)"Error in Update0",neighbour(upper)
c                  end select
c                  end if
      if(state(upper)==0) then
          t(upper)=tau+(t(upper)-tau)*k1
          if(neighheap(1,c3)>1) then
              if(z(neighheap(1,c3)/2).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)/2
                  z(neighheap(1,c3)/2)=.true.
              end if
          end if
      else
          t(upper)=tau+(t(upper)-tau)*k3
          if(neighheap(1,c3).le.no2) then
              if(z(neighheap(1,c3)).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)
                  z(neighheap(1,c3))=.true.
              end if
          end if
      end if
      end subroutine update0

      subroutine heap
          do while(ch)

c      Keep correcting the sub-heaps until they all satify the heap
c      condition.

              ch=.false.

c      Establish the heap condition for the topmost element

              if(z(1).eqv..true.) then
              tau=t(y(2,1))
              a=t(y(3,1))
              b=t(y(1,1))
              if(b>tau) then
                  if(a<tau) then
                      nu=y(1,1)
                      y(1,1)=y(3,1)
                      y(3,1)=nu
                      nu=y(1,2)
                      y(1,2)=y(3,2)
                      y(3,2)=nu
                      nu=y(1,3)
                      y(1,3)=y(3,3)
                      y(3,3)=nu
                      nu=y(1,4)
                      y(1,4)=y(3,4)
                      y(3,4)=nu
                      nu=y(1,5)
                      y(1,5)=y(3,5)
                      y(3,5)=nu
                      upper=neighheap(3,1)
                      left=neighheap(3,2)
                      right=neighheap(3,3)
                      lower=neighheap(3,4)
                      neighheap(neighheap(1,1),4)=3
                      neighheap(neighheap(1,2),3)=3
                      neighheap(neighheap(1,3),2)=3
                      neighheap(neighheap(1,4),1)=3
                      neighheap(upper,4)=1
                      neighheap(left,3)=1
                      neighheap(right,2)=1
                      neighheap(lower,1)=1
                      upper=neighheap(1,1)
                      left=neighheap(1,2)
                      right=neighheap(1,3)
                      lower=neighheap(1,4)
                      neighheap(1,1)=neighheap(3,1)
                      neighheap(1,2)=neighheap(3,2)
                      neighheap(1,3)=neighheap(3,3)
                      neighheap(1,4)=neighheap(3,4)
                      neighheap(3,1)=upper
                      neighheap(3,2)=left
                      neighheap(3,3)=right
                      neighheap(3,4)=lower
                      ch=.true.

c      Since a change was made, the heap condition was not satisfied at
c      the beginning of this iteration, so run it again.

                      z(3)=.true.
                  else
                      nu=y(1,1)
                      y(1,1)=y(2,1)
                      y(2,1)=nu
                      nu=y(1,2)
                      y(1,2)=y(2,2)
                      y(2,2)=nu
                      nu=y(1,3)
                      y(1,3)=y(2,3)
                      y(2,3)=nu
                      nu=y(1,4)
                      y(1,4)=y(2,4)
                      y(2,4)=nu
                      nu=y(1,5)
                      y(1,5)=y(2,5)
                      y(2,5)=nu
                      upper=neighheap(2,1)
                      left=neighheap(2,2)
                      right=neighheap(2,3)
                      lower=neighheap(2,4)
                      neighheap(neighheap(1,1),4)=2
                      neighheap(neighheap(1,2),3)=2
                      neighheap(neighheap(1,3),2)=2
                      neighheap(neighheap(1,4),1)=2
                      neighheap(upper,4)=1
                      neighheap(left,3)=1
                      neighheap(right,2)=1
                      neighheap(lower,1)=1
                      upper=neighheap(1,1)
                      left=neighheap(1,2)
                      right=neighheap(1,3)
                      lower=neighheap(1,4)
                      neighheap(1,1)=neighheap(2,1)
                      neighheap(1,2)=neighheap(2,2)
                      neighheap(1,3)=neighheap(2,3)
                      neighheap(1,4)=neighheap(2,4)
                      neighheap(2,1)=upper
                      neighheap(2,2)=left
                      neighheap(2,3)=right
                      neighheap(2,4)=lower
                      ch=.true.
                      z(2)=.true.
                  end if
              else
                  if(a<b) then
                      nu=y(1,1)
                      y(1,1)=y(3,1)
                      y(3,1)=nu
                      nu=y(1,2)
                      y(1,2)=y(3,2)
                      y(3,2)=nu
                      nu=y(1,3)
                      y(1,3)=y(3,3)
                      y(3,3)=nu
                      nu=y(1,4)
                      y(1,4)=y(3,4)
                      y(3,4)=nu
                      nu=y(1,5)
                      y(1,5)=y(3,5)
                      y(3,5)=nu
                      upper=neighheap(3,1)
                      left=neighheap(3,2)
                      right=neighheap(3,3)
                      lower=neighheap(3,4)
                      neighheap(neighheap(1,1),4)=3
                      neighheap(neighheap(1,2),3)=3
                      neighheap(neighheap(1,3),2)=3
                      neighheap(neighheap(1,4),1)=3
                      neighheap(upper,4)=1
                      neighheap(left,3)=1
                      neighheap(right,2)=1
                      neighheap(lower,1)=1
                      upper=neighheap(1,1)
                      left=neighheap(1,2)
                      right=neighheap(1,3)
                      lower=neighheap(1,4)
                      neighheap(1,1)=neighheap(3,1)
                      neighheap(1,2)=neighheap(3,2)
                      neighheap(1,3)=neighheap(3,3)
                      neighheap(1,4)=neighheap(3,4)
                      neighheap(3,1)=upper
                      neighheap(3,2)=left
                      neighheap(3,3)=right
                      neighheap(3,4)=lower
                      ch=.true.
                      z(3)=.true.

c      Since nothing had to be changed, the heap condition may have been
c      satisfied at the beginning of this iteration.

                  end if
              end if
              z(1)=.false.
              end if

c      Establish the heap condition for each of the other sub-heaps,
c      except perhaps the last.

              do l=2,rno2
c              do l=2,ne/2
              if(z(l).eqv..true.) then
              gam=2*l
              tau=t(y(gam,1))
              a=t(y(gam+1,1))
              b=t(y(l,1))
              if(b>tau) then
                  if(a<tau) then
                      nu=y(l,1)
                      y(l,1)=y(gam+1,1)
                      y(gam+1,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam+1,2)
                      y(gam+1,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam+1,3)
                      y(gam+1,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam+1,4)
                      y(gam+1,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam+1,5)
                      y(gam+1,5)=nu
                      upper=neighheap(gam+1,1)
                      left=neighheap(gam+1,2)
                      right=neighheap(gam+1,3)
                      lower=neighheap(gam+1,4)
                      neighheap(neighheap(l,1),4)=gam+1
                      neighheap(neighheap(l,2),3)=gam+1
                      neighheap(neighheap(l,3),2)=gam+1
                      neighheap(neighheap(l,4),1)=gam+1
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(l,3)=neighheap(gam+1,3)
                      neighheap(l,4)=neighheap(gam+1,4)
                      neighheap(gam+1,1)=upper
                      neighheap(gam+1,2)=left
                      neighheap(gam+1,3)=right
                      neighheap(gam+1,4)=lower
                      ch=.true.
                      z(l/2)=.true.
                      z(gam+1)=.true.
                  else
                      nu=y(l,1)
                      y(l,1)=y(gam,1)
                      y(gam,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam,2)
                      y(gam,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam,3)
                      y(gam,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam,4)
                      y(gam,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam,5)
                      y(gam,5)=nu
                      upper=neighheap(gam,1)
                      left=neighheap(gam,2)
                      right=neighheap(gam,3)
                      lower=neighheap(gam,4)
                      neighheap(neighheap(l,1),4)=gam
                      neighheap(neighheap(l,2),3)=gam
                      neighheap(neighheap(l,3),2)=gam
                      neighheap(neighheap(l,4),1)=gam
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam,1)
                      neighheap(l,2)=neighheap(gam,2)
                      neighheap(l,3)=neighheap(gam,3)
                      neighheap(l,4)=neighheap(gam,4)
                      neighheap(gam,1)=upper
                      neighheap(gam,2)=left
                      neighheap(gam,3)=right
                      neighheap(gam,4)=lower
                      ch=.true.
                      z(l/2)=.true.
                      z(gam)=.true.
                  end if
              else
                  if(a<b) then
                      nu=y(l,1)
                      y(l,1)=y(gam+1,1)
                      y(gam+1,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam+1,2)
                      y(gam+1,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam+1,3)
                      y(gam+1,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam+1,4)
                      y(gam+1,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam+1,5)
                      y(gam+1,5)=nu
                      upper=neighheap(gam+1,1)
                      left=neighheap(gam+1,2)
                      right=neighheap(gam+1,3)
                      lower=neighheap(gam+1,4)
                      neighheap(neighheap(l,1),4)=gam+1
                      neighheap(neighheap(l,2),3)=gam+1
                      neighheap(neighheap(l,3),2)=gam+1
                      neighheap(neighheap(l,4),1)=gam+1
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(l,3)=neighheap(gam+1,3)
                      neighheap(l,4)=neighheap(gam+1,4)
                      neighheap(gam+1,1)=upper
                      neighheap(gam+1,2)=left
                      neighheap(gam+1,3)=right
                      neighheap(gam+1,4)=lower
                      ch=.true.
                      z(l/2)=.true.
                      z(gam+1)=.true.
                  end if
              end if
              z(l)=.false.
              end if
              end do

c      Establish the heap condition for the last sub-heap, if not
c      already done.

              if(rno2-no2==-1) then
              if(z(no2).eqv..true.) then
                  if(t(y(n,1))<t(y(no2,1))) then
                      nu=y(no2,1)
                      y(no2,1)=y(n,1)
                      y(n,1)=nu
                      nu=y(no2,2)
                      y(no2,2)=y(n,2)
                      y(n,2)=nu
                      nu=y(no2,3)
                      y(no2,3)=y(n,3)
                      y(n,3)=nu
                      nu=y(no2,4)
                      y(no2,4)=y(n,4)
                      y(n,4)=nu
                      nu=y(no2,5)
                      y(no2,5)=y(n,5)
                      y(n,5)=nu
                      upper=neighheap(n,1)
                      left=neighheap(n,2)
                      right=neighheap(n,3)
                      lower=neighheap(n,4)
                      neighheap(neighheap(no2,1),4)=n
                      neighheap(neighheap(no2,2),3)=n
                      neighheap(neighheap(no2,3),2)=n
                      neighheap(neighheap(no2,4),1)=n
                      neighheap(upper,4)=no2
                      neighheap(left,3)=no2
                      neighheap(right,2)=no2
                      neighheap(lower,1)=no2
                      upper=neighheap(no2,1)
                      left=neighheap(no2,2)
                      right=neighheap(no2,3)
                      lower=neighheap(no2,4)
                      neighheap(no2,1)=neighheap(n,1)
                      neighheap(no2,2)=neighheap(n,2)
                      neighheap(no2,3)=neighheap(n,3)
                      neighheap(no2,4)=neighheap(n,4)
                      neighheap(n,1)=upper
                      neighheap(n,2)=left
                      neighheap(n,3)=right
                      neighheap(n,4)=lower
                      z(no2/2)=.true.
                      ch=.true.
                  end if
                  z(no2)=.false.
              end if
              end if
          end do
      end subroutine heap

      subroutine heap1
      do while(c>0)

c      Keep correcting the sub-heaps until they all satify the heap
c      condition.

          do c3=1,c
              cl1(c3)=cl(c3)
          end do
          c2=c
          c=0
          do c3=1,c2
              l=cl1(c3)
              if(l.le.no2) then
              gam=2*l
              ch=gam+1>n
              tau=t(y(gam,1))
              b=t(y(l,1))
              if(b>tau) then
                  if(ch) then
                      nu=y(l,1)
                      y(l,1)=y(gam,1)
                      y(gam,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam,2)
                      y(gam,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam,3)
                      y(gam,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam,4)
                      y(gam,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam,5)
                      y(gam,5)=nu
                      upper=neighheap(gam,1)
                      left=neighheap(gam,2)
                      right=neighheap(gam,3)
                      lower=neighheap(gam,4)
                      neighheap(neighheap(l,1),4)=gam
                      neighheap(neighheap(l,2),3)=gam
                      neighheap(neighheap(l,3),2)=gam
                      neighheap(neighheap(l,4),1)=gam
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam,1)
                      neighheap(l,2)=neighheap(gam,2)
                      neighheap(l,3)=neighheap(gam,3)
                      neighheap(l,4)=neighheap(gam,4)
                      neighheap(gam,1)=upper
                      neighheap(gam,2)=left
                      neighheap(gam,3)=right
                      neighheap(gam,4)=lower
                      if(l.ne.1) then
                      if(z(l/2).eqv..false.) then
                          c=c+1
                          cl(c)=l/2
                          z(l/2)=.true.
                      end if
                      end if
                      if(z(gam).eqv..false.) then
                          c=c+1
                          cl(c)=gam
                          z(gam)=.true.
                      end if
                  else
                  a=t(y(gam+1,1))
                  if(a<tau) then
                      nu=y(l,1)
                      y(l,1)=y(gam+1,1)
                      y(gam+1,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam+1,2)
                      y(gam+1,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam+1,3)
                      y(gam+1,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam+1,4)
                      y(gam+1,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam+1,5)
                      y(gam+1,5)=nu
                      upper=neighheap(gam+1,1)
                      left=neighheap(gam+1,2)
                      right=neighheap(gam+1,3)
                      lower=neighheap(gam+1,4)
                      neighheap(neighheap(l,1),4)=gam+1
                      neighheap(neighheap(l,2),3)=gam+1
                      neighheap(neighheap(l,3),2)=gam+1
                      neighheap(neighheap(l,4),1)=gam+1
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(l,3)=neighheap(gam+1,3)
                      neighheap(l,4)=neighheap(gam+1,4)
                      neighheap(gam+1,1)=upper
                      neighheap(gam+1,2)=left
                      neighheap(gam+1,3)=right
                      neighheap(gam+1,4)=lower
                      if(l.ne.1) then
                      if(z(l/2).eqv..false.) then
                          c=c+1
                          cl(c)=l/2
                          z(l/2)=.true.
                      end if
                      end if
                      if(z(gam+1).eqv..false.) then
                          c=c+1
                          cl(c)=gam+1
                          z(gam+1)=.true.
                      end if
                  else
                      nu=y(l,1)
                      y(l,1)=y(gam,1)
                      y(gam,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam,2)
                      y(gam,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam,3)
                      y(gam,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam,4)
                      y(gam,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam,5)
                      y(gam,5)=nu
                      upper=neighheap(gam,1)
                      left=neighheap(gam,2)
                      right=neighheap(gam,3)
                      lower=neighheap(gam,4)
                      neighheap(neighheap(l,1),4)=gam
                      neighheap(neighheap(l,2),3)=gam
                      neighheap(neighheap(l,3),2)=gam
                      neighheap(neighheap(l,4),1)=gam
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam,1)
                      neighheap(l,2)=neighheap(gam,2)
                      neighheap(l,3)=neighheap(gam,3)
                      neighheap(l,4)=neighheap(gam,4)
                      neighheap(gam,1)=upper
                      neighheap(gam,2)=left
                      neighheap(gam,3)=right
                      neighheap(gam,4)=lower
                      if(l.ne.1) then
                      if(z(l/2).eqv..false.) then
                          c=c+1
                          cl(c)=l/2
                          z(l/2)=.true.
                      end if
                      end if
                      if(z(gam).eqv..false.) then
                          c=c+1
                          cl(c)=gam
                          z(gam)=.true.
                      end if
                  end if
                end if
              else
                  if(.not.ch) then
                  a=t(y(gam+1,1))
                  if(a<b) then
                      nu=y(l,1)
                      y(l,1)=y(gam+1,1)
                      y(gam+1,1)=nu
                      nu=y(l,2)
                      y(l,2)=y(gam+1,2)
                      y(gam+1,2)=nu
                      nu=y(l,3)
                      y(l,3)=y(gam+1,3)
                      y(gam+1,3)=nu
                      nu=y(l,4)
                      y(l,4)=y(gam+1,4)
                      y(gam+1,4)=nu
                      nu=y(l,5)
                      y(l,5)=y(gam+1,5)
                      y(gam+1,5)=nu
                      upper=neighheap(gam+1,1)
                      left=neighheap(gam+1,2)
                      right=neighheap(gam+1,3)
                      lower=neighheap(gam+1,4)
                      neighheap(neighheap(l,1),4)=gam+1
                      neighheap(neighheap(l,2),3)=gam+1
                      neighheap(neighheap(l,3),2)=gam+1
                      neighheap(neighheap(l,4),1)=gam+1
                      neighheap(upper,4)=l
                      neighheap(left,3)=l
                      neighheap(right,2)=l
                      neighheap(lower,1)=l
                      upper=neighheap(l,1)
                      left=neighheap(l,2)
                      right=neighheap(l,3)
                      lower=neighheap(l,4)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(l,3)=neighheap(gam+1,3)
                      neighheap(l,4)=neighheap(gam+1,4)
                      neighheap(gam+1,1)=upper
                      neighheap(gam+1,2)=left
                      neighheap(gam+1,3)=right
                      neighheap(gam+1,4)=lower
                      if(l.ne.1) then
                      if(z(l/2).eqv..false.) then
                          c=c+1
                          cl(c)=l/2
                          z(l/2)=.true.
                      end if
                      end if
                      if(z(gam+1).eqv..false.) then
                          c=c+1
                          cl(c)=gam+1
                          z(gam+1)=.true.
                      end if
                  end if
                  end if
              end if
              end if
              z(l)=.false.
          end do

c      Establish the heap condition for the last sub-heap, if not
c      already done.
      end do
      end subroutine heap1
      end Program
