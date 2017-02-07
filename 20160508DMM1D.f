      Program DMM1D
      implicit none
      integer i,j,s,n,te,l,gam,nu,c1,no2,rno2,c,nevents,infin,nmsp1
      integer lower,sp1,sm1,nms,s2,nm1,nm2sp1,nmsp2,upper,left,right
      integer sp(4),ne,c2,c3,no2p1,kl,mf,pmf
c                           ,c4
      parameter(s=10,n=s,no2=n/2,nmsp1=(n-s)+1,sp1=s+1)
      parameter(nm1=n-1,nm2sp1=nmsp1-s,nmsp2=(n-s)+2,s2=2*s,nms=n-s)
      parameter(sm1=s-1,no2p1=no2+1,nevents=1000000)
      integer y(n,3),state(n),neighbour(n),neighheap(n,2),cl(n)
      integer cl1(n)
c                   ,cl2(n)
      logical z(n)
c                 ,z1(n)
      real*8 t(n),tau,K,ti1,t2,tto1,tto2,a,b,d(1),Infinity,to
      real*8 ti2,t3,tto3,k0,k1,k2,k0s,k2s,ptau
      parameter(Infinity=Huge(1.d0))
      logical ch
      rno2=no2-(1-mod(n,2))
      open(unit=10,file='DMM1DOutput')
      open(unit=11,file="DMM1DWT")
c      close(11)
      do kl=1,1
      K=(kl-1)*(10.d0/300)
      te=1
      tto1=0.d0
      tto2=0.d0
      tto3=0.d0
      c1=0
      c2=0
      c3=0
c      c4=0
c      nevents=Huge(c1)
c             Huge(c1)
      k0s=dexp(-2*K)
      k0=dexp(-K)
      k1=1.d0
      k2=dexp(K)
      k2s=dexp(2*K)
      do gam=1,n
          call neighbours(gam,left,right)
          y(gam,1)=gam
          y(gam,2)=left
          y(gam,3)=right
          neighheap(gam,1)=left
          neighheap(gam,2)=right
c     tr was used to keep track of the Time Remaining when a unit's
c     state was "frozen"; it is not needed for the DMM with finite K.
c          tr(gam)=-1.d0
c          cl2(gam)=0
c          z1(gam)=.false.
      end do
c      call cpu_time(ti1)
      do i=1,te

c     Initialize the state of the system.  For the AB model, I choose 1
c     to represent the state A and 0 to represent the state B.

          j=0
c          c=0
c          z=.false.
c          cl=0
          neighbour=0
          mf=-n
          ptau=0.d0
          do gam=1,n
              call neighbours(gam,left,right)
              if((.5d0+kiss64()/18446744073709551616.d0).ge.0.5d0) then
	          state(gam)=1
                  mf=mf+2
                  neighbour(left)=neighbour(left)+1
                  neighbour(right)=neighbour(right)+1
              else
                  state(gam)=0
              end if
          end do
          do gam=1,n
              select case(neighbour(gam)+3*state(gam))
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
     &                       )/k2
                  case(4)
                      t(gam)=-dlog(.5d0+kiss64()/18446744073709551616.d0
     &                       )/k1
                  case(5)
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
c      write(*,*)j
c      write(*,*)'before',y(:,1)
          call heap1
c      write(*,*)'after',y(:,1),t(y(:,1))
c      call cpu_time(t3)
c      tto2=tto2+t3-ti2
c      c=0
c          do while(t(y(1,1))<1000000.d0.and.j.le.nevents)
          do j=1,nevents
c           j=j+1

c     Generate a new value for the smallest

              gam=y(1,1)
c      call cpu_time(ti1)
c          tau=minval(t)
c      call cpu_time(t2)
c      write(*,*)j,t(gam),tau
              tau=t(gam)
c      tto2=tto2+t3-ti2
c      tto3=tto3+t2-ti1
      to=tau
              left=y(1,2)
              right=y(1,3)
              if(state(gam)==1) then
                  state(gam)=0
                  call wtd(mf,2,tau,ptau,n)
                  neighbour(left)=neighbour(left)-1
                  neighbour(right)=neighbour(right)-1
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
                  call update1(left)
                  call update1(right)
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
                  call wtd(mf,-2,tau,ptau,n)
                  neighbour(left)=neighbour(left)+1
                  neighbour(right)=neighbour(right)+1
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
                  call update0(left)
                  call update0(right)
c                  t(gam)=tau-dlog(.5d0+kiss64()/18446744073709551616
c     &            .d0)
                  select case(neighbour(gam))
                      case(0)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k2+tau
                      case(1)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k1+tau
                      case(2)
                          t(gam)=-dlog(.5d0+kiss64()/1844674407370955161
     &                           6.d0)/k0+tau
                      case default
                          write(*,*)'Error in Min Update',neighbour(gam)
                  end select
              end if
              if(z(1).eqv..false.) then
                  c=c+1
                  cl(c)=1
                  z(1)=.true.
              end if
              ch=.true.
              
c     Update the heap
c      write(*,*)j,t(gam)
c      do l=1,no2
c          z(l)=.true.
c          cl(l)=no2p1-l
c      end do
c      call cpu_time(ti2)
          call heap1
c      call cpu_time(t3)
c      write(*,*)'after',y(:,1),t(y(:,1)),state,tau
c      cl=0
c      z=.false.
c      tto2=tto2+t3-ti2
c      c=0
          end do
      tto1=tto1+to
      c1=c1+j
      end do
c      call cpu_time(t2)
c      write(*,*) c1,tto1,tto2,t(y(1,1))
c                       ,tto3               ,minval(t)
      write(10,*)K,sum(state)
      end do
      close(11)
      close(10)
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
      subroutine neighbours(lq,left,right)
      integer, intent(in) :: lq
      integer, intent(out) :: left,right
      nu=0
      if(lq==1) nu=nu+1
      if(lq==n) nu=nu+2
      select case (nu)
          case (0)
              left=lq-1
              right=lq+1
          case (1)
              left=n
              right=2
          case (2)
              left=nm1
              right=1
          case default
              write(*,*)"Error in Neighbours"
      end select
      end subroutine neighbours

      subroutine update1(upper)
      integer upper
      if(upper==left) then
          c3=1
      else
          c3=2
      end if
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
          t(upper)=tau+(t(upper)-tau)*k2
          if(neighheap(1,c3).le.no2) then
              if(z(neighheap(1,c3)).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)
                  z(neighheap(1,c3))=.true.
              end if
          end if
      else
          t(upper)=tau+(t(upper)-tau)*k0
          if(neighheap(1,c3)>1) then
              if(z(neighheap(1,c3)/2).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)/2
                  z(neighheap(1,c3)/2)=.true.
              end if
          end if
      end if
      end subroutine update1

      subroutine update0(upper)
      integer upper
      if(upper==left) then
          c3=1
      else
          c3=2
      end if
      if(state(upper)==0) then
          t(upper)=tau+(t(upper)-tau)*k0
          if(neighheap(1,c3)>1) then
              if(z(neighheap(1,c3)/2).eqv..false.) then
                  c=c+1
                  cl(c)=neighheap(1,c3)/2
                  z(neighheap(1,c3)/2)=.true.
              end if
          end if
      else
          t(upper)=tau+(t(upper)-tau)*k2
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
              c=0
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
                      left=neighheap(3,1)
                      right=neighheap(3,2)
                      neighheap(neighheap(1,1),2)=3
                      neighheap(neighheap(1,2),1)=3
                      neighheap(left,2)=1
                      neighheap(right,1)=1
                      left=neighheap(1,1)
                      right=neighheap(1,2)
                      neighheap(1,1)=neighheap(3,1)
                      neighheap(1,2)=neighheap(3,2)
                      neighheap(3,1)=left
                      neighheap(3,2)=right
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
                      left=neighheap(2,1)
                      right=neighheap(2,2)
                      neighheap(neighheap(1,1),2)=2
                      neighheap(neighheap(1,2),1)=2
                      neighheap(left,2)=1
                      neighheap(right,1)=1
                      left=neighheap(1,1)
                      right=neighheap(1,2)
                      neighheap(1,1)=neighheap(2,1)
                      neighheap(1,2)=neighheap(2,2)
                      neighheap(2,1)=left
                      neighheap(2,2)=right
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
                      left=neighheap(3,1)
                      right=neighheap(3,2)
                      neighheap(neighheap(1,1),2)=3
                      neighheap(neighheap(1,2),1)=3
                      neighheap(left,2)=1
                      neighheap(right,1)=1
                      left=neighheap(1,1)
                      right=neighheap(1,2)
                      neighheap(1,1)=neighheap(3,1)
                      neighheap(1,2)=neighheap(3,2)
                      neighheap(3,1)=left
                      neighheap(3,2)=right
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
                      left=neighheap(gam+1,1)
                      right=neighheap(gam+1,2)
                      neighheap(neighheap(l,1),2)=gam+1
                      neighheap(neighheap(l,2),1)=gam+1
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(gam+1,1)=left
                      neighheap(gam+1,2)=right
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
                      left=neighheap(gam,1)
                      right=neighheap(gam,2)
                      neighheap(neighheap(l,1),2)=gam
                      neighheap(neighheap(l,2),1)=gam
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam,1)
                      neighheap(l,2)=neighheap(gam,2)
                      neighheap(gam,1)=left
                      neighheap(gam,2)=right
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
                      left=neighheap(gam+1,1)
                      right=neighheap(gam+1,2)
                      neighheap(neighheap(l,1),2)=gam+1
                      neighheap(neighheap(l,2),1)=gam+1
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(gam+1,1)=left
                      neighheap(gam+1,2)=right
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
                      left=neighheap(n,1)
                      right=neighheap(n,2)
                      neighheap(neighheap(no2,1),2)=n
                      neighheap(neighheap(no2,2),1)=n
                      neighheap(left,2)=no2
                      neighheap(right,1)=no2
                      left=neighheap(no2,1)
                      right=neighheap(no2,2)
                      neighheap(no2,1)=neighheap(n,1)
                      neighheap(no2,2)=neighheap(n,2)
                      neighheap(n,1)=left
                      neighheap(n,2)=right
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
                      left=neighheap(gam,1)
                      right=neighheap(gam,2)
                      neighheap(neighheap(l,1),2)=gam
                      neighheap(neighheap(l,2),1)=gam
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam,1)
                      neighheap(l,2)=neighheap(gam,2)
                      neighheap(gam,1)=left
                      neighheap(gam,2)=right
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
                      left=neighheap(gam+1,1)
                      right=neighheap(gam+1,2)
                      neighheap(neighheap(l,1),2)=gam+1
                      neighheap(neighheap(l,2),1)=gam+1
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(gam+1,1)=left
                      neighheap(gam+1,2)=right
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
                      left=neighheap(gam,1)
                      right=neighheap(gam,2)
                      neighheap(neighheap(l,1),2)=gam
                      neighheap(neighheap(l,2),1)=gam
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam,1)
                      neighheap(l,2)=neighheap(gam,2)
                      neighheap(gam,1)=left
                      neighheap(gam,2)=right
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
                      left=neighheap(gam+1,1)
                      right=neighheap(gam+1,2)
                      neighheap(neighheap(l,1),2)=gam+1
                      neighheap(neighheap(l,2),1)=gam+1
                      neighheap(left,2)=l
                      neighheap(right,1)=l
                      left=neighheap(l,1)
                      right=neighheap(l,2)
                      neighheap(l,1)=neighheap(gam+1,1)
                      neighheap(l,2)=neighheap(gam+1,2)
                      neighheap(gam+1,1)=left
                      neighheap(gam+1,2)=right
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
      subroutine wtd(mf,dmf,tau,ptau,n)
      integer, intent(inout) :: mf
      integer, intent(in) :: dmf,n
      real*8, intent(in) :: tau
      integer pmf
      real*8, intent(inout) :: ptau
      pmf=mf
      mf=mf+dmf
      if(mod(n,2)==0) then
       if(pmf==0.and.mf.ne.0) then
c        open(unit=11,file="DMM1DWT",status='OLD')
        write(11,*) tau-ptau
c        close(11)
       end if
       if(mf==0.and.pmf.ne.0) ptau=tau
      else
       if(pmf*mf<0) then
c        open(unit=11,file="DMM1DWT",status='OLD')
        write(11,*) tau-ptau
c        close(11)
        ptau=tau
       end if
      end if
      end subroutine wtd
      end Program
