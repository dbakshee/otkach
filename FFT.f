program    FFT
!        use params
!  implicit none

   implicit double precision (a-h,o-z)
   character*(*) fname
   character*(*) FILEgt
   complex sum

   parameter (fname='w1Vr2_b_n_dos_400.txt')
   parameter (FILEgt='FFT_w1Vr2_400a.dat')
!        parameter (constBT=0.05*0.0862) ! meV:50mK
!        double precision  b(1:402), rho(1:402, 1:101), DoS(1:402, 1:101), F(1:501, 1:101)
!        double precision  b(1:402), rho(1:402, 1:201), DoS(1:402, 1:201), F(1:1001, 1:201)
   double precision  b(1:402), rho(1:402, 1:401), DoS(1:402, 1:401), F(1:1001, 1:401)
   open(3,file=fname,status='unknown',iostat=ierr)
   if (ierr.ne.0) then
      write (STDERR,*) 'Error opening file', fname, '; ierr=', ierr
      stop 1
   end if
!	write(*,*) 'File G(Ef)'
   Nn=400 ! 200 ! 100
   Nb=401

   Nq=1001 ! 501
   DoSav=0
   do ib=1,Nb
      do in=1,Nn
!        do ib=1,Nb
         read(3,*) b(ib),rho(ib,in),DoS(ib,in)
         if(ib.eq.1) DoSav=DoSav+DoS(ib,in)
!
!        DoS(ib,in)=z
      enddo
   enddo

   close(3)
   DoSav=DoSav/Nn
   open(4,file=FILEgt,status='replace',&
   &iostat=ierr)
   if(ierr.ne.0) then
      write(*,*) 'Error writing file ',FILEgt
      Stop
   endif
   hq=0.04
   do iq=1,Nq
      q=hq*(iq-1)
      do in=1,Nn
         sum=cmplx(0.,0.)
         do ib=50,Nb-1
            hb=b(ib+1)-b(ib)
!           sum = sum + hb*(0.5*( DoS(ib,in) + DoS(ib+1,in) )-DoSav)*exp( -cmplx( 0.0, q/b(ib) ) )/b(ib)/b(ib)
            sum = sum + hb*0.5*( (DoS(ib,in)-DoSav)*exp(-cmplx( 0.0, q/b(ib)))/b(ib)/b(ib) + (DoS(ib+1,in)-DoSav)*exp(-cmplx( 0.0, q/b(ib+1)))/b(ib+1)/b(ib+1) )
         enddo
         F(iq,in) = abs(sum)
      enddo
   enddo
!        do in=1,Nn
!        do iq=1,Nq
!           write(4,'(1x,g17.7$)') F(iq,in)
!        end do
!       write(4,'(1x)')
!        end do


   do in=1,Nn
      do iq=1,Nq
         write(4,*) rho(1, in), hq*(iq-1), F(iq,in)
!!           write(4,*) hq*(iq-1), rho(ib, in), F(iq,in)
      enddo
   enddo
   close(4)
   stop
end
