        real*8 u(1000), ut(1000)
        real*8 t,h,nu,c,time,sum
        n=1000
        t=0.00005
        h=0.01
        nu=1
            open(10,file='test.dat')
        do i=1,n
        x=h*(i-1-n/2)
        u(i)=exp(-x**2/(0.5)**2)
        enddo

        c=t*nu/h**2
        time=0

        mt=10
        mtt=2000

        do j=1,mt
                do jj=1,mtt
                    ut(1)=u(1)+c*(u(2)-2*u(1)+u(2))
                    do i=2,n-1
                    ut(i)=u(i)+c*(u(i-1)-2*u(i)+u(i+1))
                    enddo
                    i=n
                    ut(n)=u(n)+c*(u(n-1)-2*u(n)+u(n-1))
                time=time+t
                    do i=1,n
                    u(i)=ut(i)
                    enddo
                enddo
            sum=0
                do i=1,n-1
                sum=sum+(h/2)*(u(i)+u(i+1))
                enddo
        
            write(*,*)time,sum
            
                do i=1,n
                write(10,'(2f15.7)')h*(i-1-n/2),u(i)
                enddo
            write(10,*)'       '
            write(10,*)'       '
        
        enddo
        
        READ(*,*)
        end
