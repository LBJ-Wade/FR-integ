!===========================================================================================
!A dvode driver, for solving perturbation evolution in WANG Jun & WANG Hao 2012.
!
!hwang.phy@gmail.com
!
! v1.1 2012-Dec-02:
!       1. fiducial values changed, especially, n changed from 0.1 to 0.0001;
!       2. relative errors are modified
! v1.0 2012-Jul-15
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! Independent variable is A (eta), the conformal time
! lower case 'a' is the scale factor
!Y(1) = a
!Y(2) = a'
!Y(3) = a''
!Y(4) = delta
!Y(5) = delta'
!
!The system consists of the following 5 rate
!equations:
!1, YDOT(1) = a' = Y(2)
!2, YDOT(2) = a'' = Y(3)
!3, YDOT(3) = a''' = ...
!4, YDOT(4) = delta' = Y(5)
!5, YDOT(5) = delta'' = ...
!
!All length/time quantities are in unit of Gyr (or Gyr^-1);
!All density quantities are in unit of 8 \pi G * ['normal' density],
!but density appears in formulas only in form of dimensionless \Omega,
!making its actual unit irrelevant;
!
!The initial conditions are *not* trivial, neither (Y)s and A, nor
!(Y')s; 
!Initial value of Y(5) (delta') is set to a' at the same A (eta);
!this is argued to be due to similar behavior of $delta(a)$ in 
!f(R) theory and in LCDM, which seems to be the case in the figures
!shown in the paper. 
!
!Problems to be solved:
!1. A clever way to pick the wanted k (used to plot the growth history)
!===========================================================================================        
    MODULE wjnew 
    
        integer, parameter          :: NEQ = 5
        integer, parameter          :: WP = KIND(1.0D0)
        real(kind=wp), parameter :: kmsmpc2gyr = 0.0010227300227 !km s Mpc^-1 in Gyr^-1
        real(kind=wp), parameter :: mpc2gyr = 299792.458 * kmsmpc2gyr !Mpc^-1 in Gyr^-1
        real(kind=wp), parameter :: H0 = 75*kmsmpc2gyr
        real(kind=wp), parameter :: Omega = 0.25
        real(kind=wp), parameter :: q0 = 1.5*Omega - 1

        integer, parameter :: NPARA = 3  !NPARA paras to investigate: n, alpha, lambda
        character(10), dimension(NPARA) :: namepara = ['n_', 'a_', 'l_']
        integer, parameter :: enumn = 1
        integer, parameter :: enuma = 2
        integer, parameter :: enuml = 3

        real(kind=wp), dimension(NPARA) :: curpara
        real(kind=wp) :: k

        integer, parameter :: NTEST = 4  !NTEST values for each para for illustration

    CONTAINS

      SUBROUTINE ptnew(NEQ,A,Y,YDOT)
        implicit none
        integer, intent(in) :: NEQ
        real(kind=wp), intent(in), dimension(NEQ) :: Y
        real(kind=wp), intent(out), dimension(NEQ) :: YDOT
        real(kind=wp), intent(in) :: A
        real(kind=wp) :: n, alpha, lambda

        n = curpara(enumn)
        alpha = curpara(enuma)
        lambda = curpara(enuml)

        YDOT(1) = Y(2)

        YDOT(2) = Y(3)

        YDOT(3) = Y(3)*(6**n*alpha*(Y(1)**9)*(Y(3)/Y(1)**3)**n*&
        (n*(-4+3*n)*Y(2)**2+(-1+n)*Y(1)*Y(3))+6*Y(3)*(-Y(1)**6*&
        Y(2)**2+H0**2*(Y(1)**7-360*lambda*Y(2)**2*Y(3)-&
        36*lambda*Y(1)*Y(3)**2)*Omega))/(Y(1)*Y(2)*(6**n*(-1+n)*n&
        *alpha*Y(1)**9*(Y(3)/Y(1)**3)**n-432*lambda*H0**2*Y(3)**2*&
        Omega))

        YDOT(4) = Y(5)

        YDOT(5) = Y(1)**(-1)*(-Y(2)*Y(5)+6*H0**2*Y(3)*Y(4)*Omega*(&
        2*(Y(1)**6+36*lambda*Y(3)**2)/(6**n*n*alpha*Y(1)**9*&
        (Y(3)/Y(1)**3)**n+6*Y(3)*(Y(1)**6-72*lambda*H0**2*Y(3)*&
        Omega)) - Y(3)*(Y(1)**6+72*k**2*lambda*Y(1)*Y(3)+&
        36*lambda*Y(3)**2)/(6**n*n*alpha*Y(1)**9*(Y(3)/Y(1)**3&
        )**n*(k**2*(-1+n)*Y(1)+2*Y(3))+12*Y(3)**2*(Y(1)**6&
        -36*lambda*H0**2*(k**2*Y(1)+2*Y(3))*Omega))))

        RETURN
      END SUBROUTINE ptnew

END MODULE wjnew


PROGRAM propagate_wjnew

        use DVODE_F90_M
        use wjnew

        implicit none
        real(kind=wp) :: RTOL, A, AOUT, RSTATS
        integer ITASK, ISTATE, ISTATS, IOUT, IERROR, I
        dimension RSTATS(22), ISTATS(31)
        type (VODE_OPTS) :: OPTIONS

        integer ia, nia, np, nt
        real(kind=wp), dimension(NEQ) :: Y
        real(kind=wp) A_ini, A_fin
        real(kind=wp), dimension(NEQ) :: ATOL 
        
        !------------------------------------------------
        !The 'reference curves' in the paper, LCDM and non-coupling, are produced by letting 
        !(n = 0.0001, lambda = 0.0001) and (n = 0.1, lambda  = 0.0001), respectively
        !------------------------------------------------
        real(kind=wp), dimension(NPARA) :: para_fid, para_fid_norm = [0.00001, -4.3, -0.000001] !n, alpha, lambda
        real(kind=wp), dimension(NPARA) :: uconv_vec = [1._wp, H0**2, 1e-6_wp]
        real(kind=wp), dimension(NTEST, NPARA) :: para_mat_norm, para_mat !para matrix in sensible units
        real(kind=wp), dimension(NPARA, NPARA) :: uconv_mat = 0 !unit-converting matrix
        character(50) :: paravalue !containing name of para being investigated now

        !set ks, in base-10 logarithm value
        real(kind=wp) :: logk, dlogk = 0.1_wp, midlogk = -1.4_wp, maxlogk = 0._wp
        
        real(kind=wp) :: k_in_hmpc
        integer :: idxlogk, NLOGK
        NLOGK = 2 * (maxlogk - midlogk)/dlogk

        para_fid = para_fid_norm * uconv_vec !contain fiducial para values
        !fortran arrange array elements in a column-wise way
        para_mat_norm = reshape([ 0.01_wp,  0.02_wp, 0.10_wp,  0.20_wp, &! n
                                 -0.43_wp, -1.00_wp, -4.3_wp, -9.0_wp, &! alpha, in [H0^2], fid 4.3 (-2Lambda = -6*OL*H0^2)
                                 -0.01_wp, -0.02_wp, -0.2_wp, -0.4_wp],&! lambda, in [e^-6], [lambda] = [R^-2] = [t^4]
                                 shape(para_mat_norm))! mainly used in output file name

        forall(i=1 : NPARA) uconv_mat(i, i) = uconv_vec(i) !construct a diagonal matrix

        para_mat = matmul(para_mat_norm, uconv_mat) !parameter matrix used in calculation

        A_ini = 2.07_wp 
        A_fin = 46._wp

        !-------------------------------------------------------------
        nia = 1000
        do np = 1, NPARA
                curpara = para_fid
                do nt = 1, NTEST
                        curpara(np) = para_mat(nt, np)
!print *, np, nt
!print *, curpara/uconv_vec
                        !get output file name
                        write(paravalue, '(f15.3)') para_mat_norm(nt, np)

                        open(unit=6, file=trim('dat/Hist_')//trim(namepara(np))//trim(adjustl(paravalue))//'.dat')
                        WRITE (6,'(7a15)') 'k', 'eta', 'a', "a'", "a''", 'delta', "delta'"

                        open(unit=7, file=trim('dat/Tran_')//trim(namepara(np))//trim(adjustl(paravalue))//'.dat')
                        WRITE (7,'(6a15)') 'k', 'a', "a'", "a''", 'delta', "delta'"
                        
                        logk = 2 * midlogk - maxlogk !each (np, nt) reset logk to minimum value
                        do idxlogk = 1, NLOGK+1
                                logk = logk + dlogk
                                k_in_hmpc = 10._wp**(logk)
                                k = k_in_hmpc * (H0/100./kmsmpc2gyr) * (mpc2gyr)
                                !set the initial value for integration
                                IERROR = 0
                                ITASK = 1
                                ISTATE = 1
                                A = A_ini 
                                !-------------------------------------------------------------
                                !In determining these initial values, the evolution of Y(1:3),
                                !i.e., a, a', a'', are solved by integrating from eta = 40 back to 
                                !eta = 0.001; then insert the initial values of Y(5) according to
                                !delta' = a'; 
                                !
                                !-------------------------------------------------------------
                                Y(1) = 0.01_wp !-> Y(1) \equiv a, choose z_ini=100
                                Y(2) = 0.007816_wp !-> Y(2) \equiv a'
                                Y(3) = 0.00000981_wp !-> Y(3) = a''

                                !ODEs of delta are in fact decoupled from those above, so 
                                !the initial value of Y(4) should not matter much except for 
                                !an offset
                                !here Y(4)_ini is chosen to make delta_m today approximately
                                !0.9, i.e., consistent with observations on the 10Mpc-scale
                                Y(4) = 0._wp !->Y(4) = delta_m + c 
                                Y(5) = Y(2) !->Y(5) = delta_m' = a'

                                RTOL = 1.D-5
                        !      RTOL(2) = 1.e-5
                        !      RTOL(3) = 1.e-5
                        !      RTOL(4) = 1.e-3
                        !      RTOL(5) = 1.e-3
                        !        ATOL = 0
                                ATOL(1) = 1.D-7 !a, abserr matches relerr
                                ATOL(2) = 1.D-7 !a'
                                ATOL(3) = 1.D-14 !a'', a''' is very small
                                ATOL(4) = 1.D-7 !delta
                                ATOL(5) = 1.D-14 !delta', delta'' is very small
                                OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, RELERR=RTOL)

                                do ia = 1, nia
                                        AOUT = A_ini + (A_fin-A_ini) * (ia*1.0_wp/nia)
                                        CALL DVODE_F90(ptnew,NEQ,Y,A,AOUT,ITASK,ISTATE,OPTIONS)!,J_FCN=JEX)
                                                
                                        if (abs(logk - midlogk) .lt. dlogk/2.)  WRITE (6,'(5f12.7)') &
                                                k_in_hmpc, AOUT, Y(1), Y(4), Y(5) 
        !                                if (abs(logk - midlogk) .lt. dlogk/2.)  WRITE (6,'(7f15.10)') &
        !                                        k_in_hmpc, AOUT, Y(1), Y(2), Y(3), Y(4), Y(5) 
                                        CALL GET_STATS(RSTATS,ISTATS)
                                enddo
                                WRITE(7, '(3f12.7)') k_in_hmpc, Y(1), Y(4)
        !                        WRITE(7, '(6f15.10)') k_in_hmpc, Y(1), Y(2), Y(3), Y(4), Y(5)
                        enddo
                        
                        close(6)
                        close(7)
                enddo
        enddo

90000 FORMAT (/'  No. steps =',I4,'   No. f-s =',I4,'  No. J-s =',I4, &
        '   No. LU-s =',I4/'  No. nonlinear iterations =', &
        I4/'  No. nonlinear convergence failures =', &
        I4/'  No. error test failures =',I4/)
90001 FORMAT (/' An error occurred.')
90002 FORMAT (/' No errors occurred.')
90003 FORMAT (' At t =',D12.4,'   y =',3D14.6)
90004 FORMAT (///' Error halt: ISTATE =',I3)
        STOP

    END PROGRAM propagate_wjnew
