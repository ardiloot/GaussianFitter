
subroutine gaussian2d(xs, ys, n, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi, res)
    ! Input variables
    integer, intent(in) :: n
    double precision, intent(in) :: xs(n), ys(n)
    double precision, intent(in) :: amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi
    double precision, intent(out) :: res(n)
        
    ! Local variables
    double precision :: x, y, xP, yP
    
    ! Code
    do i = 1,n
        x = xs(i) - x0
        y = ys(i) - y0
        xP = cosPhi * x - sinPhi * y
        yP = sinPhi * x + cosPhi * y
        res(i) = amp * exp(-((xP) ** 2.0) / (2.0 * sigmaX ** 2.0) - ((yP) ** 2.0) / (2.0 * sigmaY ** 2.0)) + offset
    enddo    
end
    
subroutine gaussian2d_sym(xs, ys, n, amp, x0, y0, sigma, offset, res)
    ! Input variables
    integer, intent(in) :: n
    double precision, intent(in) :: xs(n), ys(n)
    double precision, intent(in) :: amp, x0, y0, sigma, offset
    double precision, intent(out) :: res(n)
        
    ! Local variables
    double precision :: x, y
    
    ! Code
    do i = 1,n
        x = xs(i) - x0
        y = ys(i) - y0
        res(i) = amp * exp(-((y) ** 2.0 + (x) ** 2.0) / (2.0 * sigma ** 2.0)) + offset
    enddo    
end

subroutine fit_gaussian_sigma(xs, ys, values, m, initial, n, cosPhi, sinPhi, ftol, &
        xtol, gtol, maxfev, factor, solution, fvec, info, nfev, fjac, ipvt, qtf)
    
    ! Input variables
    integer, intent(in) :: m, n
    double precision, intent(in) :: xs(m), ys(m), values(m), initial(n)
    double precision, intent(in) :: cosPhi, sinPhi, ftol, xtol, gtol, factor
    integer, intent(in) :: maxfev
    double precision, intent(out) :: fvec(m), fjac(m, n), qtf(n), solution(n)
    integer, intent(out) :: info, nfev, ipvt(n)
    
    ! Local variables
    integer :: mode, nprint, ldfjac
    double precision :: epsfcn, diag(n), wa1(n), wa2(n), wa3(n), wa4(m)
    
    ! Copy initial to solution array
    do i = 1,n
        solution(i) = initial(i)
    enddo
    
    ! Init params
    epsfcn = 0.0d0
    mode = 1
    nprint = 0
    ldfjac = m
    
    ! Fit
    call lmdif(fit_func, m, n, solution, fvec, ftol, xtol, gtol, & 
        maxfev, epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, &
        ipvt, qtf, wa1, wa2, wa3, wa4)
    
    ! Contains
    contains
        subroutine fit_func(m, n, x, fvec, iflag)
            ! Input variables
            integer, intent(in) :: m, n
            integer, intent(inout) :: iflag
            double precision, intent(in) :: x(n)
            double precision, intent(inout) :: fvec(m)
            double precision, external :: enorm
            ! x = (amp, x0, y0, sigmaX, sigmaY, offset)
            
            ! Save gaussian values to fvec
            call gaussian2d(xs, ys, m, x(1), x(2), x(3), x(4), x(5), x(6), cosPhi, sinPhi, fvec)
            
            ! Subtract values
            do i = 1, m
                fvec(i) = fvec(i) - values(i)
            enddo
        end
end
    
subroutine fit_gaussian_sigma_sym(xs, ys, values, m, initial, n, ftol, &
        xtol, gtol, maxfev, factor, solution, fvec, info, nfev, fjac, ipvt, qtf)
    
    ! Input variables
    integer, intent(in) :: m, n
    double precision, intent(in) :: xs(m), ys(m), values(m), initial(n)
    double precision, intent(in) :: ftol, xtol, gtol, factor
    integer, intent(in) :: maxfev
    double precision, intent(out) :: fvec(m), fjac(m, n), qtf(n), solution(n)
    integer, intent(out) :: info, nfev, ipvt(n)
    
    ! Local variables
    integer :: mode, nprint, ldfjac
    double precision :: epsfcn, diag(n), wa1(n), wa2(n), wa3(n), wa4(m)
    
    ! Copy initial to solution array
    do i = 1,n
        solution(i) = initial(i)
    enddo
    
    ! Init params
    epsfcn = 0.0d0
    mode = 1
    nprint = 0
    ldfjac = m
    
    ! Fit
    call lmdif(fit_func, m, n, solution, fvec, ftol, xtol, gtol, & 
        maxfev, epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, &
        ipvt, qtf, wa1, wa2, wa3, wa4)
    
    ! Contains
    contains
        subroutine fit_func(m, n, x, fvec, iflag)
            ! Input variables
            integer, intent(in) :: m, n
            integer, intent(inout) :: iflag
            double precision, intent(in) :: x(n)
            double precision, intent(inout) :: fvec(m)
            double precision, external :: enorm
            ! x = (amp, x0, y0, sigma, offset)
            
            ! Save gaussian values to fvec
            call gaussian2d_sym(xs, ys, m, x(1), x(2), x(3), x(4), x(5), fvec)
            
            ! Subtract values
            do i = 1, m
                fvec(i) = fvec(i) - values(i)
            enddo
        end
end

subroutine fit_gaussian_z(xs, ys, values, m, initial, n, cosPhi, sinPhi, &
		tSigmaX, nSigmaX, cSigmaX, kSigmaX, extSigmaX, &
		tSigmaY, nSigmaY, cSigmaY, kSigmaY, extSigmaY, &
        ftol, xtol, gtol, maxfev, factor, solution, fvec, info, nfev, fjac, ipvt, qtf)
    
    ! Input variables
    integer, intent(in) :: m, n
    double precision, intent(in) :: xs(m), ys(m), values(m), initial(n)
    double precision, intent(in) :: cosPhi, sinPhi, ftol, xtol, gtol, factor
    integer, intent(in) :: maxfev
    double precision, intent(out) :: fvec(m), fjac(m, n), qtf(n), solution(n)
    integer, intent(out) :: info, nfev, ipvt(n)
    integer, intent(in) :: nSigmaX, kSigmaX, extSigmaX
    double precision, intent(in) :: tSigmaX(nSigmaX), cSigmaX(nSigmaX)
    integer, intent(in) :: nSigmaY, kSigmaY, extSigmaY
    double precision, intent(in) :: tSigmaY(nSigmaY), cSigmaY(nSigmaY)
    
    ! Local variables
    integer :: mode, nprint, ldfjac
    double precision :: epsfcn, diag(n), wa1(n), wa2(n), wa3(n), wa4(m)
    
    ! Copy initial to solution array
    do i = 1,n
        solution(i) = initial(i)
    enddo
    
    ! Init params
    epsfcn = 0.0d0
    mode = 1
    nprint = 0
    ldfjac = m
    
    ! Fit
    call lmdif(fit_func, m, n, solution, fvec, ftol, xtol, gtol, & 
        maxfev, epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, &
        ipvt, qtf, wa1, wa2, wa3, wa4)
    
    ! Contains
    contains
        subroutine fit_func(m, n, x, fvec, iflag)
            ! Input variables
            integer, intent(in) :: m, n
            integer, intent(inout) :: iflag
            double precision, intent(in) :: x(n)
            double precision, intent(inout) :: fvec(m)
            double precision, external :: enorm
            ! x = (amp, x0, y0, z0, offset)
            
            ! Local variabled
            double precision :: sigmaXTmp(1), sigmaYTmp(1), zTmp(1)
            integer :: ier
            
            ! Get sigmas by z
            zTmp(1) = x(4)
            call splev(tSigmaX, nSigmaX, cSigmaX, kSigmaX, zTmp, sigmaXTmp, 1, 0, ier)
            call splev(tSigmaY, nSigmaY, cSigmaY, kSigmaY, zTmp, sigmaYTmp, 1, 0, ier)
            
            ! Save gaussian values to fvec
            ! gaussian2d(xs, ys, n, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi, res)
            call gaussian2d(xs, ys, m, x(1), x(2), x(3), sigmaXTmp(1), sigmaYTmp(1), x(5), cosPhi, sinPhi, fvec)
            
            ! Subtract values
            do i = 1, m
                fvec(i) = fvec(i) - values(i)
            enddo
        end
    end

    
subroutine fit_gaussian_z_sym(xs, ys, values, m, initial, n, &
		tSigma, nSigma, cSigma, kSigma, extSigma, &
        ftol, xtol, gtol, maxfev, factor, solution, fvec, info, nfev, fjac, ipvt, qtf)
    
    ! Input variables
    integer, intent(in) :: m, n
    double precision, intent(in) :: xs(m), ys(m), values(m), initial(n)
    double precision, intent(in) :: ftol, xtol, gtol, factor
    integer, intent(in) :: maxfev
    double precision, intent(out) :: fvec(m), fjac(m, n), qtf(n), solution(n)
    integer, intent(out) :: info, nfev, ipvt(n)
    integer, intent(in) :: nSigma, kSigma, extSigma
    double precision, intent(in) :: tSigma(nSigma), cSigma(nSigma)

    
    ! Local variables
    integer :: mode, nprint, ldfjac
    double precision :: epsfcn, diag(n), wa1(n), wa2(n), wa3(n), wa4(m)
    
    ! Copy initial to solution array
    do i = 1,n
        solution(i) = initial(i)
    enddo
    
    ! Init params
    epsfcn = 0.0d0
    mode = 1
    nprint = 0
    ldfjac = m
    
    ! Fit
    call lmdif(fit_func, m, n, solution, fvec, ftol, xtol, gtol, & 
        maxfev, epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, &
        ipvt, qtf, wa1, wa2, wa3, wa4)
    
    ! Contains
    contains
        subroutine fit_func(m, n, x, fvec, iflag)
            ! Input variables
            integer, intent(in) :: m, n
            integer, intent(inout) :: iflag
            double precision, intent(in) :: x(n)
            double precision, intent(inout) :: fvec(m)
            double precision, external :: enorm
            ! x = (amp, x0, y0, z0, offset)
            
            ! Local variabled
            double precision :: sigmaTmp(1), zTmp(1)
            integer :: ier
            
            ! Get sigmas by z
            zTmp(1) = x(4)
            call splev(tSigma, nSigma, cSigma, kSigma, zTmp, sigmaTmp, 1, 0, ier)
            
            ! Save gaussian values to fvec
            call gaussian2d_sym(xs, ys, m, x(1), x(2), x(3), sigmaTmp(1), x(5),fvec)
            
            ! Subtract values
            do i = 1, m
                fvec(i) = fvec(i) - values(i)
            enddo
        end
end
