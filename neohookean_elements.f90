
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

!
! Element routines for the Mooney-Rivlin incompressible elastic solid
!

module neohookean_elements_m

  use tfem_elem_m

  implicit none


contains


! Internal element routine for the elastic model 

  subroutine neohookean_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use neohookean_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: funcnr, vfuncnr
    integer :: i, j, ip, i1, i2, i3, is, ie, N, M, k, q, z, l

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_elastic_up ( mesh, coefficients, elgrp )
      call set_elastic_model ( coefficients )

!     allocate arrays

      allocate ( wg(ninti), Jwg(ninti), fpg(ninti), fg(ninti,ndim) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( x0(nodalp,ndim), u(ndf*ndim) )
      allocate ( F0(ninti,ndim,ndim) )
      allocate ( F0inv(ninti,ndim,ndim), detF0(ninti) )
      allocate ( Fc(ninti,ndim,ndim) )
      allocate ( Fcinv(ninti,ndim,ndim), detFc(ninti) )
      allocate ( Jv(ninti) )
      allocate ( psi(ninti,ndfp) )
      allocate ( xg(ninti,ndim), xg0(ninti,ndim) )
      allocate ( tmp(ndf,ndim) )
      allocate ( dphi(ninti,ndf,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( B(ninti,ndim,ndim), Binv(ninti,ndim,ndim) )
      allocate ( tau(ninti,ndim,ndim), pg(ninti), p(ndfp) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( work(ninti,ndim,ndim,ndim), work1(ninti,ndim) )
      allocate ( work2(ninti,ndim), work3(ninti,ndim,ndim,ndim) )
      allocate ( work4(ninti), work6(ninti,ndf,ndim), work5(ndf,ndim)  )
      allocate ( Smat(ndf,ndf,ndim,ndim) )
      allocate ( Lu(ndfp,ndf), Lv(ndfp,ndf), Lw(ndfp,ndf) ) 
      allocate ( F_T(ninti,ndim,ndim), Finv_T(ninti,ndim,ndim),F_iso(ninti,ndim,ndim))
      allocate ( B_iso_tr(ninti), F_iso_T(ninti,ndim,ndim))
      allocate ( B_iso(ninti,ndim,ndim), B_iso_d(ninti,ndim,ndim))
      allocate ( I2(ninti,ndim,ndim), I4(ninti,ndim,ndim,ndim,ndim), &
      I4RT(ninti,ndim,ndim,ndim,ndim), II(ninti,ndim,ndim,ndim,ndim), &
      I4S(ninti,ndim,ndim,ndim,ndim),dphidx0(ninti,ndf,ndim))
      allocate ( M41(ninti,ndim,ndim,ndim,ndim), M42(ninti,ndim,ndim,ndim,ndim) ,&
      M4(ninti,ndim,ndim,ndim,ndim), B41(ninti,ndim,ndim,ndim,ndim), &
      B42(ninti,ndim,ndim,ndim,ndim), B43(ninti,ndim,ndim,ndim,ndim), &
      B44(ninti,ndim,ndim,ndim,ndim), B45(ninti,ndim,ndim,ndim,ndim), &
      B4(ninti,ndim,ndim,ndim,ndim), K1(ninti,ndim,ndim,ndim,ndim), &
      K2(ninti,ndim,ndim,ndim), K3(ninti,ndim,ndim))         

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if ( coefficients%i(41) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
      end if

    end if

!   coordinates in the reference configuration

    call get_coordinates ( mesh, elgrp, elem, x0 )

    call isoparametric_deformation ( x0, dphi, F0, F0inv, detF0 )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error in elastic_elem:', &
        ' axi-symmetrical elements not yet implemented'
      stop 
    end if

!   coordinates in the current configuration

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=(/physqdisp/), layer=layer )

!   NOTE: the following only works for isoparametric elements

    x = reshape( u, (/ndf,ndim/) ) + x0

    call isoparametric_deformation ( x, dphi, Fc, Fcinv, detFc )

    call isoparametric_coordinates ( x, phi, xg )
    call isoparametric_coordinates ( x0, phi, xg0 )

    Jwg = detF0 * wg

    call shape_derivative ( dphi, F0inv, dphidx0 )
    call shape_derivative ( dphi, Finv, dphidx )

    if ( any ( coefficients%i(41) == [1,2] ) ) then
!     global shape function for pressure
      call set_elastic_shape_function_global ( shapefuncp, x0, xg0, psi, &
        coefficients%i(41)==1 )
    end if
!! 
!   initialize the identity tensors used
!!
    do ip = 1, ninti
      do i = 1, ndim
        I2(ip,i,i) = 1
      enddo
    enddo
    
    do ip = 1, ninti
      do i = 1, ndim
        do j = 1, ndim
        I4(ip,i,j,j,i) = 1
        enddo
      enddo
    enddo
        
    do i = 1, ndim
      do j = 1, ndim
        I4RT(:,:,:,i,j) = I4(:,:,:,j,i)
      enddo
    enddo

    do i = 1, ndim
      do j = 1, ndim
        do k = 1, ndim
          do l = 1, ndim
            II(:,i,j,k,l) = I2(:,i,j)*I2(:,k,l)
          enddo
        enddo
      enddo
    enddo
     
    I4S = 0.5*(I4+I4RT)
    
!   current deformation and stresses 
!   use reference element as intermediate body

    Jv = detFc / detF0

    do ip = 1, ninti
      F(ip,:,:)         = matmul( Fc(ip,:,:), F0inv(ip,:,:) )
      F_T(ip,:,:)       = transpose(F(ip,:,:))
      Finv(ip,:,:)      = matmul( F0(ip,:,:), Fcinv(ip,:,:) )
      F_inv_T(ip,:,:)   = transpose(Finv(ip,:,:))
      F_iso(ip,:,:)     = Jv**(1/3)*transpose(F(ip,:,:))
      F_iso_T(ip,:,:)   = transpose(F_iso(ip,:,:))
      B(ip,:,:)         = matmul(F(ip,:,:),transpose(F(ip,:,:)))
      B_iso(ip,:,:)     = matmul(F_iso(ip,:,:),transpose(F_iso(ip,:,:)))
    end do
    
    B_iso_tr = 0
    do ip = 1, ninti
      do i = 1, ndim
        B_iso_tr(ip)    = B_iso_tr(ip) + B_iso(ip,ndim,ndim)
      enddo
    enddo
       
    do ip = 1, ninti
      B_iso_d(ip,:,:)   = B_iso - ((1/3)*B_iso_tr(ip)*I2)
    enddo  
     
    tau                 = kappa*(Jv-1)*I2+G*B_d_iso

    
!   current pressures

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, p, &
      physq=(/physqpress/), layer=layer )

    do ip = 1, ninti
      pg(ip) = dot_product( p, psi(ip,:) )
    end do
    print*, ninti
!   pointers in unknowns

    i1 = ndf
    i2 = 2*ndf
    i3 = 3*ndf

    is = ndim*ndf
    ie = ndim*ndf+ndfp

!   compute vector for Newton-Raphson iteration

    if ( vector ) then

!     - (nabla v)^T:tau

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tau(ip,:,:) )
      end do
      do j = 1, ndim
        work5(:,j) = - matmul ( Jwg, work6(:,:,j) )
      end do
      elemvec(1:is) = reshape ( work5, (/ ndf*ndim /) )

!     - ( J - 1 ) / J

      work4 = ( Jv - 1 ) / Jv * Jwg

      do N = 1, ndfp
        elemvec(is+N) = - sum ( psi(:,N) * work4 )
      end do

!     body force

      vfuncnr = coefficients%i(14)

      if ( vfuncnr > 0 ) then

        do ip = 1, ninti
          fg(ip,:) = coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) )
        end do

        do j = 1, ndim
          do N = 1, ndf
            tmp(N,j) = sum ( fg(:,j) * phi(:,N) * Jwg )
          end do
        end do

        elemvec(1:is) = elemvec(1:is) + reshape ( tmp, (/ndim*ndf/) )

      end if

!     fictitious body force for the incompressibility
  
      funcnr = coefficients%i(12)

      if ( funcnr > 0 ) then

        do ip = 1, ninti
          fpg(ip) = coefficients%func ( funcnr, xg(ip,:) )
        end do

        do i = 1, ndfp
          elemvec(is+i) = elemvec(is+i) + sum ( fpg * psi(:,i) * Jwg )
        end do

      end if
  
    end if

!   compute matrix for Newton-Raphson iteration

    if ( matrix ) then

      do M = 1, ndf

        work = 0   ! work is equal to A^M_ijk

!      -F^{-T}*I4RT 
        do i = 1, ndim
          do j = 1, ndim
            do k = 1, ndim
              do l = 1, ndim
                do q =1, ndim
                  M41(:,i,j,k,l) = M41(:,i,j,k,l)+(-F_inv_T(:,i,q)*I4RT(:,q,j,k,l))
                end do
              enddo
            enddo
          enddo
        enddo

!       -F^{-T}*I4RT*tau

        do i = 1, ndim
          do j = 1, ndim
            do k = 1, ndim
              do l = 1, ndim
                do q =1, ndim
                  M42(:,i,j,k,l) = M42(:,i,j,k,l)+(M41(:,i,j,k,q)*tau(:,q,l))
               end do
              enddo
            enddo
          enddo
        enddo

!        -F^{-T}*I4RT*tau*F^{-1}

        do i = 1, ndim
          do j = 1, ndim
            do k = 1, ndim
              do l = 1, ndim
                do q = 1, ndim
                  M4(:,i,j,k,l) = M4(:,i,j,k,l)+(M42(:,i,j,k,q)*F_inv(:,q,l))
                enddo
              enddo
            enddo
          enddo
        enddo

!       Now comes the model dependency

!       4I-((1/3)*II)

        B41 = I4-((1/3)*II)

!       2*I4S*F_iso_T
        
        do i = 1, ndim
          do j = 1, ndim
            do k = 1, ndim
              do l = 1, ndim
                do q =1, ndim
                  B42(:,i,j,k,l) = B42(:,i,j,k,l)+((2*I4S(:,i,j,k,q))*F_iso_T(:,q,l))
               end do
              enddo
            enddo
          enddo
        enddo

!       (-(1/3)*J^(-1/3)*F^T*F^-T)+J^1/3*I4

        do i = 1, ndim
          do j = 1, ndim
            do k = 1, ndim
              do l = 1, ndim
              B43(:,i,j,k,l) = ((-1/3)*(Jv**(-1/3)))*F_T(:,i,j)*F_inv_T(:,k,l)+&
                               (Jv**(-1/3))*I4(:,i,j,k,l)
              enddo
            enddo
          enddo
        enddo
        

!   F^{-1}*4I-((1/3)*II)

        do i = 1, ndim
          do j = 1, ndim
            do k = 1, ndim
              do l = 1, ndim
                do q =1, ndim
                  B44(:,i,j,k,l) = B44(:,i,j,k,l)+(F_inv(:,i,q)*B41(:,q,j,k,l))
                end do
              enddo
            enddo
          enddo
        enddo

!  F^{-1}*4I-((1/3)*II):(2*I4S*F_iso_T)

        do i = 1 , ndim
          do j = 1, ndim      
            do q = 1, ndim
              do z = 1, ndim
                do k = 1, ndim
                   do l = 1, ndim
                     B45(:,i,j,k,l) = B45(:,i,j,k,l)+B44(i,j,z,q)*B42(q,z,k,l)
                   enddo
                enddo
              enddo
            enddo
          enddo
        enddo

! F^{-1}*4I-((1/3)*II):(2*I4S*F_iso_T):(-(1/3)*J^(-1/3)*F^T*F^-T)+J^1/3*I4

        do i = 1 , ndim
          do j = 1, ndim      
            do q = 1, ndim
              do z = 1, ndim
                do k = 1, ndim
                   do l = 1, ndim
                     B4(:,i,j,k,l) = B4(:,i,j,k,l)+B45(i,j,z,q)*B43(q,z,k,l)
                   enddo
                enddo
              enddo
            enddo
          enddo
        enddo



!       (4B+M4):(nabla v)
       K1 = B4+M4

        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              do k = 1, ndim
                do q = 1, ndim
                  K2(:,i,j,k) = K2(:i,j,k)+ K1(:,i,j,k,q)*dphidx(:,M,q)
                enddo
              enddo
            enddo
          enddo
        enddo


!       (nabla v)^T:(4B+M4):(nabla v)

        do N = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              do q = 1, ndim
                  K3(:,i,j) = K3(:i,j)+ dphidx0(:,N,q)* K2(:,q,i,j)              
              enddo
            enddo
          enddo
        enddo
        
        do N = 1, ndf
          do M = 1 ,ndf
            do i = 1, ndim
              do j = 1, ndim
                Smat(N,M,i,j) = K(:,i,j) * Jwg
              enddo
            enddo
          enddo
        enddo

       
!     pressure part

      do i = 1, ndfp
        do j = 1, ndf
          Lu(i,j) = sum ( psi(:,i) * dphidx(:,j,1) * Jwg )
          Lv(i,j) = sum ( psi(:,i) * dphidx(:,j,2) * Jwg )
        end do
      end do

      if ( coorsys == 2 ) then
        do i = 1, ndfp
          do j = 1, ndf
            Lw(i,j) = sum ( psi(:,i) * dphidx(:,j,3) * Jwg )
          end do
        end do
      end if

!     fill element matrix

      if ( coorsys == 2 ) then

!       3D

        elemmat(    1:i1,    1:i1 ) = Smat(:,:,1,1)
        elemmat(    1:i1, i1+1:i2 ) = Smat(:,:,1,2)
        elemmat(    1:i1, i2+1:i3 ) = Smat(:,:,1,3)
        elemmat( i1+1:i2,    1:i1 ) = Smat(:,:,2,1)
        elemmat( i1+1:i2, i1+1:i2 ) = Smat(:,:,2,2)
        elemmat( i1+1:i2, i2+1:i3 ) = Smat(:,:,2,3)
        elemmat( i2+1:i3,    1:i1 ) = Smat(:,:,3,1)
        elemmat( i2+1:i3, i1+1:i2 ) = Smat(:,:,3,2)
        elemmat( i2+1:i3, i2+1:i3 ) = Smat(:,:,3,3)
        elemmat( is+1:ie,    1:i1 ) = Lu
        elemmat( is+1:ie, i1+1:i2 ) = Lv
        elemmat( is+1:ie, i2+1:i3 ) = Lw
        elemmat(    1:i1, is+1:ie ) = -transpose(Lu)
        elemmat( i1+1:i2, is+1:ie ) = -transpose(Lv)
        elemmat( i2+1:i3, is+1:ie ) = -transpose(Lw)
        elemmat( is+1:ie, is+1:ie)  = 0

      else

!       2D and axisymmetric

        elemmat(    1:i1,    1:i1 ) = Smat(:,:,1,1)
        elemmat(    1:i1, i1+1:i2 ) = Smat(:,:,1,2)
        elemmat( i1+1:i2,    1:i1 ) = Smat(:,:,2,1)
        elemmat( i1+1:i2, i1+1:i2 ) = Smat(:,:,2,2)
        elemmat( is+1:ie,    1:i1 ) = Lu
        elemmat( is+1:ie, i1+1:i2 ) = Lv
        elemmat(    1:i1, is+1:ie ) = -transpose(Lu)
        elemmat( i1+1:i2, is+1:ie ) = -transpose(Lv)
        elemmat( is+1:ie, is+1:ie)  = 0

      end if

    end if 

    if ( last ) then

!     last element in this group
   
      deallocate ( wg, Jwg, fpg, fg )
      deallocate ( xig, phi, x )
      deallocate ( x0, u )
      deallocate ( F0 )
      deallocate ( F0inv, detF0 )
      deallocate ( Fc )
      deallocate ( Fcinv, detFc )
      deallocate ( Jv )
      deallocate ( psi )
      deallocate ( xg, xg0 )
      deallocate ( tmp )
      deallocate ( dphi )
      deallocate ( F, Finv )
      deallocate ( B, Binv )
      deallocate ( tau, pg, p )
      deallocate ( dphidx )
      deallocate ( work, work1 )
      deallocate ( work2, work3 )
      deallocate ( work4, work6, work5 )
      deallocate ( Smat )
      deallocate ( Lu, Lv, Lw ) 
      deallocate ( F_inv_T, F_iso, F_iso_T, B_iso, B_iso_tr, B_iso_tr)
      deallocate ( I2, I4, I4RT, II, I4S, dphidx0)
      deallocate ( M41, M42, M4, B41, B42, B43, B44, B45, B4, K2, K1, K)

    end if

  end subroutine neohookean_elem


! compute pressures in all nodes

  subroutine elastic_pressure ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    if ( first ) then

!     first element in this group
   
!     set globals

      call set_globals_elastic_up ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( psi(nodalp,ndfp), x0(nodalp,ndim) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( p(ndfp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      if ( coefficients%i(41) == 0 ) then
        call set_shape_function ( shapefuncp, xrnod, psi )
      end if

    end if

    if ( any ( coefficients%i(41) == [1,2] ) ) then
!     global shape function for pressure
      call get_coordinates ( mesh, elgrp, elem, x0 )
      call set_elastic_shape_function_global ( shapefuncp, x0, x0, psi, &
        coefficients%i(41)==1 )
    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, p, &
      physq=(/physqpress/), layer=layer )

    elemvec = matmul ( psi, p )

    elemwts = 1

    if ( last ) then

!     last element in this group
   
      deallocate ( psi, x0, xrnod, p )

    end if

  end subroutine elastic_pressure


! compute elastic energy in all nodes

  subroutine elastic_derive ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: i, ip


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_elastic_up ( mesh, coefficients, elgrp )
      call set_elastic_model ( coefficients )

!     allocate arrays

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( x0(nodalp,ndim), u(ndf*ndim) )
      allocate ( F0(nodalp,ndim,ndim) )
      allocate ( F0inv(nodalp,ndim,ndim), detF0(nodalp) )
      allocate ( Fc(nodalp,ndim,ndim) )
      allocate ( Fcinv(nodalp,ndim,ndim), detFc(nodalp) )
      allocate ( dphi(nodalp,ndf,ndim) )
      allocate ( F(nodalp,ndim,ndim), Finv(nodalp,ndim,ndim) )
      allocate ( B(nodalp,ndim,ndim), Binv(nodalp,ndim,ndim) )
      allocate ( dphidx(nodalp,ndf,ndim), tr(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

!   coordinates in the reference configuration

    call get_coordinates ( mesh, elgrp, elem, x0 )

    call isoparametric_deformation ( x0, dphi, F0, F0inv, detF0 )

!   coordinates in the current configuration

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=(/physqdisp/), layer=layer )

!   NOTE: the following only works for isoparametric elements

    x = reshape( u, (/ndf,ndim/) ) + x0

    call isoparametric_deformation ( x, dphi, Fc, Fcinv, detFc )

    call shape_derivative ( dphi, Fcinv, dphidx )

!   current deformation and stresses 
!   use reference element as intermediate body

    do ip = 1, nodalp
      F(ip,:,:) = matmul( Fc(ip,:,:), F0inv(ip,:,:) )
      B(ip,:,:) = matmul( F(ip,:,:), transpose(F(ip,:,:)) )
    end do

    if ( model == 2 ) then
!     Mooney-Rivlin
      do ip = 1, nodalp
        Finv(ip,:,:) = matmul( F0(ip,:,:), Fcinv(ip,:,:) )
        Binv(ip,:,:) = matmul( transpose(Finv(ip,:,:)), Finv(ip,:,:) )
      end do
    end if

    select case ( coefficients%i(13) )
    case(1) ! energy

!     energy = C1 * ( I_B - 3 ) + C2 * ( II_B - 3 )
!     note that II_B=I_Binv for J=1

      tr = 0
      do i = 1, ndim
        tr = tr + B(:,i,i) - 1
      end do
    
      elemvec = C1*tr

      if ( model == 2 ) then
!       Mooney-Rivlin
        tr = 0
        do i = 1, ndim
          tr = tr + Binv(:,i,i) - 1
        end do
        elemvec = elemvec + C2*tr
      end if

    case(2) ! relative volume change

      elemvec = detFc / detF0

    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( xrnod, phi, x )
      deallocate ( x0, u )
      deallocate ( F0 )
      deallocate ( F0inv, detF0 )
      deallocate ( Fc )
      deallocate ( Fcinv, detFc )
      deallocate ( dphi )
      deallocate ( F, Finv )
      deallocate ( B, Binv )
      deallocate ( dphidx, tr )

    end if

  end subroutine elastic_derive


! set elastic model

  subroutine set_elastic_model ( coefficients )

    use elastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

    model = coefficients%i(3)

    C1 = coefficients%r(1)
    C2 = coefficients%r(2)

  end subroutine set_elastic_model 


! set global parameters (internal element)

  subroutine set_globals_elastic_up ( mesh, coefficients, elgrp )

    use elastic_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp


!   check size of coefficients 

    call check ( coefficients, 'set_globals_elastic_up', ncoefi=100, &
      ncoefr=50, indexarray=(/23,41/), minimum=(/0,0/), maximum=(/1,2/) )
   
    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    physqdisp = coefficients%i(6)
    physqpress = coefficients%i(7)
    layer = coefficients%i(38)
    if ( ndim == 2 ) then
      coorsys = coefficients%i(23)
    else
      coorsys = 2
    end if
    nsides = mesh%element(elgrp)%numsides
    nodalpb = mesh%element(elgrp)%sidnumnod
    globalshape = mesh%element(elgrp)%globalshape

!   set number of degrees of freedom displacement

    intpol = coefficients%i(1)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(4)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_elastic_up', ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/4(a/))') 'Error in set_globals_elastic_up:', &
    ' Number of degrees of freedom of the displacement shape function (ndf) ', &
    ' is different from the number of nodal points in the element (nodalp) ', &
    ' This possibility (ndf /= nodalp) is not available.'
      stop 
    end if

!   set number of degrees of freedom pressure

    intpolp = coefficients%i(2)

    shapefuncp%globalshape = globalshape
    shapefuncp%interpolation = intpolp
    shapefuncp%numbering = 'regular'
    shapefuncp%p = coefficients%i(5)
    shapefuncp%spec_eval = 'gauss'
!   needed for spectral quads/hexahedra only to call Pp_at_GLL routine:
    shapefuncp%intrule = shapefunc%p + 1

    call set_ndf ( shapefuncp, 'set_globals_elastic_up', ndf=ndfp )

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_elastic_up:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop 
    end if

    if ( coefficients%i(42) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(10) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(15) )
    else
      intrule = coefficients%i(10)
      intrule2 = coefficients%i(15)
    end if
    nsubint = get_coefficient ( coefficients, index=32, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_elastic_up:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape 
      stop 
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_elastic_up:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape, ' and inttype = ', inttype
      stop 
    end if

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%intrule2 = intrule2
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_elastic_up:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_elastic_up


! set global parameters (boundary element)

  subroutine set_globals_elastic_up_boun ( mesh, coefficients, curve, surface, &
    volume, ndimr, geometry )

    use elastic_globals_m
    use set_optional_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in), optional :: curve, surface, volume
    integer, intent(in), optional :: ndimr, geometry

    integer :: lcurve, lsurface, lvolume


    physqdisp = coefficients%i(6)
    physqpress = coefficients%i(7)
    layer = coefficients%i(38)

!   traditional interface (legacy)

    lcurve = set_optional ( variable=curve, default=0 )
    lsurface = set_optional ( variable=surface, default=0 )
    lvolume = set_optional ( variable=volume, default=0 )

!   interface with dimension of reference space of geometry

    if ( present(ndimr) .and. present(geometry) ) then
      select case (ndimr)
        case(1); lcurve = geometry
        case(2); lsurface = geometry
        case(3); lvolume = geometry
      end select
    end if

    if ( lcurve > 0 ) then
      ndim = mesh%curves(lcurve)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
      else
        coorsys = coefficients%i(23)
      endif
      globalshape = mesh%curves(lcurve)%element%globalshape
      nodalp = mesh%curves(lcurve)%elnumnod
    else if ( lsurface > 0 ) then
      ndim = mesh%surfaces(lsurface)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
      else
        coorsys = coefficients%i(23)
      endif
      globalshape = mesh%surfaces(lsurface)%element%globalshape
      nodalp = mesh%surfaces(lsurface)%elnumnod
    else if ( lvolume > 0 ) then
      ndim = mesh%volumes(lvolume)%ndim
      coorsys = 2 ! always 3D
      globalshape = mesh%volumes(lvolume)%element%globalshape
      nodalp = mesh%volumes(lvolume)%elnumnod
    end if

!   set number of degrees of freedom displacement

    intpol = coefficients%i(1)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(4)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_elastic_up_boun', ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_elastic_vp_boun:', &
        ' Number of degrees of freedom of the shape function (ndf) ', &
        ' is different from the number of nodal points in the element', &
        ' (nodalp) ', &
        ' This possibility (ndf /= nodalp) is not available.'
      stop
    end if

!   set number of degrees of freedom pressure

    intpolp = coefficients%i(2)

    shapefuncp%globalshape = globalshape
    shapefuncp%interpolation = intpolp
    shapefuncp%numbering = 'regular'
    shapefuncp%p = coefficients%i(5)
    shapefuncp%spec_eval = 'gauss'
!   needed for spectral quads only to call Pp_at_GLL routine:
    shapefuncp%intrule = shapefunc%p + 1

    call set_ndf ( shapefuncp, 'set_globals_elastic_up_boun', ndf=ndfp )

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_elastic_up_boun:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(42) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(11) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(16) )
    else
      intrule = coefficients%i(11)
      intrule2 = coefficients%i(16)
    end if
    nsubint = get_coefficient ( coefficients, index=33, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_elastic_up_boun:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape 
      stop 
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_elastic_up_boun:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape, ' and inttype = ', inttype
      stop 
    end if

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%intrule2 = intrule2
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_elastic_up_boun:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_elastic_up_boun


! set global parameters (boundary Langrangian multiplier)

  subroutine set_globals_elastic_l_boun ( coefficients )

    use elastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   Lagrange multiplier interpolation

!   set number of degrees of freedom

    intpol = coefficients%i(1)

!   set number of degrees of freedom

    call set_ndf ( shapefunc, name_of_routine='set_globals_elastic_l_boun', &
      ndfl=ndfl, intpoll=intpoll )

    shapefuncl = shapefunc
    shapefuncl%interpolation = intpoll

  end subroutine set_globals_elastic_l_boun


! set shapefunction based on scaled global coordinates

  subroutine set_elastic_shape_function_global ( shapefunc, x, xg, phi, &
    scaling )

    type(shapefunc_t), intent(in) :: shapefunc

!   coordinates of the nodes of the element
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x 

!   Global coordinates where phi must be computed:
!   xg(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xg 

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi 

!   perform scaling of the shape function
    logical, intent(in) :: scaling

    integer :: ndim, ninti, i
    real(dp), allocatable, dimension(:) :: xmin, xmax
    real(dp), allocatable, dimension(:,:) :: xg_scale

    if ( scaling ) then

      ndim = size(x,2)
      ninti = size(xg,1)

      allocate ( xmin(ndim), xmax(ndim), xg_scale(ninti,ndim) )

      xmin = maxval ( x, dim=1 )  ! maximum in all coordinate directions
      xmax = minval ( x, dim=1 )  ! minimum in all coordinate directions

!     scale global coordinates:
!       - coordinates relative to the midpoint of the element (xmax+xmin)/2
!       - scale with half the extend xmax - xmin 
      do i = 1, ndim
        xg_scale(:,i) = &
                  ( 2 * xg(:,i) - xmax(i) - xmin(i) ) / ( xmax(i) - xmin(i) ) 
      end do

      call set_shape_function_global ( shapefunc, xg_scale, phi ) 

      deallocate ( xmin, xmax, xg_scale )

    else

      call set_shape_function_global ( shapefunc, xg, phi ) 

    end if

  end subroutine set_elastic_shape_function_global

end module neohookean_elements_m

