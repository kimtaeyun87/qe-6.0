!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_nsdiago (ik, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the superposition of atomic wavefunctions
  ! for k-point "ik" - output in "wfcatom"
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : omega, tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,      ONLY : xk, igk_k, ngk
  USE wvfct,      ONLY : npwx
  USE us,         ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE mp_bands,   ONLY : inter_bgrp_comm, set_bgrp_indices
  USE mp,         ONLY : mp_sum
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_U, Hubbard_l
  USE lsda_mod,   ONLY : nspin, isk, current_spin
  !
  implicit none
  !
  integer, intent(in) :: ik
  complex(DP), intent(out) :: wfcatom (npwx, npol, natomwfc)
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm, npw
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  complex(DP), allocatable :: sk (:), aux(:)
  complex(DP) :: kphase, lphase
  real(DP) :: arg, px, ux, vx, wx
  integer :: ig_start, ig_end

  call start_clock ('atomic_wfc_nsdiago')

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  do nt = 1, ntyp
     lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  enddo
  !
  nwfcm = MAXVAL ( upf(1:ntyp)%nwfc )
  npw = ngk(ik)
  !
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
     iig = igk_k (ig,ik)
     gk (1,ig) = xk(1, ik) + g(1,iig)
     gk (2,ig) = xk(2, ik) + g(2,iig)
     gk (3,ig) = xk(3, ik) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)

  ! from now to the end of the routine the ig loops are distributed across bgrp
  call set_bgrp_indices(npw,ig_start,ig_end)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = ig_start, ig_end
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do nt = 1, ntyp
     do nb = 1, upf(nt)%nwfc
        if ( upf(nt)%oc (nb) >= 0.d0) then
           do ig = ig_start, ig_end
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           enddo
        endif
     enddo
  enddo

  deallocate (qg, gk)
  allocate ( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  do na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     !
     !     sk is the structure factor
     !
     do ig = ig_start, ig_end
        iig = igk_k (ig,ik)
        sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                           eigts2 (mill (2,iig), na) * &
                           eigts3 (mill (3,iig), na)
     enddo
     !
     nt = ityp (na)
     do nb = 1, upf(nt)%nwfc
        if (upf(nt)%oc(nb) >= 0.d0) then
           l = upf(nt)%lchi(nb)
           lphase = (0.d0,1.d0)**l
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           if ( Hubbard_U(nt) /= 0.0d0 ) then
              call atomic_wfc_nsdiago___ ( )
           else
              call atomic_wfc___ ( )
           end if
           !
        end if
        !
     end do
     !
  end do

  if (n_starting_wfc /= natomwfc) call errore ('atomic_wfc_nsdiago', &
       'internal error: some wfcs were lost ', 1)

  deallocate(aux, sk, chiq, ylm)

  ! collect results across bgrp
  call mp_sum(wfcatom, inter_bgrp_comm)

  call stop_clock ('atomic_wfc_nsdiago')
  return

CONTAINS

   SUBROUTINE atomic_wfc___( )
   !
   ! ... LSDA or nonmagnetic case
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      if (n_starting_wfc > natomwfc) call errore &
         ('atomic_wfc___', 'internal error: too many wfcs', 1)
      !
      do ig = ig_start, ig_end
         wfcatom (ig, 1, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, nb, nt)
      ENDDO
      !
   END DO
   !
   END SUBROUTINE atomic_wfc___

   SUBROUTINE atomic_wfc_nsdiago___( )
   !
   ! ... LSDA or nonmagnetic case with diagonalization of ns
   !
   complex(DP)        :: f(Hubbard_lmax,Hubbard_lmax),vet(Hubbard_lmax,Hubbard_lmax)
   real(DP)           :: lambda(Hubbard_lmax)
   integer            :: ind, ind1, m1, m2, ldim
   !
   ldim = 2*l + 1
   !
   if ( ldim > Hubbard_lmax ) &
        call errore ('atomic_wfc_nsdiago___', 'Hubbard_lmax is too small', 1)
   !
   if ( nspin == 1 ) then
      current_spin = 1
   else if ( nspin == 2 ) then
      current_spin = isk ( ik )
   else
      call errore ('atomic_wfc_nsdiago___',' called in the wrong case ',1)
   end if
   !
   do m1 = 1, ldim
      do m2 = 1, ldim
         f (m1, m2) = rho%ns (m1, m2, current_spin, na)
      end do
   end do
   !
   call cdiagh(ldim, f, Hubbard_lmax, lambda, vet)
   !
   do m = 1, 2*l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      if (n_starting_wfc > natomwfc) call errore &
         ('atomic_wfc_nsdiago___', 'internal error: too many wfcs', 1)
      !
      aux(:) = (0.0d0, 0.0d0)
      do ind = 1, ldim
         ind1 = l**2 + ind
         aux(:) = aux(:) + vet(ind,m)*ylm(:,ind1)
      end do
      !
      do ig = ig_start, ig_end
         wfcatom (ig, 1, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, nb, nt)
      end do
      !
   end do
   !
   END SUBROUTINE atomic_wfc_nsdiago___

END SUBROUTINE atomic_wfc_nsdiago
