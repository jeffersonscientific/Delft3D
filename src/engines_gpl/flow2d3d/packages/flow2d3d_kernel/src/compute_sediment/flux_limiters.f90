module flux_limiters
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2024.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
    private
    
    public no_fluxlim
    public fluxlim_koren
    public fluxlim_mc
    public fluxlim_minmod
    
    contains
    
    function no_fluxlim(slp1, slp2) result (phi)
    use precision
    implicit none
    !
    real(fp), intent(in) :: slp1
    real(fp), intent(in) :: slp2
    real(fp)             :: phi
    !
    phi = 0.0_fp
    end function no_fluxlim
    
    function fluxlim_koren(slp1, slp2) result (phi)
    ! Koren ... based on Exnereq_newscheme_v03.pdf by Mart Borsboom
    use precision
    implicit none
    !
    real(fp), intent(in) :: slp1
    real(fp), intent(in) :: slp2
    real(fp)             :: phi
    !
    real(fp)             :: r     !< ratio of slp1 and slp2
    real(fp)             :: slp2_ !< local (limited) copy of slp2
    !
    if (abs(slp2) < 0.1_fp * abs(slp1)) then
        slp2_ = sign(0.1_fp * abs(slp1), slp2)
    else
        slp2_ = slp2
    endif
    r = slp1 / slp2_
    phi = max(0.0_fp, &
        &     min(2.0_fp * r, &
        &         min((2.0_fp + r) / 3.0_fp, &
        &             2.0_fp) &
        &        ) &
        &    )
    end function fluxlim_koren
    
    function fluxlim_mc(slp1, slp2) result (phi)
    ! monotized central (van Leer, 1977)
    use precision
    implicit none
    !
    real(fp), intent(in) :: slp1
    real(fp), intent(in) :: slp2
    real(fp)             :: phi
    !
    real(fp)             :: r     !< ratio of slp1 and slp2
    real(fp)             :: slp2_ !< local (limited) copy of slp2
    !
    if (abs(slp2) < 0.1_fp * abs(slp1)) then
        slp2_ = sign(0.1_fp * abs(slp1), slp2)
    else
        slp2_ = slp2
    endif
    r = slp1 / slp2_
    phi = max(0.0_fp, &
        &     min(2.0_fp * r, &
        &         min(0.5_fp * (1.0_fp + r), &
        &             2.0_fp) &
        &        ) &
        &    )
    end function fluxlim_mc
    
    function fluxlim_minmod(slp1, slp2) result (phi)
    ! minmod (Roe, 1986)
    use precision
    implicit none
    !
    real(fp), intent(in) :: slp1
    real(fp), intent(in) :: slp2
    real(fp)             :: phi
    !
    real(fp)             :: r     !< ratio of slp1 and slp2
    real(fp)             :: slp2_ !< local (limited) copy of slp2
    !
    if (abs(slp2) < 0.1_fp * abs(slp1)) then
        slp2_ = sign(0.1_fp * abs(slp1), slp2)
    else
        slp2_ = slp2
    endif
    r = slp1 / slp2_
    phi = max(0.0_fp, &
        &     min(r, &
        &         1.0_fp) &
        &    )
    end function fluxlim_minmod

end module flux_limiters
