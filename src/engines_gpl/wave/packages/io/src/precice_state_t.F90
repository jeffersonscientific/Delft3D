module m_precice_state_t
   use, intrinsic :: iso_c_binding, only: c_int, c_char
   implicit none(type, external)
   private
   type, public :: precice_state_t
#if defined(HAS_PRECICE_FM_WAVE_COUPLING)
      integer(kind=c_int), dimension(:), allocatable :: vertex_ids
      character(kind=c_char, len=4) :: component_name = "wave"
      character(kind=c_char, len=10) :: swan_mesh_name = "wave_nodes"
      character(kind=c_char, len=10) :: bed_levels_name = "bed_levels"
      character(kind=c_char, len=2) :: hs_name = "hs"
      character(kind=c_char, len=3) :: dir_name = "dir"
      character(kind=c_char, len=4) :: dirc_name = "dirc"
      character(kind=c_char, len=4) :: dirs_name = "dirs"
      character(kind=c_char, len=6) :: period_name = "period"
      character(kind=c_char, len=5) :: depth_name = "depth"
      character(kind=c_char, len=2) :: fx_name = "fx"
      character(kind=c_char, len=2) :: fy_name = "fy"
      character(kind=c_char, len=7) :: wsbodyu_name = "wsbodyu"
      character(kind=c_char, len=7) :: wsbodyv_name = "wsbodyv"
      character(kind=c_char, len=2) :: mx_name = "mx"
      character(kind=c_char, len=2) :: my_name = "my"
      character(kind=c_char, len=7) :: dissip1_name = "dissip1"
      character(kind=c_char, len=7) :: dissip2_name = "dissip2"
      character(kind=c_char, len=7) :: dissip3_name = "dissip3"
      character(kind=c_char, len=7) :: dissip4_name = "dissip4"
      character(kind=c_char, len=4) :: ubot_name = "ubot"
      character(kind=c_char, len=5) :: steep_name = "steep"
      character(kind=c_char, len=4) :: wlen_name = "wlen"
      character(kind=c_char, len=1) :: u_name = "u"
      character(kind=c_char, len=1) :: v_name = "v"
      character(kind=c_char, len=4) :: dspr_name = "dspr"
      character(kind=c_char, len=5) :: rleak_name = "rleak"
      character(kind=c_char, len=2) :: qb_name = "qb"
      character(kind=c_char, len=1) :: x_name = "x"
      character(kind=c_char, len=1) :: y_name = "y"
      character(kind=c_char, len=3) :: rtp_name = "rtp"
      character(kind=c_char, len=4) :: hrms_name = "hrms"
      character(kind=c_char, len=2) :: tp_name = "tp"
      character(kind=c_char, len=4) :: pdir_name = "pdir"
      character(kind=c_char, len=5) :: windu_name = "windu"
      character(kind=c_char, len=5) :: windv_name = "windv"
      character(kind=c_char, len=3) :: tps_name = "tps"
      character(kind=c_char, len=4) :: tm02_name = "tm02"
      character(kind=c_char, len=5) :: tmm10_name = "tmm10"
      character(kind=c_char, len=6) :: dhsign_name = "dhsign"
      character(kind=c_char, len=6) :: drtm01_name = "drtm01"
      character(kind=c_char, len=5) :: setup_name = "setup"
#endif
   end type precice_state_t
end module m_precice_state_t
