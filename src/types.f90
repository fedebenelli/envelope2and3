module constants
   implicit none

   integer, parameter :: wp = kind(0.0q1)
end module constants

module dtypes
   use constants

   implicit none

   private
   public :: envelope
   public :: cross
   public :: kfcross

   type :: envelope
      real(8), allocatable :: z(:)
      real(8), allocatable :: t(:)
      real(8), allocatable :: p(:)
      real(8), allocatable :: logk(:, :)
      type(critical_point), allocatable :: critical_points
   end type envelope

   type :: critical_point
      real(8) :: t
      real(8) :: p
   end type critical_point

   type :: cross
      real(wp) :: x
      real(wp) :: y
      integer :: i
      integer :: j
   end type cross

contains

   function kfcross(i, t_values, logk, target_t)
      !! Estimate the Kvalues of an envelope by interpolation from a near point.
      integer, intent(in) :: i
      real(8), allocatable, intent(in) :: t_values(:) !! Envelope's temperatures
      real(8), allocatable, intent(in) :: logk(:, :) !! Envelope's kvalues
      real(8), intent(in) :: target_t !! Target temperature where to interpolate
      real(8), allocatable :: kfcross(:) !! Kvalues at desired point

      kfcross = (logk(i, :) - logk(i - 1, :)) &
                /&
                (t_values(i) - t_values(i - 1)) &
                * (target_t - t_values(i - 1)) &
                + logk(i - 1, :)
   end function

end module dtypes
