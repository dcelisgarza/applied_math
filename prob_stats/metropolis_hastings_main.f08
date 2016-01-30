program main
  use numbers
  use distributions
  use metropolis_hastings_mod
  implicit none

  !open(1,file = 'sine.dat')
  !open(2,file = 'sine_cos.dat')
  open(3,file = 'exp_decay_cos.dat')
  call circular_met_hast(dist = radial_sine_dist, radial = 2._dp,&
                         angular = [0._dp,tau], auto_corr = 5,&
                         nbr_samples = int(1e6), file_unit=1)

 call circular_met_hast(dist = radial_sine_cos_dist, radial = 2._dp,&
                        angular = [0._dp,tau], auto_corr = 5,&
                        nbr_samples = int(1e6), file_unit=2)

 call circular_met_hast(dist = exp_decay_cos, radial = 2._dp,&
                       angular = [0._dp,tau], auto_corr = 5,&
                       nbr_samples = int(1e6), file_unit=3)

end program main
