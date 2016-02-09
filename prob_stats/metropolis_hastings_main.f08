program main
  use nrtype
  use distributions
  use metropolis_hastings_mod
  implicit none

  !open(1,file = 'sine.dat')
  !open(2,file = 'sine_cos.dat')
  !open(3,file = 'exp_decay_cos.dat')
  open(4,file = 'bivar_gauss.dat')
  !call circular_met_hast(dist = radial_sine_dist, radial = 2._dp,&
  !                       angular = [0._dp,tau], auto_corr = 5,&
  !                       nbr_samples = int(1e6), file_unit = 1)

 !call circular_met_hast(dist = radial_sine_cos_dist, radial = 2._dp,&
 !                       angular = [0._dp,tau], auto_corr = 5,&
 !                       nbr_samples = int(1e6), file_unit = 2)

 !call circular_met_hast(dist = radial_exp_decay_cos, radial = 2._dp,&
 !                      angular = [0._dp,tau], auto_corr = 5,&
 !                      nbr_samples = int(1e6), file_unit = 3)

 call n_dim_met_hast(dist = bivariate_gaussian, n_dim = 2, &
                     aux = [-.5_dp,-5._dp,7._dp,6._dp,0.5_dp], mutagen = [-2._dp,2._dp], &
                     auto_corr = 5, nbr_samples = int(1e6), file_unit = 4)

 ! Equivalence.
 ! x = x(1)
 ! y = x(2)
 ! rho = aux(1)
 ! mu_x = aux(2)
 ! mu_y = aux(3)
 ! sig_x = aux(4)
 ! sig_y = aux(5)

end program main
