DOCUMENTATION FOR FILE RANDGEN.F
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

randgen.f is a collection of FORTRAN functions and subroutines
to produce pseudo-random numbers from a variety of distributions. There
is some doubt as to the adequacy of the routines provided in the NAG
library; hence the need for something like this. The (alleged) superiority
of the routines provided here is due to the use of a super-wizz-o
algorithm, due to Marsaglia & Zaman (1991), to produce random numbers
uniformly distributed on the interval (0,1). This algorithm has a period of
2^(1376) iterations (cf. 2^(57) for the NAG generator) and is suitable for
use when very long runs of data are required.

All of the functions and subroutines provided here have names
beginning with ZBQL. While this wasn't an entirely random choice (does
anybody remember the French & Saunders `Countdown' sketch?) it was felt that
it was unlikely to clash with anything already in existence. In addition,
there is one COMMON block, called ZBQL0001, and one BLOCK DATA unit called
ZBQLBD01 - so avoid these in programs if you possibly can ;-)

As far as possible, the specification for these subroutines is
exactly the same as that for the corresponding NAG routines, so converting
your existing programs should be a doddle. Nonetheless, a full specification
may be found below. Bug reports, abuse, fan mail etc. should be directed at
Richard Chandler (email: richard@stats.ucl.ac.uk).

The authors of these routines are as follows:

REC	Richard Chandler (email as above)
PJN	Paul Northrop (email: paul@stats.ucl.ac.uk)

Please feel free to use and adapt these routines for your own use. We
do, however, request that if the routines are used in work which is later
presented in public, or published, the source of the code (ie. us!) is
acknowledged. If for no other reason, this will allow other people to test
the routines to see if there's anything horribly wrong with them that will
invalidate your results!


REFERENCES:

Johnk, M.D. (1964) Erzegeugung von Betarerteilten und Gammaverteilten
Zuffallszahlen. Metrika, 8, 5-15.
Fishman, G.S. (1978) Principles of Discrete Event Simulation. Wiley,
New York.
Marsaglia, G. & Zaman, A. (1991) A new class of random number generators.
Annals of Applied Probability, vol. 1, pp.462-480.
Press, W.H., Teukolsky, S.A., Vetterling, W.T. and Flannery, B.P. (1992)
Numerical Recipes in FORTRAN - the Art of Scientific Computing
(second edition). CUP.
Ripley, B.D. (1987) Stochastic Simulation. Wiley, New York.
Wakefield, J.C., Smith, A.F.M. & Gelfand, A.E. (1991) Efficient computation
of random variates via the ratio-of-uniforms method. Statistics
and Computing, 1, 129-133.


  **********************************


    THE ROUTINES
    ^^^^^^^^^^^^

1)      SUBROUTINE                              ZBQLINI(SEED)

Author: REC

Replaces NAG routines G05CBF and G05CCF.

Parameters:     SEED            (integer, input)

Purpose: Sets the seed for the random number generator to a
repeatable or non-repeatable initial state.

Description: The Marsaglia-Zaman algorithm requires a seed array
rather than a single seed. If SEED is equal to zero, the first element in
this array is assigned a value obtained from the system clock. Otherwise,
the first element is taken to be the value of SEED, modulo 424967291 (just
in case there are any wise guys out there . . .). Thus setting SEED = 0
allows a non-repeatable initial state. The remaining elements in the seed
array are set using the usual congruential-type thing that most random
number generators use. The default values of the seed array are contained
within the BLOCK DATA unit ZBQLBD01.

Remarks: a temporary file called zbql1234.tmp is opened when the
seed is set from the system clock. This should be erased on exit. The
default file number is 801; however, if this is already open, the routine
searches till it finds an available file handle. NOTE: this only works on
systems for which the compiler recognizes the CALL SYSTEM command (NOT
part of the FORTRAN77 standard, but available on most Unix machines).
On compilers without this facility, some lines form this routine should
be commented out. These lines are clearly indicated in the source code.


2)      DOUBLE PRECISION FUNCTION               ZBQLU01(DUMMY)

Author: REC

Replaces NAG routine G05CAF.

Parameters:     DUMMY           (double precision, input)

Purpose: returns a pseudo-random number taken from a uniform
distribution between 0 & 1.

Description: uses the Marsaglia-Zaman algorithm. The argument DUMMY
is not used, but is required by FORTRAN syntax.

Remarks: None.


3)      DOUBLE PRECISION FUNCTION               ZBQLUAB(A,B)

Author: REC

Replaces NAG function G05DAF.

Parameters:     A,B             (both double precision, input)

Purpose: returns a pseudo-random number taken from a uniform
distribution between A & B.

Description: just shifts & scales the output from ZBQLU01. If A is
greater than B, a uniform random number between B and A is returned. If A
and B are equal, a warning message is generated and the function returns A.

Remarks: None.

Other routines called: ZBQLU01


4)      DOUBLE PRECISION FUNCTION               ZBQLEXP(MU)

Author: REC

Replaces NAG function G05DBF.

Parameters: MU          (double precision, input)

Purpose: returns a pseudo-random number taken from an exponential
distribution with mean MU.

Description: uses the fact that if U has a uniform (0,1)
distribution, -log(U)*MU has the required distribution. LAMBDA must be
positive - otherwise an error message is generated and the function returns
zero.

Remarks: none.

Other routines called: ZBQLU01


5)      DOUBLE PRECISION FUNCTION               ZBQLNOR(MU,SIGMA)

Author: REC

Replaces NAG function G05DDF.

Parameters: MU,SIGMA            (both double precision, input)

Purpose: returns a pseudo-random number taken from a normal
distribution with mean MU and standard deviation |SIGMA|.

Description: uses the Box-Muller algorithm (see Ripley(1987) in the
reference list above for a description), which is exact.

Remarks: if SIGMA is negative, its absolute value is used.

Other routines called: ZBQLU01, ZBQLEXP


6)      INTEGER FUNCTION                        ZBQLBIN(N,P)

Author: REC

No single NAG function does this . . .

Parameters: N                   (integer, input)
    P                   (double precision, input)

Purpose: returns a pseudo-random number taken from a binomial
distribution with parameters N (number of trials) and P (probability of
success at each trial).

Description: If NP and N(1-P) are both large, we guess the
value of ZBQLBIN: suppose we guess it to be K, then we can simulate the
Kth out of N order statistics for U(0,1) variates, which is Beta(K,N-K).
If this is less than P then there are at least K successes out of
N trials; otherwise there are at most K-1. We continue in this way till
our order statistic counting can be done by directly generating a binomial
with a small mean. See Fishman (1978) for the general idea.

When NP is small, geometric random variates are generated until their
sum exceeds N. If P lies outside the interval [0,1], or N<1, an
error message is generated and the function returns zero.

Remarks: none.

Other routines called: ZBQLGEO,ZBQLBET1


7)      INTEGER FUNCTION                        ZBQLGEO(P)

Author: REC

No corresponding NAG function.

Parameters: P                   (double precision, input)

Purpose: returns a pseudo-random number taken from a geometric
distribution with parameter P (ie. mean 1/P).

Description: For large P (greater than 0.9), counts Bernoulli
trials with success probability P until first success is obtained.
Otherwise, repeated Bernoulli trials are expensive and it's worth the effort
of a couple of logs, so use the inverse cumulative method: U is generated
from U(0,1), then the required quantity is 1 + INT( log(U)/log(1-P) ). If
P lies outside the interval [0,1], an error message is generated and the
function returns zero.

Remarks: none.

Other routines called: ZBQLU01


8)      INTEGER FUNCTION                        ZBQLPOI(MU)

Author: REC

No single NAG function does this . . .

Parameters: MU                  (double precision, input)

Purpose: returns a pseudo-random number taken from a Poisson
distribution with mean MU.

Description: When MU is less than 10, generates uniform random
variables until their product falls below exp(-MU), and counts the number
of trials required. When MU>10 this is expensive and we `guess' the value
in a similar way to the algorithm for binomial generation - if we guess the
value is K, then simulate time to Kth event in a Poisson Process using
a Gamma variate Y. If this is greater than 1 then the actual number of
events in (0,1] is binomial (K-1,1/Y); otherwise, we know that there are
at least K events and the remainder is a Poisson variate with mean
MU*(1-Y). For huge MU, even this could be slow so we use a rejection method
as described in Press et al. I've found that this can be slightly inaccurate,
probably due to things blowing up in the calculations - this is why I only
use it when absolutely necessary. MU must be positive, otherwise an error
message is generated and the function returns zero.

Remarks: none.

Other routines called: ZBQLGAM, ZBQLBIN, ZBQLU01, ZBQLLG


9)      DOUBLE PRECISION FUNCTION               ZBQLGAM(G,H)

Author: PJN

Replaces NAG function G05DGF.

Parameters: G,H            (both double precision, input)

Purpose: returns a pseudo-random number taken from a gamma
distribution with mean G/H and variance G/(H^2).  N.B. The functional part
of the Gamma density used in this case is x^(G-1) exp(-Hx).

Description: for 1 < G < 10 uses the ratio-of-uniforms method (see Ripley(1987),
page 65, in the reference list above for a description). This is a rejection
method and so the probability of acceptance has been increased using the
combined strategies of relocation of the density by mode-1 (i.e. G-2) and
the use of a generalised ratio-of-uniforms method with r=1/2 (see Wakefield,
Gelfand and Smith(1991) in the reference list above).  This function
calls ZBQLU01.  For G < 1 uses Ahrens and Dieter (1974) (see  Ripley (1987),
page 88).  For G > 10 uses Cheng and Feast (1979) (see Ripley (1987),
page 90).

Remarks: the only restrictions on G and H are that they are both
positive.  Simulation from a Chi-squared(v) distribution can be achieved
by noting that Chi-squared(v) = Gamma(v/2,1/2).

Other routines called: ZBQLU01


10)	DOUBLE PRECISION FUNCTION		ZBQLBET1(NU1,NU2)

Author: REC

    Replaces NAG function G05DLF

Parameters: NU1, NU2		(both double precision, input)

Purpose: returns a pseudo-random number taken from a beta
distribution of the first kind, with degrees of freedom NU1 and NU2.

Description: this uses the fact that, if X1 and X2 are independent
random variables, gamma distributed with parameters (NU1,1) and (NU2,1)
respectively, then Y = X1/(X1+X2) has the required distribution. X1 and
X2 are obtained from ZBQLGAM. The exception is when both NU1 and NU2 are
close to zero, when this method has a tendency to return 0/(0+0), even
with double precision arithmetic. In this case a computationally stable
version of Johnks' method (Johnk 1964) is used: we generate U1 and U2 from a
Uniform (0,1) distribution and then set

Y = 1/(1 + ( (U2**(1/NU2)) * (U1**(-1/NU1)) ) )

if (U1**(1/NU1)) + (U2**(1/NU2)) <= 1. The product can be evaluated stably
using logs and then transforming back.

Remarks: NU1 and NU2 must both be strictly positive.

Other routines called: ZBQLGAM


11)	DOUBLE PRECISION FUNCTION		ZBQLWEI(A,B)

Author: REC

    Replaces NAG function G05DPF

Parameters: A,B			(both double precision, input)

Purpose: returns a pseudo-random number taken from a Weibull
distribution with shape parameter A and location parameter B i.e.
density is f(x) = ( A/(B**A) ) * x**(A-1) * EXP( -(x/B)**A )

Description: Uses the fact that if X is exponential then
(X/B)**A has the required distribution. Exponential calculation is done
within this routine to avoid some small extra cost of a call to ZBQLEXP.

Remarks: A and B must both be strictly positive.

Other routines called: ZBQLU01


11)	INTEGER FUNCTION		ZBQLNB(R,P)

Author: REC

    No single NAG function does this ...

Parameters: R,P			(both double precision, input)

Purpose: returns a pseudo-random number taken from a negative
binomial distribution with parameters R and P. Note that it is not
restricted to integer values of R. The form of the distribution is *not*
the `number of trials till rth success' version. The results are in the
range 0,1,2,... . The probability mass function used here has the form

      (x+r-1)!  x     r
P(X=x) =       -------- p (1-p)
            x!(r-1)!

as in Fishman (1978). To get the number of trials until the 4th success in a
sequence of Bernoulli trials with success probability 0.6, you'd call
ZBQLNB(4,0.4) and then add 4 to the result.

Description: The routine uses the fact that the negative binomial
can be regarded as a Gamma mixture of Poissons. See Fishman (1978) for
details.

Remarks: R must be strictly positive, and P must lie in (0,1).

Other routines called: ZBQLGAM, ZBQLPOI.

12)	DOUBLE PRECISION FUNCTION		ZBQLPAR(A,B)

Author: REC

    No single NAG function does this ...

Parameters: A,B			(both double precision, input)

Purpose: returns a pseudo-random number taken from a Pareto
distribution with parameters A and B. The density is

       a     -(a+1)
f(x) = a b (b+x)		for x > 0.

The Pareto distribution is usually defined slightly differently, over
the range x > b. If this is required, call ZBQLPAR and then add B to
the result.

Description: The routine uses the inverse cumulative method,
which is very straightforward as the quantile function is available
in closed form. A nd B must be strictly positive; otherwise an error
message is generated and the function returns zero.

Remarks: A and B must be strictly positive.

Other routines called: ZBQLU01.



AUXILIARY ROUTINES
==================

1)	DOUBLE PRECISION FUNCTION		ZBQLLG(X)

Author: REC

Returns the logarithm of the gamma function evaluated at X, using the
series expansion thingy in Press et al (1992). It works for all real numbers
(except negative integers where the thing has poles in any case).
