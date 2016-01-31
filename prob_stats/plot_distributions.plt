set ticslevel 0
set isosample 100
splot sin(sqrt(x**2 + y**2))**2./sqrt(x**2 + y**2)**3.
splot sin(sqrt(x**2 + y**2))**2.*cos(sqrt(x**2 + y**2))**2./sqrt(x**2 + y**2)**3.
splot cos(sqrt(x**2 + y**2))**2.*exp(-sqrt(x**2 + y**2))
