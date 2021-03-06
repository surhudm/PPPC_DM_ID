import numpy as np
from scipy.integrate import quad


def mu(x):
    return np.log(1.+x) - x/(1.+x)


gee = 4.2994e-9
Omegam = 0.315
rhocrit = 3.E4/(8.*np.pi*gee)
rhobar = rhocrit * Omegam

# Read milky way mass
Milky_way_mass = np.double(input("Mass of the milky way (Mvir) in Msun/h? \n"))
print "Milky way mass: ", Milky_way_mass

# Virial radius: Virial contrast for Planck is 104.2 wrt critical, see Bryan and
# Norman 92
Milky_way_radius = (Milky_way_mass/(4./3.*np.pi*rhocrit*104.2))**(1./3.)
print "Milky way radius: ", Milky_way_radius

# Concentration parameter
Milky_way_conc = 10.0 ** (1.025-0.097*np.log10(Milky_way_mass/1E12))
print "Milky way concentration: ", Milky_way_conc

# Scale radius
Milky_way_rs = Milky_way_radius/Milky_way_conc
print "Milky way rs: ", Milky_way_rs

# Rho_s
Milky_way_rhos = Milky_way_mass / \
    (4. * np.pi * Milky_way_rs**3 * mu(Milky_way_conc))
print "Milky way rhos: ", Milky_way_rhos

# Rsol (in hinv Mpc)
Rsol = 8.3e-3 * 0.7


def rho_iso(r, Milky_way_rhos, Milky_way_rs, Milky_way_radius):
    if r > Milky_way_radius:
        return 0
    return Milky_way_rhos/(r/Milky_way_rs)**2


def rho_nfw(r, Milky_way_rhos, Milky_way_rs, Milky_way_radius):
    if r > Milky_way_radius:
        return 0
    return Milky_way_rhos/(r/Milky_way_rs)/(1.+r/Milky_way_rs)**2


def Jfactor(
        theta,
        Milky_way_rhos,
        Milky_way_rs,
        Milky_way_radius,
        rhofunction,
        Rsol,
        rhopow=1):
    Rhosol = rhofunction(Rsol, Milky_way_rhos, Milky_way_rs, Milky_way_radius)

    def integ_func(
            s,
            theta,
            Milky_way_rhos,
            Milky_way_rs,
            Rsol,
            Rhosol,
            rhofunction,
            Milky_way_radius):
        r = (Rsol**2 + s**2 - 2 * Rsol * s * np.cos(theta))**0.5
        return 1. / Rsol * (rhofunction(r,
                                        Milky_way_rhos,
                                        Milky_way_rs,
                                        Milky_way_radius) / Rhosol) ** rhopow

    res, err = quad(
        integ_func, 0, np.inf, args=(
            theta, Milky_way_rhos, Milky_way_rs, Rsol, Rhosol, rhofunction, Milky_way_radius))

    return res

rho = rho_nfw
theta = np.arange(0.01, 181.0, 5.0) * np.pi/180.0
Jdecay = theta * 0.0
Jannihilate = theta * 0.0
for i in range(theta.size):
    Jdecay[i] = Jfactor(
        theta[i],
        Milky_way_rhos,
        Milky_way_rs,
        Milky_way_radius,
        rho,
        Rsol)
    Jannihilate[i] = Jfactor(
        theta[i],
        Milky_way_rhos,
        Milky_way_rs,
        Milky_way_radius,
        rho,
        Rsol,
        rhopow=2)

import pylab as pl
ax = pl.subplot(221)
ax.set_yscale("log")
ax.plot(theta * 180./np.pi, Jdecay, color="c")
ax.plot(theta * 180./np.pi, Jannihilate, color="m")
pl.savefig("Decay.pdf")


# Now integrate over the angles
def integrate_angle(spl, thtmin, thtmax):
    def func(x, spl, norm):
        return np.sin(x) * spl(x) * norm
    norm = 1./(np.cos(thtmin)-np.cos(thtmax))
    res, err = quad(func, thtmin, thtmax, args=(spl, norm))
    return res


from scipy.interpolate import interp1d

# Let us slightly cheat on the lower limit. Is this reasonable?
theta[0] = 0.0
Jdecay_spl = interp1d(theta, Jdecay, kind="cubic")
Jannihilate_spl = interp1d(theta, Jannihilate, kind="cubic")

print integrate_angle(Jdecay_spl, 0.0, np.pi)
print integrate_angle(Jannihilate_spl, 0.0, np.pi)
