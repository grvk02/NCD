from astropy import units as u
from poliastro.plotting._base import BaseOrbitPlotter

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
r = [-6045, -3490, 2500] * u.km
v = [-3.457, 6.618, 2.533] * u.km / u.s

orb = Orbit.from_vectors(Earth, r, v)
orb.plot()