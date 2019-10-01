#!/usr/bin/env python3
import myokit
import matplotlib.pyplot as plt

m, p, _ = myokit.load('example')
s = myokit.Simulation(m, p)

d = s.run(100)
d = s.run(100, log=d)

s = myokit.Simulation(m, p)

d = s.run(100)

plt.figure()
plt.plot(d.time(), d['membrane.V'])

d = s.run(100, log=d)

