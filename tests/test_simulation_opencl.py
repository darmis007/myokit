#!/usr/bin/env python
#
# Tests the OpenCL simulation classes
#
# This file is part of Myokit
#  Copyright 2011-2018 Maastricht University, University of Oxford
#  Licensed under the GNU General Public License v3.0
#  See: http://myokit.org
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import os
import unittest
import numpy as np

import myokit

from shared import OpenCL_FOUND, DIR_DATA


# Show simulation output
debug = False


class SimulationOpenCL1dTest(unittest.TestCase):
    """
    Tests the OpenCL simulation in 1d mode.
    """
    '''
    def test_set_state(self):
        if not OpenCL_FOUND:
            print('OpenCL support not found, skipping test.')
            return

        m, p, x = myokit.load('example')
        n = 10
        s = myokit.SimulationOpenCL(m, p, n)
        sm = m.state()
        ss = [s.state(x) for x in range(n)]
        for si in ss:
            self.assertEqual(sm, si)

        # Test setting a single, global state
        sx = [0.0] * 8
        self.assertNotEqual(sm, sx)
        s.set_state(sx)
        for i in range(n):
            self.assertEqual(sx, s.state(i))
        self.assertEqual(sx * n, s.state())
        s.set_state(sm)
        self.assertEqual(sm * n, s.state())
        # Test setting a single state
        j = 1
        s.set_state(sx, j)
        for i in range(n):
            if i == j:
                self.assertEqual(s.state(i), sx)
            else:
                self.assertEqual(s.state(i), sm)

    def test_sim(self):
        if not OpenCL_FOUND:
            print('OpenCL support not found, skipping test.')
            return

        m, p, x = myokit.load('example')
        s = myokit.SimulationOpenCL(m, p, 20)

        # Run, log state and intermediary variable (separate logging code!)
        d = s.run(1, log=['engine.time', 'membrane.V', 'ina.INa'])
        self.assertIn('engine.time', d)
        self.assertIn('0.membrane.V', d)
        self.assertIn('19.membrane.V', d)
        self.assertIn('0.ina.INa', d)
        self.assertIn('19.ina.INa', d)
        self.assertEqual(len(d), 41)


@unittest.skipIf(not OpenCL_FOUND, 'OpenCL not found on this system.')
class SimulationOpenCL2dTest(unittest.TestCase):
    """
    Tests the OpenCL simulation in 2d mode.
    """
    def test_sim(self):
        m, p, x = myokit.load('example')
        n = (8, 8)
        s = myokit.SimulationOpenCL(m, p, n)
        s.set_paced_cells(4, 4)

        # Run, log state and intermediary variable (separate logging code!)
        d = s.run(1, log=['engine.time', 'membrane.V', 'ina.INa'])
        self.assertEqual(len(d), 129)
        self.assertIn('engine.time', d)
        self.assertIn('0.0.membrane.V', d)
        self.assertIn('7.7.membrane.V', d)
        self.assertIn('0.0.ina.INa', d)
        self.assertIn('7.7.ina.INa', d)


@unittest.skipIf(not OpenCL_FOUND, 'OpenCL not found on this system.')
class FiberTissueSimulationTest(unittest.TestCase):
    """
    Tests the fiber-tissue simulation.
    """
    def test_simple(self):
        # Load models
        mf = os.path.join(DIR_DATA, 'dn-1985-normalised.mmt')
        mt = os.path.join(DIR_DATA, 'lr-1991.mmt')
        mf = myokit.load_model(mf)
        mt = myokit.load_model(mt)

        # Run times
        run = .1

        # Create pacing protocol
        p = myokit.pacing.blocktrain(1000, 2.0, offset=.01)

        # Fiber/Tissue sizes
        nfx = 8
        nfy = 4
        ntx = 8
        nty = 6

        # Create simulation
        s = myokit.FiberTissueSimulation(
            mf,
            mt,
            p,
            ncells_fiber=(nfx, nfy),
            ncells_tissue=(ntx, nty),
            nx_paced=10,
            g_fiber=(235, 100),
            g_tissue=(9, 5),
            g_fiber_tissue=9
        )
        s.set_step_size(0.0012)
        # Set up logging
        logf = [
            'engine.time',
            'membrane.V',
            'isi.isiCa',
        ]
        logt = [
            'membrane.V',
            'ica.Ca_i',
            'ica.ICa',
        ]
        # Run simulation
        with myokit.PyCapture():
            logf, logt = s.run(run, logf=logf, logt=logt, log_interval=0.01)

        self.assertEqual(len(logf), 1 + 2 * nfx * nfy)
        self.assertIn('engine.time', logf)
        self.assertIn('0.0.membrane.V', logf)
        self.assertIn(str(nfx - 1) + '.' + str(nfy - 1) + '.membrane.V', logf)
        self.assertIn('0.0.isi.isiCa', logf)
        self.assertIn(str(nfx - 1) + '.' + str(nfy - 1) + '.isi.isiCa', logf)

        self.assertEqual(len(logt), 3 * ntx * nty)
        self.assertIn('0.0.membrane.V', logt)
        self.assertIn(str(ntx - 1) + '.' + str(nty - 1) + '.membrane.V', logt)
        self.assertIn('0.0.ica.Ca_i', logt)
        self.assertIn(str(ntx - 1) + '.' + str(nty - 1) + '.ica.Ca_i', logt)
        self.assertIn('0.0.ica.ICa', logt)
        self.assertIn(str(ntx - 1) + '.' + str(nty - 1) + '.ica.ICa', logt)
    '''

    def test_against_cvode(self):
        # Compare the SimulationOpenCL output with CVODE output

        # Load model and event with appropriate level/duration
        m, p, _ = myokit.load(os.path.join(DIR_DATA, 'lr-1991.mmt'))
        e = p.head()

        #
        # Make protocol to compare t=0 event implementations
        #
        dt = 0.1    # Note: this is too big to be very accurate
        tmax = 1300
        p = myokit.Protocol()
        p.schedule(level=e.level(), duration=e.duration(), period=600, start=0)
        logvars = ['engine.time', 'membrane.V', 'engine.pace']

        s1 = myokit.SimulationOpenCL(
            m, p, ncells=1, rl=True, precision=myokit.DOUBLE_PRECISION)
        s1.set_step_size(dt)
        d1 = s1.run(tmax, logvars, log_interval=dt).npview()

        s2 = myokit.Simulation(m, p)
        s2.set_tolerance(1e-8, 1e-8)
        d2 = s2.run(tmax, logvars, log_interval=dt).npview()

        # Check implementation of logging point selection
        e0 = np.max(np.abs(d1.time() - d2.time()))

        # Check implementation of pacing
        r1 = d1['engine.pace'] - d2['engine.pace']
        e1 = np.sum(r1**2)

        # Check membrane potential (will have some error!)
        # Using MRMS from Marsh, Ziaratgahi, Spiteri 2012
        r2 = d1['membrane.V', 0] - d2['membrane.V']
        r2 /= (1 + np.abs(d2['membrane.V']))
        e2 = np.sqrt(np.sum(r2**2) / len(r2))

        if debug:
            import matplotlib.pyplot as plt
            print('Event at t=0')

            plt.figure()
            plt.suptitle('Time points')
            plt.plot(d1.time(), label='Euler')
            plt.plot(d1.time(), label='CVODE')
            plt.legend()
            print(d1.time()[:7])
            print(d2.time()[:7])
            print(d1.time()[-7:])
            print(d2.time()[-7:])
            print(e0)

            plt.figure()
            plt.suptitle('Pacing signals')
            plt.subplot(2, 1, 1)
            plt.plot(d1.time(), d1['engine.pace'], label='OpenCL')
            plt.plot(d2.time(), d2['engine.pace'], label='CVODE')
            plt.legend()
            plt.subplot(2, 1, 2)
            plt.plot(d1.time(), r1)
            print(e1)

            plt.figure()
            plt.suptitle('Membrane potential')
            plt.subplot(2, 1, 1)
            plt.plot(d1.time(), d1['membrane.V', 0], label='OpenCL')
            plt.plot(d2.time(), d2['membrane.V'], label='CVODE')
            plt.legend()
            plt.subplot(2, 1, 2)
            plt.plot(d1.time(), r2)
            print(e2)

            plt.show()

        self.assertLess(e0, 1e-14)
        self.assertLess(e1, 1e-14)
        self.assertLess(e2, 0.1)   # Note: The step size is really too big here

        #
        # Make protocol to compare with event at step-size multiple
        #
        dt = 0.1   # Note: this is too big to be very accurate
        tmax = 1300
        p = myokit.Protocol()
        p.schedule(
            level=e.level(),
            duration=e.duration(),
            period=600,
            start=10)

        logvars = ['engine.time', 'membrane.V', 'engine.pace']

        s1.reset()
        s1.set_protocol(p)
        s1.set_step_size(dt)
        d1 = s1.run(tmax, logvars, log_interval=dt).npview()

        s2.reset()
        s2.set_protocol(p)
        s2.set_tolerance(1e-8, 1e-8)
        d2 = s2.run(tmax, logvars, log_interval=dt).npview()

        # Check implementation of logging point selection
        e0 = np.max(np.abs(d1.time() - d2.time()))

        # Check implementation of pacing
        r1 = d1['engine.pace', 0] - d2['engine.pace']
        e1 = np.sum(r1**2)

        # Check membrane potential (will have some error!)
        # Using MRMS from Marsh, Ziaratgahi, Spiteri 2012
        r2 = d1['membrane.V', 0] - d2['membrane.V']
        r2 /= (1 + np.abs(d2['membrane.V']))
        e2 = np.sqrt(np.sum(r2**2) / len(r2))

        if debug:
            import matplotlib.pyplot as plt
            print('Event at t=0')

            plt.figure()
            plt.suptitle('Time points')
            plt.plot(d1.time(), label='OpenCL')
            plt.plot(d1.time(), label='CVODE')
            plt.legend()
            print(d1.time()[:7])
            print(d2.time()[:7])
            print(d1.time()[-7:])
            print(d2.time()[-7:])
            print(e0)

            plt.figure()
            plt.suptitle('Pacing signals')
            plt.subplot(2, 1, 1)
            plt.plot(d1.time(), d1['engine.pace', 0], label='OpenCL')
            plt.plot(d2.time(), d2['engine.pace'], label='CVODE')
            plt.legend()
            plt.subplot(2, 1, 2)
            plt.plot(d1.time(), r1)
            print(e1)

            plt.figure()
            plt.suptitle('Membrane potential')
            plt.subplot(2, 1, 1)
            plt.plot(d1.time(), d1['membrane.V', 0], label='OpenCL')
            plt.plot(d2.time(), d2['membrane.V'], label='CVODE')
            plt.legend()
            plt.subplot(2, 1, 2)
            plt.plot(d1.time(), r2)
            print(e2)

            plt.show()

        self.assertLess(e0, 1e-14)
        self.assertLess(e1, 1e-14)
        self.assertLess(e2, 0.1)   # Note: The step size is really too big here

        #
        # Make protocol to compare with event NOT at step-size multiple
        #
        dt = 0.1   # Note: this is too big to be very accurate
        tmax = 1300
        p = myokit.Protocol()
        p.schedule(
            level=e.level(),
            duration=e.duration(),
            period=600,
            start=1.05)

        logvars = ['engine.time', 'membrane.V', 'engine.pace']

        s1.reset()
        s1.set_protocol(p)
        s1.set_step_size(dt)
        d1 = s1.run(tmax, logvars, log_interval=dt).npview()

        s2.reset()
        s2.set_protocol(p)
        s2.set_tolerance(1e-8, 1e-8)
        d2 = s2.run(tmax, logvars, log_interval=dt).npview()

        # Check implementation of logging point selection
        e0 = np.max(np.abs(d1.time() - d2.time()))

        # Check implementation of pacing
        r1 = d1['engine.pace', 0] - d2['engine.pace']
        e1 = np.sum(r1**2)

        # Check membrane potential (will have some error!)
        # Using MRMS from Marsh, Ziaratgahi, Spiteri 2012
        r2 = d1['membrane.V', 0] - d2['membrane.V']
        r2 /= (1 + np.abs(d2['membrane.V']))
        e2 = np.sqrt(np.sum(r2**2) / len(r2))

        if debug:
            import matplotlib.pyplot as plt
            print('Event at t=0')

            plt.figure()
            plt.suptitle('Time points')
            plt.plot(d1.time(), label='Euler')
            plt.plot(d1.time(), label='CVODE')
            plt.legend()
            print(d1.time()[:7])
            print(d2.time()[:7])
            print(d1.time()[-7:])
            print(d2.time()[-7:])
            print(e0)

            plt.figure()
            plt.suptitle('Pacing signals')
            plt.subplot(2, 1, 1)
            plt.plot(d1.time(), d1['engine.pace', 0], label='Euler')
            plt.plot(d2.time(), d2['engine.pace'], label='CVODE')
            plt.legend()
            plt.subplot(2, 1, 2)
            plt.plot(d1.time(), r1)
            print(e1)

            plt.figure()
            plt.suptitle('Membrane potential')
            plt.subplot(2, 1, 1)
            plt.plot(d1.time(), d1['membrane.V', 0], label='Euler')
            plt.plot(d2.time(), d2['membrane.V'], label='CVODE')
            plt.legend()
            plt.subplot(2, 1, 2)
            plt.plot(d1.time(), r2)
            print(e2)

            plt.show()

        self.assertLess(e0, 1e-14)
        self.assertLess(e1, 1e-14)
        self.assertLess(e2, 0.2)   # Note: The step size is really too big here

        #
        # Test again, with event at step multiple, but now without Rush-Larsen
        #
        dt = 0.01   # Note: this is too big to be very accurate
        tmax = 200
        p = myokit.Protocol()
        p.schedule(
            level=e.level(),
            duration=e.duration(),
            period=600,
            start=10)

        logvars = ['engine.time', 'membrane.V', 'engine.pace']
        s1 = myokit.SimulationOpenCL(
            m, p, ncells=1, rl=False, precision=myokit.DOUBLE_PRECISION)
        s1.set_step_size(dt)
        d1 = s1.run(tmax, logvars, log_interval=dt).npview()
        s2.reset()
        s2.set_protocol(p)
        s2.set_tolerance(1e-8, 1e-8)
        d2 = s2.run(tmax, logvars, log_interval=dt).npview()

        # Check membrane potential (will have some error!)
        # Using MRMS from Marsh, Ziaratgahi, Spiteri 2012
        r2 = d1['membrane.V', 0] - d2['membrane.V']
        r2 /= (1 + np.abs(d2['membrane.V']))
        e2 = np.sqrt(np.sum(r2**2) / len(r2))
        self.assertLess(e2, 0.05)  # Note: The step size is really too big here



if __name__ == '__main__':
    import sys
    if '-v' in sys.argv:
        print('Running in debug/verbose mode')
        debug = True
    else:
        print('Add -v for more debug output')
    unittest.main()
