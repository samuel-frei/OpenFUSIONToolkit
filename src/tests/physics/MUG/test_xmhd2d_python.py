#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: the Python interface to the 2D xMHD model (OpenFUSIONToolkit.MUG.XMHD2D)
#
# Each case runs in its own process (one OFT_env per process) on a small periodic box with a
# smooth, non-uniform state. The Jacobian is checked against central differences of the
# residual, an advance against the step equation it solves, and restart files by round trip.
import os
import sys
import numpy as np
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing

MESH_GROUPS = {
    'mesh_options': {'meshname': "'box'", 'cad_type': '92', 'nlevels': '1', 'nbase': '1',
                     'grid_order': '1', 'fix_boundary': 'T'},
    'cube_options': {'mesh_type': '1', 'ref_per': 'T,F,F', 'ni': '4,3,0', 'rscale': '2.0,1.0,0.0',
                     'shift': '-1.0,-0.5,0.0'},
}
# The preconditioner the harris_sheet runs use: GMRES around block-Jacobi with local LU
XML = """<oft>
  <xmhd2d>
    <pre type="gmres">
      <its>2</its>
      <nrits>2</nrits>
      <pre type="block_jacobi">
        <nlocal>-1</nlocal>
        <solver type="lu"></solver>
      </pre>
    </pre>
  </xmhd2d>
</oft>
"""
PARAMS = dict(chi=50.0, eta=20.0, nu=100.0, gamma=5.0/3.0, D_diff=1.0, den_scale=1.e19)


def _state(sim, rng, flow=10.0, flux=1.e-5, spread=0.01):
    '''A positive density and temperature with random flow and flux (global numbering)

    The defaults are gentle enough for the Newton-GMRES solve to converge tightly; the
    Jacobian check uses a rougher state, since it needs no solve and a strong flow makes
    every velocity-dependent term large.
    '''
    f = sim.nfield
    blocks = [1.0 + spread*rng.random(f), flow*rng.standard_normal(f), np.zeros(f),
              flow*rng.standard_normal(f), 50.0*(1.0 + spread*rng.random(f)),
              flux*rng.standard_normal(f), np.zeros(f)]
    return np.concatenate(blocks)


def _xmhd2d_case(order, mp_q):
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.MUG import XMHD2D
        env = OFT_env(nthreads=2, quiet=True)
        with open('oft_in_xmhd2d_pytest.xml', 'w') as fid:
            fid.write(XML)
        sim = XMHD2D(env, MESH_GROUPS, order, xml_file='oft_in_xmhd2d_pytest.xml')
        sim.set_params(**PARAMS, nl_tol=1.e-10)
        rng = np.random.default_rng(0)
        u = _state(sim, rng)
        dt = 1.e-7
        out = {'nglobal': sim.nglobal, 'covered': bool(np.unique(sim.lg).size == sim.nglobal)}
        # Jacobian against central differences, per field block. At h = 1e-3 the truncation
        # error of an exact Jacobian is ~1e-9 and roundoff ~1e-8; a missing or wrong term
        # (the viscous density derivative was off by 1.6e-4 before it was fixed) is not
        rough = _state(sim, rng, flow=1.e3, flux=1.e-3, spread=0.1)
        jac = sim.jacobian(rough, dt)
        v = _state(sim, rng, spread=0.1) - _state(sim, rng, spread=0.1)
        jv = jac @ v
        h = 1.e-3
        fd = (sim.residual(rough + h*v, dt)[0] - sim.residual(rough - h*v, dt)[0])/(2*h)
        out['jac_err'] = max(float(np.linalg.norm(a - b)/np.linalg.norm(b))
                             for a, b in zip(np.split(jv, 7), np.split(fd, 7)) if np.linalg.norm(b) > 0)
        # The mass term is linear, and its Jacobian is the dt = 0 assembly
        mass = sim.jacobian(u, 0.0)
        out['mass_err'] = float(np.linalg.norm(mass @ u - sim.residual(u, 0.0)[0])/np.linalg.norm(mass @ u))
        # An advance solves F_dt(u_{n+1}) = F_0(u_n) at every step
        res = sim.advance(u, dt, 3, store=True)
        traj = res['trajectory']
        out['failed'] = res['failed_step']
        out['step_err'] = max(float(np.linalg.norm(sim.residual(traj[n + 1], dt)[0] - sim.residual(traj[n], 0.0)[0])
                                    / np.linalg.norm(sim.residual(traj[n], 0.0)[0])) for n in range(3))
        out['final_matches'] = bool(np.array_equal(traj[-1], res['u']))
        # Restart round trip
        sim.save_restart('xmhd2d_pytest.rst', traj[-1], t=3*dt, dt=dt)
        out['restart_err'] = float(np.abs(sim.load_restart('xmhd2d_pytest.rst') - traj[-1]).max())
        os.remove('xmhd2d_pytest.rst')
        os.remove('oft_in_xmhd2d_pytest.xml')
        mp_q.put(out)
    except Exception as err:
        mp_q.put({'error': repr(err)})
        raise


@pytest.mark.parametrize('order', (2, 3))
def test_xmhd2d_python(order):
    out = oft_testing.mp_run(_xmhd2d_case, (order,), test_dir)
    assert out is not None and 'error' not in out, out
    assert out['covered']
    assert out['jac_err'] < 1.e-6
    assert out['mass_err'] < 1.e-12
    assert out['failed'] == 0
    assert out['step_err'] < 1.e-8
    assert out['final_matches']
    assert out['restart_err'] == 0.0
