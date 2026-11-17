#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python interface for the 2D MUG xMHD model (`xmhd_2d`)

Runs the model in-process, so a loop that needs many short advances, residuals or
Jacobians (an optimisation, an adjoint, a Newton-Krylov continuation) pays for the
mesh and finite element setup once instead of once per executable launch.

The state is the composite 7-field vector `(n, v_x, v_y, v_z, T, psi, b_y)`. In
serial OFT keeps periodic seam copies of a degree of freedom in the local vector,
so every method works in *global* numbering: the same numbering `mat_save` writes
(`lg`) and the `nonmodal` package reads. A local vector is `u_g[lg]`; a local
matrix becomes the global operator `P^T J P` with `P[i, lg[i]] = 1`.

Only one simulation may exist per process, because `xmhd_2d` keeps its finite
element spaces in module variables.

@ingroup doxy_oft_python
'''
import ctypes
import os
import re
import numpy
from scipy import sparse
from .._interface import *

## @cond
# xmhd2d_setup(xmhd2d_ptr,mesh_ptr,order,xml_ptr,nlocal,nglobal,xml_pre,error_str)
xmhd2d_setup = ctypes_subroutine(oftpy_lib.xmhd2d_setup,
    [c_void_ptr_ptr, c_void_p, c_int, c_void_p, c_int_ptr, c_int_ptr, c_bool_ptr, c_char_p])
# xmhd2d_set_params(xmhd2d_ptr,chi,eta,nu,gamma,d_diff,den_scale,b0,lin_tol,nl_tol,error_str)
xmhd2d_set_params = ctypes_subroutine(oftpy_lib.xmhd2d_set_params,
    [c_void_p, c_double, c_double, c_double, c_double, c_double, c_double,
     ctypes_numpy_array(float64,1), c_double, c_double, c_char_p])
# xmhd2d_get_lg(xmhd2d_ptr,lg,error_str)
xmhd2d_get_lg = ctypes_subroutine(oftpy_lib.xmhd2d_get_lg,
    [c_void_p, ctypes_numpy_array(int32,1), c_char_p])
# xmhd2d_state(xmhd2d_ptr,vals,set,error_str)
xmhd2d_state = ctypes_subroutine(oftpy_lib.xmhd2d_state,
    [c_void_p, ctypes_numpy_array(float64,1), c_bool, c_char_p])
# xmhd2d_restart(xmhd2d_ptr,filename,vals,t,dt,save,error_str)
xmhd2d_restart = ctypes_subroutine(oftpy_lib.xmhd2d_restart,
    [c_void_p, c_char_p, ctypes_numpy_array(float64,1), c_double, c_double, c_bool, c_char_p])
# xmhd2d_residual(xmhd2d_ptr,u,dt,res,diag,error_str)
xmhd2d_residual = ctypes_subroutine(oftpy_lib.xmhd2d_residual,
    [c_void_p, ctypes_numpy_array(float64,1), c_double, ctypes_numpy_array(float64,1),
     ctypes_numpy_array(float64,1), c_char_p])
# xmhd2d_jacobian(xmhd2d_ptr,u,dt,nr,nnz,error_str)
xmhd2d_jacobian = ctypes_subroutine(oftpy_lib.xmhd2d_jacobian,
    [c_void_p, ctypes_numpy_array(float64,1), c_double, c_int_ptr, c_int_ptr, c_char_p])
# xmhd2d_jacobian_get(xmhd2d_ptr,kr,lc,m,error_str)
xmhd2d_jacobian_get = ctypes_subroutine(oftpy_lib.xmhd2d_jacobian_get,
    [c_void_p, ctypes_numpy_array(int32,1), ctypes_numpy_array(int32,1),
     ctypes_numpy_array(float64,1), c_char_p])
# xmhd2d_advance(xmhd2d_ptr,u,dt,nsteps,traj,diag,iters,verbose,nfailed,error_str)
xmhd2d_advance = ctypes_subroutine(oftpy_lib.xmhd2d_advance,
    [c_void_p, ctypes_numpy_array(float64,1), c_double, c_int, c_void_p,
     ctypes_numpy_array(float64,2), ctypes_numpy_array(int32,2), c_bool, c_int_ptr, c_char_p])
## @endcond

#: Field names in state-vector order
FIELD_NAMES = ('n', 'velx', 'vely', 'velz', 'T', 'psi', 'by')


def read_namelists(path):
    '''! Read every `&group ... /` of a Fortran namelist file as raw strings

    Values are kept as written (`'test'`, `T`, `6.4,3.2,0.0`), so groups can be
    copied verbatim into an OFT input file.

    @param path Input deck
    @result `{group: {key: value}}`, keys lower case
    '''
    with open(path) as fid:
        text = fid.read()
    groups = {}
    for name, body in re.findall(r'&(\w+)(.*?)\n\s*/', text, flags=re.S):
        entries = {}
        # Split on `key=` so values may themselves contain commas
        parts = re.split(r'(\b\w+)\s*=', body)
        for key, value in zip(parts[1::2], parts[2::2]):
            entries[key.lower()] = value.strip().rstrip(',').strip()
        groups[name.lower()] = entries
    return groups


def _float(value):
    return float(value.replace('d', 'e').replace('D', 'E'))


class XMHD2D():
    '''! In-process 2D MUG xMHD simulation

    Typical use, from a `harris_sheet` deck:

        env = OFT_env(nthreads=8)
        sim = XMHD2D.from_deck(env, 'oft.in', xml_file='oft_in.xml')
        u = sim.load_restart('xmhd2d_00000.rst')
        out = sim.advance(u, dt=1e-7, nsteps=100, store=True)
    '''

    def __init__(self, oft_env, mesh_groups, order, xml_file=None):
        '''! Build the mesh from namelist groups, then the simulation on it

        @param oft_env OFT runtime environment (@ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        @param mesh_groups `{group: {key: value}}` defining the mesh (`mesh_options` and,
        for generated meshes, `cube_options`), written verbatim into the input file
        @param order Lagrange element order
        @param xml_file XML file with an `<xmhd2d><pre>` preconditioner definition (optional)
        '''
        self._oft_env = oft_env
        for group, entries in mesh_groups.items():
            self._oft_env.oft_in_groups[group] = dict(entries)
        self._oft_env.update_oft_in()
        # ndim=0 builds the mesh from the input file rather than from arrays
        rfake = numpy.ones((1, 1), dtype=numpy.float64)
        lcfake = numpy.ones((1, 1), dtype=numpy.int32)
        regfake = numpy.ones((1,), dtype=numpy.int32)
        nregs = c_int()
        self._mesh_ptr = c_void_p()
        oft_setup_smesh(0, 0, rfake, 0, 0, lcfake, regfake, ctypes.byref(nregs), ctypes.byref(self._mesh_ptr))
        xml_ptr = c_void_p()
        if xml_file is not None:
            if not os.path.isfile(xml_file):
                raise FileNotFoundError(xml_file)
            oftpy_load_xml(self._oft_env.path2c(os.path.abspath(xml_file)), ctypes.byref(xml_ptr))
            if not xml_ptr:
                raise ValueError('Unable to parse XML file "{0}"'.format(xml_file))
        self._ptr = c_void_p()
        nlocal = c_int()
        nglobal = c_int()
        xml_pre = c_bool()
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_setup(ctypes.byref(self._ptr), self._mesh_ptr, order, xml_ptr,
                     ctypes.byref(nlocal), ctypes.byref(nglobal), ctypes.byref(xml_pre), error_string)
        self._check(error_string)
        ## Preconditioner of the Newton steps' GMRES: 'xml' (from `xml_file`) or 'diagonal'
        self.preconditioner = 'xml' if xml_pre.value else 'diagonal'
        if xml_file is not None and not xml_pre.value:
            raise ValueError('"{0}" has no <xmhd2d><pre> preconditioner definition'.format(xml_file))
        ## Local size of the state vector
        self.nlocal = nlocal.value
        ## Global size of the state vector
        self.nglobal = nglobal.value
        ## Element order
        self.order = order
        lg = numpy.zeros((self.nlocal,), dtype=numpy.int32)
        xmhd2d_get_lg(self._ptr, lg, error_string)
        self._check(error_string)
        ## Global index of each local entry
        self.lg = lg.astype(numpy.int64)
        if self.lg.min() < 0 or self.lg.max() >= self.nglobal:
            raise ValueError('local-to-global map out of range')
        if numpy.unique(self.lg).size != self.nglobal:
            raise ValueError('local-to-global map does not cover the global vector')
        ## Local-to-global scatter `P` (`nlocal x nglobal`)
        self.scatter = sparse.csr_array(
            (numpy.ones(self.nlocal), (numpy.arange(self.nlocal), self.lg)),
            shape=(self.nlocal, self.nglobal))
        # One representative local entry per global index
        self._rep = numpy.zeros((self.nglobal,), dtype=numpy.int64)
        self._rep[self.lg[::-1]] = numpy.arange(self.nlocal)[::-1]
        ## Number of degrees of freedom per field (global)
        self.nfield = self.nglobal // len(FIELD_NAMES)
        ## Physical parameters last set
        self.params = {}

    @classmethod
    def from_deck(cls, oft_env, deck, xml_file=None, options_group='harris_sheet_options',
                  lin_tol=-1.0, nl_tol=-1.0):
        '''! Build a simulation from a `harris_sheet` input deck

        Reads the mesh from `mesh_options`/`cube_options` and the order and transport
        coefficients from `options_group`, so the result matches what the
        `harris_sheet` executable would run with the same deck.
        '''
        groups = read_namelists(deck)
        mesh_groups = {name: groups[name] for name in ('mesh_options', 'cube_options') if name in groups}
        if 'mesh_options' not in mesh_groups:
            raise ValueError('"{0}" has no mesh_options group'.format(deck))
        if options_group not in groups:
            raise ValueError('"{0}" has no {1} group'.format(deck, options_group))
        opts = groups[options_group]
        sim = cls(oft_env, mesh_groups, int(opts['order']), xml_file=xml_file)
        sim.set_params(chi=_float(opts['chi']), eta=_float(opts['eta']), nu=_float(opts['nu']),
                       gamma=_float(opts['gamma']), D_diff=_float(opts['d_diff']),
                       den_scale=_float(opts['den_scale']), lin_tol=lin_tol, nl_tol=nl_tol)
        ## Deck the simulation was built from
        sim.deck = os.path.abspath(deck)
        ## Deck timestep
        sim.deck_dt = _float(opts['dt'])
        return sim

    def _check(self, error_string):
        if error_string.value != b'':
            raise Exception(error_string.value.decode())

    def set_params(self, chi, eta, nu, gamma, D_diff, den_scale, B0=(0.0, 0.0, 0.0),
                   lin_tol=-1.0, nl_tol=-1.0):
        '''! Set transport coefficients (uniform over regions) and solver tolerances

        Tolerances <= 0 keep OFT's defaults (linear 1e-8, Newton 1e-5, absolute).
        '''
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_set_params(self._ptr, chi, eta, nu, gamma, D_diff, den_scale,
                          numpy.ascontiguousarray(B0, dtype=numpy.float64), lin_tol, nl_tol, error_string)
        self._check(error_string)
        self.params = dict(chi=chi, eta=eta, nu=nu, gamma=gamma, D_diff=D_diff,
                           den_scale=den_scale, B0=tuple(B0), lin_tol=lin_tol, nl_tol=nl_tol)

    def to_global(self, u_local):
        '''! Global vector from a local one (seam copies must agree)'''
        return numpy.ascontiguousarray(numpy.asarray(u_local, dtype=numpy.float64)[self._rep])

    def to_local(self, u_global):
        '''! Local vector from a global one'''
        return numpy.ascontiguousarray(numpy.asarray(u_global, dtype=numpy.float64)[self.lg])

    def fields(self, u_global):
        '''! Split a global vector into `{name: block}`'''
        return dict(zip(FIELD_NAMES, numpy.split(numpy.asarray(u_global), len(FIELD_NAMES))))

    def load_restart(self, filename):
        '''! Read a restart file written by `xmhd_2d`; returns the global state'''
        if not os.path.isfile(filename):
            raise FileNotFoundError(filename)
        vals = numpy.zeros((self.nlocal,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_restart(self._ptr, self._oft_env.path2c(filename), vals, 0.0, 0.0, False, error_string)
        self._check(error_string)
        return self.to_global(vals)

    def save_restart(self, filename, u_global, t=0.0, dt=0.0):
        '''! Write a global state as an `xmhd_2d` restart file'''
        vals = self.to_local(u_global)
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_restart(self._ptr, self._oft_env.path2c(filename), vals, t, dt, True, error_string)
        self._check(error_string)

    def residual(self, u_global, dt):
        '''! Global time-advance function `F_dt(u)` and OFT's diagnostics of `u`

        A backward-Euler step solves `F_dt(u_{n+1}) = F_0(u_n)`; `F_0` is the mass term.
        Diagnostics: magnetic energy [J*2*mu0], kinetic energy [J], integrated density,
        integrated temperature, area.
        '''
        res = numpy.zeros((self.nlocal,), dtype=numpy.float64)
        diag = numpy.zeros((5,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_residual(self._ptr, self.to_local(u_global), dt, res, diag, error_string)
        self._check(error_string)
        return self.to_global(res), diag

    def jacobian(self, u_global, dt):
        '''! Global Jacobian of `F_dt` at `u`, as `P^T J_local P` (CSR)'''
        nr = c_int()
        nnz = c_int()
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_jacobian(self._ptr, self.to_local(u_global), dt, ctypes.byref(nr), ctypes.byref(nnz), error_string)
        self._check(error_string)
        kr = numpy.zeros((nr.value + 1,), dtype=numpy.int32)
        lc = numpy.zeros((nnz.value,), dtype=numpy.int32)
        vals = numpy.zeros((nnz.value,), dtype=numpy.float64)
        xmhd2d_jacobian_get(self._ptr, kr, lc, vals, error_string)
        self._check(error_string)
        local = sparse.csr_array((vals, lc, kr), shape=(nr.value, nr.value))
        return (self.scatter.T @ local @ self.scatter).tocsr()

    def advance(self, u_global, dt, nsteps, store=False, verbose=False):
        '''! Advance a global state by `nsteps` backward-Euler steps of fixed size `dt`

        Uses OFT's own Newton solver and preconditioner, as `run_simulation` does, but
        never adapts the timestep and writes nothing to disk. `verbose` prints OFT's
        Newton and GMRES progress.

        @result dict with `u` (final global state), `trajectory` (`(nsteps+1, nglobal)`
        global states if `store`, else None), `energy` (`(nsteps+1, 2)` magnetic
        [J*2*mu0] and kinetic [J]), `iterations` (`(nsteps, 2)` Newton and linear),
        `failed_step` (0 on success)
        '''
        vals = self.to_local(u_global)
        diag = numpy.zeros((nsteps + 1, 2), dtype=numpy.float64)
        iters = numpy.zeros((nsteps, 2), dtype=numpy.int32)
        traj = None
        traj_ptr = c_void_p()
        if store:
            traj = numpy.zeros((nsteps + 1, self.nlocal), dtype=numpy.float64)
            traj_ptr = traj.ctypes.data_as(c_void_p)
        nfailed = c_int()
        error_string = self._oft_env.get_c_errorbuff()
        xmhd2d_advance(self._ptr, vals, dt, nsteps, traj_ptr, diag, iters, verbose, ctypes.byref(nfailed), error_string)
        self._check(error_string)
        out = {
            'u': self.to_global(vals),
            'trajectory': None if traj is None else numpy.ascontiguousarray(traj[:, self._rep]),
            'energy': diag,
            'iterations': iters,
            'failed_step': nfailed.value,
        }
        return out
