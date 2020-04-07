import os
import ctypes as ct
import numpy as np

def array_ctype(ndim, dtype=np.float64, flags="F"):
    return [ct.POINTER(ct.c_int)]*ndim + [np.ctypeslib.ndpointer(ndim=ndim, dtype=dtype, flags=flags)]

def array_arg(a):
    arr = a
    return (*(ct.c_int(s) for s in arr.shape), arr)

class HMx:
    libname = "libhmx_wrapper.so"
    module_name = "HMx_wrapper"

    def __init__(self):
        self.load_lib()

        self.constants = HMxConstants(self.lib)

    def load_lib(self, path=None):
        if path is None:
            path = os.path.dirname(__file__)
        libpath = os.path.abspath(os.path.join(path, self.libname))
        self.lib = ct.CDLL(libpath)

    def get_function(self, name, c_bind=True):
        if c_bind:
            return getattr(self.lib, name)
        else:
            return getattr(self.lib, f"__{self.module_name}_MOD_{name}")

    def _run_hmx(self, omm, omb, omv, h, ns, sigma8, mnu, w, wa,
                       halo_model_mode, Theat, eta0, As,
                       fields=None,
                       k=None, a=None,
                       k_lin=None, a_lin=None,
                       pk_lin=None,
                    #    q=None,
                    #    Tcold=None,
                       verbose=True):
        f = self.get_function("run_HMx")
        f.restype = ct.c_int
        f.argtypes = [ct.POINTER(ct.c_double),     # omm
                      ct.POINTER(ct.c_double),     # omb
                      ct.POINTER(ct.c_double),     # omv
                      ct.POINTER(ct.c_double),     # mnu
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # ns
                      ct.POINTER(ct.c_double),     # sigma8
                      ct.POINTER(ct.c_double),     # w
                      ct.POINTER(ct.c_double),     # wa
                      ct.POINTER(ct.c_int),        # halo_model_mode
                      ct.POINTER(ct.c_double),     # Theat
                      ct.POINTER(ct.c_double),     # eta0
                      ct.POINTER(ct.c_double),     # As
                      *array_ctype(ndim=1, dtype=np.int), # fields
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # a
                      *array_ctype(ndim=1, dtype=np.float64), # k_lin
                      *array_ctype(ndim=1, dtype=np.float64), # a_lin
                      *array_ctype(ndim=2, dtype=np.float64), # Pk_lin
                    #   *array_ctype(ndim=1, dtype=np.float64), # q
                    #   *array_ctype(ndim=2, dtype=np.float64), # Tcold
                      *array_ctype(ndim=4, dtype=np.float64), # Pk_hmx
                      ct.POINTER(ct.c_bool),        # verbose
                      ]      

        Pk_hmx = np.zeros((len(fields), len(fields), len(k), len(a)), dtype=np.float64, order="F")

        # q = q or k
        # Tcold = Tcold or np.ones((len(k), len(a)), order="F", dtype=np.float64)

        if len(a) > 1 and a[0] > a[1]:
            raise ValueError("Scale factor needs to be increasing.")

        status = f(ct.c_double(omm), ct.c_double(omb), ct.c_double(omv), ct.c_double(mnu), 
          ct.c_double(h), ct.c_double(ns), ct.c_double(sigma8), ct.c_double(w), ct.c_double(wa),
          ct.c_int(halo_model_mode), ct.c_double(Theat), ct.c_double(eta0), ct.c_double(As),
          *array_arg(fields),
          *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
          *array_arg(np.ascontiguousarray(a, dtype=np.float64)),
          *array_arg(np.ascontiguousarray(k_lin, dtype=np.float64)),
          *array_arg(np.ascontiguousarray(a_lin, dtype=np.float64)),
          *array_arg(np.asfortranarray(pk_lin, dtype=np.float64)),
        #   *array_arg(np.ascontiguousarray(q, dtype=np.float64)),
        #   *array_arg(np.asfortranarray(Tcold, dtype=np.float64)),
          *array_arg(Pk_hmx),
          ct.c_bool(verbose)
          )
        if status != 0:
            raise RuntimeError("HMx failed.")
        return Pk_hmx

class HMxConstants:
    def __init__(self, lib):
        for c in ["HMCode2016", "HMCode2020", 
                  "field_dmonly", "field_matter", "field_cdm", "field_gas", "field_stars", "field_electron_pressure"]:
            setattr(self, c, ct.c_int.in_dll(lib, f"constant_{c.lower()}").value)

if __name__ == "__main__":
    import camb

    hmx = HMx()

    h = 0.7
    omc = 0.25
    omb = 0.048
    mnu = 0.12
    w = -1.0
    wa = 0.0
    ns = 0.97
    As = 2.1e-9

    halo_model_mode = hmx.constants.HMCode2016
    Theat = 10**7.8
    A = 2.6
    eta0 = 0.98 - 0.12*As

    fields = np.array([hmx.constants.field_dmonly])

    # Get linear power spectrum
    camb.set_feedback_level(2)
    p = camb.CAMBparams(WantTransfer=True, NonLinearModel=camb.nonlinear.Halofit(halofit_version="mead", HMCode_A_baryon=A, HMCode_eta_baryon=eta0))
    p.set_cosmology(H0=h*100, omch2=omc*h**2, ombh2=omb*h**2, mnu=mnu)
    p.set_dark_energy(w=w)
    p.set_initial_power(camb.InitialPowerLaw(As=As, ns=ns))
    
    a_lin = np.linspace(0.1, 1.0, 16, endpoint=True)
    p.set_matter_power(redshifts=1/a_lin-1, kmax=10.0, nonlinear=False)

    r = camb.get_results(p)
    sigma8 = r.get_sigma8()[-1]
    k_lin, z_lin, pofk_lin_camb = r.get_matter_power_spectrum(minkh=1e-3, maxkh=2e2, npoints=128)

    pofk_lin_camb = pofk_lin_camb.T
    pofk_lin_camb = pofk_lin_camb[:,::-1]

    z_lin = np.array(z_lin[::-1])
    a_lin = 1/(1+z_lin)
    
    omv = r.omega_de + r.get_Omega("photon") + r.get_Omega("neutrino")
    omm = p.omegam

    # Ranges to get from HMx
    k = np.logspace(-3, 1, 100)
    a = np.linspace(0.2, 1.0, 5)

    Pk_HMx = hmx._run_hmx(omm, omb, omv, h, ns, sigma8, mnu, w, wa,
                          halo_model_mode, Theat, eta0, A,
                          fields=fields,
                          k=k,
                          a=a,
                          k_lin=k_lin, a_lin=a_lin,
                          pk_lin=pofk_lin_camb,
                          verbose=True)

    z = 1/a - 1
    p.set_matter_power(redshifts=z, kmax=k[-1], nonlinear=True)
    r = camb.get_results(p)
    Pk_nl_CAMB_interpolator = r.get_matter_power_interpolator()
    pofk_nonlin_camb = Pk_nl_CAMB_interpolator.P(z[::-1], k, grid=True)

    pofk_nonlin_camb = pofk_nonlin_camb.T
    pofk_nonlin_camb = pofk_nonlin_camb[:,::-1]

    import matplotlib.colorbar
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0, right=0.95)
    cmap = plt.get_cmap("inferno_r")

    cb_ax = matplotlib.colorbar.make_axes(ax)
    norm = matplotlib.colors.Normalize(vmin=a[0], vmax=a[-1])
    cb1 = matplotlib.colorbar.ColorbarBase(cb_ax[0], cmap=cmap,
                                    norm=norm, **cb_ax[1])
    cb1.set_label('a')

    for i in range(len(a_lin)):
        ax[0].loglog(k_lin, pofk_lin_camb[:,i], ls=":", c=cmap(a_lin[i]), label="Linear" if i == 0 else None)
    for i in range(len(a)):
        ax[0].loglog(k, pofk_nonlin_camb[:,i], ls="--", c=cmap(a[i]), label="HMCode CAMB" if i == 0 else None)
        ax[0].loglog(k, Pk_HMx[0,0,:,i], ls="-", c=cmap(a[i]), label="HMCode HMx" if i == 0 else None)

        ax[1].semilogx(k, Pk_HMx[0,0,:,i]/pofk_nonlin_camb[:,i]-1, c=cmap(a[i]))

    ax[0].legend(frameon=False)
    ax[0].set_ylabel("$P(k)$ [Mpc$^3$ $h^{-3}$]")
    ax[1].set_ylabel("Frac. diff. HMCode")
    ax[1].set_xlabel("$k$ [$h$ Mpc$^{-1}$]")

    ax[0].set_title("HMCode vs HMx (ihm = 1)")
    fig.savefig("HMCode_test_CAMB_vs_HMx.png", dpi=300)
