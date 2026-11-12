"""
Pemodelan dan Simulasi Sistem Dinamik: Pengisian dan Pengosongan Kapasitor (Rangkaian RC)

Model:  RC * dVc/dt = -Vc + Vin(t)
Input:  Vin(t) = 0 V (t < 2 s), 5 V (2 s <= t < 12 s), 0 V (t >= 12 s)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ---------------------------------------------------------------
# 1. Parameter rangkaian
# ---------------------------------------------------------------
R = 10e3          # resistor (ohm)  = 10 kOhm
C = 220e-6        # kapasitor (farad) = 220 uF
tau = R * C       # konstanta waktu (s) = 2.2 s
Vs = 5.0          # tegangan sumber saat saklar ON (V)
T_ON, T_OFF = 2.0, 12.0   # saklar ON pada 2 s, OFF pada 12 s
V0 = 0.0          # kapasitor awalnya kosong


def vin(t):
    """Tegangan input: sumber 5 V aktif antara T_ON dan T_OFF."""
    return Vs if T_ON <= t < T_OFF else 0.0


# ---------------------------------------------------------------
# 2. Model ODE untuk odeint
# ---------------------------------------------------------------
def rc_model(vc, t, R, C):
    return (-vc + vin(t)) / (R * C)


# ---------------------------------------------------------------
# 3. Simulasi dengan ODEINT (per langkah 0.05 s agar step input terbaca)
# ---------------------------------------------------------------
t = np.linspace(0, 25, 501)
vc = np.zeros_like(t)
vc[0] = V0
for k in range(1, len(t)):
    # input dianggap konstan pada tiap interval [t[k-1], t[k]]
    t_mid = 0.5 * (t[k - 1] + t[k])
    u = vin(t_mid)
    z = odeint(lambda y, tt: (-y + u) / (R * C), vc[k - 1], [t[k - 1], t[k]])
    vc[k] = z[-1, 0]

vin_arr = np.array([vin(0.5 * (t[max(k - 1, 0)] + t[k])) for k in range(len(t))])
vin_arr[0] = 0.0
i_arr = (vin_arr - vc) / R * 1e3   # arus dalam mA

# ---------------------------------------------------------------
# 4. Solusi analitik sebagai pembanding
# ---------------------------------------------------------------
def analitik(tt):
    v_off = Vs * (1 - np.exp(-(T_OFF - T_ON) / tau))   # tegangan saat saklar OFF
    out = np.zeros_like(tt)
    for n, x in enumerate(tt):
        if x < T_ON:
            out[n] = V0
        elif x < T_OFF:
            out[n] = Vs + (V0 - Vs) * np.exp(-(x - T_ON) / tau)
        else:
            out[n] = v_off * np.exp(-(x - T_OFF) / tau)
    return out


vc_an = analitik(t)
galat_maks = np.max(np.abs(vc - vc_an))
print(f"tau = {tau:.2f} s | galat maksimum odeint vs analitik = {galat_maks:.2e} V")

# ---------------------------------------------------------------
# 5. Simulasi dengan GEKKO (opsional, perlu: pip install gekko)
# ---------------------------------------------------------------
try:
    from gekko import GEKKO
    m = GEKKO(remote=False)
    tg = np.linspace(0, 25, 2501)              # grid lebih rapat (dt = 0.01 s)
    m.time = tg
    # GEKKO memakai nilai input di akhir interval, sehingga input digeser satu langkah
    ustep = np.array([vin(x - (tg[1] - tg[0])) for x in tg])
    u_g = m.Param(value=ustep)
    v_g = m.Var(value=V0)
    m.Equation(R * C * v_g.dt() == -v_g + u_g)
    m.options.IMODE = 4
    m.solve(disp=False)
    vc_gekko_full = np.array(list(v_g.value), dtype=float)
    vc_gekko = vc_gekko_full[::5]               # ambil tiap 5 titik agar sama dengan grid odeint
    print(f"galat maksimum GEKKO vs analitik = {np.max(np.abs(vc_gekko - vc_an)):.2e} V")
except Exception as e:
    vc_gekko = None
    print("GEKKO dilewati:", e)

# ---------------------------------------------------------------
# 6. Grafik 1: respon tegangan kapasitor dan input
# ---------------------------------------------------------------
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax[0].plot(t, vc, 'b-', label='Vc(t) odeint')
ax[0].plot(t, vc_an, 'k:', label='Vc(t) analitik')
if vc_gekko is not None:
    ax[0].plot(t, vc_gekko, 'g--', lw=1, label='Vc(t) GEKKO')
ax[0].axhline(0.632 * Vs, color='gray', ls='--', lw=0.8)
ax[0].set_ylabel('Tegangan kapasitor Vc (V)')
ax[0].legend(frameon=False)
ax[0].grid(alpha=0.3)
ax[1].plot(t, vin_arr, 'r-')
ax[1].set_ylabel('Input Vin (V)')
ax[1].set_xlabel('t (s)')
ax[1].grid(alpha=0.3)
plt.tight_layout()
plt.savefig('grafik1_respon_tegangan.png', dpi=150)

# ---------------------------------------------------------------
# 7. Grafik 2: arus pengisian dan pengosongan
# ---------------------------------------------------------------
plt.figure(figsize=(8, 3.5))
plt.plot(t, i_arr, 'm-')
plt.axhline(0, color='k', lw=0.6)
plt.xlabel('t (s)')
plt.ylabel('Arus i(t) (mA)')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('grafik2_arus.png', dpi=150)

# ---------------------------------------------------------------
# 8. Grafik 3: pengaruh variasi R (variasi konstanta waktu)
# ---------------------------------------------------------------
plt.figure(figsize=(8, 4))
for R_var, sty in [(4.7e3, 'r-'), (10e3, 'b--'), (22e3, 'g:')]:
    tau_v = R_var * C
    v_var = np.zeros_like(t)
    for k in range(1, len(t)):
        u = vin(0.5 * (t[k - 1] + t[k]))
        z = odeint(lambda y, tt: (-y + u) / tau_v, v_var[k - 1], [t[k - 1], t[k]])
        v_var[k] = z[-1, 0]
    plt.plot(t, v_var, sty, label=f'R = {R_var/1e3:.1f} k$\\Omega$ ($\\tau$ = {tau_v:.2f} s)')
plt.xlabel('t (s)')
plt.ylabel('Vc (V)')
plt.legend(frameon=False)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('grafik3_variasi_R.png', dpi=150)

# ---------------------------------------------------------------
# 9. Tabel perhitungan manual: Euler (h = 1 s) vs analitik vs odeint
# ---------------------------------------------------------------
print("\nTabel perhitungan manual (Euler, h = 1 s):")
print(f"{'t':>4} {'Vin':>5} {'Vc Euler':>10} {'Vc analitik':>12} {'Vc odeint':>10}")
v_e = V0
h = 1.0
for n in range(0, 17):
    tn = n * h
    u_n = vin(tn)
    idx = int(round(tn / (t[1] - t[0])))
    print(f"{tn:4.0f} {u_n:5.1f} {v_e:10.4f} {analitik(np.array([tn]))[0]:12.4f} {vc[idx]:10.4f}")
    v_e = v_e + h * (-v_e + u_n) / tau
