# Pemodelan dan Simulasi Rangkaian RC (Pengisian dan Pengosongan Kapasitor)

Project mata kuliah Dynamic System Modelling. Model: `R*C*dVc/dt = -Vc + Vin(t)` dengan R = 10 kOhm, C = 220 uF (tau = 2,2 s) dan input step 0 -> 5 V (t = 2 s) -> 0 V (t = 12 s).

## Isi
- `rc_simulasi.ipynb` : notebook Google Colab (odeint, GEKKO, analitik, Euler)
- `rc_simulasi.py` : versi skrip Python
- `grafik*.png`, `diagram_rc.png` : hasil grafik dan diagram

## Menjalankan
```
pip install numpy scipy matplotlib gekko
python rc_simulasi.py
```
Atau buka `rc_simulasi.ipynb` di Google Colab lalu pilih Runtime > Run all.
