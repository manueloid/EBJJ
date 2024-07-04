using EBJJ, Test, QuadGK
N = 50
Λ0 = 10
U, Ωf = 2Λ0 / N, 0.1
t0, tf = 0.001 / U, 0.04 / U
tfs = range(t0, tf, length=102)
c = ControlFull(N, Ωf, U, t0, 2, 2:2:2)
q = ConstantQuantity(c)

corrs = corrections(corrections(c))
