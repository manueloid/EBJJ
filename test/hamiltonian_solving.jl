N = 10
Λ0 = 2.5
U, Ωf = 2Λ0 / N, 0.1
t0, tf = 0.0002 / U, 0.005 / U
tfs = range(t0, tf, length=10)
c = ControlFull(N, Ωf, U, t0, 2, 2:2:2)
q = ConstantQuantity(c)
