using EBJJ
N = 400
Λ0 = 10
U, Ωf = 2Λ0 / N, 0.1
t0, tf = 0.0001 / U, 0.005 / U
tfs = range(t0, tf, length=10)
c = ControlFull(N, Ωf, U, tf, 2, 2:2:2)
cs = ControlSTA(c)
q = ConstantQuantity(c)
corrections(corrections(c))
