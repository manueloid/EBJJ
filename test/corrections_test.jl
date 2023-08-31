#=
# Testing the corrections

In other tests I have alredy proved that the simplificatons that I made are correct.
Here I am going to test how the corrections are calculated.

The plan is to calculate the corrections in two different ways and compare the results.
The first way is to calculate the corrections by evaluating them using the full wave function, something of the form $ \chi_n^{*}(t,z) \hat{H} \chi_0(t,z) $, where $ \hat{H} $ is any of the operators I need to use and the function $ \chi_n $ are the full STA wave function with no simplificatons made.

The second one is to evaluate the corrections using the simplified wave function, something of the form $ \phi_n(t) \theta_n(t,z)^{*} \hat{H} \theta_0(t,z) $, where I gathered the parts that are only dependent on the time variable $t$.

=#

