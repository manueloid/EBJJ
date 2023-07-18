#=
# STA Wave function 

Here I will implement the STA wave function that I will use to calculate the eSTA corrections.
All the calculations have been implemented in another file and I am not going to repeat them here.
In writing this code, I will assume that the calculations are correct.

There is a lot to unpack here, so let us start from the part of the wave function that is not dependent on the position variable $z$.

### 1. Normalization factor

This is probably the easiest part of the wave function. The normalization factor is given by
$$
\begin{equation}
     \left(\frac{1}{\pi}\sqrt{\frac{8J_{0}N}{U}}\right)^{1/4} % normalisation factor 
     \tag{1}
\begin{equation}
$$
In this case then I need to extract all the parameters from the `Control` type and then calculate the normalization factor.
This is pretty straightforward so I can use an one-liner function.
=#

normalisation(c::Control) = (1/pi*sqrt(8*c.J0*c.N/c.U))^(1/4)

#=
### 2. Time dependend phase

This part is the one that is purely imaginary and it is given by
$$
\begin{equation}
	\exp\left\{-i(n + 1/2)\ \int_{0}^{t}d \tau \frac{J_{0}}{d(\tau)^{2}}\right\} % time dependent phase factor
    \tag{2}
\end{equation}
$$

I had this problem before, where it looks like the integral part of the phase factor creates some issues when calculating the integral.
=#

