#=
# Fourier checks
In this set of tests I am going to numerically check if the subsitutions I made to turn the STA wave function from momentum to position representation are correct.

My plan is to find make as much simplifications as possible, but in order to do that, I need to know if the substitutions I made are correct.

Plans:
I have already checked that the STA wf in monentum representation is correct, so I will use that as a reference.
I will hence perform a numerical Fourier transform on the STA wf in momentum representation and then I will compare it with the analytical solution that I got from the internet.

First, I will define the analytical solution of the Fourier transform of the product of a Gaussian function of the form $ e^{-p^2/2\alpha}$ and a hermite polynomial of the form $ \mathcal_n(\beta p) $

It is given by:
$$ i^n \sqrt{\alpha} \exp\left\{\frac { -\alpha x^2}{2} \right\} \gamma^n \mathcal{H}_n\left( \frac{\alpha\beta}{\gamma}x\right)$$
where
$ \gamma = \sqrt{2 \alpha \beta^2 - 1} $

I also will define the Fourier transform of a function depending on the variabl $p$ as
    $$ F^{-1}[f(p)](z) =\frac{1}{\sqrt{2\pi h}} \int_{-\infty}^{\infty} f(p) \exp\left\{\frac{ipz}{h}\right\} dp $$

For the moment I will only focus on the product of the Gaussian and the hermite polynomial with general $\alpha$ and $\beta$, and only later I will try to map the two parameters with the respective functions in the argument of theSTA wave function.

I want to compare three quantities:
- the analytical solution
- the numerical solution obtained by using the Fourier transform on the full STA wave function in monentum representation
- the numerical solution obtained by using the Fourier transform on only the spatial part of the STA wave function in monentum representation and then multiplying the time dependent one

In the following then I will assume that the STA wave function in momentum representation is of the form

=#

