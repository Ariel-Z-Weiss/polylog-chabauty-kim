# Bounding the number of polylogarithmic Chabauty-Kim functions

Let $S$ be a finite set of primes and let $`{X = \mathbb{P}^1\setminus\{0,1,\infty\}}`$ over $`{\mathbb Z_S}`$. This SageMath code is described in the paper <i><a>Polylogarithmic motivic Chabauty-Kim for $`{\mathbb{P}^1 \setminus \{ 0,1,\infty \}}`$: the geometric step via resultants</a></i>, by David Jarossay, David T.-B. G. Lilienfeldt, Francesco M. Saettone, Ariel Weiss, and Sa'ar Zehavi.

It is used in the paper to give a provable upper bound on the dimension of the kernel of the map
```math
\theta^{\#}_{d,v} \colon \mathcal{O}(U_S)_{\leq d}[\mathrm{log}^{\mathfrak{u}}, \mathrm{Li}_1^{\mathfrak{u}}, \mathrm{Li}_2^{\mathfrak{u}}, \ldots, \mathrm{Li}_d^{\mathfrak{u}}]_v \to \mathcal{O}(U_S)_{\leq d}[\Phi, d]_v
```
It can also compute the exact dimension of this space, with high probability.

<h1>Using the code</h1>

<h3>To compute a provable upper bound (and likely lower bound) of the dimension of ${\ker(\theta^{\#}_{d, v})}$ for $|S| = n$</h3>

```
T = ThetaSharpOperator(d, n)
upper_bound = T.upper_bound_on_dimension_of_kernel(d, v)
```

<h3>To compute the image of ${\theta^{\#}(\mathrm{log}^{\mathfrak u})}$ or ${\theta^{\#}(\mathrm{Li}_i^{\mathfrak u})}$ in the dual PBWL basis for $|S| = n$</h3>

```
T = ThetaSharpOperator(i, n) #or ThetaSharpOperator(d, n) for any d >= i
Li_i = T.Li(i)
log = T.log()
```

The output represents the letters $\sigma_3, \sigma_5, \sigma_7\ldots$ as $A, B, C, \ldots$ and $\tau_{p_1}, \tau_{p_2}, \tau_{p_3},\ldots$ for $`{S = \{p_1, p_2,p_3,\ldots\}}`$ as $a, b, c,\ldots$.

The outputs can be manipulated like polynomials:

```
T = ThetaSharpOperator(4, 2)

T.log() * T.Li(1) + T.Li(3)
# 3/2*Sa^2*phi0t0*phi1t0 + (Sa*Sb + Sab)*phi0t1*phi1t0 + (2*Sa*Sb - Sab)*phi0t0*phi1t1 + 3/2*Sb^2*phi0t1*phi1t1

T.log() * T.Li(3) + T.Li(1)
# 1/6*Sa^4*phi0t0^3*phi1t0 + (1/6*Sa^3*Sb + (Sa*Sab - Saab)*Sa)*phi0t0^2*phi0t1*phi1t0 + (Sa*Sabb + (Sa*Sab - Saab)*Sb)*phi0t0*phi0t1^2*phi1t0 + Sabb*Sb*phi0t1^3*phi1t0 + (1/2*(Sa^2*Sb - 2*Sa*Sab + 2*Saab)*Sa)*phi0t0^3*phi1t1 + (1/2*(Sa*Sb^2 - 2*Sabb)*Sa + 1/2*(Sa^2*Sb - 2*Sa*Sab + 2*Saab)*Sb)*phi0t0^2*phi0t1*phi1t1 + (1/6*Sa*Sb^3 + 1/2*(Sa*Sb^2 - 2*Sabb)*Sb)*phi0t0*phi0t1^2*phi1t1 + 1/6*Sb^4*phi0t1^3*phi1t1 + SA*Sa*phi0t0*phisigma3 + SA*Sb*phi0t1*phisigma3 + Sa*phi1t0 + Sb*phi1t1
````

<h3>To reproduce the results of the paper</h3>

**To verify when $|S| = 2$ that the kernel of $`{\theta^{\#}_{17, 17}}`$ is trivial and the kernel of $`{\theta^{\#}_{17, 18}}`$ is at most one-dimensional** 

```
# Process took 66 hours and 432GB of memory on The Ohio State University's HPC cluster.
compute_upper_bound_on_dimension_of_kernel(depth=17, parallel=True, nprocesses=95)
```

**To verify when $|S| = 2$ that the kernel of $`{\theta^{\#}_{6, 18}}`$ is at most one-dimensional**
```
# Process takes < 10 seconds
compute_upper_bound_on_dimension_of_kernel(depth=6, parallel=False)
```
