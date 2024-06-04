# Bounding the number of polylogarithmic Chabauty-Kim functions

Let $S$ be a finite set of primes and let $`{X = \mathbb{P}^1\setminus\{0,1,\infty\}}`$ over $`{\mathbb Z_S}`$. This SageMath code is described in the paper "<a>Polylogarithmic motivic Chabauty--Kim for $\mathbb{P}^1 \setminus \{ 0,1,\infty \}$: \\ the geometric step via resultants</a>", by David Jarossay, David T.-B. G. Lilienfeldt, Francesco M. Saettone, Ariel Weiss, and Sa'ar Zehavi.

It is used in the paper to give a provable upper bound on the dimension of the kernel of the map
```math
\theta^{\#}_{d,v} \colon \mathcal{O}(U_S)_{\leq d}[\mathrm{log}^{\mathfrak{u}}, \mathrm{Li}_1^{\mathfrak{u}}, \mathrm{Li_2}^{\mathfrak{u}}, \ldots, \mathrm{Li_d}^{\mathfrak{u}}]_v \to \mathcal{O}(U_S)_{\leq d}[\Phi, d]_v
```
It can also compute the exact dimension of this space, with high probability.

<h1>Using the code</h1>

<h3>To compute a provable upper bound (and likely lower bound) of the dimension of ${\theta^{\#}_{d, v}}$ for $|S| = n$</h3>

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

<h3>To reproduce the results of the paper</h3>

**To verify when $|S| = 2$ that the kernel of $`{\theta^{\#}_{17, 17}}`$ is trivial and the kernel of $`{\theta^{\#}_{17, 18}}`$ is at most one-dimensional** 

```
# Process took 65 hours and 450GB of memory on The Ohio State University's HPC cluster.
compute_upper_bound_on_dimension_of_kernel(depth=17, parallel=True, nprocesses=95)
```

**To verify when $|S| = 2$ that the kernel of $`{\theta^{\#}_{6, 18}}`$ is at most one-dimensional**
```
# Process takes < 10 seconds
compute_upper_bound_on_dimension_of_kernel(depth=6, parallel=False)
```
