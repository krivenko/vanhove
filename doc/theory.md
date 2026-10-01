# The mathematics underneath

A density of states is not necessarily a smooth curve. It can have delta-peaks, cusps and
logarithmic divergences at the van Hove points, and it has a smooth background in
between. This crate makes a clear distinction between qualitatively different contributions
to a spectrum: it handles discrete resonances and singular pieces in closed form, and
calls the quadrature only for the smooth part that can be efficiently integrated.

## The spectral function splitting

A spectral function is a weighted sum of discrete resonances and continuous
contributions. The discrete part $D(\omega)$ is exact - a list of resonance
positions $\varepsilon_i$ and their respective weights $w_i$.
Each continuous contribution is split again, into a regular part $R_c(\omega)$ and
one singular part $S_{c,p}(\omega)$ for every frequency $\Omega_{c,p}$ where the
contribution stops being smooth:

$$
A(\omega) = \sum_i w_i \delta(\omega - \varepsilon_i)
    + \sum_c w_c \Big[ R_c(\omega) + \sum_p S_{c,p}(\omega) \Big].
$$

The index $c$ is dropped below wherever a single contribution is in view, leaving
$R(\omega)$, $S_p(\omega)$ and $\Omega_p$.

$R_c(\omega)$ must be **smooth between consecutive singular points $\Omega_{c,p}$**.
A singular point still must be declared as such even if $A(\omega)$ stays finite there:
a band edge where the density of states merely goes to zero with infinite slope counts,
and so does a point where two bands touch. Failing to declare a singular point does not
necessarily break the algorithms, but it can degrade the accuracy of the underlying
quadrature method severely.

<figure>
<svg viewBox="0 0 620 414" role="img" aria-label="The Lieb lattice density of states, with the flat band drawn as an arrow, split into that discrete resonance, a smooth regular part, and the three singular parts, drawn separately, at the two logarithmic van Hove points and the band touching point." style="max-width:100%;height:auto">
<defs><marker id="vh-spike" viewBox="0 0 8 8" refX="8" refY="4" markerWidth="7" markerHeight="7" orient="auto"><polygon points="0,0 8,4 0,8" fill="#E8537F"/></marker><marker id="vh-spike-a" viewBox="0 0 8 8" refX="8" refY="4" markerWidth="7" markerHeight="7" orient="auto"><polygon points="0,0 8,4 0,8" fill="currentColor"/></marker></defs>
<line x1="58" y1="78" x2="598" y2="78" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="51.0" text-anchor="end" font-family="monospace" font-size="12" fill="currentColor">A(ω)</text>
<polyline points="0,62.0 3,62.0 7,62.0 10,62.0 14,62.0 15,62.0 15,38.7 17,38.6 20,38.3 24,37.9 27,37.6 30,37.2 34,36.7 37,36.3 41,35.8 44,35.3 47,34.8 51,34.1 54,33.4 57,32.8 61,31.7 64,30.9 68,29.5 71,28.3 74,26.8 78,24.2 81,21.6 84,17.8 90,8.3 95,18.3 98,23.6 101,27.3 105,30.9 108,33.1 111,34.9 115,37.0 118,38.4 122,40.0 125,41.1 128,42.1 132,43.3 135,44.1 138,44.9 142,45.9 145,46.5 149,47.4 152,48.0 155,48.6 159,49.3 162,49.8 165,50.3 169,50.9 172,51.4 176,52.0 179,52.4 182,52.8 186,53.4 189,53.7 192,54.1 196,54.6 199,55.0 203,55.4 206,55.7 209,56.1 213,56.5 216,56.8 219,57.1 223,57.6 226,57.9 230,58.3 233,58.5 236,58.8 240,59.2 243,59.5 246,59.8 250,60.2 253,60.4 257,60.8 260,61.1 263,61.4 267,61.7 270,62.0 273,61.7 277,61.4 280,61.1 284,60.7 287,60.4 290,60.2 294,59.8 297,59.5 300,59.2 304,58.8 307,58.5 311,58.2 314,57.9 317,57.6 321,57.1 324,56.8 327,56.5 331,56.1 334,55.7 338,55.3 341,55.0 344,54.6 348,54.1 351,53.7 354,53.4 358,52.8 361,52.4 365,51.8 368,51.4 371,50.9 375,50.3 378,49.8 381,49.3 385,48.6 388,48.0 392,47.2 395,46.5 398,45.9 402,44.9 405,44.1 408,43.3 412,42.1 415,41.1 419,39.6 422,38.4 425,37.0 429,34.9 432,33.1 435,30.9 439,27.3 442,23.6 446,15.8 450,8.3 456,17.8 459,21.6 462,24.2 466,26.8 469,28.3 473,29.9 476,30.9 479,31.7 483,32.8 486,33.4 489,34.1 493,34.8 496,35.3 500,35.9 503,36.3 506,36.7 510,37.2 513,37.6 516,37.9 520,38.3 523,38.6 525,38.7 525,62.0 527,62.0 530,62.0 533,62.0 537,62.0 540,62.0" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linejoin="round" transform="translate(58,16)"/>
<line x1="328" y1="78" x2="328" y2="22" stroke="currentColor" stroke-width="2" marker-end="url(#vh-spike-a)"/>
<line x1="58" y1="170" x2="598" y2="170" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="143.0" text-anchor="end" font-family="monospace" font-size="12" fill="#E8537F">D(ω)</text>
<line x1="328" y1="170" x2="328" y2="114" stroke="#E8537F" stroke-width="2" marker-end="url(#vh-spike)"/>
<line x1="58" y1="262" x2="598" y2="262" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="235.0" text-anchor="end" font-family="monospace" font-size="12" fill="#159C8F">R(ω)</text>
<polyline points="0,62.0 3,62.0 7,62.0 10,62.0 14,62.0 15,62.0 15,62.0 17,61.8 20,61.4 24,60.9 27,60.6 30,60.4 34,60.1 37,59.9 41,59.8 44,59.8 47,59.9 51,60.1 54,60.3 57,60.6 61,61.2 64,61.7 68,62.7 71,63.6 74,64.6 78,66.4 81,68.0 84,70.1 88,73.9 90,77.3 91,79.3 95,83.4 98,85.2 101,86.5 105,87.7 108,88.2 111,88.5 115,88.7 118,88.6 122,88.4 125,88.1 128,87.7 132,87.2 135,86.7 138,86.1 142,85.3 145,84.7 149,83.8 152,83.1 155,82.4 159,81.4 162,80.7 165,80.0 169,79.0 172,78.2 176,77.2 179,76.4 182,75.7 186,74.7 189,74.0 192,73.3 196,72.3 199,71.6 203,70.7 206,70.0 209,69.4 213,68.6 216,68.0 219,67.4 223,66.7 226,66.1 230,65.5 233,65.0 236,64.6 240,64.0 243,63.7 246,63.3 250,62.9 253,62.7 257,62.4 260,62.2 263,62.1 267,62.0 270,62.0 273,62.0 277,62.1 280,62.2 284,62.5 287,62.7 290,62.9 294,63.3 297,63.7 300,64.0 304,64.6 307,65.0 311,65.6 314,66.1 317,66.7 321,67.4 324,68.0 327,68.6 331,69.4 334,70.0 338,70.9 341,71.6 344,72.3 348,73.3 351,74.0 354,74.7 358,75.7 361,76.4 365,77.4 368,78.2 371,79.0 375,80.0 378,80.7 381,81.4 385,82.4 388,83.1 392,84.0 395,84.7 398,85.3 402,86.1 405,86.7 408,87.2 412,87.7 415,88.1 419,88.5 422,88.6 425,88.7 429,88.5 432,88.2 435,87.7 439,86.5 442,85.2 446,82.6 449,79.3 450,77.3 452,73.9 456,70.1 459,68.0 462,66.4 466,64.6 469,63.6 473,62.4 476,61.7 479,61.2 483,60.6 486,60.3 489,60.1 493,59.9 496,59.8 500,59.9 503,59.9 506,60.1 510,60.4 513,60.6 516,60.9 520,61.4 523,61.8 525,62.0 525,62.0 527,62.0 530,62.0 533,62.0 537,62.0 540,62.0" fill="none" stroke="#159C8F" stroke-width="1.6" stroke-linejoin="round" transform="translate(58,200)"/>
<line x1="58" y1="354" x2="598" y2="354" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="327.0" text-anchor="end" font-family="monospace" font-size="12" fill="#B5810C">S<tspan baseline-shift = "sub">p</tspan>(ω)</text>
<polyline points="0,62.0 3,62.0 7,62.0 10,62.0 14,62.0 15,62.0 15,53.2 17,53.0 20,52.6 24,52.0 27,51.5 30,51.0 34,50.4 37,49.8 41,49.0 44,48.4 47,47.7 51,46.8 54,46.0 57,45.1 61,43.8 64,42.7 68,41.0 71,39.6 74,37.9 78,35.0 81,32.1 84,28.1 88,17.1 90,9.4 91,10.2 95,26.3 98,31.0 101,34.1 105,37.2 108,39.0 111,40.6 115,42.3 118,43.4 122,44.8 125,45.7 128,46.5 132,47.5 135,48.2 138,48.8 142,49.6 145,50.2 149,50.9 152,51.4 155,51.8 159,52.4 162,52.9 165,53.3 169,53.8 172,54.2 176,54.6 179,55.0 182,55.3 186,55.7 189,56.0 192,56.3 196,56.7 199,57.0 203,57.4 206,57.6 209,57.9 213,58.2 216,58.4 219,58.7 223,59.0 226,59.2 230,59.5 233,59.7 236,59.9 240,60.2 243,60.4 246,60.6 250,60.8 253,61.0 257,61.3 260,61.4 263,61.6 267,61.8 270,62.0 273,62.2 277,62.4 280,62.5 284,62.7 287,62.9 290,63.1 294,63.2 297,63.4 300,63.5 304,63.7 307,63.9 311,64.0 314,64.2 317,64.3 321,64.5 324,64.6 327,64.7 331,64.9 334,65.0 338,65.2 341,65.3 344,65.4 348,65.6 351,65.7 354,65.8 358,66.0 361,66.1 365,66.2 368,66.3 371,66.4 375,66.6 378,66.7 381,66.8 385,66.9 388,67.0 392,67.2 395,67.3 398,67.4 402,67.5 405,67.6 408,67.7 412,67.8 415,67.9 419,68.0 422,68.1 425,68.2 429,68.3 432,68.4 435,68.5 439,68.6 442,68.7 446,68.8 449,68.9 450,68.9 452,69.0 456,69.1 459,69.2 462,69.2 466,69.3 469,69.4 473,69.5 476,69.6 479,69.7 483,69.8 486,69.9 489,69.9 493,70.0 496,70.1 500,70.2 503,70.3 506,70.4 510,70.4 513,70.5 516,70.6 520,70.7 523,70.8 525,70.8 525,62.0 527,62.0 530,62.0 533,62.0 537,62.0 540,62.0" fill="none" stroke="#B5810C" stroke-width="1.4" stroke-linejoin="round" opacity="1" transform="translate(58,292)"/>
<polyline points="0,62.0 3,62.0 7,62.0 10,62.0 14,62.0 15,62.0 15,39.9 17,40.0 20,40.2 24,40.6 27,40.9 30,41.1 34,41.5 37,41.7 41,42.1 44,42.3 47,42.6 51,42.9 54,43.2 57,43.5 61,43.8 64,44.1 68,44.4 71,44.7 74,44.9 78,45.3 81,45.6 84,45.8 88,46.2 90,46.3 91,46.4 95,46.8 98,47.0 101,47.3 105,47.6 108,47.9 111,48.2 115,48.5 118,48.8 122,49.1 125,49.4 128,49.6 132,50.0 135,50.3 138,50.5 142,50.9 145,51.1 149,51.5 152,51.7 155,52.0 159,52.3 162,52.6 165,52.9 169,53.2 172,53.5 176,53.8 179,54.1 182,54.3 186,54.7 189,55.0 192,55.2 196,55.6 199,55.8 203,56.2 206,56.4 209,56.7 213,57.0 216,57.3 219,57.6 223,57.9 226,58.2 230,58.5 233,58.8 236,59.0 240,59.4 243,59.7 246,59.9 250,60.3 253,60.5 257,60.9 260,61.1 263,61.4 267,61.7 270,62.0 273,61.7 277,61.4 280,61.1 284,60.8 287,60.5 290,60.3 294,59.9 297,59.7 300,59.4 304,59.0 307,58.8 311,58.4 314,58.2 317,57.9 321,57.6 324,57.3 327,57.0 331,56.7 334,56.4 338,56.1 341,55.8 344,55.6 348,55.2 351,55.0 354,54.7 358,54.3 361,54.1 365,53.7 368,53.5 371,53.2 375,52.9 378,52.6 381,52.3 385,52.0 388,51.7 392,51.4 395,51.1 398,50.9 402,50.5 405,50.3 408,50.0 412,49.6 415,49.4 419,49.0 422,48.8 425,48.5 429,48.2 432,47.9 435,47.6 439,47.3 442,47.0 446,46.7 449,46.4 450,46.3 452,46.2 456,45.8 459,45.6 462,45.3 466,44.9 469,44.7 473,44.3 476,44.1 479,43.8 483,43.5 486,43.2 489,42.9 493,42.6 496,42.3 500,42.0 503,41.7 506,41.5 510,41.1 513,40.9 516,40.6 520,40.2 523,40.0 525,39.9 525,62.0 527,62.0 530,62.0 533,62.0 537,62.0 540,62.0" fill="none" stroke="#B5810C" stroke-width="1.4" stroke-dasharray="6 3" stroke-linejoin="round" opacity=".85" transform="translate(58,292)"/>
<polyline points="0,62.0 3,62.0 7,62.0 10,62.0 14,62.0 15,62.0 15,70.8 17,70.8 20,70.7 24,70.6 27,70.5 30,70.4 34,70.4 37,70.3 41,70.2 44,70.1 47,70.0 51,69.9 54,69.9 57,69.8 61,69.7 64,69.6 68,69.5 71,69.4 74,69.3 78,69.2 81,69.2 84,69.1 88,69.0 90,68.9 91,68.9 95,68.8 98,68.7 101,68.6 105,68.5 108,68.4 111,68.3 115,68.2 118,68.1 122,68.0 125,67.9 128,67.8 132,67.7 135,67.6 138,67.5 142,67.4 145,67.3 149,67.1 152,67.0 155,66.9 159,66.8 162,66.7 165,66.6 169,66.4 172,66.3 176,66.2 179,66.1 182,66.0 186,65.8 189,65.7 192,65.6 196,65.4 199,65.3 203,65.2 206,65.0 209,64.9 213,64.7 216,64.6 219,64.5 223,64.3 226,64.2 230,64.0 233,63.9 236,63.7 240,63.5 243,63.4 246,63.2 250,63.1 253,62.9 257,62.7 260,62.5 263,62.4 267,62.2 270,62.0 273,61.8 277,61.6 280,61.4 284,61.2 287,61.0 290,60.8 294,60.6 297,60.4 300,60.2 304,59.9 307,59.7 311,59.4 314,59.2 317,59.0 321,58.7 324,58.4 327,58.2 331,57.9 334,57.6 338,57.3 341,57.0 344,56.7 348,56.3 351,56.0 354,55.7 358,55.3 361,55.0 365,54.5 368,54.2 371,53.8 375,53.3 378,52.9 381,52.4 385,51.8 388,51.4 392,50.7 395,50.2 398,49.6 402,48.8 405,48.2 408,47.5 412,46.5 415,45.7 419,44.5 422,43.4 425,42.3 429,40.6 432,39.0 435,37.2 439,34.1 442,31.0 446,24.0 449,10.2 450,9.4 452,17.1 456,28.1 459,32.1 462,35.0 466,37.9 469,39.6 473,41.5 476,42.7 479,43.8 483,45.1 486,46.0 489,46.8 493,47.7 496,48.4 500,49.2 503,49.8 506,50.4 510,51.0 513,51.5 516,52.0 520,52.6 523,53.0 525,53.2 525,62.0 527,62.0 530,62.0 533,62.0 537,62.0 540,62.0" fill="none" stroke="#B5810C" stroke-width="1.4" stroke-dasharray="1.5 3" stroke-linejoin="round" opacity=".85" transform="translate(58,292)"/>
<text x="34" y="97.0" text-anchor="middle" font-size="16" fill="currentColor" opacity=".65">=</text>
<text x="34" y="189.0" text-anchor="middle" font-size="16" fill="currentColor" opacity=".65">+</text>
<text x="34" y="281.0" text-anchor="middle" font-size="16" fill="currentColor" opacity=".65">+</text>
<line x1="148.0" y1="16" x2="148.0" y2="366" stroke="#B5810C" stroke-width="1" stroke-dasharray="2 4" opacity=".5"/>
<text x="148.0" y="383" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">ϵ-2t</text>
<line x1="328.0" y1="16" x2="328.0" y2="366" stroke="#B5810C" stroke-width="1" stroke-dasharray="2 4" opacity=".5"/>
<text x="328.0" y="383" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">ϵ</text>
<line x1="508.0" y1="16" x2="508.0" y2="366" stroke="#B5810C" stroke-width="1" stroke-dasharray="2 4" opacity=".5"/>
<text x="508.0" y="383" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">ϵ+2t</text>
<text x="73.4" y="383" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">ϵ-2&#8730;2t</text>
<text x="582.6" y="383" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">ϵ+2&#8730;2t</text>
<text x="598" y="399" text-anchor="end" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">ω</text>
</svg>
<figcaption>

**Density of states of the Lieb lattice, and its split into contributions.** The flat
band gives rise to a single discrete resonance holding a third of the spectral weight,
drawn as an arrow in $A(\omega)$ and again in $D(\omega)$.
The dispersive bands contribute two logarithmic van Hove points at $\epsilon \pm 2t$ and a
kink where they touch at $\epsilon$. What is left over is the regular part $R(\omega)$,
drawn here on a scale magnified about eight times; it is smooth, and it is the only part
the quadrature method ever sees.

The three singular parts $S_p(\omega)$ are drawn separately: solid for the one at
$\Omega_0 = \epsilon - 2t$, dashed for $\Omega_1 = \epsilon$, dotted for
$\Omega_2 = \epsilon + 2t$. Each is specified in a closed form over the whole band.
Only their sum is constrained by $A(\omega)$; how it is distributed among them is a choice
the model makes.

</figcaption>
</figure>

## Asymptotic forms expressible by a singular part

Each singular part $S_p(\omega)$ is a sum of terms written as a function of the distance
$\Delta$ to the singular point $\Omega_p$. One term is characterized by an exponent $r_k$,
and two polynomials in the logarithm of that distance:

$$
S_p(\omega) = \sum_k |\Delta|^{r_k} P^{\pm}\_k(\ln|\Delta|), \qquad
\Delta = \omega - \Omega_p.
$$

The exponents are restricted to $r_k > -1$ as required by the integrability of a spectral
function: a singular point always has finite weight. One of the two polynomials $P^\pm_k$
is selected depending on the sign of $\Delta$, and the two may differ, as they must for a
band edge, a one-sided cusp, or a jump sitting underneath a divergence.

When a model constructs a [`Singularity`] structure to define a contribution $S_p(\omega)$,
it provides an energy scale $s$ and monomial coefficients of $u^{r}\ln^{m}u$,
where $u = |\Delta|/s$ is a dimensionless variable. This keeps the coefficients free of
fractional powers of energy and the logarithms dimensionless. The scale is, however, not
stored in [`Singularity`]. Instead the term is rewritten in $|\Delta|$ on construction,
which requires a rescaling and a shift: $u^r = s^{-r}|\Delta|^r$ multiplies each
polynomial by $s^{-r}$, and $\ln u = \ln|\Delta| - \ln s$ displaces its argument.

One monomial $|\Delta|^r \ln^m|\Delta|$ of a term integrates in closed form: integrating
by parts once takes one power of the logarithm off, so the $m$-th power is removed in $m$
steps:

$$
I_m = \int_0^l x^r \ln^m x\ dx =
\frac{l^{r+1}\ln^m l - m I_{m-1}}{r+1}, \qquad
I_0 = \frac{l^{r+1}}{r+1},
$$

and the precondition $r > -1$ makes every one of them converge at the singular point
$x = 0$. This closed form is crucial for the implementation of
[`SpectralFunction::integrate()`].

## `integrate()`'s algorithm

Integrating $A(\omega)$ against a test function $f(\omega)$ - a Fermi function,
a moment-generating $\omega^n$, a Green's function kernel $1/(z-\omega)$ - amounts to
reducing the discrete part to a sum, and decomposing each continuous contribution
as follows:

$$
\begin{aligned}
\int A(\omega) f(\omega)\ d\omega &=
    \sum_i w_i f(\varepsilon_i) +\\\\&+
    \sum_c w_c \Big[\int R_c(\omega) f(\omega)\ d\omega
    + \sum_p \int S_{c,p}(\omega)[f(\omega) - f(\Omega_{c,p})] d\omega
    + \sum_p f(\Omega_{c,p}) \int S_{c,p}(\omega)\ d\omega \Big].
\end{aligned}
$$

The middle term of the second line regularizes the integral.
Passing $S_p(\omega) f(\omega)$ straight to a quadrature would
render convergence of the Gauss-Kronrod method slow due to a singular point within
the integration domain. Subtracting the value of $f$ at the singular point of $S_p$
effectively eliminates the singularity from the integrand: near $\Omega_p$
the difference $f(\omega) - f(\Omega_p)$ vanishes linearly or faster,
so the product behaves as $|\Delta|^R P^\pm(\ln|\Delta|)$ with $R \geq r + 1 > 0$, which
is bounded and goes to zero.

The third term compensates the subtraction made in the second one. It is
$\int S_{c,p}(\omega)\ d\omega$ over the support of the contribution, and is computable
in closed form as shown in the previous section.

<figure>
<svg viewBox="0 0 560 150" role="img" aria-label="Left: the product of a divergent singular part with a test function runs off the top of the panel at the singular point. Right: subtracting the value of the test function at that point leaves an integrand that vanishes there instead." style="max-width:100%;height:auto">
<defs><marker id="vh-ar" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto"><polygon points="0,0 8,4 0,8" fill="currentColor"/></marker></defs>
<g transform="translate(40,22)">
<line x1="0" y1="72" x2="150" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<line x1="0" y1="0" x2="0" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<polyline points="1.8,0.0 4.3,18.7 6.7,29.5 9.2,35.7 11.7,39.7 14.1,42.7 16.6,45.0 19.1,46.8 21.6,48.3 24.0,49.5 26.5,50.6 29.0,51.5 31.4,52.3 33.9,53.1 36.4,53.7 38.9,54.3 41.3,54.9 43.8,55.3 46.3,55.8 48.7,56.2 51.2,56.6 53.7,57.0 56.1,57.3 58.6,57.6 61.1,57.9 63.6,58.2 66.0,58.4 68.5,58.7 71.0,58.9 73.4,59.1 75.9,59.3 78.4,59.5 80.8,59.7 83.3,59.9 85.8,60.1 88.2,60.3 90.7,60.4 93.2,60.6 95.7,60.7 98.1,60.9 100.6,61.0 103.1,61.1 105.5,61.3 108.0,61.4 110.5,61.5 113.0,61.6 115.4,61.7 117.9,61.8 120.4,62.0 122.8,62.1 125.3,62.2 127.8,62.2 130.2,62.3 132.7,62.4 135.2,62.5 137.7,62.6 140.1,62.7 142.6,62.8 145.1,62.8 147.5,62.9 150.0,63.0" fill="none" stroke="#B5810C" stroke-width="1.7"/>
<text x="0" y="-8" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".7">&#937;</text>
<text x="75" y="96" text-anchor="middle" font-family="monospace" font-size="11" fill="#B5810C">S<tspan baseline-shift = "sub">p</tspan>(ω) &#183; f(ω)</text>
<text x="75" y="110" text-anchor="middle" font-family="monospace" font-size="10" fill="currentColor" opacity=".6">unbounded</text>
</g>
<g transform="translate(232,58)" color="currentColor"><line x1="0" y1="0" x2="46" y2="0" stroke="currentColor" stroke-width="1.2" opacity=".6" marker-end="url(#vh-ar)"/><text x="23" y="-9" text-anchor="middle" font-family="monospace" font-size="9.5" fill="currentColor" opacity=".7">&#8722; f(&#937;)</text></g>
<g transform="translate(340,22)">
<line x1="0" y1="72" x2="150" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<line x1="0" y1="0" x2="0" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<polyline points="1.8,64.5 4.3,60.4 6.7,57.5 9.2,55.0 11.7,52.9 14.1,50.9 16.6,49.2 19.1,47.5 21.6,46.0 24.0,44.6 26.5,43.2 29.0,41.9 31.4,40.6 33.9,39.4 36.4,38.2 38.9,37.1 41.3,36.0 43.8,35.0 46.3,33.9 48.7,32.9 51.2,31.9 53.7,31.0 56.1,30.0 58.6,29.1 61.1,28.2 63.6,27.4 66.0,26.5 68.5,25.7 71.0,24.8 73.4,24.0 75.9,23.2 78.4,22.4 80.8,21.7 83.3,20.9 85.8,20.1 88.2,19.4 90.7,18.7 93.2,18.0 95.7,17.2 98.1,16.5 100.6,15.8 103.1,15.2 105.5,14.5 108.0,13.8 110.5,13.2 113.0,12.5 115.4,11.8 117.9,11.2 120.4,10.6 122.8,9.9 125.3,9.3 127.8,8.7 130.2,8.1 132.7,7.5 135.2,6.9 137.7,6.3 140.1,5.7 142.6,5.1 145.1,4.6 147.5,4.0 150.0,3.4" fill="none" stroke="#159C8F" stroke-width="1.7"/>
<circle cx="1.8" cy="64.5" r="2.6" fill="#159C8F"/>
<text x="0" y="-8" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".7">&#937;</text>
<text x="75" y="96" text-anchor="middle" font-family="monospace" font-size="11" fill="#159C8F">S<tspan baseline-shift = "sub">p</tspan>(ω) &#183; [f(ω) &#8722; f(Ω)]</text>
<text x="75" y="110" text-anchor="middle" font-family="monospace" font-size="10" fill="currentColor" opacity=".6">vanishes as |Δ|<tspan baseline-shift = "super">r+1</tspan></text>
</g>
</svg>
<figcaption>

**Subtracting the test function at the singular point.** Both panels show one term of
the integrand near $\Omega_p$, drawn for $r = -1/2$, the band-edge exponent of a linear
chain. On the left $S_p(\omega) f(\omega)$ diverges at $\Omega_p$, where Gauss-Kronrod
converges slowly. On the right $f(\omega) - f(\Omega_p)$ vanishes at least linearly at
$\Omega_p$, so the product tends to zero at least as fast as $|\Delta|^{r+1}$. The
subtracted constant is restored by the third term of the decomposition above, which is
evaluated in closed form.

</figcaption>
</figure>

## Convolution of two spectral functions

[`SpectralFunction::conv()`] performs convolution of two spectral functions,
$A \ast B = \int A(\nu) B(\omega-\nu)\ d\nu$. Writing each factor as a discrete part
plus a continuous one, the convolution has four terms:

$$
A \ast B = (D_A + C_A) \ast (D_B + C_B) = D_A \ast D_B + D_A \ast C_B + C_A \ast D_B
    + C_A \ast C_B.
$$

The $D_A \ast D_B$ term is treated exactly: two resonances convolve into a resonance at
the sum of their positions, with the product of their weights. A resonance convolved with
a continuous part displaces the whole band to the resonance's position and scales it by
the weight, also exact, since [`ContinuousSF`] requires every model to produce a displaced
copy of itself. The first three terms therefore require no quadrature.

The fourth term calls for numerical work, and is currently implemented only where $C_A$
and $C_B$ both have bounded support. It splits into four blocks of its own. Writing each
continuous part out as $R + \sum_p S_p$,

$$
C_A \ast C_B = R_A \ast R_B + \sum_p S_p \ast R_B + \sum_q R_A \ast S_q
             + \sum_{p,q} S_p \ast S_q,
$$

and these are not equally hard. The first has nothing singular in it at all. The middle
two have one divergence each, against a smooth factor, which the subtraction
described in the previous section removes. Only the last can put two divergences at one
and the same $\nu$, and that one is never passed to a quadrature: it has a closed form,
which is the subject of the subsections below.

The singular structure of the result, the positions of its singular points and the
asymptotic form at each, is *derived* algebraically from the factors, i.e. computed in
closed form and recorded as a singular part. It comes from the last three blocks: every
pair of frequencies at which the two factors stop being smooth, including the ends of
their supports, contributes a singularity.

In principle, $R_A \ast R_B$ can produce a kink going as $|\Delta|^1$ or milder due to a
collision of the support ends of the factors. This contribution to the singular structure
of the result is skipped. A kink that mild is a polynomial on either side of the point and
the regular part reproduces it exactly, while deriving it would count the same kink twice
wherever two pairs of support ends land on one frequency. For example, a convolution of
two sharp-edged boxes of equal width would produce the $|\Delta|$ kink at the centre
twice: once where the first box's upper edge meets the second box's lower edge, and again
with the two boxes exchanged.

Subtracting the singular structure from all four blocks leaves a smooth function,
stored as a Chebyshev interpolation ([`InterpolatedSF`]).

### Where the convolution stops being smooth

The convolution integral runs over $\nu$, with the first factor evaluated at $\nu$
and the second at $\omega - \nu$. It loses smoothness in $\omega$ where both
factors have a *feature* at one and the same $\nu$ - where $\nu - \Omega_1 = 0$ and
$\omega - \nu - \Omega_2 = 0$ hold together. Eliminating $\nu$,

$$
\omega = \Omega_1 + \Omega_2,
$$

so the singular points of the result are the pairwise sums. A *feature* here is a
singular point or an end of a support, where a spectral function abruptly turns to
zero. An end of a support is equivalent to a singular term with $r = 0$ and two
zero-degree polynomials $P^+ \neq P^-$ (a step function).

### Three stretches of the integration interval

This section computes a singular-against-singular contribution to the convolution,

$$
(S_p \ast S_q)(\omega) = \int S_p(\nu)\\, S_q(\omega - \nu)\ d\nu,
$$

the integral running over the frequencies where both factors are within their supports.
Near a singular point of the convolution, write $\Delta = \omega - (\Omega_1 + \Omega_2)$
and take one asymptotic term from each factor. Measuring from the first singular point,
$x = \nu - \Omega_1$, the second factor is met at $\Delta - x$, so a pair of terms
contributes

$$
\int c_1^{\pm}\\,|x|^{r_1}\ln^{m_1}\\!|x|\\;\\;
      c_2^{\pm}\\,|\Delta - x|^{r_2}\ln^{m_2}\\!|\Delta - x|\ dx.
$$

The sign of $x$ selects one coefficient from the pair $c_1^\pm$, and similarly the sign
of $\Delta - x$ selects from the pair $c_2^\pm$.
The integrand therefore changes form in two places: at
$\nu = \Omega_1$, that is $x = 0$, where the first factor changes side, and at
$\nu = \omega - \Omega_2$, that is $x = \Delta$, where the second does. Which cut comes
first follows from the sign of $\Delta$. Of the four possible sign combinations only
three can be realized at a given $\Delta$.

<figure>
<svg viewBox="0 -14 580 192" role="img" aria-label="The integration axis is cut at Ω_1 and at ω-Ω_2. For positive Δ the cuts sit in that order and the middle stretch takes both factors from above; for negative Δ the cuts swap and the middle stretch takes both from below. The two outer stretches keep their sides either way." style="max-width:100%;height:auto">
<g transform="translate(46,30)">
<text x="-30" y="4" font-family="monospace" font-size="11" fill="#B5810C">&#916;&gt;0</text>
<line x1="0" y1="0" x2="480" y2="0" stroke="currentColor" stroke-width="1" opacity=".3"/>
<line x1="170" y1="-13" x2="170" y2="13" stroke="currentColor" stroke-width="1.5"/>
<line x1="290" y1="-13" x2="290" y2="13" stroke="currentColor" stroke-width="1.5"/>
<text x="170" y="-20" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor">&#937;&#8321;</text>
<text x="290" y="-20" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor">&#969;&#8722;&#937;&#8322;</text>
<line x1="170" y1="24" x2="290" y2="24" stroke="#B5810C" stroke-width="1"/>
<line x1="170" y1="20" x2="170" y2="28" stroke="#B5810C" stroke-width="1"/>
<line x1="290" y1="20" x2="290" y2="28" stroke="#B5810C" stroke-width="1"/>
<text x="230" y="38" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">|&#916;|</text>
<text x="85" y="-30" text-anchor="middle" font-family="monospace" font-size="11" fill="currentColor" opacity=".7">c&#8321;&#8315; c&#8322;&#8314;</text>
<text x="230" y="-30" text-anchor="middle" font-family="monospace" font-size="11" fill="#B5810C">c&#8321;&#8314; c&#8322;&#8314;</text>
<text x="385" y="-30" text-anchor="middle" font-family="monospace" font-size="11" fill="currentColor" opacity=".7">c&#8321;&#8314; c&#8322;&#8315;</text>
</g>
<g transform="translate(46,128)">
<text x="-30" y="4" font-family="monospace" font-size="11" fill="#B5810C">&#916;&lt;0</text>
<line x1="0" y1="0" x2="480" y2="0" stroke="currentColor" stroke-width="1" opacity=".3"/>
<line x1="170" y1="-13" x2="170" y2="13" stroke="currentColor" stroke-width="1.5"/>
<line x1="290" y1="-13" x2="290" y2="13" stroke="currentColor" stroke-width="1.5"/>
<text x="170" y="-20" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor">&#969;&#8722;&#937;&#8322;</text>
<text x="290" y="-20" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor">&#937;&#8321;</text>
<line x1="170" y1="24" x2="290" y2="24" stroke="#B5810C" stroke-width="1"/>
<line x1="170" y1="20" x2="170" y2="28" stroke="#B5810C" stroke-width="1"/>
<line x1="290" y1="20" x2="290" y2="28" stroke="#B5810C" stroke-width="1"/>
<text x="230" y="38" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">|&#916;|</text>
<text x="85" y="-30" text-anchor="middle" font-family="monospace" font-size="11" fill="currentColor" opacity=".7">c&#8321;&#8315; c&#8322;&#8314;</text>
<text x="230" y="-30" text-anchor="middle" font-family="monospace" font-size="11" fill="#B5810C">c&#8321;&#8315; c&#8322;&#8315;</text>
<text x="385" y="-30" text-anchor="middle" font-family="monospace" font-size="11" fill="currentColor" opacity=".7">c&#8321;&#8314; c&#8322;&#8315;</text>
<text x="480" y="20" text-anchor="end" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">&#957;</text>
</g>
<line x1="276" y1="42" x2="276" y2="92" stroke="#B5810C" stroke-width="1" stroke-dasharray="3 3" opacity=".55"/>
<text x="290" y="86" font-family="monospace" font-size="9.5" fill="#B5810C">only the middle swaps</text>
</svg>
<figcaption>

**The three stretches.** The middle stretch has length $|\Delta|$ and takes both factors
from the same side, the sign of $\Delta$ deciding which. The outer two run out to where
the support ends and keep their sides for either sign of $\Delta$.

</figcaption>
</figure>

The three stretches contribute to a single asymptotic term of the result,

$$
|\Delta|^{\rho} P^{\pm}(\ln|\Delta|), \quad \rho = r_1 + r_2 + 1,
$$

which has the same shape as the terms of $S_p$ and $S_q$. The two exponents
alone fix $\rho$, so a pair of terms produces one term of the result whatever powers of
the logarithm the two hold.

A pair of powers $(m_1, m_2)$ contributes to the polynomials

$$
P^+ \mathrel{+}= c_1^+ c_2^+ X^{\mathrm{mid}}
    + c_1^- c_2^+ X^{\mathrm{lo}}
    + c_1^+ c_2^- X^{\mathrm{hi}},
$$
$$
P^- \mathrel{+}= c_1^- c_2^- X^{\mathrm{mid}}
    + c_1^+ c_2^- X^{\mathrm{lo}}
    + c_1^- c_2^+ X^{\mathrm{hi}}.
$$

Each $X$ is itself a polynomial in $\ln|\Delta|$, of degree $m_1 + m_2$, or one higher
where $\rho$ is a whole number.

One substitution reduces the integral to a dimensionless form on all three stretches.
Measuring the integration variable in units of the distance between the two singular
points, $x = \Delta v$, takes that distance out of the integral:
$$
|\Delta|^{\rho}\int |v|^{r_1}\\,|1-v|^{r_2}\\;
    \big[\ln|\Delta| + \ln|v|\big]^{m_1}\\,
    \big[\ln|\Delta| + \ln|1-v|\big]^{m_2}\ dv.
$$

The cuts now sit at $v = 0$ and $v = 1$, and the three stretches are $v < 0$,
$0 < v < 1$ and $v > 1$. Expanding the two brackets binomially leaves $|\Delta|^{\rho}$
times a polynomial in $\ln|\Delta|$ whose coefficients are integrals of
$|v|^{r_1}|1-v|^{r_2}$ against $\ln^j|v|\\,\ln^k|1-v|$. Those are beta functions
differentiated with respect to their *exponents*, $\partial_a^j \partial_b^k$, with $j$
counting the logarithm powers of one factor and $k$ those of the other. No standard library
supplies such derivatives, so the crate computes them from scratch.

The integration runs only where both factors are inside their own supports $\mathcal{S}\_1$,
$\mathcal{S}\_2$, that is where $\Omega_1 + \Delta v \in \mathcal{S}\_1$ and $\Omega_2 +
\Delta(1-v) \in \mathcal{S}\_2$.
Writing $\mathcal{R}$ for the set of $v$ meeting both conditions, the three stretches are

$$
\mathcal{R} \cap (-\infty, 0), \quad
\mathcal{R} \cap (0, 1), \quad
\mathcal{R} \cap (1, \infty).
$$

An outer stretch therefore stops at the nearer support end instead of running to infinity,
the middle one can begin after $v = 0$ or end before $v = 1$, and an empty
intersection means that a stretch does not contribute at all.

#### The middle stretch

Over $0 < v < 1$ the two factors are met from the same side of their respective
singular points. Writing $v_1, v_2 \in [0, 1]$ for the ends of the stretch,
the integral takes the form

$$
\int_{v_1}^{v_2} v^{r_1}(1-v)^{r_2} \ln^j v \\, \ln^k(1-v)\\, dv
    = \partial_a^j \partial_b^k \left[B_{v_2}(a, b) - B_{v_1}(a, b)\right],
\qquad a = r_1+1,\quad b = r_2+1.
$$

Both parameters of the beta function are positive, $r_{1,2} > -1$, so this stretch never
meets a pole. A stretch running the whole way between the points, $v_1 = 0$ and
$v_2 = 1$, results in the complete beta function $B(a,b)$. Only such a stretch
contributes to $P^{\pm}(\ln|\Delta|)$ - both its ends are singular points that collide
in the limit $\Delta\to 0$ and produce a singular point in the result.

A stretch cut short by the end of a support keeps the incomplete beta function
$B_z(a,b)$ at the cut, whose argument $x/\Delta$ moves with $\Delta$ because the support
end sits at a fixed $x$. Since

$$
|\Delta|^{\rho} B_{x/\Delta}(a, b) = \int_0^{x} s^{a-1}(\Delta-s)^{b-1}\\, ds
$$

is analytic in $\Delta$ past $\Delta = x$, such a stretch contributes only to the regular
part of the convolution.

The partial derivatives $\partial_a^j\partial_b^k B_z(a,b)$ are evaluated using the
series representation,

$$
B_z(a,b) = \sum_{n\ge0} \frac{(1-b)\_n}{n!}\\,\frac{z^{a+n}}{a+n},
$$

differentiated term by term. Its two factors depend on one parameter each, so no Leibniz
rule couples them: $\partial_a$ reaches only $z^{a+n}/(a+n)$, where $a > 0$ keeps every
denominator away from zero, and $\partial_b$ reaches only the Pochhammer symbol.
The $k$-th derivative of that symbol is computed from the recurrence
$P^{(k)}\_{n+1} = \frac{1}{n+1}[(1-b+n)\\,P^{(k)}\_n - k\\,P^{(k-1)}\_n]$ with
$P^{(0)}\_n = \frac{(1-b)\_n}{n!}$.

The series converges slowly as $z$ approaches $1$, and a second series in the reflected
argument $y = 1 - z$ is used there instead. Thanks to the reflection symmetry of the
integral representation, one has
$B_z(a, b) = B_1(b, a) - B_y(b, a) = B_{1/2}(a, b) + B_{1/2}(b, a) - B_y(b, a)$,
which gives

$$
B_z(a,b) = B_{1/2}(a,b)
    + \sum_{n\ge0}\frac{(1-a)\_n}{n!}\\,\frac{(1/2)^{\\,b+n} - y^{\\,b+n}}{b+n}.
$$

The cut is placed at $\tfrac12$ because both the series for $B_{1/2}(a,b)$ and the series
with the difference $(1/2)^{\\,b+n} - y^{\\,b+n}$ in the numerator converge quickly for
$0 \leq y \leq 1/2$.

#### The outer stretches

An outer stretch runs away from its singular point, and it stops where the
support $\mathcal{R}$ does, a finite distance $L$ away. Measured outwards from that
point, $u = -v$ below $v = 0$ and $u = v - 1$ above $v = 1$, its near end is some
$U_1 \ge 0$ and its far end $U_2 = L/\Delta$, and it gives

$$
\int_{U_1}^{U_2} (1+u)^{a-1} u^{b-1} \ln^j(1+u)\\, \ln^k u\\, du
    = \partial_a^j \partial_b^k
        \left[B^{\mathrm{out}}\_{U_2}(a, b) - B^{\mathrm{out}}\_{U_1}(a, b)\right],
$$

where

$$
B^{\mathrm{out}}\_U(a, b) = \int_0^U (1+u)^{a-1} u^{b-1}\\, du.
$$

The two outer stretches are described by the same formula with $r_1$ and $r_2$
exchanged: $a = r_2+1$, $b = r_1+1$ below the lower cut, and the reverse above the upper
one.

Both parameters are positive, each being some $r+1$ with $r > -1$, and the exponent
the pair produces enters only through the combination $\rho = a + b - 1$. At a finite $U$
the integral is therefore entire in $\rho$: the integrand is potentially unbounded only
at $u = 0$, where $b > 0$ keeps it integrable.

In order to evaluate the integral, its range is cut at $u = 1$, each part applying the
substitution that turns the integrand into a series in a variable no larger than
$\tfrac12$, whatever $U$ is. Below the cut the substitution $t = u/(1+u)$ gives
$t^{b-1}(1-t)^{-a-b}$ over $[0, w]$ with $w = \min(U/(1+U), \tfrac12)$, and the factor
$(1-t)^{-a-b}$ is expanded in powers of $t$.
Above the cut another substitution $y = 1/(1+u)$ turns the integrand into
$(1-y)^{b-1}y^{-\rho-1}$ over a range reaching no higher than
$y = \tfrac12$, and one expands $(1-y)^{b-1}$ instead.

The two logarithms the derivatives put under the integral follow each substitution.
Below the cut $\ln^j(1+u)\\,\ln^k u$ reads $(-\ln(1-t))^j(\ln t - \ln(1-t))^k$: the
$\ln(1-t)$ powers join the expansion, every logarithm of $t$ comes from the second
factor alone, and the $n$-th term of the expansion is proportional to

$$
\int_0^w t^{b+n-1}\ln^M t\\,dt, \qquad M \le k.
$$
This integral is finite because $b + n > 0$ and is known in the closed form.
Above the cut it reads $(-\ln y)^j(\ln(1-y) - \ln y)^k$, where both factors give
logarithms of $y$, and the $i$-th term of that expansion leaves

$$
\int_{1/(1+U)}^{1/2} y^{\\,i-\rho-1}\ln^M\\! y\\; dy, \qquad M \le j+k.
$$

Writing $B^{\mathrm{out}}\_U(a, b) = B(b, -\rho) - T_U(a, b)$ splits the stretch
contribution into two pieces:

$$
|\Delta|^{\rho} B^{\mathrm{out}}\_U(a, b)
  = \underbrace{|\Delta|^{\rho}\\,B(b, -\rho)}\_{\text{singular}}
  \\;-\\; \underbrace{|\Delta|^{\rho}\\,T_U(a, b)}\_{\text{analytic in }\Delta}.
$$

The first of them is the $U\to\infty$ limit of the integral for $\rho < 0$.
At $\rho \ge 0$ there is no such limit, the integrand is going as
$u^{\rho-1}$ at large $u$ and $\rho = 0$ already makes it diverge logarithmically.
$B(b, -\rho)$ must then be understood as the analytical continuation of the integral
from the $\rho<0$ region. Either way it is finite wherever $\Gamma(-\rho)$ is,
which is everywhere except for an integer $\rho \ge 0$.
The second term $|\Delta|^{\rho}\\,T_U(a, b)$ is a remainder analytic in $\Delta$.

Only a stretch with $U_1 = 0$ contributes to the asymptotic term
$|\Delta|^\rho P^\pm(\ln|\Delta|)$. There $B^{\mathrm{out}}\_{U_1}(a, b)$ vanishes, the
stretch is the split above at $U = U_2$, and the singular piece
$|\Delta|^{\rho}\\,B(b, -\rho)$ is what reaches $P^{\pm}$. With $U_1 > 0$ the complete
beta cancels between the two ends, leaving
$|\Delta|^{\rho}(T_{U_1}(a, b) - T_{U_2}(a, b))$, which contributes to the regular part
alone.

#### When the resulting exponent is a whole number

The derivation of the contribution to the singular structure from the outer stretches
breaks down at a whole-number resulting exponent $\rho$. The coefficient
of the $|\Delta|^{\rho}$ family is the value $B(b, -\rho)$, which has a simple pole
at every integer $\rho \ge 0$, unless $b$ is itself a positive integer no greater
than $\rho$.

That pole is not a defect of the formula but of the singular/analytic split it
belongs to. Let $\rho$ approach an integer $n$. Expanding $(1+u)^{a-1}$ in powers of
$1/u$ and integrating out to $U_2 = L/|\Delta|$,

$$
|\Delta|^{\rho}B^{\mathrm{out}}\_{U_2}(a, b) = B(b, -\rho)\\,|\Delta|^{\rho}
    + \sum_{i\ge0}\binom{a-1}{i}\frac{L^{\rho-i}}{\rho-i}\\,|\Delta|^{i},
$$

where two coefficients diverge at $\rho = n$: the complete beta function,
and the $i = n$ term of the sum. At $\rho = n-\delta$ their polar parts combine as

$$
\frac{\binom{a-1}{n}}{\delta}|\Delta|^{n-\delta}
    - \frac{\binom{a-1}{n}}{\delta}L^{-\delta}|\Delta|^{n}
    = \frac{\binom{a-1}{n}}{\delta}|\Delta|^{n}
        \left[|\Delta|^{-\delta} - L^{-\delta}\right]
    \\;\xrightarrow[\\,\delta\to0\\,]{}\\;
        \binom{a-1}{n}\\,|\Delta|^{n}\ln\frac{L}{|\Delta|}.
$$

$|\Delta|^{\rho}B^{\mathrm{out}}\_{U_2}(a, b)$ is continuous in $\rho$ while neither
of the two coefficients is, and their collision produces the logarithm.

The crate tests $\rho$ against a band of fixed width,
$10^{-4}$, around the nearest integer and, inside the band, uses a formula written for a
whole $\rho$. A finite width band is used because the generic formula is already
unusable near an integer $\rho$. In the small $\delta$ limit $B(b, -\rho)$ grows
as $1/\delta$, and so does $T_U(a, b)$. Their difference, the stretch itself, stays
finite, but forming it as a sum of the singular contribution and the regular part in
floating point arithmetic costs about $\log_{10}(1/|\delta|)$ digits of accuracy.

Inside the band the polynomial in $\ln|\Delta|$ is given degree $m_1 + m_2 + 1$
instead of $m_1 + m_2$.

An outer stretch cannot use $B(b, -\rho)$, so the pole is divided out of it once,
$B(b, -\rho) = F(\delta)/\delta$, where
$F(\delta) = \sum_{l\ge0} f_l\delta^l$. The cut at $u = 1$ gives $F$ in closed form:

$$
F(\delta) = 2^{-\delta}\left[\frac{(1-b)\_n}{n!}
    + \delta\sum_{i\ne n}\frac{(1-b)\_i}{i!}\frac{2^{\\,n-i}}{i-n+\delta}\right]
    + \delta\\, B^{\mathrm{out}}\_1(n{+}1{-}b,\\, b).
$$

The bracket is the expansion above the cut,
$\sum_{i\ge0}\frac{(1-b)\_i}{i!}\frac{2^{\rho-i}}{i-\rho}$, whose $i = n$
denominator is
$\delta$ itself and so loses the factor $\delta$ the other terms keep.
$\delta\\, B^{\mathrm{out}}\_1$ is the integral below the cut.
The crate tabulates $\partial_b^{j} f_l$, the derivatives of the
coefficients of $F$.

The stretch's own parameters are $a = r_2+1$ and $b = r_1+1$, while the table is indexed
by $b$ and $\delta$, the pole living in $\delta$ alone. Reading it therefore means a
change of variables, $\partial_a \mapsto -\partial_\delta$ and
$\partial_b \mapsto \partial_b - \partial_\delta$, with the $\partial_b$ on the right
taken at fixed $\delta$.

Differentiating $B = F(\delta)/\delta$ brings in negative powers of $\delta$, though none
of them is ever stored: the $1/\delta$ stays explicit outside the table. An entry
reached with the $k$-th derivative in $\delta$ runs down to
$\delta^{-(k+1)}$. Writing it as $\sum_{m\ge-(k+1)} c_m\delta^{m}$, the result
keeps the coefficient of $\delta^{0}$ in its product with
$|\Delta|^{\rho} = |\Delta|^{n}e^{-\delta\ln|\Delta|}$, the polar parts cancelling against
the remainder as above:

$$
\left[\delta^{0}\right]\left(\sum_{m\ge-(k+1)} c_m\delta^{m}
    \\;e^{-\delta\ln|\Delta|}\right)
    = \sum_{p=0}^{k+1} c_{-p}\\,\frac{(-\ln|\Delta|)^{p}}{p!}.
$$

Every negative power therefore contributes a logarithm, $c_{-p}$ giving
$\ln^{p}|\Delta|$, while $c_0$ gives none. Since $c_{-p}$ comes from $f_{k+1-p}$, the
required coefficients are $f_0$ through $f_{k+1}$, which is why the table of
$\partial_b^{j}f_l$ is the right thing to store.

Of those, $f_0$ alone is elementary, and it is the one that matters most: it supplies the
deepest pole $c_{-(k+1)}$, so it sets the top power of $\ln|\Delta|$, the power the
whole-number case adds. It is

$$
f_0 = \frac{(1-b)\_n}{n!} = \binom{n-b}{n},
$$

a polynomial in $b$, so nothing of the pole of $B(b, -\rho)$ survives into that
logarithm's coefficient. It vanishes identically wherever $b$ is a positive integer no
greater than $n$, and there the extra logarithm is absent altogether.

The coefficient of $|\Delta|^{n}$ with no logarithm beside it is settled separately. For
$n \ge 1$ it is set to zero and left to the regular part. It is analytic, so it belongs
there in any case, and since $|\Delta|^{n}$ vanishes at $\Delta = 0$, so does every other
term the pair derives. For $n = 0$ the term does not vanish, and there the stretches
keep their own $\ln^0$ entries. The constant the pair leaves at $\Omega_1 + \Omega_2$
is added to those entries where $\rho$ is exactly zero, and stands as a term of its own
otherwise.

### The regular part of the convolution

Nothing above computes $S_p \ast S_q$ by quadrature. The quadrature gets
$R_A \ast R_B$, which has no divergence at all, and the two $S \ast R$ terms,
which have one each. Those are subtracted the way
[`SpectralFunction::integrate()`] subtracts one, with the other factor evaluated at
the singular point:

$$
\int_{\mathcal{S}\_p} S_p(\nu)\\, R_B(\omega-\nu)\\,d\nu
    = \int_{\mathcal{S}\_p} S_p(\nu)
        \left[R_B(\omega-\nu) - R_B(\omega-\Omega_p)\right] d\nu
    + R_B(\omega-\Omega_p) \int_{\mathcal{S}\_p} S_p.
$$

The regular part of the result is then what remains once the derived structure is
subtracted:

$$
R(\omega) = (C_A \ast C_B)(\omega) - \sum_p S^{\mathrm{derived}}\_p(\omega).
$$

A convolution $S_p \ast S_q$ contributes to $R$ its subleading singular structure,
and anything with $\rho \geq 2$ is discarded as smooth enough.

### Is the $S_p(\omega)$ ansatz closed under convolution?

Convolving two terms of the ansatz produces terms of the same form: the exponent goes
to $r_1 + r_2 + 1$ and the logarithmic degree to at most $m_1 + m_2 + 1$. The
closed form of $S_p(\omega)$ is therefore enough not only for the models in this crate
but for anything they generate under repeated convolution.

The exponents also form an additive semigroup. Writing $\sigma = r + 1 > 0$ for
the margin a single term keeps against non-integrability, convolution simply gives
$\sigma = \sigma_1 + \sigma_2$. Every pre-defined lattice model has a half-integer
$\sigma$ - chain edges at $\tfrac12$, logarithms and band-edge constants at $1$,
semicircle and Bethe edges at $\tfrac32$, Dirac and band-touching points at $2$ - so
every convolution of lattice models lands on $\sigma \in \tfrac12\mathbb{Z}$. Since
$\sigma$ is strictly positive and additive, after finitely many convolutions every term
passes the fixed exponent ceiling $r_\text{max}=2$. Therefore, the truncated algebra
is finite.
