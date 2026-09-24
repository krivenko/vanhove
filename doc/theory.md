# The mathematics underneath

A density of states is not necessarily a smooth curve. It can have delta-peaks, cusps and
logarithmic divergences at the van Hove points, and it has a smooth remainder in
between. Everything in this crate follows from taking that split seriously:
writing the singular pieces in closed form, and handing the quadrature only what
it can actually integrate.

## The splitting

A spectral function is a weighted sum of discrete resonances and continuous
contributions. The discrete part is exact - positions and weights, nothing to
approximate. Each continuous contribution is split again, into a regular part and
one singular part for every frequency where it stops being smooth:

$$
A(\omega) = \sum_i w_i \delta(\omega - \varepsilon_i)
    + \sum_c w_c \Big[ R_c(\omega) + \sum_p S_p(\omega) \Big].
$$

The contract on $R_c$ is the load-bearing one: it must be **smooth between
consecutive singular points**. Every frequency where it stops being smooth has to
be declared as a singular point, whether or not $A(\omega)$ diverges there. A band
edge where the density of states merely goes to zero with infinite slope counts,
and so does a point where two bands touch. Falling short costs convergence rather
than correctness: the quadrature and the interpolation both slow to an algebraic
crawl, but they still converge on the right answer.

<figure>
<svg viewBox="0 0 620 400" role="img" aria-label="The Lieb lattice density of states split into a discrete resonance, a smooth regular part, and three singular parts at the two logarithmic van Hove points and the band touching point." style="max-width:100%;height:auto">
<defs><marker id="vh-spike" viewBox="0 0 8 8" refX="8" refY="4" markerWidth="7" markerHeight="7" orient="auto"><polygon points="0,0 8,4 0,8" fill="#E8537F"/></marker></defs>
<line x1="58" y1="78" x2="598" y2="78" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="51.0" text-anchor="end" font-family="monospace" font-size="12" fill="currentColor">A(w)</text>
<polyline points="0,62.0 7,62.0 13,62.0 20,62.0 20,38.3 27,37.6 34,36.8 41,35.9 47,34.8 54,33.4 61,31.8 68,29.7 74,26.7 81,21.6 88,8.3 95,17.1 101,27.6 108,33.1 115,36.9 122,39.8 128,42.2 135,44.1 142,45.8 148,47.3 155,48.6 162,49.8 169,50.9 176,51.9 182,52.9 189,53.7 196,54.6 202,55.4 209,56.1 216,56.8 223,57.5 230,58.2 236,58.9 243,59.5 250,60.1 256,60.8 263,61.4 270,62.0 277,61.4 284,60.8 290,60.1 297,59.5 304,58.9 311,58.2 317,57.5 324,56.8 331,56.1 338,55.4 344,54.6 351,53.7 358,52.9 364,51.9 371,50.9 378,49.8 385,48.6 392,47.3 398,45.8 405,44.1 412,42.2 418,39.8 425,36.9 432,33.1 439,27.6 446,17.1 452,8.3 459,21.6 466,26.7 472,29.7 479,31.8 486,33.4 493,34.8 499,35.9 506,36.8 513,37.6 520,38.3 520,62.0 526,62.0 533,62.0 540,62.0" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linejoin="round" transform="translate(58,16)"/>
<line x1="58" y1="170" x2="598" y2="170" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="143.0" text-anchor="end" font-family="monospace" font-size="12" fill="#E8537F">D(w)</text>
<line x1="328" y1="170" x2="328" y2="114" stroke="#E8537F" stroke-width="2" marker-end="url(#vh-spike)"/>
<line x1="58" y1="262" x2="598" y2="262" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="235.0" text-anchor="end" font-family="monospace" font-size="12" fill="#159C8F">R(w)</text>
<polyline points="0,8.9 7,8.9 13,8.9 20,8.9 20,7.8 27,6.4 34,5.5 41,5.1 47,5.1 54,5.8 61,7.3 68,9.7 74,13.7 81,19.7 88,29.9 95,46.7 101,53.1 108,56.0 115,56.9 122,56.5 128,55.1 135,53.4 142,50.9 148,48.3 155,45.4 162,42.5 169,39.4 176,36.5 182,33.4 189,30.3 196,27.5 202,24.8 209,22.1 216,19.7 223,17.3 230,15.3 236,13.5 243,12.0 250,10.6 256,9.7 263,9.1 270,8.9 277,9.1 284,9.7 290,10.6 297,12.0 304,13.5 311,15.3 317,17.3 324,19.7 331,22.1 338,24.8 344,27.5 351,30.3 358,33.4 364,36.5 371,39.4 378,42.5 385,45.4 392,48.3 398,50.9 405,53.4 412,55.1 418,56.5 425,56.9 432,56.0 439,53.1 446,46.7 452,29.9 459,19.7 466,13.7 472,9.7 479,7.3 486,5.8 493,5.1 499,5.1 506,5.5 513,6.4 520,7.8 520,8.9 526,8.9 533,8.9 540,8.9" fill="none" stroke="#159C8F" stroke-width="1.6" stroke-linejoin="round" transform="translate(58,200)"/>
<line x1="58" y1="354" x2="598" y2="354" stroke="currentColor" stroke-width="1" opacity=".25"/>
<text x="46" y="327.0" text-anchor="end" font-family="monospace" font-size="12" fill="#B5810C">S(w)</text>
<polyline points="0,62.0 7,62.0 13,62.0 20,62.0 20,39.5 27,38.9 34,38.2 41,37.3 47,36.3 54,35.0 61,33.4 68,31.2 74,28.0 81,22.8 88,9.4 95,16.7 101,26.3 108,31.3 115,34.9 122,37.7 128,40.0 135,42.0 142,43.8 148,45.4 155,46.8 162,48.2 169,49.4 176,50.6 182,51.7 189,52.7 196,53.7 202,54.6 209,55.5 216,56.4 223,57.2 230,58.0 236,58.7 243,59.4 250,60.1 256,60.8 263,61.4 270,62.0 277,61.4 284,60.8 290,60.1 297,59.4 304,58.7 311,58.0 317,57.2 324,56.4 331,55.5 338,54.6 344,53.7 351,52.7 358,51.7 364,50.6 371,49.4 378,48.2 385,46.8 392,45.4 398,43.8 405,42.0 412,40.0 418,37.7 425,34.9 432,31.3 439,26.3 446,16.7 452,9.4 459,22.8 466,28.0 472,31.2 479,33.4 486,35.0 493,36.3 499,37.3 506,38.2 513,38.9 520,39.5 520,62.0 526,62.0 533,62.0 540,62.0" fill="none" stroke="#B5810C" stroke-width="1.6" stroke-linejoin="round" transform="translate(58,292)"/>
<text x="34" y="97.0" text-anchor="middle" font-size="16" fill="currentColor" opacity=".65">=</text>
<text x="34" y="189.0" text-anchor="middle" font-size="16" fill="currentColor" opacity=".65">+</text>
<text x="34" y="281.0" text-anchor="middle" font-size="16" fill="currentColor" opacity=".65">+</text>
<line x1="148.0" y1="16" x2="148.0" y2="354" stroke="#B5810C" stroke-width="1" stroke-dasharray="2 4" opacity=".5"/>
<text x="148.0" y="371" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">e-2t</text>
<line x1="328.0" y1="16" x2="328.0" y2="354" stroke="#B5810C" stroke-width="1" stroke-dasharray="2 4" opacity=".5"/>
<text x="328.0" y="371" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">e</text>
<line x1="508.0" y1="16" x2="508.0" y2="354" stroke="#B5810C" stroke-width="1" stroke-dasharray="2 4" opacity=".5"/>
<text x="508.0" y="371" text-anchor="middle" font-family="monospace" font-size="10.5" fill="#B5810C">e+2t</text>
<text x="73.4" y="371" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">e-2&#8730;2t</text>
<text x="582.6" y="371" text-anchor="middle" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">e+2&#8730;2t</text>
<text x="598" y="387" text-anchor="end" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".6">w</text>
</svg>
<figcaption>

**The Lieb lattice, actually split.** The flat band is a single resonance carrying
a third of the weight; the dispersive bands contribute two logarithmic van Hove
points at $\epsilon \pm 2t$ and a kink where they touch at $\epsilon$. What is
left over - the regular part, drawn here on a scale magnified about twenty times -
is smooth, and it is the only piece the quadrature ever sees.

</figcaption>
</figure>

## What a singular part may say

Each singular part is a sum of terms, written in the distance to the singular
point measured in units of a scale $s$:

$$
S_p(\omega) = \sum_k c_k^{\pm}\, u^{r_k} \ln^{m_k} u, \qquad
u = \frac{|\omega - \Omega_p|}{s}.
$$

Three things about that form earn their keep. The exponent is restricted to
$r > -1$, which is exactly integrability: a singular point may diverge, but it
always carries finite weight. The coefficients $c^{\pm}$ may differ on the two
sides, which is what lets a band edge, a one-sided cusp, or a jump sitting
underneath a divergence be written down at all. And the scale $s$ keeps the
coefficients free of fractional powers of the bandwidth and the logarithms
dimensionless.

The sides may differ in a bare constant too, which looks at first like a jump in
$A(\omega)$ and is not one. $S_p$ is a piece of the splitting rather than the
spectral function, and what $A$ does at $\Omega_p$ is decided by every piece
there at once. A convolution really does produce such a term, and the section on
whole-number exponents is where it comes from.

The whole family integrates in closed form. Integrating by parts once takes one
logarithm off, so $m$ of them come off in $m$ steps:

$$
I_m = \frac{l^{r+1}\ln^m l - m I_{m-1}}{r+1}, \qquad
I_0 = \frac{l^{r+1}}{r+1},
$$

and $r > -1$ is what makes every one of them converge at $u = 0$. That closed form
is not a convenience. It is what the next section spends.

## What `integrate()` does

Integrating $A(\omega)$ against a test function $f$ - a Fermi function, a moment,
a Green's function kernel - the discrete part is a sum, and each continuous
contribution is taken apart:

$$
\int A f = \sum_i w_i f(\varepsilon_i) + \sum_c w_c \Big[
    \int R f
    + \sum_p \int S_p\,[f - f(\Omega_p)]
    + \sum_p f(\Omega_p) \int S_p \Big].
$$

The middle term is the trick. Handing $\int S_p f$ straight to a quadrature puts a
divergence under it, and Gauss-Kronrod converges too slowly there to notice that
it has not converged. Subtracting the value of $f$ at the singular point costs
nothing - the third term puts it back - and it changes the character of the
integrand completely. Near $\Omega_p$ the difference $f - f(\Omega_p)$ vanishes
linearly, so the product behaves as $u^{r+1}$ with $r + 1 > 0$: bounded, and going
to zero.

The third term is the one that costs nothing at all. It is $\int S_p$ over the
support, which the closed form above supplies exactly.

<figure>
<svg viewBox="0 0 560 150" role="img" aria-label="Left: the product of a divergent singular part with a test function runs off the top of the panel at the singular point. Right: subtracting the value of the test function at that point leaves an integrand that vanishes there instead." style="max-width:100%;height:auto">
<defs><marker id="vh-ar" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto"><polygon points="0,0 8,4 0,8" fill="currentColor"/></marker></defs>
<g transform="translate(40,22)">
<line x1="0" y1="72" x2="150" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<line x1="0" y1="0" x2="0" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<polyline points="1.8,0.0 4.3,18.7 6.7,29.5 9.2,35.7 11.7,39.7 14.1,42.7 16.6,45.0 19.1,46.8 21.6,48.3 24.0,49.5 26.5,50.6 29.0,51.5 31.4,52.3 33.9,53.1 36.4,53.7 38.9,54.3 41.3,54.9 43.8,55.3 46.3,55.8 48.7,56.2 51.2,56.6 53.7,57.0 56.1,57.3 58.6,57.6 61.1,57.9 63.6,58.2 66.0,58.4 68.5,58.7 71.0,58.9 73.4,59.1 75.9,59.3 78.4,59.5 80.8,59.7 83.3,59.9 85.8,60.1 88.2,60.3 90.7,60.4 93.2,60.6 95.7,60.7 98.1,60.9 100.6,61.0 103.1,61.1 105.5,61.3 108.0,61.4 110.5,61.5 113.0,61.6 115.4,61.7 117.9,61.8 120.4,62.0 122.8,62.1 125.3,62.2 127.8,62.2 130.2,62.3 132.7,62.4 135.2,62.5 137.7,62.6 140.1,62.7 142.6,62.8 145.1,62.8 147.5,62.9 150.0,63.0" fill="none" stroke="#B5810C" stroke-width="1.7"/>
<text x="0" y="-8" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".7">&#937;</text>
<text x="75" y="96" text-anchor="middle" font-family="monospace" font-size="11" fill="#B5810C">S_p(w) &#183; f(w)</text>
<text x="75" y="110" text-anchor="middle" font-family="monospace" font-size="10" fill="currentColor" opacity=".6">unbounded</text>
</g>
<g transform="translate(232,58)" color="currentColor"><line x1="0" y1="0" x2="46" y2="0" stroke="currentColor" stroke-width="1.2" opacity=".6" marker-end="url(#vh-ar)"/><text x="23" y="-9" text-anchor="middle" font-family="monospace" font-size="9.5" fill="currentColor" opacity=".7">&#8722; f(&#937;)</text></g>
<g transform="translate(340,22)">
<line x1="0" y1="72" x2="150" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<line x1="0" y1="0" x2="0" y2="72" stroke="currentColor" stroke-width="1" opacity=".3"/>
<polyline points="1.8,64.5 4.3,60.4 6.7,57.5 9.2,55.0 11.7,52.9 14.1,50.9 16.6,49.2 19.1,47.5 21.6,46.0 24.0,44.6 26.5,43.2 29.0,41.9 31.4,40.6 33.9,39.4 36.4,38.2 38.9,37.1 41.3,36.0 43.8,35.0 46.3,33.9 48.7,32.9 51.2,31.9 53.7,31.0 56.1,30.0 58.6,29.1 61.1,28.2 63.6,27.4 66.0,26.5 68.5,25.7 71.0,24.8 73.4,24.0 75.9,23.2 78.4,22.4 80.8,21.7 83.3,20.9 85.8,20.1 88.2,19.4 90.7,18.7 93.2,18.0 95.7,17.2 98.1,16.5 100.6,15.8 103.1,15.2 105.5,14.5 108.0,13.8 110.5,13.2 113.0,12.5 115.4,11.8 117.9,11.2 120.4,10.6 122.8,9.9 125.3,9.3 127.8,8.7 130.2,8.1 132.7,7.5 135.2,6.9 137.7,6.3 140.1,5.7 142.6,5.1 145.1,4.6 147.5,4.0 150.0,3.4" fill="none" stroke="#159C8F" stroke-width="1.7"/>
<circle cx="1.8" cy="64.5" r="2.6" fill="#159C8F"/>
<text x="0" y="-8" font-family="monospace" font-size="10.5" fill="currentColor" opacity=".7">&#937;</text>
<text x="75" y="96" text-anchor="middle" font-family="monospace" font-size="11" fill="#159C8F">S_p(w) &#183; [f(w) &#8722; f(&#937;)]</text>
<text x="75" y="110" text-anchor="middle" font-family="monospace" font-size="10" fill="currentColor" opacity=".6">vanishes as u^(r+1)</text>
</g>
</svg>
<figcaption>

**Why the subtraction is free.** Drawn for an inverse square root, the worst case
a linear chain offers. On the left the quadrature is asked to integrate through a
divergence; on the right it is handed something that goes to zero. The constant
taken out comes back through $\int S_p$, which is known exactly, so nothing is
approximated twice.

</figcaption>
</figure>

## Convolving two spectral functions

Writing each operand as a discrete part plus a continuous one, the convolution has
four terms, and they are not equally hard:

$$
(D_A + C_A) \ast (D_B + C_B) = D_A \ast D_B + D_A \ast C_B + C_A \ast D_B
    + C_A \ast C_B.
$$

Two resonances convolve into a resonance at the sum of their positions carrying
the product of their weights - exact, and closed under the operation. A resonance
against a band displaces the whole band to the resonance's position and scales it
by the weight - also exact, since every model knows how to shift itself. Those
three terms are bookkeeping.

The fourth is the real problem, and it is four problems. Writing each continuous
part out as $R + \sum_p S_p$,

$$
C_A \ast C_B = R_A \ast R_B + \sum_p S_p \ast R_B + \sum_q R_A \ast S_q
             + \sum_{p,q} S_p \ast S_q,
$$

and these are not equally hard either. The first has nothing singular in it at
all. The middle two carry one divergence each, against a factor that has none, and
a single subtraction disposes of it - unconditionally, because $R$ is at worst
Lipschitz, so anchoring it at the singular point always leaves a bounded
integrand. Only the last can put two divergences at one and the same $\nu$, and
that one is never handed to a quadrature: it has a closed form, which is the
subject of most of what follows.

What the quadrature does get is then interpolated, and the interpolation is only
as good as the singular structure it is told to expand around. That structure is
*derived*, not supplied.

### Where it stops being smooth

The convolution integral runs over $\nu$, with the first factor evaluated at $\nu$
and the second at $\omega - \nu$. It loses smoothness in $\omega$ where *both*
factors are singular at one and the same $\nu$ - where $\nu - \Omega_1 = 0$ and
$\omega - \nu - \Omega_2 = 0$ hold together. Eliminating $\nu$,

$$
\omega = \Omega_1 + \Omega_2,
$$

so the singular points of the result are the pairwise sums. A *feature* here is a
singular point **or an end of a support**: where a factor stops contributing at
all is as much a loss of smoothness as where it diverges. Written about a feature,
a factor is a handful of singular terms plus the value everything else takes
there, and a band edge is that constant on one side and nothing on the other. So a
band edge pairs with a singular point the same way two singular points pair with
each other, and the band edges of the result stop being a special case.

Every such pairing is derived but one - a constant against a constant, which is
$R \ast R$ and is left alone for reasons taken up at the end.

## Three stretches

Near such a point, write $\Delta = \omega - \Omega_1 - \Omega_2$ and take one local
term from each factor. The integration axis is cut in two places: at
$\nu = \Omega_1$, where the first factor changes side, and at
$\nu = \omega - \Omega_2$, where the second does. Which cut comes first follows the
sign of $\Delta$. Of the four possible sign combinations only three are ever live.

<figure>
<svg viewBox="0 0 580 214" role="img" aria-label="The integration axis is cut at Omega-one and at omega minus Omega-two. For positive Delta the cuts sit in that order and the middle stretch takes both factors from above; for negative Delta the cuts swap and the middle stretch takes both from below. The two outer stretches keep their sides either way." style="max-width:100%;height:auto">
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
<text x="290" y="72" font-family="monospace" font-size="9.5" fill="#B5810C">only the middle swaps</text>
</svg>
<figcaption>

**The three stretches.** The middle one has length $|\Delta|$ and takes both
factors from the same side, swapping with the sign of $\Delta$; the outer two run
out to where the support ends and keep their sides throughout. Each is a Beta
function - the middle a complete one, the outer two incomplete, cut at the length
they actually have.

</figcaption>
</figure>

Adding the three, the coefficient on the upper side is

$$
C^+_k = c_1^+ c_2^+ X^{\mathrm{mid}}_k
      + c_1^- c_2^+ X^{\mathrm{lo}}_k
      + c_1^+ c_2^- X^{\mathrm{hi}}_k,
$$

with $C^-_k$ the same with every sign flipped, and the term itself
$C^{\pm}_k |\Delta|^{\rho} \ln^k|\Delta|$, where $\rho = r_1 + r_2 + 1$. Scaling the length
of a stretch out of the integral turns every logarithm into $\ln|\Delta|$ plus the
logarithm of something of order one, and expanding those binomials leaves a
*polynomial* in $\ln|\Delta|$ of degree $m_1 + m_2$. Its coefficients are Beta
integrals,

$$
\int_0^1 t^{r_1}(1-t)^{r_2} \ln^a t \, \ln^b(1-t)\, dt
    = \partial_\alpha^a \partial_\beta^b B(\alpha, \beta)
$$

over the middle stretch, and

$$
\int_0^\infty u^{r_1}(1+u)^{r_2} \ln^a u \, \ln^b(1+u)\, du
    = (\partial_\alpha - \partial_\gamma)^a (-\partial_\gamma)^b B(\alpha, \gamma)
$$

over an outer one, at $\alpha = r_1 + 1$, $\beta = r_2 + 1$ and $\gamma = -\rho$.

Written that way the outer integral runs to infinity, but the stretch it stands
for does not: it stops where the support does, a finite distance $L$ from the
singular point. Keeping that limit is the whole of what follows. Substituting
$s = \Delta u$ and then $u = t/(1-t)$ carries the outer stretch to an **incomplete**
Beta function,

$$
\int_0^{L} s^{r_1}(\Delta+s)^{r_2}\,ds = \Delta^{\rho}\,B_z(a, b),
\qquad a = r_1+1,\quad b = -\rho,\quad z = \frac{L}{\Delta+L},
$$

and the logarithms come out of the same two substitutions:
$\ln s = \ln\Delta + \ln t - \ln(1-t)$ and $\ln(\Delta+s) = \ln\Delta - \ln(1-t)$,
so expanding both binomially leaves nothing but entries of the derivative table
$\partial_a^j\partial_b^k B_z(a,b)$. The middle stretch is the same story with
$x = \Delta t$ alone, and keeps both parameters positive.

Reflecting, $B_z(a,b) = B(a,b) - B_{1-z}(b,a)$ with $1-z = \Delta/(\Delta+L)$, splits
the stretch into exactly two pieces:

$$
\Delta^{\rho} B_z(a,b) = \underbrace{\Delta^{\rho} B(a,b)}_{\text{singular}}
  \;-\; \underbrace{(\Delta+L)^{\rho} \sum_{n\ge0} \frac{(1-a)_n}{n!\,(b+n)}
        \left(\frac{\Delta}{\Delta+L}\right)^{n}}_{\text{analytic in }\Delta}
$$

- the prefactor collapsing because $\Delta^{\rho}(1-z)^{b} = (\Delta+L)^{\rho}$. The
first term is what an infinite limit would have given on its own. The second is
what taking it to infinity throws away, and it is analytic, which is why throwing
it away is legitimate *for the asymptotics*. It is not legitimate for the value,
and it is not thrown away here.

### Why not through $\ln B$

The obvious way to differentiate a Beta function is through its logarithm, which
turns derivatives into polygammas. It fails here. The outer stretches land on
$\alpha + \gamma = -r_2$, a pole of $\Gamma$ whenever the other exponent is a whole
number - and a constant or a bare logarithm is exactly $r_2 = 0$, the commonest
case there is. There $B$ vanishes while its derivatives do not, so through $\ln B$
the two meet as $0 \cdot \infty$. Taking the derivatives through

$$
B(\alpha, \gamma) = \Gamma(\alpha)\,\Gamma(\gamma)\cdot\frac{1}{\Gamma(\alpha+\gamma)}
$$

instead makes the pole a zero of an entire function, and everything stays finite by
construction.

## Computing the incomplete Beta

Nothing in the Rust ecosystem offers $\partial_a^j \partial_b^k B_z(a,b)$, and
nothing offers $B_z(a,b)$ for $b \leq 0$ at all - the usual routines are written
for the regularized $I_x(a,b)$ with both parameters positive. It is written from
scratch here, out of the series

$$
B_z(a,b) = \sum_{n\ge0} \frac{(1-b)_n}{n!}\,\frac{z^{a+n}}{a+n},
$$

differentiated term by term. Its two factors depend on one parameter each, so
there is no Leibniz rule to apply between them: $\partial_a$ reaches only
$z^{a+n}/(a+n)$, where $a > 0$ keeps every denominator away from zero, and
$\partial_b$ reaches only the Pochhammer.

The Pochhammer is carried by the recurrence
$C^{(k)}_{n+1} = [(1-b+n)\,C^{(k)}_n - k\,C^{(k-1)}_n]/(n+1)$ on
$C_n = (1-b)_n/n!$, rather than through $\ln (1-b)_n$ and differences of
polygammas. That divides the factorial out before it can overflow, and it stays
right where a factor vanishes: $b$ a positive integer terminates the series
without terminating its derivatives.

The series converges for $z < 1$ and any real $b$, but only geometrically, so past
$z = \tfrac12$ the reflection $B_z(a,b) = B(a,b) - B_{1-z}(b,a)$ turns it into the
same series in $1-z$ with the parameters exchanged. That is the branch a
convolution spends its time on: $1 - z = \Delta/(\Delta+L)$, so the argument
approaches one exactly as the frequency approaches a singular point, and the
series converges fastest where the accuracy matters most.

## When the exponent is a whole number

The reflection has a catch, and it is the interesting case rather than an awkward
one. $b = -\rho$, so $b$ is a non-positive integer exactly when $\rho$ is a whole
number $n$ - and there both halves are infinite. $\Gamma(b)$ has a pole, and so
does the one term of the series whose denominator $b + n$ vanishes.

They cancel, but not symmetrically, and what survives is the point. Write
$\rho = n - \epsilon$. The reflected half carries $\Delta^{\rho}(1-z)^{b+n}$, which
is $(\Delta+L)^{\rho}$ up to $e^{-\epsilon\ln(\Delta+L)}$ - no $\ln\Delta$ anywhere
in it. The first half carries $\Delta^{\rho} = \Delta^{n}e^{-\epsilon\ln\Delta}$
against a $1/\epsilon$, and that does leave one:

$$
\Delta^{n}\left[\frac{P}{\epsilon} - P\ln\Delta + \ldots\right].
$$

So the poles cancel between the halves while the logarithm comes from $\Delta^{\rho}
B(a,b)$ alone. **A convolution gains a logarithm exactly where the generic formula
loses one to a pole of $\Gamma(-\rho)$**, with

$$
X^{\mathrm{lo}}_1 = -\binom{r_2}{n}, \qquad X^{\mathrm{hi}}_1 = -\binom{r_1}{n}.
$$

This is where the physics is. Two inverse square roots - the band edges of a
linear chain - meet at $\rho = 0$, and the logarithm that comes out is the van Hove
singularity at the centre of the square lattice band. The library is never told
this; it falls out of the two binomials.

And it explains an absence. Two locally constant factors are $r_1 = r_2 = 0$,
$n = 1$, where both binomials ask for more than they have and vanish. Two smooth
constants convolve into no logarithm at all - which is why a van Hove saddle in
three dimensions is a square-root cusp and not a peak.

Beside the logarithm sits a term in $\Delta^n$ with no logarithm on it, and it is
the one term in this whole business the generic formula cannot supply. Its value
depends on how far each stretch actually runs - on $\ln L$, which an infinite
limit has thrown away - so it is the pole's own finite part, and what survives
against it is

$$
\int_0^L s^{-1}\ln^M\\!s\\,ds \;\longrightarrow\; \frac{\ln^{M+1}L}{M+1}.
$$

For $n \geq 1$ the term vanishes at the point along with everything else here, and
it goes to the regular part. For $n = 0$ it does not vanish: it *is* the constant
the pair leaves behind, and it is derived.

That constant is where the two sides part company. The two outer stretches swap
their lengths along with their coefficients and so come to the same number either
way, but the middle stretch does not: it takes both factors from above below the
point and both from below above it, and at $n = 0$ it no longer has zero length to
hide behind. So

$$
C^{\pm} = c_1^{\pm}c_2^{\pm} X^{\mathrm{mid}}_0 + (\text{the same both ways}),
$$

and the constant an [`AsymptTerm`](crate::singularity::AsymptTerm) may carry per
side is there for exactly this.

## What is left over

Nothing above computes $S_p \ast S_q$ by quadrature. What the quadrature does get
is $R_A \ast R_B$, which has no divergence at all, and the two $S \ast R$ terms,
which have one each. Those are subtracted the way
[`SpectralFunction::integrate()`] subtracts one, with the other factor anchored at
the singular point:

$$
\int_{\mathcal{S}_p} S_p(\nu) B(\omega-\nu)\,d\nu
    = \int_{\mathcal{S}_p} S_p(\nu)\left[B(\omega-\nu) - B(\omega-\Omega_p)\right] d\nu
    + B(\omega-\Omega_p) \int_{\mathcal{S}_p} S_p,
$$

the first integrand bounded because the bracket vanishes at least linearly at
$\Omega_p$ while $|\nu-\Omega_p|^{r}$ blows up more slowly than a pole, and the
second in closed form. Handing over the bare product instead would leave an
inverse square root sitting on the end of a panel, where Gauss-Kronrod converges
too slowly to notice that it is not converging.

The regular part of the result is then what remains once the derived structure is
taken away:

$$
R(\omega) = (C_A \ast C_B)(\omega) - \sum_p S^{\mathrm{derived}}_p(\omega).
$$

What lands in $R$ from a pair is therefore its subleading structure and nothing
else: only the leading local form is derived, and anything with $\rho \geq 2$ is
discarded as smooth enough. The constant is not among the leftovers. Every pair of
singular parts derives the one it leaves at $\Omega_p + \Omega_q$, so its share of
$R$ goes to *zero* there rather than to a step the two panels would straddle. What
is left at the point is $R_A \ast R_B$ and the two $S \ast R$ terms, which are
continuous across it anyway.

Those leftovers are visible. A linear chain against itself is the square lattice,
so its $R$ is known exactly, and near the band centre it behaves like this:

| $d$ | $R(+d)$ | $R(-d)$ | $\dfrac{R(d) - R(d/10)}{d^2\ln d - (d/10)^2\ln(d/10)}$ |
| --- | --- | --- | --- |
| $10^{-2}$ | $-0.21614760404361$ | $-0.21614760404361$ | $-0.001098$ |
| $10^{-3}$ | $-0.21614810201956$ | $-0.21614810201956$ | $-0.000995$ |
| $10^{-4}$ | $-0.21614810880399$ | $-0.21614810880399$ | $-0.000966$ |

Two things to read off it. $R$ tends to a finite limit rather than to zero, all of
which comes from the terms a quadrature computes: the $S \ast S$ pairs colliding at
the band centre have already taken their own share away. And the approach to that
limit goes as $d^2\ln d$ - a ratio steady to two digits over three decades - which is
exactly the subleading term the exponent ceiling discards, $C^1$ with a second
derivative that diverges logarithmically, and mild enough that the fit still reaches
$8 \times 10^{-13}$.

The last column differences two samples rather than measuring each against $R(0^+)$.
The limit is not known independently, and estimating it costs more accuracy than the
departure from it has to spare: at $d = 10^{-4}$ that departure is a part in $10^{10}$
of $R$ itself, where differencing cancels the unknown and leaves the law standing.

Wherever a genuine divergence sits, $R$ is still not a number *at* the point: the
convolution and the term subtracted from it are both infinite there, and only the
limit exists. That is why the Chebyshev nodes lie strictly inside each panel. A
panel boundary sits on every singular point, and the fit never asks what $R$ is at
one - it asks either side, and the derived constant is what makes the two answers
agree.

## Does it close?

Convolving two terms of the family produces terms of the family: the exponent goes
to $r_1 + r_2 + 1$ and the logarithmic degree to at most $m_1 + m_2 + 1$. So the
closed form above is not just enough for the models in this crate - it is enough
for anything they generate under repeated convolution.

Better, the exponents form an additive semigroup. Writing $\sigma = r + 1 > 0$ for
the margin a single term keeps against non-integrability, convolution is simply
$\sigma = \sigma_1 + \sigma_2$. Every seed in the lattice models has $\sigma$ a
half-integer - chain edges at $\tfrac12$, logarithms and band-edge constants at $1$,
semicircle and Bethe edges at $\tfrac32$, Dirac and band-touching points at $2$ - so
every convolution of lattice models lands on $r \in \tfrac12\mathbb{Z}$. And since
$\sigma$ is strictly positive and additive,
after finitely many convolutions every term passes any fixed exponent ceiling.
Discarding the smooth ones is not a heuristic; the truncated algebra really is
finite.

### What comes out

Two checks where the answer is known independently. A linear chain against itself
is the square lattice, since the dispersion is a sum of two independent
one-dimensional bands; a chain against a square lattice is the simple cubic band,
whose van Hove coefficients are known in closed form.

| Quantity | Known | Derived agrees to |
| -------- | ----- | ----------------- |
| square lattice peak, chain $\ast$ chain | $-1/2\pi^2 t$ | the last bit |
| chain $\ast$ chain across the band, against [`models::square()`] | - | $8 \times 10^{-12}$ |
| simple cubic band edge | $1/4\pi^2 t^{3/2}$ | $5 \times 10^{-8}$ |
| simple cubic saddle | $-3/4\pi^2 t^{3/2}$ | $2 \times 10^{-8}$ |
| simple cubic $\mu_2, \mu_4, \mu_6$ | $6t^2,\ 90t^4,\ 1860t^6$ | $2 \times 10^{-11}$ |

The peak is exact because the derived term *is* $c\ln|\Delta|$, so reading $c$ off two
samples returns it unchanged. The two three-dimensional points are read off a square
root instead, and there the eight digits are the reading rather than the coefficient:
a square root sitting on a constant takes two samples to separate, and that is what
costs them.

None of these coefficients is written anywhere in the library.
[`models::simple_cubic()`] is one line - it convolves a chain with a square
lattice. Everything else is a consequence of the algebra above.

### The one pairing left alone

A constant against a constant is not derived. That product is $R_A \ast R_B$,
which leaves a kink going as $|\Delta|$ or milder - and $|\Delta|$ is a polynomial
either side of the panel boundary it sits on, which is where every singular point
of the result puts one. The interpolation fits it exactly. There is nothing to
gain.

There is something to lose. Two band edges of equal width put two pairs at the
same frequency, one for each way of matching them up, and the local form at a band
edge is not local: it is the value everything else takes there, which for a box is
the whole box. Each pair would then compute the entire convolution rather than a
piece of it, and the two added together would claim the kink twice. No choice of
coefficient repairs that - the ambiguity is in the pairing, not the arithmetic.

Excluding the product removes the case rather than patching it. Two sharp-edged
boxes have no singular parts at all, so no pair forms, nothing is derived at any
of the three frequencies, and their triangle comes out of the regular part exactly
- the fit reaching $3 \times 10^{-15}$, which is rounding rather than a limit.
Anything with a singular part on one side or the other is derived as usual, which
is how the simple cubic gets a coefficient at a band edge where the square lattice
merely stops.
