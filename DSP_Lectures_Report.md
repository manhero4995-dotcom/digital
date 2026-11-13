# Digital Signal Processing — Simple Summary Report

**Coverage:** Lectures 1, 3, 4, 5, 6, and 6-1  
**Lecturer:** Dr. Rehab Mahmoud  
**Report date:** 13 November 2026

> **How to read this report:** Each part starts with a **simple idea**, then gives the **key formulas**, and ends with an **exam tip**. Lectures 4 and 5 repeat and extend the Z-transform material from Lecture 3, so they are summarized together.

<div class="cover-box">
<p class="cover-title">Report at a glance</p>
<ul>
<li><strong>Lecture 1 — Signals:</strong> continuous vs discrete signals, impulse, step, exponentials, and Fourier series.</li>
<li><strong>Lectures 3–5 — Z-transform:</strong> definition, ROC, poles and zeros, stability, standard pairs, and worked examples.</li>
<li><strong>Lecture 6 — Inverse Z-transform:</strong> inspection, partial fractions, and power series expansion.</li>
<li><strong>Lecture 6-1 — Properties:</strong> time shift, exponential multiplication, differentiation, time reversal, convolution, and the initial-value theorem.</li>
</ul>
<p class="cover-rule"><strong>Golden rule:</strong> a Z-transform is incomplete without its ROC.</p>
</div>

<div class="page-break"></div>

[TOC]

---

## 1. The Big Picture

Digital Signal Processing (DSP) is the study of signals that exist at discrete instants of time, and the systems that process them.

The course moves in this order:

1. **Signals** — continuous and discrete signals, basic signals, Fourier series.
2. **Z-transform** — a tool that turns a discrete-time sequence into a function of a complex variable `z`.
3. **Inverse Z-transform** — going back from `X(z)` to `x[n]`.
4. **LTI systems and stability** — using poles and the ROC to judge system behaviour.
5. **Filter design (FIR/IIR)** and the **DFT** — later topics.

The Z-transform is the centre of Lectures 3–6. Think of it as the discrete-time cousin of the Laplace transform.

---

## 2. Lecture 1 — Signals and Fourier Series

### 2.1 Continuous-time vs discrete-time signals

| | Continuous-time (C.T.) | Discrete-time (D.T.) |
|---|---|---|
| Defined | At every instant of time | Only at discrete instants |
| Symbol | `t` | `n` |
| Notation | `x(t)` — parentheses | `x[n]` — square brackets |
| Examples | Video, audio (naturally continuous) | Stock market, temperature samples |

Most real signals are continuous, but computers and digital signal processors need discrete values. So continuous signals are **sampled** to become discrete.

**Simple idea:** continuous time is a smooth film; discrete time is a sequence of photographs.

### 2.2 Basic discrete-time signals

**Unit impulse**

```
δ[n] = 1,  n = 0
       0,  n ≠ 0
```

It is one single dot at the origin. It is the simplest building block of discrete signals.

**Unit step**

```
u[n] = 1,  n ≥ 0
       0,  n < 0
```

It "switches on" at `n = 0` and stays on forever.

**Complex exponential / sinusoid**

```
x[n] = C αⁿ
```

- If `α` is real, the sequence grows or decays.
- If `α` is complex, the real and imaginary parts are sinusoids.
- If `|α| = 1`, the sequence is a pure sinusoid.
- If `|α| < 1`, it is a decaying sinusoid.
- If `|α| > 1`, it is a growing sinusoid.

Euler's relation links exponentials and sinusoids:

```
e^(jω₀t) = cos(ω₀t) + j sin(ω₀t)

cos(ω₀t) = (e^(jω₀t) + e^(−jω₀t)) / 2
sin(ω₀t) = (e^(jω₀t) − e^(−jω₀t)) / (2j)
```

**Simple idea:** a sinusoid is just two complex exponentials spinning in opposite directions.

### 2.3 Transformations of the independent variable

1. **Time shifting**
   - `x[n − n₀]` → shifted **right** (delayed).
   - `x[n + n₀]` → shifted **left** (advanced).

2. **Time reversal**
   - `x[−n]` → mirror image about `n = 0`.

3. **Time scaling**
   - `x[2n]` → plays twice as fast (in discrete time this is downsampling).
   - `x[n/2]` → plays half as fast (upsampling/interpolation).

**Example from the lecture:** for `x[n] = u[n]`, the signals `x[n−1]`, `x[n+1]`, and `x[n+3]` are the same step, just delayed or advanced.

**Exam tip:** for `x[n − n₀]`, the graph moves in the **positive** direction when `n₀` is positive. Do not reverse this.

### 2.4 Fourier series of periodic signals

**Continuous time**

```
x(t) = Σ_{k=−∞}^{∞} a_k e^(jkω₀t)

a_k = (1/T) ∫_T x(t) e^(−jkω₀t) dt
```

**Discrete time**

```
x[n] = Σ_{k=<N>} a_k e^(jk(2π/N)n)

a_k = (1/N) Σ_{n=<N>} x[n] e^(−jk(2π/N)n)
```

The notation `k = <N>` means the sum is over one period of `N` terms.

**Worked example 1:** `x(t) = sin(ω₀t)`

Using Euler:

```
sin(ω₀t) = (1/2j)e^(jω₀t) − (1/2j)e^(−jω₀t)
```

So:

```
a₁ = 1/(2j),   a₋₁ = −1/(2j),   all other a_k = 0
```

**Worked example 2:**

```
x(t) = 1 + sin(ω₀t) + 2cos(ω₀t) + cos(2ω₀t + π/4)
```

The coefficients are:

```
a₀ = 1
a₁ = 1 + 1/(2j) = 1 − 0.5j
a₋₁ = 1 − 1/(2j) = 1 + 0.5j
a₂ = (1/2)e^(jπ/4) = (√2/4)(1 + j)
a₋₂ = (1/2)e^(−jπ/4) = (√2/4)(1 − j)
```

**Simple idea:** Fourier series answers the question "how much of each sinusoid is inside this periodic signal?"

---

## 3. Lectures 3–5 — The Z-Transform

### 3.1 Definition

The **bilateral (two-sided) Z-transform** of a sequence `x[n]` is:

```
X(z) = Σ_{n=−∞}^{∞} x[n] z^(−n)
```

For causal signals that are zero for `n < 0`, the one-sided version is often used:

```
X(z) = Σ_{n=0}^{∞} x[n] z^(−n)
```

Here `z` is a complex variable:

```
z = r e^(jω)
```

Substituting:

```
X(re^(jω)) = Σ x[n] (re^(jω))^(−n)
           = Σ (x[n] r^(−n)) e^(−jωn)
```

So the Z-transform is the **DTFT of the sequence after multiplying it by `r^(−n)`**.

- If `r = 1`, then `z = e^(jω)` and the Z-transform becomes the **DTFT**.
- The circle `|z| = 1` is called the **unit circle**.

**Simple idea:** the Z-transform is a frequency lens with an extra "growth/decay" knob `r`. The unit circle is the special setting where it behaves like the Fourier transform.

### 3.2 Why the Z-transform is useful

- It solves **difference equations**.
- It evaluates the response of **LTI systems**.
- It helps design and implement **digital filters**.
- It handles signals for which the DTFT does not converge, such as `u[n]`.
- It generalizes the DTFT and is often more convenient in analytical work.

**Advantage over the DTFT:** the DTFT needs absolute summability, `Σ |x[n]| < ∞`. Many useful signals fail this test. The Z-transform can still handle them because the factor `r^(−n)` can make the sum converge.

### 3.3 Region of Convergence (ROC)

The **ROC** is the set of all values of `z` for which the Z-transform sum converges.

**Simple idea:** the ROC is the "valid region" of the transform. Outside it, the infinite sum blows up, so the formula is not valid there.

Important facts:

- The ROC **cannot contain any poles**.
- The ROC must be a **connected region**.
- For a right-sided sequence, the ROC is **outside** the outermost pole.
- For a left-sided sequence, the ROC is **inside** the innermost pole.
- For a two-sided sequence, the ROC is usually an **annulus** (a ring) between two poles.
- For a finite-length sequence, the ROC is the **whole z-plane**, except possibly `z = 0` or `z = ∞`.

The Fourier transform exists (and the system has a frequency response) **if and only if the ROC includes the unit circle**.

**Geometric series used everywhere:**

```
Σ_{n=0}^{∞} aⁿ = 1 / (1 − a),   |a| < 1
```

For the Z-transform, the same idea gives the condition `|a z^(−1)| < 1`, or `|z| > |a|`.

### 3.4 Poles and zeros

For a rational Z-transform:

```
X(z) = P(z) / Q(z)
```

- **Zeros** are the values of `z` where `X(z) = 0` (roots of the numerator).
- **Poles** are the values of `z` where `X(z) = ∞` (roots of the denominator).

On the pole-zero plot:

- a zero is drawn as `○`
- a pole is drawn as `×`

**Exam tip:** the ROC is always separated by poles. It never jumps over a pole.

### 3.5 Stability and causality

For an LTI system with impulse response `h[n]` and transfer function `H(z)`:

- The system is **stable** if the ROC of `H(z)` includes the **unit circle**.
- A **causal** system has a right-sided impulse response, so its ROC extends **outward** from the outermost pole.
- A causal system is **stable** if and only if **all its poles lie inside the unit circle**.

Why? If all poles are inside the unit circle, the ROC is `|z| > r_max` with `r_max < 1`, so the ROC includes `|z| = 1`.

**Lecture example:** suppose the poles are at `z = 1/2` and `z = 2`.

- The ROC `1/2 < |z| < 2` includes the unit circle → the system is **stable**, but it is **not causal** (two-sided).
- The ROC `|z| > 2` is causal, but it does **not** include the unit circle → the system is **unstable**.
- Therefore, this system cannot be both causal and stable.

### 3.6 Common Z-transform pairs

| Sequence `x[n]` | Transform `X(z)` | ROC |
|---|---|---|
| `δ[n]` | `1` | all `z` |
| `u[n]` | `1 / (1 − z⁻¹)` | `abs(z) > 1` |
| `−u[−n−1]` | `1 / (1 − z⁻¹)` | `abs(z) < 1` |
| `δ[n−m]` | `z⁻ᵐ` | all `z` except 0 or ∞ |
| `aⁿu[n]` | `1 / (1 − a z⁻¹)` | `abs(z) > abs(a)` |
| `−aⁿu[−n−1]` | `1 / (1 − a z⁻¹)` | `abs(z) < abs(a)` |
| `n aⁿu[n]` | `a z⁻¹ / (1 − a z⁻¹)²` | `abs(z) > abs(a)` |
| `−n aⁿu[−n−1]` | `a z⁻¹ / (1 − a z⁻¹)²` | `abs(z) < abs(a)` |
| `cos(ω₀n)u[n]` | `(1 − cosω₀ z⁻¹) / (1 − 2cosω₀ z⁻¹ + z⁻²)` | `abs(z) > 1` |
| `sin(ω₀n)u[n]` | `sinω₀ z⁻¹ / (1 − 2cosω₀ z⁻¹ + z⁻²)` | `abs(z) > 1` |
| `rⁿcos(ω₀n)u[n]` | `(1 − r cosω₀ z⁻¹) / (1 − 2r cosω₀ z⁻¹ + r²z⁻²)` | `abs(z) > r` |
| `rⁿsin(ω₀n)u[n]` | `r sinω₀ z⁻¹ / (1 − 2r cosω₀ z⁻¹ + r²z⁻²)` | `abs(z) > r` |
| `aⁿ, 0 ≤ n ≤ N−1` | `(1 − aᴺz⁻ᴺ) / (1 − a z⁻¹)` | `abs(z) > 0` |

In the table, `abs(z)` means `|z|`.

### 3.7 Worked examples

#### Example A — Right-sided exponential

Find the Z-transform of:

```
x[n] = aⁿu[n]
```

Using the definition:

```
X(z) = Σ_{n=0}^{∞} aⁿ z^(−n)
     = Σ_{n=0}^{∞} (a z⁻¹)ⁿ
     = 1 / (1 − a z⁻¹),   |a z⁻¹| < 1
```

So:

```
x[n] = aⁿu[n]  ↔  X(z) = 1 / (1 − a z⁻¹) = z / (z − a),   |z| > |a|
```

**Check:** the pole is at `z = a`; the ROC is outside that pole. Correct for a right-sided sequence.

#### Example B — Left-sided exponential

Now:

```
x[n] = −aⁿu[−n−1]
```

The same algebra gives:

```
X(z) = 1 / (1 − a z⁻¹),   |z| < |a|
```

This is the **same algebraic expression** as Example A, but the ROC is different. Therefore the inverse is different:

```
−aⁿu[−n−1]  ↔  1 / (1 − a z⁻¹),   |z| < |a|
```

**Very important:** `X(z)` alone is not enough. You must also know the ROC.

#### Example C — Sum of two exponentials

```
x[n] = (1/2)ⁿu[n] + (−1/3)ⁿu[n]
```

By linearity:

```
X(z) = 1/(1 − (1/2)z⁻¹) + 1/(1 + (1/3)z⁻¹)
```

Combining:

```
X(z) = [2 − (1/6)z⁻¹] / [(1 − (1/2)z⁻¹)(1 + (1/3)z⁻¹)]
```

The ROC is the overlap:

```
|z| > 1/2
```

because `1/2` is the larger pole magnitude.

#### Example D — Two-sided sequence

```
x[n] = (−1/3)ⁿu[n] − (1/2)ⁿu[−n−1]
```

The first term is right-sided, so its ROC is `|z| > 1/3`.  
The second term is left-sided, so its ROC is `|z| < 1/2`.

The overall ROC is the ring:

```
1/3 < |z| < 1/2
```

The Z-transform is:

```
X(z) = 1/(1 + (1/3)z⁻¹) + 1/(1 − (1/2)z⁻¹)
```

#### Example E — Finite-length signal

```
x[n] = δ[n] + δ[n−5]
```

```
X(z) = 1 + z⁻⁵
```

The ROC is the whole z-plane except `z = 0` (because of the `z⁻⁵` term):

```
|z| > 0
```

#### Example F — Finite geometric sequence

```
x[n] = aⁿ,   0 ≤ n ≤ N−1
```

```
X(z) = (1 − aᴺz⁻ᴺ) / (1 − a z⁻¹)
     = (zᴺ − aᴺ) / [z^(N−1)(z − a)]
```

The zeros satisfy `zᴺ = aᴺ`; there is a pole at `z = 0` and a pole at `z = a`. The pole at `z = a` is cancelled by one of the zeros, so the only remaining pole is at `z = 0`.

ROC: `|z| > 0`.

---

## 4. Lecture 6 — Inverse Z-Transform

The inverse problem is: given `X(z)` and its ROC, find `x[n]`.

The lecture presents three informal methods:

1. **Inspection** — recognise a standard pair.
2. **Partial fraction expansion** — split a complicated rational function into simple table entries.
3. **Power series expansion** — read the coefficients of `z⁻ⁿ` directly.

### 4.1 Inspection method

This works when `X(z)` is already in a known form.

**Example 1:**

```
X(z) = 1 / (1 − (1/2)z⁻¹),   |z| > 1/2
```

This is the standard right-sided pair:

```
x[n] = (1/2)ⁿu[n]
```

**Example 2:**

```
X(z) = 1 / (1 − (1/2)z⁻¹),   |z| < 1/2
```

Same expression, but now the ROC is inside the pole, so the sequence is left-sided:

```
x[n] = −(1/2)ⁿu[−n−1]
```

**Exam tip:** always look at the ROC before choosing the right-sided or left-sided version.

### 4.2 Partial fraction expansion

Assume `X(z)` is a ratio of polynomials in `z⁻¹`:

```
X(z) = Σ_{k=0}^{M} b_k z⁻ᵏ / Σ_{k=0}^{N} a_k z⁻ᵏ
```

The denominator has `N` poles at `d₁, d₂, …, d_N`.

#### Case 1: `M < N` and all poles are first order

```
X(z) = Σ_{k=1}^{N} A_k / (1 − d_k z⁻¹)
```

The coefficients are found by multiplying by the factor and evaluating at the pole:

```
A_k = (1 − d_k z⁻¹) X(z) |_{z = d_k}
```

**Simple idea:** each coefficient is "the residue" left after removing that pole's factor. It tells you how much of that simple exponential is present.

**Lecture example:**

```
X(z) = 1 / [(1 − (1/4)z⁻¹)(1 − (1/2)z⁻¹)],   |z| > 1/2
```

Let:

```
X(z) = A₁/(1 − (1/4)z⁻¹) + A₂/(1 − (1/2)z⁻¹)
```

Compute:

```
A₁ = [1 / (1 − (1/2)z⁻¹)] at z = 1/4 = −1
A₂ = [1 / (1 − (1/4)z⁻¹)] at z = 1/2 = 2
```

So:

```
X(z) = −1/(1 − (1/4)z⁻¹) + 2/(1 − (1/2)z⁻¹)
```

Inverse:

```
x[n] = 2(1/2)ⁿu[n] − (1/4)ⁿu[n]
```

#### Case 2: `M ≥ N`

The numerator degree is at least the denominator degree, so first perform **long division**:

```
X(z) = Σ_{r=0}^{M−N} B_r z⁻ʳ + Σ_{k=1}^{N} A_k/(1 − d_k z⁻¹)
```

The `B_r` terms come from the long division and produce impulses `δ[n]`, `δ[n−1]`, etc.

**Lecture example:**

```
X(z) = (1 + 2z⁻¹ + z⁻²) / (1 − (3/2)z⁻¹ + (1/2)z⁻²),   |z| > 1
```

Factor:

```
X(z) = (1 + z⁻¹)² / [(1 − (1/2)z⁻¹)(1 − z⁻¹)]
```

Long division gives:

```
X(z) = 2 + (−1 + 5z⁻¹) / [(1 − (1/2)z⁻¹)(1 − z⁻¹)]
```

Partial fractions:

```
A₁ = −9,   A₂ = 8
```

So:

```
X(z) = 2 − 9/(1 − (1/2)z⁻¹) + 8/(1 − z⁻¹)
```

Inverse:

```
x[n] = 2δ[n] − 9(1/2)ⁿu[n] + 8u[n]
```

#### Case 3: Repeated poles

If a pole `d_i` has order `S`, the expansion needs terms like:

```
C₁/(1 − d_i z⁻¹) + C₂/(1 − d_i z⁻¹)² + … + C_S/(1 − d_i z⁻¹)^S
```

The coefficients are obtained from derivatives of the remaining expression. The lecture gives the general derivative formula; in practice, the idea is:

- first-order poles → one coefficient each;
- repeated pole of order `S` → `S` coefficients.

**Exam tip:** if you see `(1 − d z⁻¹)²` in the denominator, you must include **two** terms, not one.

### 4.3 Power series expansion

If:

```
X(z) = Σ_{n=−∞}^{∞} x[n]z⁻ⁿ
       = … + x[−2]z² + x[−1]z + x[0] + x[1]z⁻¹ + x[2]z⁻² + …
```

then the coefficient of `z⁻ⁿ` is exactly `x[n]`.

**Example from the lecture:**

```
X(z) = z²(1 − (1/2)z⁻¹)(1 + z⁻¹)(1 − z⁻¹)
```

Expand:

```
(1 + z⁻¹)(1 − z⁻¹) = 1 − z⁻²
```

Then:

```
X(z) = z²(1 − (1/2)z⁻¹)(1 − z⁻²)
     = z² − (1/2)z − 1 + (1/2)z⁻¹
```

Reading coefficients:

```
x[−2] = 1
x[−1] = −1/2
x[0]  = −1
x[1]  = 1/2
all other x[n] = 0
```

So:

```
x[n] = δ[n+2] − (1/2)δ[n+1] − δ[n] + (1/2)δ[n−1]
```

**Simple idea:** the Z-transform is a polynomial in `z⁻¹`; the sequence values are its coefficients.

### 4.4 Which inverse method should I use?

- The expression matches a standard pair → **inspection**.
- The expression is a proper rational function → **partial fractions**.
- The numerator degree is greater than or equal to the denominator degree → **long division first**, then partial fractions.
- The expression is already a series or a product that is easy to expand → **power series**.
- There are repeated poles → **partial fractions with repeated-pole terms**.

---

## 5. Lecture 6-1 — Z-Transform Properties

These properties let you avoid starting from the definition every time. They are especially useful for inverse transforms.

### 5.1 The eight properties

| Property | Formula | ROC / meaning |
|---|---|---|
| Linearity | `a x₁[n] + b x₂[n] ↔ aX₁(z) + bX₂(z)` | ROC contains `R₁ ∩ R₂` |
| Time shifting | `x[n − n₀] ↔ z⁻ⁿ⁰X(z)` | ROC is `R_x`, except possible changes at `z = 0` or `z = ∞` |
| Multiply by exponential | `z₀ⁿx[n] ↔ X(z/z₀)` | ROC is `abs(z₀)R_x` |
| Differentiation | `n x[n] ↔ −z dX(z)/dz` | ROC is `R_x` |
| Conjugation | `x*[n] ↔ X*(z*)` | ROC is `R_x` |
| Time reversal | `x*[−n] ↔ X*(1/z*)` | ROC is `1/R_x` |
| Convolution | `x₁[n] * x₂[n] ↔ X₁(z)X₂(z)` | ROC contains `R₁ ∩ R₂` |
| Initial-value theorem | if `x[n] = 0` for `n < 0`, then `x[0] = lim_{z→∞} X(z)` | works for causal sequences |

For a **real** sequence, the time-reversal property simplifies to:

```
x[−n] ↔ X(1/z),   ROC = 1/R_x
```

### 5.2 Property examples

#### Time shifting

Given:

```
X(z) = 1/(z − 1/4),   |z| > 1/4
```

Rewrite:

```
X(z) = z⁻¹ / (1 − (1/4)z⁻¹)
```

The factor `z⁻¹` means a shift of one sample to the right:

```
x[n] = (1/4)^(n−1) u[n−1]
```

Equivalently:

```
x[n] = −4δ[n] + 4(1/4)ⁿu[n]
```

Both answers represent the same sequence.

#### Multiplication by an exponential

Starting from `u[n] ↔ 1/(1 − z⁻¹)`, one can derive:

```
rⁿcos(ω₀n)u[n]  ↔  (1 − r cosω₀ z⁻¹) / (1 − 2r cosω₀ z⁻¹ + r²z⁻²),   |z| > r
```

This is how damped sinusoids get their Z-transforms.

#### Differentiation

From `aⁿu[n] ↔ 1/(1 − a z⁻¹)`:

```
n aⁿu[n]  ↔  a z⁻¹ / (1 − a z⁻¹)²,   |z| > |a|
```

The extra factor `n` in time becomes a derivative in the z-domain.

#### Time reversal

```
a^(−n)u[−n]  ↔  1/(1 − a z),   |z| < 1/|a|
```

This is the time-reversed version of `aⁿu[n]`.

#### Convolution

If:

```
x₁[n] = aⁿu[n],   x₂[n] = u[n],   |a| < 1
```

then:

```
Y(z) = 1/[(1 − a z⁻¹)(1 − z⁻¹)],   |z| > 1
```

Partial fractions give:

```
y[n] = (1/(1 − a)) [u[n] − a^(n+1)u[n]]
```

**Simple idea:** convolution in time becomes multiplication in the z-domain. This is one of the main reasons the Z-transform is so useful.

#### Initial-value theorem

If `x[n]` is causal (zero for `n < 0`), then:

```
x[0] = lim_{z→∞} X(z)
```

For example, if `X(z) = 1/(1 − a z⁻¹)`, then as `z → ∞`, `X(z) → 1`, so `x[0] = 1`. This matches `x[n] = aⁿu[n]`.

---

## 6. Formula Sheet

### Signals

```
δ[n] = 1 at n = 0, 0 elsewhere
u[n] = 1 for n ≥ 0, 0 for n < 0
x[n] = C αⁿ
e^(jω₀t) = cos(ω₀t) + j sin(ω₀t)
```

### Fourier series

```
CT: x(t) = Σ a_k e^(jkω₀t),   a_k = (1/T)∫_T x(t)e^(−jkω₀t)dt
DT: x[n] = Σ a_k e^(jk(2π/N)n),   a_k = (1/N)Σ x[n]e^(−jk(2π/N)n)
```

### Z-transform

```
X(z) = Σ_{n=−∞}^{∞} x[n]z⁻ⁿ
z = r e^(jω)
X(z) at z = e^(jω) → DTFT
```

### Key pairs

```
δ[n] ↔ 1
u[n] ↔ 1/(1 − z⁻¹),        |z| > 1
aⁿu[n] ↔ 1/(1 − a z⁻¹),    |z| > |a|
−aⁿu[−n−1] ↔ 1/(1 − a z⁻¹), |z| < |a|
n aⁿu[n] ↔ a z⁻¹/(1 − a z⁻¹)², |z| > |a|
```

### ROC rules

```
Right-sided  → outside the outermost pole
Left-sided   → inside the innermost pole
Two-sided    → annulus between poles
Finite-length → whole plane except maybe 0 or ∞
```

### Stability and causality

```
Stable  ⇔ ROC includes the unit circle |z| = 1
Causal  ⇔ ROC extends outward from the outermost pole
Causal + stable ⇔ all poles inside the unit circle
```

### Inverse Z-transform

```
Inspection: match a known pair
Partial fractions: X(z) = Σ A_k/(1 − d_k z⁻¹)
A_k = (1 − d_k z⁻¹)X(z) at z = d_k
M ≥ N: do long division first
Power series: coefficient of z⁻ⁿ is x[n]
```

### Properties

```
Linearity:        a x₁[n] + b x₂[n] ↔ aX₁(z) + bX₂(z)
Time shift:       x[n − n₀] ↔ z⁻ⁿ⁰X(z)
Exponential:      z₀ⁿx[n] ↔ X(z/z₀)
Differentiation:  n x[n] ↔ −z dX/dz
Conjugation:      x*[n] ↔ X*(z*)
Time reversal:    x*[−n] ↔ X*(1/z*)
Convolution:      x₁[n] * x₂[n] ↔ X₁(z)X₂(z)
Initial value:    x[0] = lim_{z→∞} X(z)   (causal x[n])
```

---

## 7. Common Mistakes to Avoid

1. **Forgetting the ROC.** The same `X(z)` can correspond to a right-sided or a left-sided sequence. The ROC decides which one.
2. **Including a pole in the ROC.** The ROC never contains a pole, because `X(z)` is infinite there.
3. **Saying "ROC is all z" for a finite-length signal.** It is the whole plane **except** `z = 0` and/or `z = ∞` in many cases.
4. **Using partial fractions before long division.** If `M ≥ N`, divide first; otherwise the expansion is not proper.
5. **Missing repeated-pole terms.** A pole of order 2 needs two terms.
6. **Mixing up stability and causality.** Stable means the ROC includes the unit circle; causal means the ROC is outside the outermost pole. You need both conditions for "causal and stable".
7. **Wrong shift direction.** `x[n − n₀]` moves right for positive `n₀`.
8. **Forgetting the `z⁻¹` factor in time-shift inverse problems.** It corresponds to a shift, not a constant.
9. **Treating the Z-transform as only a frequency transform.** It is a complex-plane transform; the unit circle is only one special slice.
10. **Dropping the `u[n]` or `u[−n−1]` in the final answer.** The ROC tells you which one belongs there.

---

## 8. One-Minute Recap

- **DSP** processes discrete signals `x[n]`.
- **Basic signals:** impulse, step, exponential, sinusoid.
- **Fourier series** decomposes periodic signals into sinusoids.
- **Z-transform:** `X(z) = Σ x[n]z⁻ⁿ`; it generalizes the DTFT.
- **ROC** is the region where the transform converges; it never contains poles.
- **Poles and zeros** describe the denominator and numerator roots.
- **Stable** ⇔ ROC includes the unit circle.
- **Causal** ⇔ ROC is outside the outermost pole.
- **Causal + stable** ⇔ all poles inside the unit circle.
- **Inverse Z** methods: inspection, partial fractions, power series.
- **Properties** (linearity, shift, exponential multiplication, differentiation, reversal, convolution) make calculations much easier.
- **Golden rule:** always pair `X(z)` with its ROC.

---

*End of report.*