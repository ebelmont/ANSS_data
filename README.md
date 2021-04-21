# Data related to the Adams-Novikov spectral sequence at p=3

This repository contains miscellaneous documentation, charts, and data related
to [Guozhen Wang's `MinimalResolution` program](https://github.com/pouiyter/MinimalResolution),
as [modified by me](https://github.com/ebelmont/MinimalResolution) to run at p = 3.
The program computes the algebraic Novikov spectral sequence converging to the
E<sub>2</sub> page of the Adams-Novikov spectral sequence for the sphere,
including multiplicative information.

<b>Key for reading the chart `anss_E2_158.pdf`:</b> This depicts the E<sub>2</sub> page
of the Adams-Novikov spectral sequence for the sphere at p=3. Blue dots denote
β<sub>1</sub>-divisible classes, brown lines denote α<sub>1</sub> multiplication, dashed
gray lines denote <α<sub>1</sub>, α<sub>1</sub>, -> brackets, and concentric circles
indicate 3-divisibility.

Higher differentials are complete through the 108 stem, and are
inferred from the homotopy groups computed in the green book
(Table A3.4). When a differential has a source or target of dimension > 1,
in general we do not claim to specify which linear combination
is involved. (Similarly, we do not attempt to indicate whether a differential
hits a generator, or 2 times the generator.)

The differentials beyond the 108 stem are just the easy differentials
propagated via β<sub>1</sub>-multiplications. Propagating other differentials,
as well as determining new families of differentials, is work in progress.

<b>Data files:</b>
The files in `data/` came from running the code
[here](https://github.com/ebelmont/MinimalResolution/tree/e8a54826088465f7206068d14af819323653e8ce) as follows:

```
./mr_st 185 90
./BPtab 185
./mr_BP 185 40
```

This computes the ANSS E<sub>2</sub> page through internal degree 185.

<b>Interpreting the data files:</b>
Information about multiplication by p and α<sub>1</sub> can be read fairly
straightforwardly from the files `185_BPAANSS_a0.txt` and `185_BPAANSS_h0.txt`,
respectively. Names like `v0^2[2-5]` denote specific generators in a minimal
resolution; they are illustrated in the chart in this repository.

It is less obvious how to interpret the information about beta multiplication;
the purpose of the rest of this readme is to document how to do that, and the
process is formalized in the python code in `beta.py`. If you just want to look
up specific multiplications, read the next section for a brief explanation of
how to run that code.

If you have issues running the code, or have other questions/ comments related
to this, I'd love to hear from you-- my email address is my github username
at northwestern.edu.


### How to use `beta.py` to view beta (and other) multiplications
I'll assume you've cloned this repository and started a python interpreter in
the top level of this repository.

```
>>> from beta import *
>>> data = ANSSData() # parses files in data/ (may be slow)
>>> data.in_deg(3,1) # [1-0] is the only generator in stem 3, ANSS filtration 1
['[1-0]']
>>> data.beta1("[1-0]") # multiply [1-0] by beta1
2[3-0]
>>> data.alpha1("[1-0]")
0
>>> data.alpha1(data.beta2("[2-2]")) # multiplications can be composed
2[5-4]
```

Available multiplications are `three`, `alpha1`, `beta1`, `beta2`, `beta33`,
`beta4`, `beta5`, and `beta63`.

To parse files from your own run of Guozhen's program, use
`ANSSData(file_prefix)` where `file_prefix` is a string such that
e.g. `file_prefix + "_BPBocSS_table.txt"` is a
valid (relative) path to your output files.


### Notes on specific files

For the rest of this readme I will be talking about the files in the `data/`
directory.

* `185_BPAANSS_table.txt` describes the algebraic Novikov spectral sequence
  (algNSS) converging to the `E_2` page of the Adams-Novikov spectral sequence
  (ANSS). Lines containing a single element name indicate that the element is a
  permanent cycle. Lines containing `y <- x | dr` indicate an algebraic Novikov
  differential `d_r(x) = y`. The degree pair deg=(a,b) means stem = a and Adams
  filtration = b. **Warning: the grading convention in this file has been changed
  from Guozhen's original!** At odd primes, the `E_2` page of the algebraic Novikov
  spectral sequence coincides with the `E_2` page of the Adams spectral sequence,
  so by taking all elements in this table, including those involved in
  differentials, one obtains a list of elements in the Adams `E_2` page.

* `185_BPBocSS_table.txt` describes the Bockstein spectral sequence from the ANSS
  `E_2` page for the mod 3 Moore space, to the ANSS `E_2` page for the sphere.
  These are different generators from the algNSS table. The single degree
  listed is Adams-Novikov internal degree (i.e., stem + ANSS filtration).

* `185_BPB2A_table.txt` provides a translation between these two sets of
  generators: a line `x -> y` indicates that `x` is the Bockstein name for the
  image of the algebraic Novikov element `y` along the natural map `q: E_2(S)
  --> E_2(S/p)`. In other words, the inverse of this mapping specifies `q`.


### Using the data files to determine beta multiplications

The files `185_BPBocSS_theta*.txt` describe multiplication on `E_2(S/3)` by an
element a<sub>i/j</sub> whose image under the boundary map `delta: E_2(S/3) -->
E_2(S)` is β<sub>i/j</sub>. This repository contains information about six
different multiplication types:

* theta2 --> β<sub>1</sub>
* theta3 --> β<sub>2</sub>
* theta4 --> β<sub>3/3</sub>
* theta5 --> β<sub>4</sub>
* theta6 --> β<sub>5</sub>
* theta7 --> β<sub>6/3</sub>

It should be easy to change the `mult_theta` function in `BP_init.cpp` and the
`thetas` function in `BP.cpp` to generate other beta's.

Thus the full procedure for determining β<sub>i/j</sub> * x for `x` in `E_2(S)`
is:

1. Apply `q: E_2(S) --> E_2(S/3)` to `x` by applying the inverse of the mapping
   described in `185_BPB2A_table.txt`.
2. Multiply by a<sub>i/j</sub> using the appropriate `185_BPBocSS_theta*.txt`
   file.
3. Apply the boundary map `delta: E_2(S/3) --> E_2(S)` as described in the "Top
   cell classes and the boundary map" section below.


### Top cell classes and the boundary map

The boundary map `δ: E_2(S/p) --> E_2(S)` arises from the long exact sequence
in ANSS `E_2` pages associated to the cofiber sequence `S --> S --> S/p`. 
Differentials in the Bockstein spectral sequence are related to a filtered
version of this map. More precisely, the Bockstein spectral sequence comes from
the filtration `E_2(S) <-- 3 E_2(S) <-- 9 E_2(S) <--...`. A Bockstein
differential `d_r(x) = 3^r y` means that `3 δ(x) = 3^r y` mod 3<sup>r+1</sup>.

The line `y <- x | dr` in the file `185_BPBocSS_table.txt` indicates `d_r(x) =
y`, and hence `δ(x) = y/3`. A permanent cycle in the Bockstein spectral
sequence (indicated in the file `185_BPBocSS_table.txt` by a line containing a
single element and a degree) is sent to zero under δ, and thus is in the image
of the natural map `q: E_2(S) --> E_2(S/p)` (these are the "bottom cell classes").
Bockstein d0's can be ignored as extraneous information about the computation;
if we are told that `d_0(x) = y`, we may treat `y = 0`, as that is true in the
Bockstein `E_1` page.

In order to specify δ in a more useful way, we need to convert `y/3` (specified
in terms of Bockstein generators) back to algNSS notation. The file
`185_BPBocSS_a0.txt` gives the multiplication-by-3 table for Bockstein
generators. Thus to calculate the algNSS name for `y/3`, we invert that mapping
and then apply the name translation in `185_BPB2A_table.txt`. (In more detail:
the `B2A` conversion is only defined for non-3-divisible names, so we must use
`185_BPBocSS_a0.txt` to divide by the highest possible power of 3, and then
multiply the converted algNSS element by all but one of those powers of 3.)
Beware that multiplication by 3 might not correspond to adding a power of `v0`
to the generator name; see the next section for more on what the names really
mean.


### More details on generator names

There are two different sets of generator names in these files. The names in
the `185_BPAANSS_*.txt` file refers to generators in the algNSS, which is built
by constructing a `BP_*BP`-comodule resolution of `BP_*` and filtering by powers of
the augmentation ideal `I = (v_0, v_1, ...)`. A basis of permanent cycles is
chosen, and the name of each generator (probably a long linear combination of
cobar terms) is abbreviated in the tables by the name of its "leading term",
where the terms in the linear combination are ranked according to an ordering
we will now describe. Write the homological-degree n part of the cobar complex
as a direct sum of copies of `BP_*`; for `v` in `BP_*`, the element `v[n-i]` in the
cobar resolution denotes `v` in the i-th copy of `BP_*`. Then the elements
`v_0^{a_0} v_1^{a_1} ... [n-i]` are ordered using the lexicographic ordering on
the tuple of exponents `(a_0, a_1, ...)`. For example, the monomial `v1^4[1-0]` in
the output files might stand for a linear combination `2 v1^4[1-0] + v2[1-0]`,
where `v1^4[1-0]` and `v2[1-0]` are elements in homological degree 1 in the cobar
complex.

The names in the `185_BPBocSS_*.txt` files represent a different set of
generators, which arise when calculating the Bockstein spectral sequence
`E_2(S/p) => E_2(S)`. A different resolution is used, but the conventions about
leading terms and lexicographic ordering are the same as above. For example,
the fact that `[8-6]` is a Bockstein permanent cycle means that there is a term
in the cobar complex `x = [8-6] + `(higher lexicographic terms such as `v2^3[8-0]`)
such that `x + p*y` for some `y` is a cycle in the cobar complex.

These conventions lead to some quirks.

* In `185_BPBocss_a0.txt`, we see lines like

  `[5-16]	->	v0^1[5-16]+v0^1v2^7[5-1]+v0^1v2^7[5-1]+ <boundaries>`

  which leads to the relation `v0^1[5-16] / 3 = [5-16] + v2^7[5-1]` (since `3 *
  v2^7[5-1] = v0^1v2^7[5-1] + <boundaries>`), not `[5-16]` as one might have
  expected by looking at the names. This is because the cycle with leading term
  `v0^1[5-16]` is not `v0` times the cycle with leading term `[5-16]`. Instead,
  the cycle represented by `[5-16]` is (up to boundaries) `[5-16] + 2v2^7[5-1]`
  whereas the cycle represented by `v0^1[5-16]` is (up to boundaries) just
  `v0^1[5-16]`.

* There are also lines like

    `v1^1[1-0]	->	v0^1v1^1[1-0]+v0^1v1^1[1-0]+...`;

  i.e., `3 * v1^1[1-0] = 2 v0^1v1^1[1-0] + ...` instead of `v0^1v1^1[1-0] +
  ...`. This is not a contradiction, since the coefficient of `v1^1[1-0]` in
  the cycle it represents does not have to be the same as the coefficient of
  `v0^1v1^1[1-0]` in the cycle that element represents.

