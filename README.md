## ModuliCells
A computational framework for constructing cellular decompositions of moduli spaces of Riemann surfaces.

# Overview
This project implements a combinatorial approach to moduli spaces
of Riemann surfaces.

The main example is the genus two moduli space $M_2$, where every
surface is hyperelliptic and can be described by six Weierstrass
points on the sphere.

Using Delaunay/Voronoi decompositions, we classify combinatorial
types and construct the corresponding polyhedral cells.

# Mathematical Background

See [needs to be filled...]

# Features

Implemented:

✓ Enumeration of genus-2 Delaunay/Voronoi types

✓ Construction of polyhedral cells

✓ Recursive face enumeration

✓ Automorphism group computation

✓ Orbifold Euler characteristic computation

# Example

For genus 2, the program reconstructs the cell complex and obtains

effective orbifold Euler characteristic:

    $$ -1/120 $$

Taking into account the universal hyperelliptic involution:

    $$ \chi_{orb}(M_2) = -1/240 $$

which agrees with the Harer-Zagier formula.

# Current limitation

The current implementation describes the unmarked moduli problem
through Weierstrass data.

Mapping class group actions invisible to the Weierstrass model
(e.g. the hyperelliptic kernel) are not represented internally
and must be accounted for separately.

# Future directions

Future extensions:

- Ribbon graph based decompositions for $M_g$, $g \geq 3$

- Mapping class group quotient structures

- Morse theoretic analysis on moduli spaces

- Comparison with Deligne-Mumford compactification