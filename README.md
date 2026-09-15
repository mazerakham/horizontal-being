# Horizontal Being

A game set in the first Heisenberg group. You can only move horizontally, and you see along horizontal geodesics. Four levels.

Play: https://mazerakham.github.io/horizontal-being/ (desktop browser with WebGL2 and a keyboard)

- `index.html` is the whole game: one file, KaTeX and three.js from cdnjs.
- Clicking the eye lays down one full turn of that geodesic as a road. From anywhere else, the eye draws the road through `HG.preimages`, which inverts the eye: every probe from you that reaches a given point.
- `node test/geometry.test.mjs` checks the geometry block shipped in `index.html` against RK4 integration of the horizontal lift and a brute-force march of the eye's hit test.
