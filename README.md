# Spherical_Path_Comparison
Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT and PmagPy.

# Work in Progress
1. IMPORTANT: Testing and debugging per-segment tests separately on angle
   difference (dif) and length dif
2. IMPORTANT: Because of complexity of the related algorithms about spherical
   surface geometry, old-version directional change for each segment (seg) was
   actually the dif between its previous seg and itself (if we think about it,
   this solution should be able to give close results to dif between the 1st seg
   and itself). This solution is also relatively fast. I've already figured out
   a solution for calculating dif between each segment and always the 1st seg.
   However, unfortunately its execution time is much more (like 3 times) than
   the previous solution. Working on the code efficiency but not very positive
   about the improvement. Maybe we need to balance complexity and time.
   * Complexity 1: Angle between two DIRECTIONAL geodesics (i.e. segments which
     are with DIRECTIONs in the order of poles' ages, but not neccessarily
     successive displacement ones); Angle between two geodesics (no constraints
     on their directions) could be 2 solutions at both intersections of the
     great circles that the two geodesics are on.
   * Complexity 2: For per-seg test on seg angle dif, it is hard to be
     implemented for the present ang dif solution (i.e. dif from the previous
     seg), so completing angle change of every seg compared to the 1st seg
     is neccessary now (achieved in v0.4.4).

# Related Algorithms

## About Segment Angle Difference

![](fig1directionalGeodesics.png?raw=true)
Figure 1: Directional change calculations, for two successive displacement
segments (directional geodesics) that describe pole wandering like
Seg<sub>1</sub><sup>'</sup> & Seg<sub>2</sub><sup>'</sup> or
Seg<sub>2</sub><sup>'</sup> & Seg<sub>3</sub><sup>'</sup>, e.g.
&Delta;&alpha;<sub>12</sub><sup>'</sup>,
&Delta;&alpha;<sub>23</sub><sup>'</sup>, and also for two separate segments like
Seg<sub>1</sub><sup>'</sup> & Seg<sub>3</sub><sup>'</sup>, e.g.
&Delta;&alpha;<sub>13</sub><sup>'</sup>: (a) Geographical; (b) Cartesian.

For example, in Figure 1(a), the three segments (directional geodesics)
Seg<sub>1</sub><sup>'</sup>, Seg<sub>2</sub><sup>'</sup> and
Seg<sub>3</sub><sup>'</sup> compose an APWP-like trajectory without pole
uncertainty shown. Figure 1(b) is an analogy of these three vectors in Cartesian
space where the related calculations are more straightforward.

Please note that &Delta;&alpha;<sub>12</sub> is not equal to
&Delta;&alpha;<sub>12</sub><sup>'</sup>, &Delta;&alpha;<sub>23</sub> not
equal to &Delta;&alpha;<sub>23</sub><sup>'</sup>, &Delta;&alpha;<sub>13</sub>
not equal to &Delta;&alpha;<sub>14</sub><sup>'</sup>.
