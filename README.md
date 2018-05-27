# Spherical_Path_Comparison
Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT and PmagPy.
Read on for more details here: https://github.com/f-i/APWP_similarity

# About the Functions
Because of complexity of the related algorithms about spherical surface geometry,
old-version directional change for each segment (seg) was actually the dif
between its previous seg and itself (if we think about it, this solution should
be able to give close results to dif between the 1st seg and itself). The
corresponding function is "spa_angpre_len_dif". This solution is also relatively
faster. I've already figured out a solution for calculating dif between each
segment and always the 1st seg, which is the function "spa_ang1st_len_dif".
However, unfortunately its execution time is much more (about 5 times) than the
previous solution. Working on the code efficiency but not very positive about
the improvement. Maybe we need to balance complexity and time. I'm also
considering to write all codes in C language for the following chapters.
* Complexity 1: Angle between two DIRECTIONAL geodesics (i.e. segments which
  are with DIRECTIONs in the order of poles' ages, but not necessarily
  successive displacement ones). Angle between two geodesics (no constraints
  on their directions) could be 2 solutions at both intersections of the
  great circles that the two geodesics are on. However, if the two geodesics
  have directions, the correct angle between them would be just one of the
  above 2 solutions.

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
&Delta;&alpha;<sub>13</sub><sup>'</sup>: (a) Geographical (orthographic
projection); (b) Cartesian.

For example, in Figure 1(a), the three segments (directional geodesics)
Seg<sub>1</sub><sup>'</sup>, Seg<sub>2</sub><sup>'</sup> and
Seg<sub>3</sub><sup>'</sup> compose an APWP-like trajectory without pole
uncertainty shown. Figure 1(b) is an analogy of these three vectors in Cartesian
space where the related calculations are more straightforward.

As we all know, azimuth &alpha;<sub>1</sub> is not equal to azimuth
&alpha;<sub>1</sub><sup>'</sup>, &alpha;<sub>2</sub> not equal to
&alpha;<sub>2</sub><sup>'</sup>, and &alpha;<sub>3</sub> not equal to
&alpha;<sub>3</sub><sup>'</sup>. In terms of directional change of each segment,
azimuth &Delta;&alpha;<sub>12</sub> is not equal to azimuth
&Delta;&alpha;<sub>12</sub><sup>'</sup>, &Delta;&alpha;<sub>23</sub> not equal
to &Delta;&alpha;<sub>23</sub><sup>'</sup>, and &Delta;&alpha;<sub>13</sub>
not equal to &Delta;&alpha;<sub>13</sub><sup>'</sup>. Although
&Delta;&alpha;<sub>13</sub> = &Delta;&alpha;<sub>12</sub> +
&Delta;&alpha;<sub>23</sub>, &Delta;&alpha;<sub>13</sub><sup>'</sup> is not
equal to the summation of &Delta;&alpha;<sub>12</sub><sup>'</sup> and
&Delta;&alpha;<sub>23</sub><sup>'</sup>. Calculating a correct
&Delta;&alpha;<sub>13</sub><sup>'</sup> is complex mainly because of the two
separate geodesics (segments) are directional like vectors. They are also in
chronological order, it does not influence the angle between these two separate
geodesics. The different situations (including the relatively simpler one shown
in Figure 1(a)) for obtaining a correct &Delta;&alpha;<sub>13</sub><sup>'</sup>
will be described in detail as follows. Proper map projections for this kind of
demonstrations are (1) Miller cylindrical, and (2) Azimuthal equidistant.
