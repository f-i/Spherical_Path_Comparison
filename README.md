# Spherical_Path_Comparison
Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT
(http://gmt.soest.hawaii.edu/) and PmagPy (https://github.com/PmagPy/PmagPy).
Read on for more details here: https://github.com/f-i/APWP_similarity

## About the Two Ideas About Segment Angle Difference
Compared with the algorithm about segment (seg) angle difference in the function
"spa_angpre_len_dif", the old-version orientation change for each seg in the
function "spa_ang1st_len_dif" was actually the difference between each seg and
always the 1st seg. However, after seeing the calculated results we realized
that the definition of the difference between each seg and the 1st seg is far
more complicated than expected in terms of spherical surface geometry,
especially for those seg far from the 1st seg. Because of this complexity,
another disadvantage of this function is its long calculating time, which is
about 5 times of the time the presently used solution (spa_angpre_len_dif)
consumes. The calculation for seg angle difference in the function
"spa_angpre_len_dif" is described in detail in the section 2.4.1 of the paper
stored in the repository: https://github.com/f-i/APWP_similarity. The
orientation change is the difference between its previous seg and itself.

[comment]: <> # Related Algorithms

[comment]: <> ## The Idea Behind the Function "spa_ang1st_len_dif"
[comment]: <>
[comment]: <> * Angle between two DIRECTIONAL geodesics (i.e. segments which are with
[comment]: <>   DIRECTIONs in the order of poles' ages, but not necessarily successive
[comment]: <>   displacement ones). Angle between two geodesics (no constraints on their
[comment]: <>   directions) could be 2 solutions at both intersections of the great circles
[comment]: <>   that the two geodesics are on. However, if the two geodesics have directions,
[comment]: <>   the correct angle between them would be just one of the above 2 solutions.
[comment]: <>
[comment]: <> ![](fig1directionalGeodesics.png?raw=true)
[comment]: <> Figure 1: Directional change calculations, for two successive displacement
[comment]: <> segments (directional geodesics) that describe pole wandering like
[comment]: <> Seg<sub>1</sub><sup>'</sup> & Seg<sub>2</sub><sup>'</sup> or
[comment]: <> Seg<sub>2</sub><sup>'</sup> & Seg<sub>3</sub><sup>'</sup>, e.g.
[comment]: <> &Delta;&alpha;<sub>12</sub><sup>'</sup>,
[comment]: <> &Delta;&alpha;<sub>23</sub><sup>'</sup>, and also for two separate segments like
[comment]: <> Seg<sub>1</sub><sup>'</sup> & Seg<sub>3</sub><sup>'</sup>, e.g.
[comment]: <> &Delta;&alpha;<sub>13</sub><sup>'</sup>: (a) Geographical (orthographic
[comment]: <> projection); (b) Cartesian.
[comment]: <>
[comment]: <> For example, in Figure 1(a), the three segments (directional geodesics)
[comment]: <> Seg<sub>1</sub><sup>'</sup>, Seg<sub>2</sub><sup>'</sup> and
[comment]: <> Seg<sub>3</sub><sup>'</sup> compose an APWP-like trajectory without pole
[comment]: <> uncertainty shown. Figure 1(b) is an analogy of these three vectors in Cartesian
[comment]: <> space where the related calculations are more straightforward.
[comment]: <>
[comment]: <> As we all know, azimuth &alpha;<sub>1</sub> is not equal to azimuth
[comment]: <> &alpha;<sub>1</sub><sup>'</sup>, &alpha;<sub>2</sub> not equal to
[comment]: <> &alpha;<sub>2</sub><sup>'</sup>, and &alpha;<sub>3</sub> not equal to
[comment]: <> &alpha;<sub>3</sub><sup>'</sup>. In terms of directional change of each segment,
[comment]: <> azimuth &Delta;&alpha;<sub>12</sub> is not equal to azimuth
[comment]: <> &Delta;&alpha;<sub>12</sub><sup>'</sup>, &Delta;&alpha;<sub>23</sub> not equal
[comment]: <> to &Delta;&alpha;<sub>23</sub><sup>'</sup>, and &Delta;&alpha;<sub>13</sub>
[comment]: <> not equal to &Delta;&alpha;<sub>13</sub><sup>'</sup>. Although
[comment]: <> &Delta;&alpha;<sub>13</sub> = &Delta;&alpha;<sub>12</sub> +
[comment]: <> &Delta;&alpha;<sub>23</sub>, &Delta;&alpha;<sub>13</sub><sup>'</sup> is not
[comment]: <> equal to the summation of &Delta;&alpha;<sub>12</sub><sup>'</sup> and
[comment]: <> &Delta;&alpha;<sub>23</sub><sup>'</sup>. Calculating a correct
[comment]: <> &Delta;&alpha;<sub>13</sub><sup>'</sup> is complex mainly because of the two
[comment]: <> separate geodesics (segments) are directional like vectors. They are also in
[comment]: <> chronological order, it does not influence the angle between these two separate
[comment]: <> geodesics. The different situations (including the relatively simpler one shown
[comment]: <> in Figure 1(a)) for obtaining a correct &Delta;&alpha;<sub>13</sub><sup>'</sup>
[comment]: <> will be described in detail as follows. Proper map projections for this kind of
[comment]: <> demonstrations are (1) Miller cylindrical, and (2) Azimuthal equidistant.
