# Spherical_Path_Comparison
Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT
(http://gmt.soest.hawaii.edu/) and PmagPy (https://github.com/PmagPy/PmagPy).
Read on for more details here: https://github.com/f-i/APWP_similarity

## How to Use
Please use demo.ipynb (openable by Jupyter Notebook https://jupyter.org/) to
see a few examples of measuring APWP similarity using the spComparison package
(stored in the "spComparison" folder).

Please refer to the README.txt file in the "data" folder for more details about
the data.

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

[//]: # (# Related Algorithms)

[//]: # (## The Idea Behind the Function "spa_ang1st_len_dif")

[//]: # (* Angle between two DIRECTIONAL geodesics i.e. segments which are with)
[//]: # (  DIRECTIONs in the order of poles' ages, but not necessarily successive displacement ones.)
[//]: # (  Angle between two great circles, i.e. geodesics with no directions,)
[//]: # (  could have 2 solutions at both intersections of the great circles)
[//]: # (  that the two geodesics are on. However, if the two geodesics have directions,)
[//]: # (  the correct angle between them would be just one of the above 2 solutions.)

<!--- ![](fig1directionalGeodesics.png?raw=true) -->
[//]: # (Figure 1: Directional change calculations, for two successive displacement)
[//]: # (segments i.e. directional geodesics, that describe pole wandering like)
[//]: # (Seg<sub>1</sub><sup>'</sup> & Seg<sub>2</sub><sup>'</sup> or)
[//]: # (Seg<sub>2</sub><sup>'</sup> & Seg<sub>3</sub><sup>'</sup>, e.g.)
[//]: # (&Delta;&alpha;<sub>12</sub><sup>'</sup>,)
[//]: # (&Delta;&alpha;<sub>23</sub><sup>'</sup>, and also for two separate segments like)
[//]: # (Seg<sub>1</sub><sup>'</sup> & Seg<sub>3</sub><sup>'</sup>, e.g.)
[//]: # (&Delta;&alpha;<sub>13</sub><sup>'</sup>: [a] Geographical [orthographicprojection];)
[//]: # ([b] Cartesian.)

[//]: # (For example, in Figure 1a, the three segments)
[//]: # (Seg<sub>1</sub><sup>'</sup>, Seg<sub>2</sub><sup>'</sup> and)
[//]: # (Seg<sub>3</sub><sup>'</sup> compose an APWP-like trajectory without pole)
[//]: # (uncertainty shown. Figure 1b is an analogy of these three vectors in Cartesian)
[//]: # (space where the related calculations are more straightforward.)

[//]: # (As we all know, azimuth &alpha;<sub>1</sub> is not equal to azimuth)
[//]: # (&alpha;<sub>1</sub><sup>'</sup>, &alpha;<sub>2</sub> not equal to)
[//]: # (&alpha;<sub>2</sub><sup>'</sup>, and &alpha;<sub>3</sub> not equal to)
[//]: # (&alpha;<sub>3</sub><sup>'</sup>. In terms of directional change of each segment,)
[//]: # (azimuth &Delta;&alpha;<sub>12</sub> is not equal to azimuth)
[//]: # (&Delta;&alpha;<sub>12</sub><sup>'</sup>, &Delta;&alpha;<sub>23</sub> not equal)
[//]: # (to &Delta;&alpha;<sub>23</sub><sup>'</sup>, and &Delta;&alpha;<sub>13</sub>)
[//]: # (not equal to &Delta;&alpha;<sub>13</sub><sup>'</sup>. Although)
[//]: # (&Delta;&alpha;<sub>13</sub> = &Delta;&alpha;<sub>12</sub> +)
[//]: # (&Delta;&alpha;<sub>23</sub>, &Delta;&alpha;<sub>13</sub><sup>'</sup> is not)
[//]: # (equal to the summation of &Delta;&alpha;<sub>12</sub><sup>'</sup> and)
[//]: # (&Delta;&alpha;<sub>23</sub><sup>'</sup>. Calculating a correct)
[//]: # (&Delta;&alpha;<sub>13</sub><sup>'</sup> is complex mainly because of the two)
[//]: # (separate geodesics, i.e. segments, are directional like vectors. They are also in)
[//]: # (chronological order, it does not influence the angle between these two separate)
[//]: # (geodesics. The different situations, including the relatively simpler one shown in Figure 1a,)
[//]: # (for obtaining a correct &Delta;&alpha;<sub>13</sub><sup>'</sup>)
[//]: # (will be described in detail as follows. Proper map projections for this kind of)
[//]: # (demonstrations are [a] Miller cylindrical, and [b] Azimuthal equidistant.)
