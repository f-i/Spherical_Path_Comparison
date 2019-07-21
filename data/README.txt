The following data
    Laurentia_RM.txt
    Laurentia_RM_NoIntp.txt: Laurentia_RM.txt without interpolated poles
    Laurentia_RM_NoErr.txt
    Laurentia_RM_NoErr_NoIntp.txt: Laurentia_RM_NoErr.txt without interpolated poles
    Laurentia_F6_RM.txt
    Baltica_StableEurope_RM.txt
    Baltica_StableEurope_F6_RM_NoIntp.txt: reproduced from original paleopoles; uninterpolated
    Baltica_StableEurope_F6_RMi.txt: reproduced from original paleopoles; interpolated same as Table 4 "RM f=0.6" columns of Torsvik et al. 2012 Earth Sci Rev
    Baltica_StableEurope_F6_RM_stp8.txt
    Baltica_StableEurope_F6_RM_stp12.txt
    Laurussia_RM.txt
    Laurussia_F6_RM.txt
    Gondwana_RM.txt
    Gondwana_F6_RM.txt
    rot_version
    |----Laurentia_RM_1deg.txt
    |----Laurentia_RM_1deg_NoIntp.txt: Laurentia_RM_1deg.txt without interpolated poles
    |----Laurentia_RM_15deg.txt
    |----Laurentia_RM_15deg_NoIntp.txt: Laurentia_RM_15deg.txt without interpolated poles
    |----Laurentia_RM_po_15deg.txt
    |----Laurentia_RM_p_15deg.txt
    |----Laurentia_RM_po_15deg_re.txt: trajII of pair e with random uncertainties
    |----Laurentia_RM_po_15deg_re_NoIntp.txt: Laurentia_RM_po_15deg_re.txt without interpolated poles
    |----Laurentia_RM_45deg.txt
    |----Laurentia_RM_45deg_NoIntp.txt: Laurentia_RM_45deg.txt without interpolated poles
    |----Laurentia_RM_NoErr_45deg.txt
    |----Laurentia_RM_NoErr_45deg_NoIntp.txt: Laurentia_RM_NoErr_45deg.txt without interpolated poles
    |----Laurentia_RM_GMAP2015BALcoord.txt
    |----Laurentia_F6_RM_GMAP2015BALcoord.txt: uninterpolated
    |----Laurentia_F6_RMi_GMAP2015BALcoord.txt: interpolated
    |----Laurentia_F6_RM_stp8_GMAP2015BALcoord.txt
    |----Laurentia_F6_RM_stp12_GMAP2015BALcoord.txt
    |----Laurussia_RM_GMAP2015AFRcoord.txt
    |----Laurussia_F6_RM_GMAP2015AFRcoord.txt
in this folder contain the following 10 columns:
dec	inc	age	dm	dp	dm_azi	k	n	possib_loest_age	possib_hiest_age
What they are:
	dec:				longitude of mean pole (if n>1), or paleopole (if n=1), or interoperated pole (if n=0)
	inc:				latitude of mean pole (if n>1), or paleopole (if n=1), or interoperated pole (if n=0)
	age:				pole age
	dm:					half of pole uncertainty major axis
	dp:					half of pole uncertainty minor axis
	dm_azi:				azimuth of pole uncertainty major axis
	k:					Fisher precision parameter
	n:					number of paleopoles that put together the pole
	possib_loest_age:	low limit of the age uncertainty
	possib_hiest_age:	high limit of the age uncertainty


The following subfolders
    Laurentia_RM
    Laurentia_RM_NoErr
    Laurentia_RM_NoErr_NoIntp
    Laurentia_F6_RM
    Baltica_StableEurope_RM
    Baltica_StableEurope_F6_RM_NoIntp
    Baltica_StableEurope_F6_RMi
    Baltica_StableEurope_F6_RM_stp8
    Baltica_StableEurope_F6_RM_stp12
    Laurussia_RM
    Laurussia_F6_RM
    Gondwana_RM
    Gondwana_F6_RM
    rot_version
    |----Laurentia_RM_1deg
    |----Laurentia_RM_1deg_NoIntp
    |----Laurentia_RM_15deg
    |----Laurentia_RM_15deg_NoIntp
    |----Laurentia_RM_po_15deg
    |----Laurentia_RM_po_15deg_NoIntp
    |----Laurentia_RM_p_15deg
    |----Laurentia_RM_p_15deg_NoIntp
    |----Laurentia_RM_45deg
    |----Laurentia_RM_45deg_NoIntp
    |----Laurentia_RM_NoErr_45deg
    |----Laurentia_RM_NoErr_45deg_NoIntp
    |----Laurentia_RM_GMAP2015BALcoord
    |----Laurentia_F6_RM_GMAP2015BALcoord
    |----Laurentia_F6_RMi_GMAP2015BALcoord
    |----Laurentia_F6_RM_stp12_GMAP2015BALcoord
    |----Laurussia_RM_GMAP2015AFRcoord
    |----Laurussia_F6_RM_GMAP2015AFRcoord
in this main folder correspond to the above listed text files, and contain the
original paleopoles for the mean poles with N>25.


All the difference-result data in the "0.result_tables" subfolder contain the
following 12 columns:
00_no	01_tstop	10_spa_pol_dif	11_spa_pol_tes	20_ang_seg_dif	21_ang_seg_tes	30_len_seg_dif	31_len_seg_tes	22_course_seg1	22_course_seg1	22_course_seg1	22_course_seg1
What they are:
	00_no:			index / order number
	01_tstop:		pole age
	10_spa_pol_dif:	spacial difference of each coeval pole pair
	11_spa_pol_tes:	1 or 0, i.e., distinguishable or indistinguishable
					spacial relationship of coeval uncertainties's
	20_ang_seg_dif:	angular difference of each coeval segment pair
	21_ang_seg_tes:	1 or 0, i.e., a significant or no significant angular difference,
					based on random sampling from pole uncertainties
	30_len_seg_dif:	length difference of each coeval segment pair
	31_len_seg_tes:	1 or 0, i.e., a significant or no significant length difference,
					based on random sampling from pole uncertainties
	22_course_seg1:	trajectory I's each segment's orientational change
	23_course_seg2:	trajectory II's each contemporary segment's orientational change
	32_len_seg1:	trajectory I's each segment's length
	33_len_seg2:	trajectory II's each contemporary segment's length
