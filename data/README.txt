The following data
    Laurentia_RM.txt
    Laurentia_RM_NoErr.txt
    Laurentia_F6_RM.txt
    Baltica_StableEurope_RM.txt
    Baltica_StableEurope_F6_RM.txt: reproduced from raw VGPs; uninterpolated
    Baltica_StableEurope_F6_RMi.txt: reproduced from raw VGPs; interpolated same as Table 4 "RM f=0.6" columns of Torsvik et al. 2012 Earth Sci Rev
    rot_version
    |----Laurentia_RM_1deg.txt
    |----Laurentia_RM_15deg.txt
    |----Laurentia_RM_po_15deg.txt
    |----Laurentia_RM_p_15deg.txt
    |----Laurentia_RM_45deg.txt
    |----Laurentia_RM_NoErr_45deg.txt
    |----Laurentia_RM_GMAP2015BALcoord.txt
    |----Laurentia_F6_RM_GMAP2015BALcoord.txt: uninterpolated
    |----Laurentia_F6_RMi_GMAP2015BALcoord.txt: interpolated
in this folder contain the following 10 columns:
dec	inc	age	dm	dp	dm_azi	k	n	possib_loest_age	possib_hiest_age

What they are:
	dec:				longitude of mean pole (if n>1), or VGP (if n=1), or interoperated pole (if n=0)
	inc:				latitude of mean pole (if n>1), or VGP (if n=1), or interoperated pole (if n=0)
	age:				pole age
	dm:					half of pole uncertainty major axis
	dp:					half of pole uncertainty minor axis
	dm_azi:				azimuth of pole uncertainty major axis
	k:					Fisher precision parameter
	n:					number of VGPs that put together the pole
	possib_loest_age:	low limit of the age uncertainty
	possib_hiest_age:	high limit of the age uncertainty


The following subfolders
    Laurentia_RM
    Laurentia_RM_NoErr
    Laurentia_F6_RM
    Baltica_StableEurope_RM
    Baltica_StableEurope_F6_RM
    Baltica_StableEurope_F6_RMi
    rot_version
    |----Laurentia_RM_1deg
    |----Laurentia_RM_15deg
    |----Laurentia_RM_po_15deg
    |----Laurentia_RM_p_15deg
    |----Laurentia_RM_45deg
    |----Laurentia_RM_NoErr_45deg
    |----Laurentia_RM_GMAP2015BALcoord
    |----Laurentia_F6_RM_GMAP2015BALcoord
    |----Laurentia_F6_RMi_GMAP2015BALcoord
in this main folder correspond to the above listed text files, and contain the
original VGPs for the mean poles with N>25.


All the data in the "0.result_tables" subfolder contain the following 8 columns:
00_no	01_tstop	10_spa_pol_dif	11_spa_pol_tes	20_ang_seg_dif	21_ang_seg_tes	30_len_seg_dif	31_len_seg_tes

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
