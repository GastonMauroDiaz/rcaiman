# rcaiman 1.0.6.9000



    minor fix for fisheye_to_equidistant() (again)

fijarse si este error no estaba en CRAN

commit 49418e2871cdaf53bc9c3aa3daaf990ae151701a
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Mon Sep 12 11:16:00 2022 -0300


    add mask_blue_sky() and remove sf from imports

commit 853b8193d77a5d0a9a4140625eabd7bf00e080fa
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Sep 9 09:35:58 2022 -0300

   
    Fix local_fuzzy_thr(), add 2 functions

commit 737ce70e7085bea2784fdc3bb88b41faf91c1794
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Tue Sep 6 17:29:07 2022 -0300

    Add examples of restricted view photography

commit fe5ad6aae00581d8d618549954d60edb1823a576
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Tue Sep 6 12:16:02 2022 -0300

    Add extract_dn() and other modifications

commit 0f4c5b7e90eb43d87f0c59ce76e53905ffcfde53
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Mon Sep 5 16:00:44 2022 -0300

    fix qtree()

commit b75e55d288a4d7523d7e5e73fcc80e800dfe84e6
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Sep 2 17:28:30 2022 -0300

    add (but not export) chessboard and qtree

commit be5306284c35c3b30b3ea3b01bb7195ff55801f2
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Thu Sep 1 11:46:06 2022 -0300

    Improve tandem working with HSP software

commit 3fbe1e3959f4caa732bf79cb8489c406e5662dc8
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Wed Aug 31 17:14:59 2022 -0300

    Make minor improvements, add colorfulness()

commit 2dfde9a77b938ea8c03c730e8d5ec25958d9c35f
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Thu Aug 25 16:48:25 2022 -0300

   
    add orientation to azimuth_image()

commit e95014a548d2cd1bbbd4138429d9191adc864907
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Aug 12 11:19:22 2022 -0300

    fix expand_noncircular()

commit 03506e5f9aa070032cf67d83d9e6d867fcb7c35f
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Wed Aug 10 11:50:16 2022 -0300

    Minor change to find_sky_pixels()

commit b8ee67ab3e2c696b706659ff27a0d2390ffc63b0
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Tue Aug 9 17:25:49 2022 -0300



    Add ootb_obia() and improve doc

commit ad96c0616f376432eeed8d7fdca9506079e1d591
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Wed Aug 3 16:26:18 2022 -0300

 
    Many adds and changes
    
export(fisheye_to_equidistant)
export(fisheye_to_pano)
(ootb_sky_reconstruction)
export(read_manual_input)
export(read_opt_sky_coef)
export(row_col_from_zenith_azimuth)
export(write_sky_points)
export(write_sun_coord)
export(zenith_azimuth_from_row_col)

commit c10537e1201e67e7a0e24c460fc517df30d04be2
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Tue Jul 5 14:20:18 2022 -0300

    Fix fit_cone*(), and many other changes
    
    Finalmente, cambie el comportamiento de fit_cone_shaped

commit 14ab5fdcdbf33d684abec924ce0bb9a0d2a027b6
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Mon Jul 4 17:13:49 2022 -0300

    Add two functions and do many changes
    
    find_sky_pixels_nonnull_criteria
    Mask_sunlit canopy
    extract_zenith_dn evolve to extract_rl()
    

commit a49ce576f03ac9e87ef0d5ba9771e569ff70a0cb
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Jul 1 10:24:37 2022 -0300

    Add obia()

commit 0d8fe7f3f9d88fe758f9c189ed740de6b836dc57
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Thu Jun 30 17:40:48 2022 -0300

    Add interpolate...() and ootb_sky_recon...()

commit 6ab950392ea48e1505b7602a76422daacc5a1116
Author: Gaston Mauro Diaz <gastonmaurodiaz@gmail.com>
Date:   Wed Jun 29 17:31:58 2022 -0300

    Add deffuzify()

commit 944f3eabbba1cb2275a861f7023f29544baea97e
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Wed Jun 29 13:00:39 2022 -0300

    Add extract_zenith_dn() and polar_qtree()

commit eeec07c63b2e1f03c86d486ad3ecc7cde7b67cc1
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Jun 24 17:40:04 2022 -0300

    Add fit_cie_sky_model()

commit ede3f2da82664ef605611e606192d600a2874754
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Jun 24 16:36:04 2022 -0300

    Remove 2 dependencies, rename 1 fun, add 1 fun
    
export(extract_sky_marks) to export(extract_sky_points)
export(extract_sun_mark) to export(extract_sun_coord)
No encontre funcion nueva


commit 4d70e1cf1d68b2f699be02e7347ac6f72e0a1b2a
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Jun 24 12:26:47 2022 -0300

    Add extract_sky_marks(), other minor changes

commit f01b4c7268edbb4cbd556dda100e73e32c860310
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Thu Jun 23 17:17:55 2022 -0300

   
    Change dependency from raster to terra

commit 0068e68afcb32caa2f906b8716c9bfa4d5e028d2 (tag: v.0.1.1)
Author: Gaston Diaz <gdiaz@ciefap.org.ar>
Date:   Fri Jan 21 15:39:41 2022 -0300


