* 320km => test decrease of Qangmin !!!






* Maintenant que "le netoyage des valeurs exteme" est retire de `deformation.py`, etre sur que tout ce qui utilise les deformations pour faire des stats applique le netoyage des valeurs exteme! En particulier le scaling (PDFs fait!)

 * Le pb, c que dans tous les cas, le tracking se base sur ce qu'il y a dans les fichiers quads `Q*.npz` qui eux n'ont pas l'info des quads qui vont disparaitre
   parce-que valeurs extremes de deformation => trouver un fix!!!!

* refaire les PDFs de base pour 10km (min area max and min for quads, take +1 -1 km2 at 10km!) pour RGPS et modele (avec la nouvelle method seeding / quads!!! pour le modele)

* control where deformations are cleaned (tiny values) when doing P90

* control `ZDiv<-zmin_div` in mkpdf.py !!!!


* check why NANs in qud generation logs about stat on quads!!!

* double-check how are really canceled quads with inconsistent time across vertices, and in particular for quad at t=t0 and the one at t=t0+dtNom !!!
  * the one at `t0+dtNom` should be checked, and the two canceled if any time fuckup at vertices !!!

* review the use of `areaR` in `lTisOK()` !

* review the `rscore` in `lQisOK()`




* add the reading and use (inside `01_selection_xy.py`) of the new array `jtBinN` saved in 1st npz file by `00_scan_batches.py`


* in TEST_trackin/SI3_only, apparently proble for: `[GetTimeSpan]: will start using model record ...` says 1 record too late!!!





* Use the nearest U and V point in ice part tracker! And not `0.5(U(i-1)+U(i))` as for now !!!

* Use `A` at `t` and not `A` at centered time `t+0.5*dt`, when estimating partial derivatives??? (Anton and Pierre say I should!)




* Fix the output file name of PDF npz files created by statistics: must include first and last stream, at least the year, and resolution !!!




* `deformation.py` rejects often because some buoys' time position is too far away from the mean time of the whole ice plate,
  this ctiterion should be applied in the 1st scrit: `01_selection_xy.py` !!! And these buoys removed!




* Commencer les PDFs plus loin du 0 que maintenant (y a encore le decrochage du shear!)
   * comme `0.0015` plutot que les `0.001` actuel??


* PDFs: ne pas utiliser les tiny values depuis le debut, pour que le PDF soit bien correctement normalise!
  * aussi: utiliser des bins croissants (ou decroissant) en taille (log ou exp) pour eviter les problemes en log-log

* Ajouter la partie sous-sampling du nuage de point qui marchait dans RGPS => ajouter a `generate_quad_mesh` !!!

* Ajouter l'option de stop a une date donnee dans `ice_part_track` !!!

* In quadrangle class, make a function that gives the nominal resolution of all the quadrangles....
  - then can be used,for example, in `deformation.py` (where it is now infered from file name...)






* pour la conversion degree-> km => utiliser uniquement fonction and util !!!

* trouver pourquoi dans `deformation.py` le `UM4` ne resemble pas au `UMc`, cela devrait seulement etre une moyenne spatiale !!!

* coder le trajectory (et son seeding initial a n'importe quelle location) direct from scratch
 
* utiliser l'interpolation en temps pour RGPS (et seeder le modele avec donnees interpolees au 1er janvier)



In recycle quadrangles: maybe add constraint on quadrangle Area??? Cause it can increas a lot in 3 days...


Retain mainly the constraint on quadrangle area in selection (be uber tolerant on H/L ration and min and max angle...)


















In `Tri2Quad()` remove the `pnam` argument and instead use the names of points defined in the class!!! (that might be strings of the point indices if explicit names given) !!!
=> so in the classes add the possibility to provide the vector of names: shape = (nP) of string










### Polygon

#### Class `Quadrangle`
- should not contain all the points from original triangles
  - `QuaPointIDs` is based on the IDs from the cloud point on which triangles were constructed => so from 0 to nT-1
  - `.QuaPointIdx` is the equivalent but in the reference of only the points used for Quads => so from 0 to nQ-1   (nQ<nT)
  I guess that should be the 2nd we use as default as in the future use of this class we won't min about all the disregarder triangle points...
  
  So `QuaPointIdx` should become `QuaPointIDs`
  and `QuaPointIDs` should be something like `QuaPointIDsTri`
  
  
