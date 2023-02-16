 * Commencer les PDFs plus loin du 0 que maintenant (y a encore le decrochage du shear!)
   * comme `0.0015` plutot que les `0.001` actuel??


* PDFs: ne pas utiliser les tiny values depuis le debut, pour que le PDF soit bien correctement normalise!
  * aussi: utiliser des bins croissants (ou decroissant) en taille (log ou exp) pour eviter les problemes en log-log

* Ajouter la partie sous-sampling du nuage de point qui marchait dans RGPS => ajouter a `generate_quad_mesh` !!!

* Ajouter l'option de stop a une date donnee dans `ice_part_track` !!!









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
  
  
