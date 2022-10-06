### Polygon

#### Class `Quadrangle`
- should not contain all the points from original triangles
  - `QuaPointIDs` is based on the IDs from the cloud point on which triangles were constructed => so from 0 to nT-1
  - `.QuaPointIdx` is the equivalent but in the reference of only the points used for Quads => so from 0 to nQ-1   (nQ<nT)
  I guess that should be the 2nd we use as default as in the future use of this class we won't min about all the disregarder triangle points...
  
  So `QuaPointIdx` should become `QuaPointIDs`
  and `QuaPointIDs` should be something like `QuaPointIDsTri`
  
  
