#include <mesh.h>


int lnofa[4][4] = {{1,2,3,0},{2,3,0,1},{3,0,1,2},{0,1,2,3}};


Mesh * msh_init()
{
  Mesh *msh = malloc(sizeof(Mesh));
  if ( ! msh ) return NULL;
  
  msh->Dim    = 0;
  msh->NbrVer = 0;
  msh->NbrTri = 0;
  msh->NbrTet = 0;
  
  msh->Ver = NULL;
  msh->Tri = NULL;  
  msh->Tet = NULL; 
 
  
  msh->bb[0] = 0.0; /* xmin  */
  msh->bb[1] = 0.0; /* xmax  */
  msh->bb[2] = 0.0; /* ymin  */
  msh->bb[3] = 0.0; /* ymax  */
  msh->bb[4] = 0.0; /* zmin  */
  msh->bb[5] = 0.0; /* zmax  */ 
  
  return msh;
  
}  
  
  
Mesh * msh_read(char *file)
{
  char   InpFil[1024];
  float  bufFlt[3];
  double bufDbl[3];
  int    i,bufTet[5],bufFac[4];
  int    FilVer, ref; 
  
  int64_t fmsh = 0;
  
  if ( ! file ) return NULL;
  
  Mesh * msh = msh_init();
    
  //--- set file name 
  strcpy(InpFil,file);
  if ( strstr(InpFil,".mesh") ) {
    if ( !(fmsh = GmfOpenMesh(InpFil,GmfRead,&FilVer,&msh->Dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".meshb");
    if ( !(fmsh = GmfOpenMesh(InpFil,GmfRead,&FilVer,&msh->Dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".mesh");
      if ( !(fmsh = GmfOpenMesh(InpFil,GmfRead,&FilVer,&msh->Dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n",InpFil,msh->Dim, FilVer);
  
  msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  msh->NbrTet = GmfStatKwd(fmsh, GmfTetrahedra);
  msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);
  
  /* allocate arrays */
  msh->Ver = calloc( (msh->NbrVer+1), sizeof(Vertex)       );
  msh->Tri = calloc( (msh->NbrTri+1), sizeof(Triangle)     );  
  msh->Tet = calloc( (msh->NbrTet+1), sizeof(Tetrahedron)  ); 
  
  
   GmfGotoKwd(fmsh, GmfVertices);
   if ( msh->Dim == 2 ) {
     if ( FilVer == GmfFloat ) {		// read 32 bits float
       for (i=1; i<=msh->NbrVer; ++i) {
         GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
         msh->Ver[i].Crd[0] = (double)bufFlt[0];
         msh->Ver[i].Crd[1] = (double)bufFlt[1];
         msh->Ver[i].Crd[2] = 0.0;
       }
     }
     else  {	// read 64 bits float
       for (i=1; i<=msh->NbrVer; ++i) {
         GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
         msh->Ver[i].Crd[0] = bufDbl[0];
         msh->Ver[i].Crd[1] = bufDbl[1];
         msh->Ver[i].Crd[2] = 0.0;
       }  
     }
   }
   else {
     if ( FilVer == GmfFloat ) {		// read 32 bits float
       for (i=1; i<=msh->NbrVer; ++i) {
         GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &bufFlt[2], &ref);
         msh->Ver[i].Crd[0] = (double)bufFlt[0];
         msh->Ver[i].Crd[1] = (double)bufFlt[1];
         msh->Ver[i].Crd[2] = (double)bufFlt[2];
       }
     }
     else  {	// read 64 bits float
       for (i=1; i<=msh->NbrVer; ++i) {
         GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &bufDbl[2], &ref);
         msh->Ver[i].Crd[0] = bufDbl[0];
         msh->Ver[i].Crd[1] = bufDbl[1];
         msh->Ver[i].Crd[2] = bufDbl[2];
       }  
     }
   }
   
  
   //--- read tetrahedra 
  GmfGotoKwd(fmsh, GmfTetrahedra);
  for (i=1; i<=msh->NbrTet; ++i) {
    GmfGetLin(fmsh, GmfTetrahedra, &bufTet[0], &bufTet[1], &bufTet[2], &bufTet[3], &bufTet[4]);
    msh->Tet[i].Ver[0]    = bufTet[0];
    msh->Tet[i].Ver[1]    = bufTet[1];
    msh->Tet[i].Ver[2]    = bufTet[2];
    msh->Tet[i].Ver[3]    = bufTet[3];
    msh->Tet[i].Ref       = bufTet[4];
  }
  
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i=1; i<=msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufFac[0], &bufFac[1], &bufFac[2], &bufFac[3]);
    msh->Tri[i].Ver[0]    = bufFac[0];
    msh->Tri[i].Ver[1]    = bufFac[1];
    msh->Tri[i].Ver[2]    = bufFac[2];
    msh->Tri[i].Ref       = bufFac[3];
  }
  
  GmfCloseMesh(fmsh);
  
  return msh;
  
}


int compar_vertex(const void *a, const void *b)
{
  Vertex *va = (Vertex *) a;
  Vertex *vb = (Vertex *) b;
  return ( vb->icrit - va->icrit );
}

int compar_triangle(const void *a, const void *b)
{
  Triangle *va = (Triangle *) a;
  Triangle *vb = (Triangle *) b;
  return ( vb->icrit - va->icrit );
}

int compar_tetrahedron(const void *a, const void *b)
{
  Tetrahedron *va = (Tetrahedron *) a;
  Tetrahedron *vb = (Tetrahedron *) b;
  return ( vb->icrit - va->icrit );
}


int msh_reorder(Mesh *msh)
{
  
  int iTet, iTri, iVer;
  
  if ( ! msh            ) return 0;
  if ( msh->NbrVer <= 0 ) return 0;
  
  /* compute bonding box */
  msh->bb[0] = msh->Ver[1].Crd[0];
  msh->bb[1] = msh->Ver[1].Crd[1];
  msh->bb[2] = msh->Ver[1].Crd[2];
  msh->bb[3] = msh->Ver[1].Crd[0];
  msh->bb[4] = msh->Ver[1].Crd[1];
  msh->bb[5] = msh->Ver[1].Crd[2];
  for(iVer=1; iVer<=msh->NbrVer; iVer++) {
    /* todo msh->bb : used to compute the Z-curve index */
    double x = msh->Ver[iVer].Crd[0];
    double y = msh->Ver[iVer].Crd[1];
    double z = msh->Ver[iVer].Crd[2];
    if(msh->bb[0] > x){msh->bb[0] = x;}
    if(msh->bb[1] > y){msh->bb[1] = y;}
    if(msh->bb[2] > z){msh->bb[2] = z;}
    if(msh->bb[3] < x){msh->bb[3] = x;}
    if(msh->bb[4] < y){msh->bb[4] = y;}
    if(msh->bb[5] < z){msh->bb[5] = z;}
  }
  
  for(iVer=1; iVer<=msh->NbrVer; iVer++) {
    msh->Ver[iVer].icrit  = rand();   /* change the randon  by Z  order */
    msh->Ver[iVer].idxNew = iVer;
    msh->Ver[iVer].idxOld = iVer;
  }
  //qsort(&msh->Ver[1],msh->NbrVer,sizeof(Vertex), compar_vertex); 

  
  /* update idxNew for vertices */
   
  
  /* re-assign triangles and tets ids */
  
  
  /* sort triangles */
  for(iTri=1; iTri<=msh->NbrTri; iTri++) {
    msh->Tri[iTri].icrit  = rand();  
  }
  qsort(&msh->Tri[1],msh->NbrTri,sizeof(Triangle), compar_triangle);
  
  /* sort tetrahedra */
  for(iTet=1; iTet<=msh->NbrTet; iTet++) {
    msh->Tet[iTet].icrit  = rand();  
  }
  qsort(&msh->Tet[1],msh->NbrTet,sizeof(Tetrahedron), compar_tetrahedron);
  
  return 1;
}

int bounding_box(Mesh *msh){
  int iVer;
  msh->bb[0] = msh->Ver[1].Crd[0];
  msh->bb[1] = msh->Ver[1].Crd[1];
  msh->bb[2] = msh->Ver[1].Crd[2];
  msh->bb[3] = msh->Ver[1].Crd[0];
  msh->bb[4] = msh->Ver[1].Crd[1];
  msh->bb[5] = msh->Ver[1].Crd[2];
  for(iVer=1; iVer<=msh->NbrVer; iVer++) {
    /* todo msh->bb : used to compute the Z-curve index */
    double x = msh->Ver[iVer].Crd[0];
    double y = msh->Ver[iVer].Crd[1];
    double z = msh->Ver[iVer].Crd[2];
    if(msh->bb[0] > x){msh->bb[0] = x;}
    if(msh->bb[1] > y){msh->bb[1] = y;}
    if(msh->bb[2] > z){msh->bb[2] = z;}
    if(msh->bb[3] < x){msh->bb[3] = x;}
    if(msh->bb[4] < y){msh->bb[4] = y;}
    if(msh->bb[5] < z){msh->bb[5] = z;}
  }
  return 1;
}

int mesh_reorder_random(Mesh* msh){
  return 1;
}

int msh_reorder_z(Mesh* msh){
  int iVer, iTri, iTet;
  int new_add[msh->NbrVer + 1];
  bounding_box(msh);
  for(iVer=1; iVer<=msh->NbrVer; iVer++) {
    msh->Ver[iVer].icrit  = get_crit(msh->Ver[iVer], msh->bb, 21);
    msh->Ver[iVer].idxNew = iVer;
    msh->Ver[iVer].idxOld = iVer;
  }
  qsort(&msh->Ver[1],msh->NbrVer,sizeof(Vertex), compar_vertex);
  for(iVer = 1; iVer <= msh->NbrVer; iVer++){
    msh->Ver[iVer].idxNew = iVer;
    new_add[msh->Ver[iVer].idxOld] = iVer;
  }

  for(iTri = 1; iTri <= msh->NbrTri; iTri++){
    msh->Tri[iTri].icrit = msh->NbrVer;
    for(int j = 0; j < 3; j++){
      msh->Tri[iTri].Ver[j] = new_add[msh->Tri[iTri].Ver[j]];
      if(msh->Tri[iTri].icrit > msh->Tri[iTri].Ver[j]){msh->Tri[iTri].icrit = msh->Tri[iTri].Ver[j];}
    }
  }
  qsort(&msh->Tri[1],msh->NbrTri,sizeof(Triangle), compar_triangle);

  for(iTet = 1; iTet <= msh->NbrTet; iTet++){
    msh->Tet[iTet].icrit = msh->NbrVer;
    for(int j = 0; j < 4; j++){
      msh->Tet[iTet].Ver[j] = new_add[msh->Tet[iTet].Ver[j]];
      if(msh->Tet[iTet].icrit > msh->Tet[iTet].Ver[j]){msh->Tet[iTet].icrit = msh->Tet[iTet].Ver[j];}
    }
  }
  qsort(&msh->Tet[1],msh->NbrTet,sizeof(Tetrahedron), compar_tetrahedron);

  return 1;
}

long long int get_crit(Vertex v, double* bb, int depth){
  if(depth == 0){return 0;}
  else{
    long long int crit = 0;
    double xmid = (bb[3] - bb[0]) / 2;
    double ymid = (bb[4] - bb[1]) / 2;
    double zmid = (bb[5] - bb[2]) / 2;
    double x = v.Crd[0];
    double y = v.Crd[1];
    double z = v.Crd[2];
    int bz = (int) z > zmid;
    int by = (int) y > ymid;
    int bx = (int) x > xmid;
    crit = (bz << 2) + ((bz ^ by) << 1) + (bz ^ bx);
    double new_bb[6] = {bb[0] + bx * xmid, bb[1] + by * ymid, bb[2] + bz * zmid, bb[3] - (1 - bx) * xmid, bb[4] - (1 - by) * ymid, bb[5] - (1 - bz) * zmid};
    return crit << (3 * (depth - 1)) + get_crit(v, new_bb, depth - 1);
  }
}

int    msh_write(Mesh *msh, char *file)
{
   int iVer, iTri, iTfr, iTet;
   int FilVer = 2;
  
   if ( ! msh  ) return 0;
   if ( ! file ) return 0;
   
   int64_t fmsh = GmfOpenMesh(file, GmfWrite, FilVer, msh->Dim);
   if ( fmsh <=  0 ) {
     printf("  ## ERROR: CANNOT CREATE FILE \n");
     return 0;
   }

   GmfSetKwd(fmsh, GmfVertices,   msh->NbrVer);
   if ( msh->Dim == 3 ) {
     for ( iVer=1; iVer<=msh->NbrVer; iVer++){
       GmfSetLin(fmsh, GmfVertices, msh->Ver[iVer].Crd[0],msh->Ver[iVer].Crd[1],msh->Ver[iVer].Crd[2],0); 
     }
   }
   else {
     for ( iVer=1; iVer<=msh->NbrVer; iVer++){
       GmfSetLin(fmsh, GmfVertices, msh->Ver[iVer].Crd[0],msh->Ver[iVer].Crd[1],0); 
     }
   }
   
   GmfSetKwd(fmsh, GmfTriangles, msh->NbrTri);
   for ( iTfr=1; iTfr<=msh->NbrTri; iTfr++)  
     GmfSetLin(fmsh, GmfTriangles, msh->Tri[iTfr].Ver[0], msh->Tri[iTfr].Ver[1], msh->Tri[iTfr].Ver[2], msh->Tri[iTfr].Ref);  
   
   GmfSetKwd(fmsh, GmfTetrahedra, msh->NbrTet);
   for ( iTet=1; iTet<=msh->NbrTet; iTet++) 
     GmfSetLin(fmsh, GmfTetrahedra, msh->Tet[iTet].Ver[0], msh->Tet[iTet].Ver[1],msh->Tet[iTet].Ver[2], msh->Tet[iTet].Ver[3], 0);
   
   GmfCloseMesh(fmsh);
   
   return 1;
   
}






int  msh_neighborsQ2(Mesh *msh)
{
  int iTet, iFac, jTet, jFac, ip1, ip2 ,ip3 , jp1, jp2, jp3;
  int count = 0;
  double time_cur = 0, time_est = 0;
  double to =  GetWallClock();
  double ti;
  
  if ( ! msh ) return 0;
  
  for(iTet=1; iTet<=msh->NbrTet; iTet++) {
    for(iFac=0; iFac<4; iFac++) {
      ip1 = msh->Tet[iTet].Ver[lnofa[iFac][0]];
      ip2 = msh->Tet[iTet].Ver[lnofa[iFac][1]];
      ip3 = msh->Tet[iTet].Ver[lnofa[iFac][2]];
      /* find the Tet different from iTet that has ip1, ip2, ip2 as vertices */
      for(jTet=1; jTet<=msh->NbrTet; jTet++) {
        if ( iTet == jTet ) continue;
        for(jFac=0; jFac<4; jFac++) {
          jp1 = msh->Tet[jTet].Ver[lnofa[jFac][0]];
          jp2 = msh->Tet[jTet].Ver[lnofa[jFac][1]];
          jp3 = msh->Tet[jTet].Ver[lnofa[jFac][2]];
          /* compare the 6 points */
          int a = (ip1 - jp1) & (ip1 - jp2) & (ip1 - jp3);
          int b = (ip2 - jp1) & (ip2 - jp2) & (ip2 - jp3);
          int c = (ip3 - jp1) & (ip3 - jp2) & (ip3 - jp3);
          if(a == 0 && b == 0 && c == 0){
            msh->Tet[iTet].Voi[lnofa[iFac][3]] = jTet;
          }
        }
      }
    }
    if(iTet%10 == 0){
      ti =  GetWallClock();
      time_cur = ti - to;
      count++;
      time_est = (time_cur / (10 * count)) * msh->NbrTet;
      printf("Temps écoulé estimé à ");
      print_time(time_cur);
      printf(", temps restant estimé à ");
      print_time(time_est - time_cur);
      printf(", temps total estimé à ");
      print_time(time_est);
      printf("\n");
    }
  }
  
  return 1;
}   

void print_time(double t){
  int h, m, s;
  s = ((int) t) % 60;
  m = (((int) t) % 3600) / 60;
  h = ((int) t) / 3600;

  printf("%d h %d m %d s", h, m, s);
}

int  msh_neighbors(Mesh *msh)
{
  int iTet, iFac, ip1, ip2 ,ip3 ;
  
  if ( ! msh ) return 0;
  
  /* initialize HashTable */
  
  for(iTet=1; iTet<=msh->NbrTet; iTet++) {
    for(iFac=0; iFac<4; iFac++) {
      ip1 = msh->Tet[iTet].Ver[lnofa[iFac][0]];
      ip2 = msh->Tet[iTet].Ver[lnofa[iFac][1]];
      ip3 = msh->Tet[iTet].Ver[lnofa[iFac][2]];
      /* compute the key : ip1+ip2+ip3   */
      /* do we have objects as that key   hash_find () */
      /*  if yes ===> look among objects and potentially update Voi */
      /*  if no  ===> add to hash table.  hash_add()   */
    }
  }
  return 1;
}



HashTable * hash_init(int SizHead, int NbrMaxObj){
  HashTable * new_table = malloc(sizeof(HashTable));
  new_table->SizHead = SizHead;
  new_table->NbrObj = 0;
  new_table->NbrMaxObj = NbrMaxObj;
  new_table->Head = malloc(SizHead * sizeof(int));
  new_table->LstObj = malloc(NbrMaxObj * sizeof(int6));
  return new_table;
}

int hash_add(HashTable *hsh, int ip1, int ip2, int ip3, int iTet){

  int key = (ip1 + ip2 + ip3) % hsh->SizHead;

  hsh->NbrObj++;
  memcpy(&(hsh->LstObj[hsh->NbrObj]), (int6) {ip1, ip2, ip3, iTet, 0, 0}, sizeof(int6));

  if(hsh->Head[key] == 0){
    hsh->Head[key] = hsh->NbrObj;
  }else{
    int6 * element = &(hsh->LstObj[hsh->Head[key]]);
    while(*element[5] != 0){
      element = &(hsh->LstObj[*element[5]]);
    }
    *element[5] = hsh->NbrObj;
  }
  return hsh->NbrObj;
}

int hash_find(HashTable * hsh, int ip1, int ip2, int ip3){
  
  int key = (ip1 + ip2 + ip3) % hsh->SizHead;

  if(hsh->Head[key] == 0){
    return 0;
  }else{
    int next = hsh->Head[key];
    int search = 1;

    while(search){
      int jp1 = hsh->LstObj[next][0];
      int jp2 = hsh->LstObj[next][1];
      int jp3 = hsh->LstObj[next][2];
      int a = (ip1 - jp1) & (ip1 - jp2) & (ip1 - jp3);
      int b = (ip2 - jp1) & (ip2 - jp2) & (ip2 - jp3);
      int c = (ip3 - jp1) & (ip3 - jp2) & (ip3 - jp3);
      if(a == 0 && b == 0 && c == 0){
        search = 0;
        return next;
      }
      next = hsh->LstObj[next][5];
      if(next == 0){
        search = 0;
        return 0;
      }
    }
  }
}
