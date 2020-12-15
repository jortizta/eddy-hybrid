#include <stdio.h>
#include <ANN/ANN.h>

using namespace std;					// make std:: accessible

extern "C"
{

  int committee(double* query, int* nQuery, double* vertex, int* nPts, int* Idx, int* k, int *iflag, int *mproc, int *imb3D)
  {

      int i,j,nq,np,nn;
      int dim=3;
      int myrank=*mproc;
      int maxPts = *nPts;
      double eps=0.0;

      np = *nPts;
      nq = *nQuery;
      nn = *k;

      static ANNpointArray dataPts;			// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;				

      if(*iflag==0)
      {
        dataPts = annAllocPts(maxPts, dim);	// allocate data points
        for (i=0; i<*nPts; i++)
        {
	  dataPts[i][0] = *(vertex+i*dim+0);
	  dataPts[i][1] = *(vertex+i*dim+1)*(*imb3D);
	  dataPts[i][2] = *(vertex+i*dim+2);
        }

        kdTree = new ANNkd_tree(		// build search structure
		dataPts,			// the data points
		*nPts,				// number of points
		dim);				// dimension of space
      }

      if(*iflag==1)
      {
        nnIdx = new ANNidx[*k];			// allocate near neigh indices
        dists = new ANNdist[*k];		// alloc

        queryPt = annAllocPt(dim);		// allocate query point
        for (i=0; i<nq; i++)
        {
          queryPt[0] = *(query+i*dim+0);
	  queryPt[1] = *(query+i*dim+1)*(*imb3D);
          queryPt[2] = *(query+i*dim+2);
	  kdTree->annkSearch(		// search
		  queryPt,		// query point
  		  nn,			// number of near neighbors
	  	  nnIdx,		// nearest neighbors (returned)
		  dists,		// distance (returned)
		  eps);						

          for(j=0; j<nn; j++)
          {
            *(Idx+j*nq+i) = nnIdx[j]+1;
          }
        }
        annDeallocPt(queryPt);

        delete [] nnIdx;
        delete [] dists;
      }


      if(*iflag==2)
      {
        delete kdTree;
        annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }


  int committee_(double* query, int* nQuery, double* vertex, int* nPts, int* Idx, int* k, int *iflag, int *mproc, int *imb3D)
  {

      int i,j,nq,np,nn;
      int dim=3;
      int myrank=*mproc;
      int maxPts = *nPts;
      double eps=0.0;

      np = *nPts;
      nq = *nQuery;
      nn = *k;

      static ANNpointArray dataPts;			// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNbd_tree*   kdTree;				

      if(*iflag==0)
      {
        dataPts = annAllocPts(maxPts, dim);	// allocate data points
        for (i=0; i<*nPts; i++)
        {
	  dataPts[i][0] = *(vertex+i*dim+0);
	  dataPts[i][1] = *(vertex+i*dim+1)*(*imb3D);
	  dataPts[i][2] = *(vertex+i*dim+2);
        }

        kdTree = new ANNbd_tree(		// build search structure
		dataPts,			// the data points
		*nPts,				// number of points
		dim);				// dimension of space
      }

      if(*iflag==1)
      {
        nnIdx = new ANNidx[*k];			// allocate near neigh indices
        dists = new ANNdist[*k];		// alloc

        queryPt = annAllocPt(dim);		// allocate query point
        for (i=0; i<nq; i++)
        {
          queryPt[0] = *(query+i*dim+0);
	  queryPt[1] = *(query+i*dim+1)*(*imb3D);
          queryPt[2] = *(query+i*dim+2);
	  kdTree->annkSearch(		// search
		  queryPt,		// query point
  		  nn,			// number of near neighbors
	  	  nnIdx,		// nearest neighbors (returned)
		  dists,		// distance (returned)
		  eps);						

          for(j=0; j<nn; j++)
          {
            *(Idx+j*nq+i) = nnIdx[j]+1;
          }
        }
        annDeallocPt(queryPt);

        delete [] nnIdx;
        delete [] dists;
      }


      if(*iflag==2)
      {
        delete kdTree;
        annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }



  int committee1_(double* query, double* vertex, int* nPts, int* Idx, int* k, int* dim1, int *iflag, int *mproc, int *imb3d)
  {

      int i,j,nq=1;
      int dim=*dim1;
      int myrank=*mproc;
      int maxPts = *nPts;
      int nn = *k;
      int np = *nPts;
      double eps=0;

      static ANNpointArray dataPts;			// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;					

      if(*iflag==0)
      {
        dataPts = annAllocPts(np,dim);	// allocate data points
        for (i=0; i<np; i++)
        {
          for (j=0; j<dim; j++)
	  {
            dataPts[i][j] = *(vertex+i*dim+j);
	  }
        }
        kdTree = new ANNkd_tree(dataPts,np,dim);
      }

      if(*iflag==1)
      {
        nnIdx = new ANNidx[nn];			// allocate near neigh indices
        dists = new ANNdist[nn];		// alloc

        queryPt = annAllocPt(dim);		// allocate query point
	queryPt[0] = *(query+0);
	queryPt[1] = *(query+1)*(*imb3d);
	queryPt[2] = *(query+2);
	kdTree->annkSearch(		// search
		  queryPt,		// query point
  		  nn,			// number of near neighbors
	  	  nnIdx,		// nearest neighbors (returned)
		  dists,		// distance (returned)
		  eps);						

        for(j=0; j<nn; j++)
        {
          *(Idx+j) = nnIdx[j]+1;
        }

        annDeallocPt(queryPt);
	delete [] nnIdx;
        delete [] dists;
      }

      if(*iflag==2)
      {
        delete kdTree;
        annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }

  int committee2_(double* query, double* vertex, int* nPts, int* Idx, int* nn, int* nnmax1, int *dim1, int *iflag, int *mproc, double *radius)
  {

      int i,j,nq,np,ntemp,nnh;
      int myrank = *mproc;
      int maxPts = *nPts;
      int nnmax = *nnmax1;
      int dim = *dim1;
      double rad = *radius;
      double eps=0;

      np = *nPts;
      nq = 1;

      static ANNpointArray dataPts;			// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;					

      if(*iflag==0)
      {

        dataPts = annAllocPts(maxPts, dim);	// allocate data points
        for (i=0; i<*nPts; i++)
        {
          for (j=0; j<dim; j++)
	  {
            dataPts[i][j] = *(vertex+i*dim+j);
	  }
        }

        kdTree = new ANNkd_tree(	// build search structure
		dataPts,		// the data points
		*nPts,			// number of points
		dim);			// dimension of space
      }

      
      if(*iflag==1)
      {
        nnIdx = new ANNidx[nnmax];		// allocate near neigh indices
        dists = new ANNdist[nnmax];		// alloc

        queryPt = annAllocPt(dim);	// allocate query point
        for (j=0; j<dim; j++)
	{
          queryPt[j] = *(query+j);
        }

        nnh = kdTree->annkFRSearch(	// search
		  queryPt,		// query point
                  rad,                  // radius
  		  nnmax,                // maximum number of near neighbors
                  nnIdx);		// array of near neighbors
        *(nn) = nnh;

        for(j=0; j<min(nnh,nnmax); j++)
        {
          *(Idx+j) = nnIdx[j]+1;
        }
        annDeallocPt(queryPt);

	delete [] nnIdx;
        delete [] dists;
      }

      
      if(*iflag==2)
      {
        delete kdTree;
        annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }


  int committee3(double* query, double* vertex, int* nPts, int* Idx, int* k, int *iflag, int *mproc, int *i3d)
  {

      int i,j;
      int dim=3;
      int myrank=*mproc;
      int maxPts = *nPts;
      double imb3d = (double)(*i3d);
      int np = *nPts;
      int nn = *k;
      double eps=0.0;

      static ANNpointArray dataPts;		// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;				

      if(*iflag==0)
      {
        dataPts = annAllocPts(np, dim);	// allocate data points
        for (i=0; i<np; i++)
        {
	  dataPts[i][0] = *(vertex+i*dim+0);
	  dataPts[i][1] = *(vertex+i*dim+1)*imb3d;
	  dataPts[i][2] = *(vertex+i*dim+2);	  
        }

        kdTree = new ANNkd_tree(dataPts,np,dim);
      }

      if(*iflag==1)
      {
        nnIdx = new ANNidx[nn];
        dists = new ANNdist[nn];
        queryPt = annAllocPt(dim);
	queryPt[0] = *(query+0);
	queryPt[1] = *(query+1)*imb3d;
	queryPt[2] = *(query+2);
	kdTree->annkSearch(queryPt,nn,nnIdx,dists,eps);						

	for(j=0; j<nn; j++)
        {
          *(Idx+j) = nnIdx[j]+1;
        }
        annDeallocPt(queryPt);

        delete [] nnIdx;
        delete [] dists;
      }

      if(*iflag==2)
      {
	//	cout << "iflag=" << *iflag << "\n";
        delete kdTree;
	annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }


  int committee3_(double* query, double* vertex, int* nPts, int* Idx, int* k, int *iflag, int *mproc, int *i3d)
  {

      int i,j;
      int dim=3;
      int myrank=*mproc;
      int maxPts = *nPts;
      double imb3d = (double)(*i3d);
      int np = *nPts;
      int nn = *k;
      double eps=0.0;

      static ANNpointArray dataPts;		// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;				

      if(*iflag==0)
      {
        dataPts = annAllocPts(np, dim);	// allocate data points
        for (i=0; i<np; i++)
        {
	  dataPts[i][0] = *(vertex+i*dim+0);
	  dataPts[i][1] = *(vertex+i*dim+1)*imb3d;
	  dataPts[i][2] = *(vertex+i*dim+2);	  
        }

        kdTree = new ANNkd_tree(dataPts,np,dim);
      }

      if(*iflag==1)
      {
        nnIdx = new ANNidx[nn];
        dists = new ANNdist[nn];
        queryPt = annAllocPt(dim);
	queryPt[0] = *(query+0);
	queryPt[1] = *(query+1)*imb3d;
	queryPt[2] = *(query+2);
	kdTree->annkSearch(queryPt,nn,nnIdx,dists,eps);						

	for(j=0; j<nn; j++)
        {
          *(Idx+j) = nnIdx[j]+1;
        }
        annDeallocPt(queryPt);

        delete [] nnIdx;
        delete [] dists;
      }

      if(*iflag==2)
      {
	//	cout << "iflag=" << *iflag << "\n";
        delete kdTree;
	annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }


  int committee4(double* query, double* vertex, int* nPts, int* Idx, int* nn, int* nnmax1, int *dim1, int *iflag, int *mproc, double *radius)
  {

      int i,j,nq,np,ntemp,nnh;
      int myrank = *mproc;
      int maxPts = *nPts;
      int nnmax = *nnmax1;
      int dim = *dim1;
      double rad = *radius;
      double eps=0;

      np = *nPts;
      nq = 1;

      static ANNpointArray dataPts;			// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;					

      if(*iflag==0)
      {

        dataPts = annAllocPts(maxPts, dim);	// allocate data points
        for (i=0; i<*nPts; i++)
        {
          for (j=0; j<dim; j++)
	  {
            dataPts[i][j] = *(vertex+i*3+j);
	  }
        }

        kdTree = new ANNkd_tree(	// build search structure
		dataPts,		// the data points
		*nPts,			// number of points
		dim);			// dimension of space
      }

      
      if(*iflag==1)
      {
        nnIdx = new ANNidx[nnmax];		// allocate near neigh indices
        dists = new ANNdist[nnmax];		// alloc

        queryPt = annAllocPt(dim);	// allocate query point
        for (j=0; j<dim; j++)
	{
          queryPt[j] = *(query+j);
        }

        nnh = kdTree->annkFRSearch(	// search
		  queryPt,		// query point
                  rad,                  // radius
  		  nnmax,                // maximum number of near neighbors
                  nnIdx);		// array of near neighbors
        *(nn) = nnh;

        for(j=0; j<min(nnh,nnmax); j++)
        {
          *(Idx+j) = nnIdx[j]+1;
        }
        annDeallocPt(queryPt);

	delete [] nnIdx;
        delete [] dists;
      }

      
      if(*iflag==2)
      {
        delete kdTree;
        annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }


  int committee4_(double* query, double* vertex, int* nPts, int* Idx, int* nn, int* nnmax1, int *dim1, int *iflag, int *mproc, double *radius)
  {

      int i,j,nq,np,ntemp,nnh;
      int myrank = *mproc;
      int maxPts = *nPts;
      int nnmax = *nnmax1;
      int dim = *dim1;
      double rad = *radius;
      double eps=0;

      np = *nPts;
      nq = 1;

      static ANNpointArray dataPts;			// data points
      ANNpoint	    queryPt;			// query point
      ANNidxArray   nnIdx;			// near neighbor indices
      ANNdistArray  dists;			// near neighbor distances
      static ANNkd_tree*   kdTree;					

      if(*iflag==0)
      {

        dataPts = annAllocPts(maxPts, dim);	// allocate data points
        for (i=0; i<*nPts; i++)
        {
          for (j=0; j<dim; j++)
	  {
            dataPts[i][j] = *(vertex+i*3+j);
	  }
        }

        kdTree = new ANNkd_tree(	// build search structure
		dataPts,		// the data points
		*nPts,			// number of points
		dim);			// dimension of space
      }

      
      if(*iflag==1)
      {
        nnIdx = new ANNidx[nnmax];		// allocate near neigh indices
        dists = new ANNdist[nnmax];		// alloc

        queryPt = annAllocPt(dim);	// allocate query point
        for (j=0; j<dim; j++)
	{
          queryPt[j] = *(query+j);
        }

        nnh = kdTree->annkFRSearch(	// search
		  queryPt,		// query point
                  rad,                  // radius
  		  nnmax,                // maximum number of near neighbors
                  nnIdx);		// array of near neighbors
        *(nn) = nnh;

        for(j=0; j<min(nnh,nnmax); j++)
        {
          *(Idx+j) = nnIdx[j]+1;
        }
        annDeallocPt(queryPt);

	delete [] nnIdx;
        delete [] dists;
      }

      
      if(*iflag==2)
      {
        delete kdTree;
        annDeallocPts(dataPts);
      }
      annClose();

      return(1);

  }

}
