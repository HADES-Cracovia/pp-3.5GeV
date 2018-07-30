#ifndef HPARTICLEPLAYERHANDLE_H
#define HPARTICLEPLAYERHANDLE_H

#include <map>
#include <vector>


class HParticlePlayerHandle {

   public:

HParticlePlayerHandle() 
{
   for (int n=1; n<15; ++n)
   {
      chooseMfromN(1, n);
      chooseMfromN(2, n);
      chooseMfromN(3, n);
      chooseMfromN(4, n);
      chooseMfromN(5, n);
      chooseMfromN(6, n);
      // here we decide that no more that 6 particles of the same type can be combined
   }
}

std::map< int, std::vector< std::vector<int> > > map1; 
std::map< int, std::vector< std::vector<int> > > map2; // select 2 particles out of int
std::map< int, std::vector< std::vector<int> > > map3; 
std::map< int, std::vector< std::vector<int> > > map4; 
std::map< int, std::vector< std::vector<int> > > map5; 
std::map< int, std::vector< std::vector<int> > > map6; 

   private:

void chooseMfromN(int m, int n)
{
   vector<vector<int> > tempVectVect;
   if (m == 1)
   {
      for (int i=0; i<n; ++i)
      {
          vector<int> temp;
          temp.push_back(i);
          tempVectVect.push_back( temp );
      }
      map1[ 10 + n ] = tempVectVect;
   }
   else if (m == 2)
   {
      for (int i=0; i<n; ++i)
#if COMB_REPETITION == 1
      for (int j=0; j<n; ++j)
      if (i!=j)
#else
      for (int j=i+1; j<n; ++j)
#endif
      {
          vector<int> temp;
          temp.push_back(i);
          temp.push_back(j);
          tempVectVect.push_back( temp );
      }
      map2[ 20 + n ] = tempVectVect;
   }
   else if (m == 3)
   {
      for (int i=0; i<n; ++i)
#if COMB_REPETITION == 1
      for (int j=0; j<n; ++j)
      for (int k=0; k<n; ++k)
      if (i!=j && i!=k && j!=k)
#else
      for (int j=i+1; j<n; ++j)
      for (int k=j+1; k<n; ++k)
#endif
      {
          vector<int> temp;
          temp.push_back(i);
          temp.push_back(j);
          temp.push_back(k);
          tempVectVect.push_back( temp );
      }
      map3[ 30 + n ] = tempVectVect;
   }
   else if (m == 4)
   {
      for (int i=0; i<n; ++i)
#if COMB_REPETITION == 1
      for (int j=0; j<n; ++j)
      for (int k=0; k<n; ++k)
      for (int l=0; l<n; ++l)
      if (i!=j && i!=k && i!=l && j!=k && j!=l && k!=l)
#else
      for (int j=i+1; j<n; ++j)
      for (int k=j+1; k<n; ++k)
      for (int l=k+1; l<n; ++l)
#endif
      {
          vector<int> temp;
          temp.push_back(i);
          temp.push_back(j);
          temp.push_back(k);
          temp.push_back(l);
          tempVectVect.push_back( temp );
      }
      map4[ 40 + n ] = tempVectVect;
   }
   else if (m == 5)
   {
      for (int i=0; i<n; ++i)
#if COMB_REPETITION == 1
      for (int j=0; j<n; ++j)
      for (int k=0; k<n; ++k)
      for (int l=0; l<n; ++l)
      for (int p=0; p<n; ++p)
      if (i!=j && i!=k && i!=l && i!=p && j!=k && j!=l && j!=p && k!=l && k!=p && l!=p)
#else
     for (int j=i+1; j<n; ++j)
     for (int k=j+1; k<n; ++k)
     for (int l=k+1; l<n; ++l)
     for (int p=l+1; p<n; ++p)
#endif
     {
         vector<int> temp;
         temp.push_back(i);
         temp.push_back(j);
         temp.push_back(k);
         temp.push_back(l);
         temp.push_back(p);
         tempVectVect.push_back( temp );
     }
     map5[ 50 + n ] = tempVectVect;
   }
   else if (m == 6)
   {
      for (int i=0; i<n; ++i)
#if COMB_REPETITION == 1
      for (int j=0; j<n; ++j)
      for (int k=0; k<n; ++k)
      for (int l=0; l<n; ++l)
      for (int p=0; p<n; ++p)
      for (int r=0; r<n; ++r)
      if (i!=j && i!=k && i!=l && i!=p && i!=r && j!=k && j!=l && j!=p && j!=r && k!=l && k!=p && k!=r && l!=p && l!=r && p!=r)
#else
     for (int j=i+1; j<n; ++j)
     for (int k=j+1; k<n; ++k)
     for (int l=k+1; l<n; ++l)
     for (int p=l+1; p<n; ++p)
     for (int r=p+1; r<n; ++r)
#endif
     {
         vector<int> temp;
         temp.push_back(i);
         temp.push_back(j);
         temp.push_back(k);
         temp.push_back(l);
         temp.push_back(p);
         temp.push_back(r);
         tempVectVect.push_back( temp );
     }
     map6[ 60 + n ] = tempVectVect;
   }
}

}; // eof class HParticlePlayerHandle

#endif // HPARTICLEPLAYERHANDLE_H
