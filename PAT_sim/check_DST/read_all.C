#include <iostream>
#include <string>
#include <fstream>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

int read_all()
{
  
  TChain* chain=new TChain();

  ifstream in("lista.list");

  if(!in)
    {
      cout << "Cannot open input file.\n";
      return 1;
    }
  else
    cout<<"File lista.list opened"<<endl;

  char line[255];

  while(in)
    {
      in.getline(line, 255);  // delim defaults to '\n'
      if(in)
	cout << line << endl;
      /*
      string line1(line);
      if(line1.substr(line1.length()-5) == ".root")
	chain->AddFile(line);
      */    
    }

  in.close();

  //return 0;
}
