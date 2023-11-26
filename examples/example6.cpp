#include "Sn.hpp"
#include "SnIrreducible.hpp"
#include <iostream>

string printAncestors(Sn::Irreducible& rho, string indenter){
  ostringstream result;
  result<<indenter<<rho.str()<<endl;
  for (unsigned i=0; i<rho.eta.size(); i++){
    result<<printAncestors(*rho.eta[i],indenter+"  ");
  }
  return result.str();
}


int main() {
  
  Sn G(5);

  for (unsigned i=0; i<G.irreducibles.size(); i++)
    cout<<printAncestors(*G.irreducibles[i],"");
      
}

