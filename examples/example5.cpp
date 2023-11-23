#include "Sn.hpp"
#include "SnIrreducible.hpp"
#include <iostream>

int main() {
  
  Sn::Sn G(5);

  Sn::Sn& Gsub=*G.subgroup->subgroup;

  for (unsigned i=0; i<Gsub.irreducibles.size(); i++)
    cout<<Gsub.irreducibles[i]->str()<<endl;
      
}
