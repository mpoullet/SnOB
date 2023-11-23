#include "Sn.hpp"
#include "SnIrreducible.hpp"
#include <iostream>

int main() {
  
  Sn::Sn G(5);

  for (unsigned i=0; i<G.irreducibles.size(); i++)
    cout<<G.irreducibles[i]->str()<<endl;
      
}
