#include "miniapp.hh"

#include <string>
#include <sstream>
#include <cstdio>

using namespace std;

int main(int argc, char* argv[])
{
  int args[9];
  string order = "sea";
  if (argc != 11)
  {
    cout << "usage: " << argv[0] << " order ng na cm fm cmz fmz up iter nt" << endl;
    cout <<"  where order = sae|sea|ase|aes|esa|eas" << endl
        << "           ng = number groups" << endl
        << "           na = number angles per octant" << endl
        << "           cm = coarse meshes along x and y" << endl
        << "           fm = fine mesh per x/y coarse mesh" << endl
        << "          cmz = coarse meshes along z" << endl
        << "          fmz = fine mesh per z coarse mesh" << endl
        << "           up = upscatter groups" << endl
        << "         iter = number of iterations (i.e., sweeps)" << endl
        << "           nt = number of threads" << endl;
    return 0;
  }
  else
  {
    int ng, na, cm, fm, cmz, fmz, up, iter, nt;
    order = argv[1];
    for (int i = 0; i < 9; ++i)
    {
      std::string s = argv[i+2];
      if (!(std::istringstream(s) >> args[i]))
      {
        cout << "invalid argument. quitting!" << endl;
        return 1;
      }
    }
    printf("order=%s, ng=%i, na=%i, cm=%i, fm=%i, cmz=%i, fmz=%i, up=%i, iter=%i, nt=%i\n",
           order.c_str(),args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8]);
  }
  Solver test(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7]);
  test.Calculate(order, args[8], 1);
  return 0;
}



