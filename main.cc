#include "miniapp.hh"

#include <string>
#include <sstream>
#include <cstdio>

using namespace std;

int main(int argc, char* argv[])
{
  // int n_eg_in, int n_a_in, int cm_xy_in, int fm_xy_in, int cm_z_in, int fm_z_in, int upscatter_in, int iter_in

  int args[12];
  std::string order;

  if (argc != 14)
  {
    cout << "use:  3D_KBA_Sweep order ng na cm_xy fm_xy cm_z fm_z up iter xbs ybs zbs nt" << endl;
    return 0;
  }
  else
  {
    order = argv[1];
    for (int i = 0; i < 12; ++i)
    {
      std::string s = argv[i+2];
      if (!(std::istringstream(s) >> args[i]))
      {
        cout << "invalid argument (i=" << i <<"). quitting!" << endl;
        return 1;
      }
      else
      {
        cout << " argument i = " << i << " is " << args[i] << endl;
      }
    }
    printf("order=%s, ng=%i, na=%i, cm_xy=%i, fm_xy=%i, cm_z=%i, fm_z=%i, up=%i, iter=%i, xbs=%i, ybs=%i, zbs=%i, nt=%i\n",
            order.c_str(),args[0],args[1],args[2],args[3],args[4],args[5],args[6],
            args[7],args[8],args[9],args[10],args[11]);
  }
  Solver test(order,args[0],args[1],args[2],args[3],args[4],args[5],args[6],
      args[7],args[8],args[9],args[10],args[11]);
  test.Calculate();
  return 0;
}



