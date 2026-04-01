//
// Created by joseph on 2026/1/3.
//
#include "specie.H"
#include "PengRobinsonGas.H"
#include "RedlichKwongGas.H"

#include <cstdlib>
#include <vector>

using namespace Foam;

template<class EOS_t>
void test_EOS(const word& EOS, const specie& sp, const realFluidProperty&prop) {

  const std::vector<scalar> T_range{73.15,493.15,913.15,1333.15,1753.15,2173.15,2593.15,3013.15,3433.15,3853.15,4273.15};
  const std::vector<scalar> p_range{1e3,1e4,1e5,1e6,1e7,1e8,1e9};

  const scalar small = 1e-10;

  const EOS_t eos{sp, prop};

  for (scalar T: T_range) {
    for (scalar p: p_range) {
      const scalar rho = eos.rho(p, T);
      if (rho<=0) {
        Warning<<EOS<<" predicted non-physical density at p = "<<p<<"[Pa], T = "<<T<<"[K]"<<endl;
        exit(1);
      }
      {
        const auto core =  eos.core(p, T);
        const scalar z = core.Z(p, T);
        // const scalar dpdT_V = core.dpdT_V(p, T, z);
        const scalar dpdV_T = core.dpdV_T(p, T, z);
        const scalar dVdp_T = core.dVdp_T(p, T, z);
        const scalar dVdT_p = core.dVdT_p(p, T, z);
        const scalar dTdp_V = core.dTdp_V(p, T, z);

        scalar error=0;
        error=dpdV_T*dVdp_T - scalar(1);
        if (std::abs(error)>small) {
          Warning<<EOS<<" violates swapping rule: dpdV_T = "<<dpdV_T<<", dVdp_T = "<<dVdp_T<<", error = "<<error<<endl;
          exit(1);
        }
        error = dpdV_T * dVdT_p * dTdp_V + 1;
        if (std::abs(error)>small) {
          Warning<<EOS<<" violates cyclic rule: dpdV_T = "<<dpdV_T<<", dVdT_p = "<<dVdT_p<<", dTdp_V = "<<dTdp_V<<", error = "<<error<<endl;
          exit(1);
        }
      }
    }
  }

}

int main(int argc, char** argv) {
  const specie sp{"H2",1,2.016};

  const realFluidProperty real_fluid_prop = []() {
    realFluidProperty prop;
    prop.Tc_=33.18;
    prop.Pc_=1292977.4;
    prop.Vc_=0.06503;
    prop.omega_=-0.22014;
    return prop;
  }();

  test_EOS<RedlichKwongGas<specie>>("RedlichKwongGas",sp,real_fluid_prop);
  test_EOS<PengRobinsonGas<specie>>("PengRobinsonGas",sp,real_fluid_prop);

  Info<<"All test passed."<<endl;
  return 0;

}