// Guards: source-record completeness and JSON emission. Assumes: synthetic CAD records.
#include "ModelCaptureRecord.hpp"
#include <iostream>
#include <limits>
using namespace cubit_mesh_export;
template<class F> void rejects(F fn,const char *name)
{
  bool rejected=false;try{fn();}catch(const std::exception &){rejected=true;}
  if(!rejected)throw std::runtime_error(std::string("Missing capture rejection: ")+name);
}
int main() try
{
  ModelCaptureRecord record;
  require_polynomial_capture_spline(false,31);
  rejects([&]{require_polynomial_capture_spline(true,31);},"rational spline");
  CapturedVolume v;v.id=7;v.name="Body";v.measure=4.1887902047863905;
  v.surfaces.push_back({11,"sphere surface",{0,0,1},{}});
  v.whole_sphere=true;v.center={1,2,3};v.radius=1;
  record.volumes.push_back(v);record.blocks.push_back({17,"Body",7,{11}});
  record.sidesets.push_back({23,"Outer\"wall",{17},{11}});
  auto bad=record;bad.blocks.clear();rejects([&]{serialize_model_capture(bad);},"missing block");
  bad=record;bad.blocks[0].volume=999;rejects([&]{serialize_model_capture(bad);},"unknown volume");
  bad=record;bad.volumes[0].center[0]=std::numeric_limits<double>::quiet_NaN();
  rejects([&]{serialize_model_capture(bad);},"nonfinite point");
  bad=record;bad.volumes[0].surfaces[0].curves.push_back(999);
  rejects([&]{serialize_model_capture(bad);},"missing curve");
  bad=record;bad.sidesets[0].blocks={999};rejects([&]{serialize_model_capture(bad);},"unknown block");
  std::cout<<serialize_model_capture(record)<<'\n';
  return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
