#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include "ExportModelCaptureCommand.hpp"
#include "ModelCaptureRecord.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include "CubitRefAccessor.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "Curve.hpp"
#include "Surface.hpp"
#include "DLIList.hpp"
#include "CubitVector.hpp"
#include "utf8_path.hpp"
#include <algorithm>
#include <cerrno>
#include <filesystem>
#include <system_error>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

namespace {
using namespace cubit_mesh_export;
std::vector<int> sorted_ids(std::vector<int> ids)
{
  std::sort(ids.begin(),ids.end());
  if(std::adjacent_find(ids.begin(),ids.end())!=ids.end() || (!ids.empty() && ids.front()<1))
    throw std::runtime_error("Duplicate or nonpositive source geometry ID");
  return ids;
}
CapturePoint point(const CubitVector &v) { return {v.x(),v.y(),v.z()}; }
CapturedCurve capture_curve(int id)
{
  RefEdge *edge=cubit_get_ref_edge(id);
  if(!edge)throw std::runtime_error("Missing source curve "+std::to_string(id));
  CapturedCurve c;c.id=id;c.type=CubitInterface::get_curve_type(id);c.length=edge->measure();
  if(c.type!="Straight curve" && c.type!="Arc curve" && c.type!="Spline curve")
    throw std::runtime_error("Unsupported capture curve "+std::to_string(id)+": "+c.type);
  for(int vid:sorted_ids(CubitInterface::get_relatives("curve",id,"vertex")))
    c.vertices.push_back({vid,CubitInterface::get_center_point("vertex",vid)});
  for(int i=0;i<=32;++i){CubitVector p;
    if(edge->position_from_fraction(double(i)/32,p)!=CUBIT_SUCCESS)
      throw std::runtime_error("Cannot sample source curve "+std::to_string(id));
    c.samples.push_back(point(p));}
  if(c.type=="Spline curve"){
    bool rational=false,reversed=false;DLIList<CubitVector> controls;DLIList<double> weights,knots;
    if(edge->get_curve_ptr()->get_spline_params(rational,c.spline_degree,controls,weights,knots,reversed)!=CUBIT_SUCCESS)
      throw std::runtime_error("Cannot read source spline degree for curve "+std::to_string(id));
  }
  return c;
}
ModelCaptureRecord capture()
{
  ModelCaptureRecord r;
  const auto volumes=sorted_ids(CubitInterface::parse_cubit_list("volume","all"));
  if(volumes.empty() || volumes.size()>10000)throw std::runtime_error("Nonempty bounded CAD volumes required");
  std::map<int,CapturedBlock> by_volume;std::set<std::string> block_names,sideset_names;
  for(int bid:sorted_ids(CubitInterface::parse_cubit_list("block","all"))){
    auto ids=sorted_ids(CubitInterface::get_block_volumes(bid));
    if(ids.size()!=1)throw std::runtime_error("Model capture requires one explicit CAD volume in block "+std::to_string(bid));
    CapturedBlock b{bid,CubitInterface::get_block_name(bid),ids[0],
                   sorted_ids(CubitInterface::get_relatives("volume",ids[0],"surface"))};
    if(b.name.empty() || !block_names.insert(b.name).second || !by_volume.emplace(b.volume,b).second)
      throw std::runtime_error("Explicit unique block names and volume ownership required");
    r.blocks.push_back(b);
  }
  if(by_volume.size()!=volumes.size())throw std::runtime_error("Every CAD volume requires exactly one source block");
  for(int vid:volumes){
    auto block=by_volume.find(vid);if(block==by_volume.end())throw std::runtime_error("Unblocked volume "+std::to_string(vid));
    CapturedVolume v;v.id=vid;v.name=block->second.name;v.measure=CubitInterface::get_volume_volume(vid);
    if(v.measure<=0)throw std::runtime_error("Nonpositive source CAD volume "+std::to_string(vid));
    for(int sid:block->second.surfaces){
      CapturedFace f;f.id=sid;f.type=CubitInterface::get_surface_type(sid);
      f.normal=CubitInterface::get_surface_normal(sid);
      f.curves=sorted_ids(CubitInterface::get_relatives("surface",sid,"curve"));
      for(int cid:f.curves)if(!r.curves.count(cid))r.curves.emplace(cid,capture_curve(cid));
      v.surfaces.push_back(f);
    }
    if(v.surfaces.size()==1 && v.surfaces[0].type=="sphere surface"){
      RefFace *face=cubit_get_ref_face(v.surfaces[0].id);
      CubitVector center,pole,origin;
      if(!face || face->get_surface_ptr()->get_sphere_params(center,v.radius,pole,origin)!=CUBIT_SUCCESS)
        throw std::runtime_error("Cannot read exact sphere parameters for volume "+std::to_string(vid));
      v.whole_sphere=true;v.center=point(center);
    }
    r.volumes.push_back(v);
  }
  for(int sid:sorted_ids(CubitInterface::parse_cubit_list("sideset","all"))){
    CapturedSideset s;s.id=sid;s.name=CubitInterface::get_exodus_entity_name("sideset",sid);
    s.surfaces=sorted_ids(CubitInterface::get_sideset_surfaces(sid));
    if(s.name.empty() || !sideset_names.insert(s.name).second || block_names.count(s.name) || s.surfaces.empty())
      throw std::runtime_error("Named nonempty CAD-face sidesets required; mesh-only sets are unsupported");
    std::set<int> covered;
    for(const auto &b:r.blocks){bool touched=false;
      for(int f:s.surfaces)if(std::binary_search(b.surfaces.begin(),b.surfaces.end(),f)){covered.insert(f);touched=true;}
      if(touched)s.blocks.push_back(b.id);}
    if(covered.size()!=s.surfaces.size())throw std::runtime_error("Sideset references uncaptured surface "+std::to_string(sid));
    r.sidesets.push_back(s);
  }
  return r;
}
void publish(const std::filesystem::path &path,const std::string &bytes)
{
  std::filesystem::path partial=path;partial+=".partial.json";
  if(std::filesystem::exists(path))throw std::runtime_error("Capture target already exists");
#ifdef _WIN32
  HANDLE file=CreateFileW(partial.c_str(),GENERIC_WRITE,0,nullptr,CREATE_NEW,FILE_ATTRIBUTE_NORMAL,nullptr);
  if(file==INVALID_HANDLE_VALUE)throw std::system_error(GetLastError(),std::system_category(),"Cannot create fresh capture partial");
  bool okay=true;std::size_t offset=0;
  while(offset<bytes.size()){
    const auto chunk=static_cast<DWORD>(std::min<std::size_t>(bytes.size()-offset,1u<<20));DWORD done=0;
    if(!WriteFile(file,bytes.data()+offset,chunk,&done,nullptr) || done!=chunk){okay=false;break;}offset+=done;
  }
  if(!FlushFileBuffers(file))okay=false;
  if(!CloseHandle(file))okay=false;
  if(!okay)throw std::runtime_error("Capture write failed; own partial retained: "+partial.u8string());
  if(!MoveFileExW(partial.c_str(),path.c_str(),MOVEFILE_WRITE_THROUGH))
    throw std::system_error(GetLastError(),std::system_category(),"Capture no-replace publish failed; own partial retained");
#else
  int fd=open(partial.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
  if(fd<0)throw std::system_error(errno,std::generic_category(),"Cannot create fresh capture partial");
  std::size_t offset=0;bool okay=true;
  while(offset<bytes.size()){ssize_t n=write(fd,bytes.data()+offset,bytes.size()-offset);
    if(n<=0){okay=false;break;}offset+=static_cast<std::size_t>(n);}
  if(fsync(fd))okay=false;if(close(fd))okay=false;
  if(!okay)throw std::runtime_error("Capture write failed; own partial retained");
  if(link(partial.c_str(),path.c_str()))throw std::system_error(errno,std::generic_category(),"Capture no-replace publish failed; own partial retained");
  if(unlink(partial.c_str()))throw std::system_error(errno,std::generic_category(),"Capture published but own partial remains");
#endif
}
} // namespace

std::vector<std::string> ExportModelCaptureCommand::get_syntax()
{ return {"export model_capture <string:label='filename',help='<filename>'> units mm"}; }
std::vector<std::string> ExportModelCaptureCommand::get_syntax_help()
{ return {"export model_capture \"filename.json\" units mm"}; }
std::vector<std::string> ExportModelCaptureCommand::get_help()
{ return {"Capture CAD topology, explicit block/sideset identity and 33 curve samples.\n"
          "The caller must explicitly declare model length units as mm. No unit inference,\n"
          "meshing, material defaults or physical boundary inference. One CAD volume per\n"
          "named block; CAD-face sidesets only; existing targets are never replaced.\n"
          "Sampled curves are source observations, not a general exact CAD conversion.\n"}; }
bool ExportModelCaptureCommand::execute(CubitCommandData &data)
{
  try {
    std::string filename;data.get_string("filename",filename);
    if(filename.empty())throw std::invalid_argument("Explicit capture filename required");
    const auto bytes=cubit_mesh_export::serialize_model_capture(capture());
    publish(u8_string_to_path(filename),bytes);
    PRINT_INFO("Neutral CAD model capture written. Physical values were not inferred.\n");
    return true;
  }catch(const std::exception &e){PRINT_ERROR("Model capture failed: %s\n",e.what());return false;}
}
