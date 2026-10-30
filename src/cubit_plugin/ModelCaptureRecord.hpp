#ifndef CUBIT_MODEL_CAPTURE_RECORD_HPP
#define CUBIT_MODEL_CAPTURE_RECORD_HPP

#include <array>
#include <cmath>
#include <iomanip>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cubit_mesh_export {
using CapturePoint = std::array<double, 3>;
struct CapturedVertex { int id; CapturePoint position; };
struct CapturedCurve {
  int id; std::string type; double length; int spline_degree = 0;
  std::vector<CapturedVertex> vertices;
  std::vector<CapturePoint> samples;
};
struct CapturedFace {
  int id; std::string type; CapturePoint normal;
  std::vector<int> curves;
};
struct CapturedVolume {
  int id; std::string name; double measure;
  std::vector<CapturedFace> surfaces;
  bool whole_sphere = false;
  CapturePoint center{}; double radius = 0;
};
struct CapturedBlock { int id; std::string name; int volume; std::vector<int> surfaces; };
struct CapturedSideset { int id; std::string name; std::vector<int> blocks, surfaces; };
struct ModelCaptureRecord {
  std::vector<CapturedVolume> volumes;
  std::map<int, CapturedCurve> curves;
  std::vector<CapturedBlock> blocks;
  std::vector<CapturedSideset> sidesets;
};

inline std::string capture_quote(const std::string &value)
{
  std::ostringstream out; out << '"';
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') out << '\\' << char(c);
    else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else out << char(c);
  }
  out << '"'; return out.str();
}
inline void capture_number(std::ostream &out, double value)
{
  if (!std::isfinite(value))
    throw std::invalid_argument("Model capture requires finite coordinates/properties");
  out << value;
}
inline void capture_point(std::ostream &out, const CapturePoint &point)
{
  out << '['; for (int i=0; i<3; ++i) { if(i) out << ','; capture_number(out,point[i]); } out << ']';
}
inline void capture_ids(std::ostream &out, const std::vector<int> &ids)
{
  out << '['; bool comma=false;
  for (int id : ids) { if(id<1) throw std::invalid_argument("Positive source IDs required");
    if(comma) out << ','; out << id; comma=true; }
  out << ']';
}
inline std::string serialize_model_capture(const ModelCaptureRecord &record)
{
  if (record.volumes.empty() || record.blocks.size()!=record.volumes.size())
    throw std::invalid_argument("Model capture requires complete one-volume blocks");
  std::ostringstream out; out.imbue(std::locale::classic()); out << std::setprecision(17);
  out << "{\"schema\":\"cubit.model-capture.v1\",\"units\":\"mm\",\"curves\":{";
  bool comma=false;
  for (const auto &entry : record.curves) {
    const auto &c=entry.second;
    if(c.id!=entry.first || c.id<1 || c.samples.size()!=33)
      throw std::invalid_argument("Curve identity and 33 explicit samples required");
    if(comma) out << ','; comma=true;
    out << capture_quote(std::to_string(c.id)) << ":{\"id\":" << c.id << ",\"type\":" << capture_quote(c.type)
        << ",\"length\":"; capture_number(out,c.length);
    if(c.type=="Spline curve") { if(c.spline_degree<1) throw std::invalid_argument("Explicit spline degree required");
      out << ",\"spline_degree\":" << c.spline_degree; }
    out << ",\"vertices\":["; bool sub=false;
    for(const auto &v:c.vertices) { if(sub)out<<',';sub=true;
      if(v.id<1)throw std::invalid_argument("Positive vertex ID required");
      out<<"{\"id\":"<<v.id<<",\"position\":";capture_point(out,v.position);out<<'}'; }
    out << "],\"samples\":[";sub=false;
    for(const auto &p:c.samples){if(sub)out<<',';sub=true;capture_point(out,p);}out<<"]}";
  }
  out << "},\"volumes\":[";comma=false;std::set<int> volumes;
  for(const auto &v:record.volumes){
    if(v.id<1 || v.name.empty() || !volumes.insert(v.id).second || v.surfaces.empty())
      throw std::invalid_argument("Unique named nonempty volumes required");
    if(comma)out<<',';comma=true;
    out<<"{\"id\":"<<v.id<<",\"name\":"<<capture_quote(v.name)<<",\"volume\":";
    capture_number(out,v.measure);out<<",\"surfaces\":[";bool sub=false;
    for(const auto &f:v.surfaces){if(sub)out<<',';sub=true;
      out<<"{\"id\":"<<f.id<<",\"type\":"<<capture_quote(f.type)<<",\"normal\":";
      capture_point(out,f.normal);out<<",\"curves\":[";bool edge=false;
      for(int id:f.curves){if(!record.curves.count(id))throw std::invalid_argument("Face references missing curve "+std::to_string(id));
        if(edge)out<<',';edge=true;out<<"{\"id\":"<<id<<'}';}out<<"]}";
    }out<<']';
    if(v.whole_sphere){if(v.surfaces.size()!=1 || v.surfaces[0].type!="sphere surface" || v.radius<=0)
        throw std::invalid_argument("Whole-sphere source declaration required");
      out<<",\"primitive\":{\"type\":\"sphere\",\"center_mm\":";capture_point(out,v.center);
      out<<",\"radius_mm\":";capture_number(out,v.radius);out<<'}';}out<<'}';
  }
  out<<"],\"labels\":{\"schema\":\"cubit.source-face-labels.v1\",\"blocks\":[";comma=false;
  std::set<int> block_ids,covered_volumes;
  for(const auto &b:record.blocks){
    if(b.id<1 || !block_ids.insert(b.id).second || !volumes.count(b.volume) || !covered_volumes.insert(b.volume).second)
      throw std::invalid_argument("Unique source block-volume coverage required");
    if(comma)out<<',';comma=true;
    out<<"{\"id\":"<<b.id<<",\"name\":"<<capture_quote(b.name)<<",\"source_volumes\":["<<b.volume
       <<"],\"source_surface_ids\":";capture_ids(out,b.surfaces);out<<'}';
  }
  out<<"],\"sidesets\":[";comma=false;std::set<int> sidesets;
  for(const auto &s:record.sidesets){
    if(s.id<1 || !sidesets.insert(s.id).second || s.name.empty() || s.blocks.empty() || s.surfaces.empty())
      throw std::invalid_argument("Unique named geometry sidesets required");
    for(int id:s.blocks)if(!block_ids.count(id))throw std::invalid_argument("Sideset references missing block");
    if(comma)out<<',';comma=true;
    out<<"{\"id\":"<<s.id<<",\"name\":"<<capture_quote(s.name)<<",\"source_blocks\":";capture_ids(out,s.blocks);
    out<<",\"source_surface_ids\":";capture_ids(out,s.surfaces);
    out<<",\"selection\":\"source_surface_subset\",\"physics\":\"selection_only\"}";
  }
  out<<"]},\"geometry_accuracy_accepted\":false,\"physical_conditions_inferred\":false}";
  return out.str();
}
} // namespace cubit_mesh_export
#endif
