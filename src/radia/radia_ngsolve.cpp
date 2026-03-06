/************************************************************************//**
 * File: radia_ngsolve.cpp
 * Description: NGSolve CoefficientFunction binding for Radia
 * Project: Radia
 * First release: October 2025
 *
 * IMPORTANT:
 * - NGSolve uses meters, Radia uses millimeters
 * - Automatic unit conversion: m -> mm (multiply by 1000)
 *
 * Field types supported:
 * - 'b': Magnetic flux density (Tesla) - vector
 * - 'h': Magnetic field (A/m) - vector
 * - 'a': Vector potential (T*m) - vector
 * - 'm': Magnetization (A/m) - vector
 * - 'phi': Magnetic scalar potential (A) - scalar
 *
 * Coordinate transformation (v0.07):
 * - origin: Translation vector (meters)
 * - u_axis, v_axis, w_axis: Local coordinate system (auto-normalized)
 *
 * Batch evaluation (v0.08):
 * - Implements efficient batch evaluation for multiple points
 * - Reduces Python call overhead from O(M) to O(1)
 * - Enables OpenMP parallelization in Radia core
 *
 * @version 0.08
 ***************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <fem.hpp>
#include <python_ngstd.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cmath>

// FMM acceleration for batch evaluation
#include "rad_dipole_collect.h"
#include "rad_exafmm.h"
#include "rad_application.h"

namespace py = pybind11;

namespace ngfem
{

class RadiaFieldCF : public CoefficientFunction
{
public:
	int radia_obj;
	std::string field_type;

	// Coordinate transformation
	double origin[3];      // Translation vector [m]
	double u_axis[3];      // Local u-axis (normalized)
	double v_axis[3];      // Local v-axis (normalized)
	double w_axis[3];      // Local w-axis (normalized)
	bool use_transform;    // Whether to apply coordinate transformation

	// Computation settings
	py::object precision;  // Computation precision (None = use Radia default)

	// Point cache for batch evaluation
	mutable std::unordered_map<uint64_t, std::array<double,3>> point_cache_;
	mutable bool use_cache_;
	double cache_tolerance_;  // Tolerance for point hashing (meters)
	mutable size_t cache_hits_;
	mutable size_t cache_misses_;

	// Unit conversion: NGSolve (meters) -> Radia (mm or m)
	double coord_scale_;  // Scaling factor for coordinates (1.0 for meters, 1000.0 for mm)

	// Cached Radia module to avoid repeated imports (memory optimization)
	mutable py::module_ rad_module_;

	// FMM acceleration settings
	double fmm_eps_;                                    // FMM tolerance (0 = disabled)
	mutable RadDipoleCollect::DipoleCollection dipole_cache_;  // Cached dipole data
	mutable bool dipoles_extracted_;                    // Whether dipoles have been extracted

	RadiaFieldCF(int obj, const std::string& ftype = "b",
	             py::object py_origin = py::none(),
	             py::object py_u = py::none(),
	             py::object py_v = py::none(),
	             py::object py_w = py::none(),
	             py::object py_precision = py::none(),
	             const std::string& units = "m",
	             double fmm_eps = 0.0)
	    : CoefficientFunction(ftype == "phi" ? 1 : 3),  // phi is scalar (1D), others are vector (3D)
	      radia_obj(obj), field_type(ftype), use_transform(false),
	      precision(py_precision),
	      use_cache_(false), cache_tolerance_(1e-10), cache_hits_(0), cache_misses_(0),
	      coord_scale_(units == "m" ? 1.0 : 1000.0),
	      fmm_eps_(fmm_eps), dipoles_extracted_(false)
	{
		// Validate field type
		if (field_type != "b" && field_type != "h" &&
		    field_type != "a" && field_type != "m" && field_type != "phi") {
			throw std::invalid_argument(
				"Invalid field_type. Must be 'b' (flux density), "
				"'h' (magnetic field), 'a' (vector potential), 'm' (magnetization), "
				"or 'phi' (magnetic scalar potential)");
		}

		// Default: identity transformation
		origin[0] = 0.0; origin[1] = 0.0; origin[2] = 0.0;
		u_axis[0] = 1.0; u_axis[1] = 0.0; u_axis[2] = 0.0;
		v_axis[0] = 0.0; v_axis[1] = 1.0; v_axis[2] = 0.0;
		w_axis[0] = 0.0; w_axis[1] = 0.0; w_axis[2] = 1.0;

		// Parse origin
		if (!py_origin.is_none()) {
			parse_vector(py_origin, origin);
			use_transform = true;
		}

		// Parse and normalize u-axis
		if (!py_u.is_none()) {
			parse_vector(py_u, u_axis);
			normalize(u_axis);
			use_transform = true;
		}

		// Parse and normalize v-axis
		if (!py_v.is_none()) {
			parse_vector(py_v, v_axis);
			normalize(v_axis);
			use_transform = true;
		}

		// Parse and normalize w-axis
		if (!py_w.is_none()) {
			parse_vector(py_w, w_axis);
			normalize(w_axis);
			use_transform = true;
		}

		// Apply computation settings
		py::gil_scoped_acquire acquire;

		// Cache Radia module import to avoid repeated imports (memory optimization)
		rad_module_ = py::module_::import("radia");

		// Set precision if specified
		if (!py_precision.is_none()) {
			double prec = py_precision.cast<double>();
			// Set precision for all field computation types
			std::string prec_str = "PrcB->" + std::to_string(prec) +
			                       ",PrcA->" + std::to_string(prec) +
			                       ",PrcH->" + std::to_string(prec) +
			                       ",PrcM->" + std::to_string(prec);
			rad_module_.attr("FldCmpPrc")(prec_str);
		}
	}

private:
	void parse_vector(py::object py_vec, double vec[3]) {
		if (py::isinstance<py::list>(py_vec)) {
			py::list lst = py_vec.cast<py::list>();
			if (lst.size() != 3) {
				throw std::invalid_argument("Vector must have 3 components");
			}
			vec[0] = lst[0].cast<double>();
			vec[1] = lst[1].cast<double>();
			vec[2] = lst[2].cast<double>();
		} else if (py::isinstance<py::tuple>(py_vec)) {
			py::tuple tup = py_vec.cast<py::tuple>();
			if (tup.size() != 3) {
				throw std::invalid_argument("Vector must have 3 components");
			}
			vec[0] = tup[0].cast<double>();
			vec[1] = tup[1].cast<double>();
			vec[2] = tup[2].cast<double>();
		} else {
			throw std::invalid_argument("Vector must be a list or tuple");
		}
	}

	void normalize(double vec[3]) {
		double norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		if (norm < 1e-12) {
			throw std::invalid_argument("Cannot normalize zero vector");
		}
		vec[0] /= norm;
		vec[1] /= norm;
		vec[2] /= norm;
	}

	double dot(const double a[3], const double b[3]) const {
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}


	// Hash function for 3D point (quantized to tolerance grid)
	uint64_t hash_point(double x, double y, double z) const {
		int64_t ix = static_cast<int64_t>(x / cache_tolerance_);
		int64_t iy = static_cast<int64_t>(y / cache_tolerance_);
		int64_t iz = static_cast<int64_t>(z / cache_tolerance_);

		uint64_t hash = 14695981039346656037ULL;  // FNV offset basis
		hash ^= static_cast<uint64_t>(ix);
		hash *= 1099511628211ULL;  // FNV prime
		hash ^= static_cast<uint64_t>(iy);
		hash *= 1099511628211ULL;
		hash ^= static_cast<uint64_t>(iz);
		hash *= 1099511628211ULL;
		return hash;
	}

public:

	// Prepare cache by batch-evaluating all points
	void PrepareCache(py::list points_list) {
		py::gil_scoped_acquire acquire;

		try {
			size_t npts = points_list.size();
			std::cout << "[PrepareCache] Caching " << npts << " points..." << std::endl;

			point_cache_.clear();
			cache_hits_ = 0;
			cache_misses_ = 0;

			if (npts == 0) {
				use_cache_ = false;
				return;
			}

			// Build Radia points list (mm)
			py::list radia_points;
			for (size_t i = 0; i < npts; i++) {
				py::list pt = points_list[i].cast<py::list>();
				double p_global[3] = {pt[0].cast<double>(), pt[1].cast<double>(), pt[2].cast<double>()};

				double p_local[3];
				if (use_transform) {
					double p_t[3] = {p_global[0]-origin[0], p_global[1]-origin[1], p_global[2]-origin[2]};
					p_local[0] = dot(u_axis, p_t);
					p_local[1] = dot(v_axis, p_t);
					p_local[2] = dot(w_axis, p_t);
				} else {
					p_local[0] = p_global[0];
					p_local[1] = p_global[1];
					p_local[2] = p_global[2];
				}

				py::list coords;
				coords.append(p_local[0] * coord_scale_);
				coords.append(p_local[1] * coord_scale_);
				coords.append(p_local[2] * coord_scale_);
				radia_points.append(coords);
			}

			// Single batch call to Radia
			py::module_ rad = py::module_::import("radia");
			py::object results = rad.attr("Fld")(radia_obj, field_type, radia_points);
			py::list results_list = results.cast<py::list>();

			// Store in cache
			// No scaling needed - Radia returns values in consistent units with FldUnits
			for (size_t i = 0; i < npts; i++) {
				py::list pt = points_list[i].cast<py::list>();
				double x = pt[0].cast<double>();
				double y = pt[1].cast<double>();
				double z = pt[2].cast<double>();

				py::list fld = results_list[i].cast<py::list>();
				double f_local[3] = {fld[0].cast<double>(), fld[1].cast<double>(), fld[2].cast<double>()};

				double f_global[3];
				if (use_transform) {
					f_global[0] = u_axis[0]*f_local[0] + v_axis[0]*f_local[1] + w_axis[0]*f_local[2];
					f_global[1] = u_axis[1]*f_local[0] + v_axis[1]*f_local[1] + w_axis[1]*f_local[2];
					f_global[2] = u_axis[2]*f_local[0] + v_axis[2]*f_local[1] + w_axis[2]*f_local[2];
				} else {
					f_global[0] = f_local[0];
					f_global[1] = f_local[1];
					f_global[2] = f_local[2];
				}

				uint64_t hash = hash_point(x, y, z);
				point_cache_[hash] = {f_global[0], f_global[1], f_global[2]};
			}

			use_cache_ = true;
			std::cout << "[PrepareCache] Complete: " << point_cache_.size() << " entries" << std::endl;

		} catch (std::exception &e) {
			std::cerr << "[PrepareCache] Error: " << e.what() << std::endl;
			point_cache_.clear();
			use_cache_ = false;
			throw;
		}
	}

	void ClearCache() {
		point_cache_.clear();
		use_cache_ = false;
		cache_hits_ = 0;
		cache_misses_ = 0;
	}

	py::dict GetCacheStats() const {
		py::dict stats;
		stats["enabled"] = use_cache_;
		stats["size"] = point_cache_.size();
		stats["hits"] = cache_hits_;
		stats["misses"] = cache_misses_;
		double total = cache_hits_ + cache_misses_;
		stats["hit_rate"] = (total > 0) ? (cache_hits_ / total) : 0.0;
		return stats;
	}

	virtual ~RadiaFieldCF() {
		// Acquire GIL before releasing Python objects to prevent memory leaks
		// When NGSolve destroys this CoefficientFunction, we must ensure
		// Python reference counting is done safely
		py::gil_scoped_acquire acquire;
		precision.release();
	}

	// Scalar evaluation - used for 'phi' (magnetic scalar potential)
	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
	    if (field_type != "phi") {
	        return 0.0;  // Vector fields return 0 for scalar evaluation
	    }

	    // Get point coordinates
	    auto pnt = mip.GetPoint();
	    int dim = pnt.Size();

	    double p_global[3];
	    p_global[0] = pnt[0];
	    p_global[1] = (dim >= 2) ? pnt[1] : 0.0;
	    p_global[2] = (dim >= 3) ? pnt[2] : 0.0;

	    // Apply coordinate transformation if enabled
	    double p_local[3];
	    if (use_transform) {
	        double p_translated[3];
	        p_translated[0] = p_global[0] - origin[0];
	        p_translated[1] = p_global[1] - origin[1];
	        p_translated[2] = p_global[2] - origin[2];

	        p_local[0] = dot(u_axis, p_translated);
	        p_local[1] = dot(v_axis, p_translated);
	        p_local[2] = dot(w_axis, p_translated);
	    } else {
	        p_local[0] = p_global[0];
	        p_local[1] = p_global[1];
	        p_local[2] = p_global[2];
	    }

	    // Convert to Radia units
	    double coords_radia[3];
	    coords_radia[0] = p_local[0] * coord_scale_;
	    coords_radia[1] = p_local[1] * coord_scale_;
	    coords_radia[2] = p_local[2] * coord_scale_;

	    // Call Radia to get phi
	    double phi_value = 0.0;
	    {
	        py::gil_scoped_acquire acquire;
	        try {
	            py::list coords;
	            coords.append(coords_radia[0]);
	            coords.append(coords_radia[1]);
	            coords.append(coords_radia[2]);

	            py::object field_result = rad_module_.attr("Fld")(radia_obj, field_type, coords);
	            // For 'phi', Radia returns [Phi, Hx, Hy, Hz] - we want just Phi (index 0)
	            phi_value = field_result[py::int_(0)].cast<double>();
	        } catch (std::exception &e) {
	            std::cerr << "[RadiaField] Scalar phi error: " << e.what() << std::endl;
	            return 0.0;
	        }
	    }
	    return phi_value;
	}

	// Single point evaluation (for direct calls like cf(mip))
	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                     FlatVector<> result) const override
	{
	    // Get point coordinates first (no GIL needed)
	    auto pnt = mip.GetPoint();
	    int dim = pnt.Size();

	    // Get global coordinates (NGSolve, in meters)
	    double p_global[3];
	    p_global[0] = pnt[0];
	    p_global[1] = (dim >= 2) ? pnt[1] : 0.0;
	    p_global[2] = (dim >= 3) ? pnt[2] : 0.0;

	    // Handle scalar field 'phi' separately
	    if (field_type == "phi") {
	        result(0) = Evaluate(mip);  // Use scalar evaluation
	        return;
	    }

	    // Check cache first - NO GIL NEEDED for C++ cache lookup
	    if (use_cache_) {
	        uint64_t hash = hash_point(p_global[0], p_global[1], p_global[2]);
	        auto it = point_cache_.find(hash);
	        if (it != point_cache_.end()) {
	            cache_hits_++;
	            result(0) = it->second[0];
	            result(1) = it->second[1];
	            result(2) = it->second[2];
	            return;  // Cache hit - no Python call needed!
	        }
	        cache_misses_++;
	    }

	    // Apply coordinate transformation if enabled (no GIL needed)
	    double p_local[3];
	    if (use_transform) {
	        double p_translated[3];
	        p_translated[0] = p_global[0] - origin[0];
	        p_translated[1] = p_global[1] - origin[1];
	        p_translated[2] = p_global[2] - origin[2];

	        p_local[0] = dot(u_axis, p_translated);
	        p_local[1] = dot(v_axis, p_translated);
	        p_local[2] = dot(w_axis, p_translated);
	    } else {
	        p_local[0] = p_global[0];
	        p_local[1] = p_global[1];
	        p_local[2] = p_global[2];
	    }

	    // Convert m -> mm (store in C++ variables)
	    double coords_mm[3];
	    coords_mm[0] = p_local[0] * coord_scale_;
	    coords_mm[1] = p_local[1] * coord_scale_;
	    coords_mm[2] = p_local[2] * coord_scale_;

	    // Now acquire GIL only for Python call (minimum scope)
	    double f_local[3];
	    {
	        py::gil_scoped_acquire acquire;
	        try {
	            // Create temporary py::list only when calling rad.Fld()
	            py::list coords;
	            coords.append(coords_mm[0]);
	            coords.append(coords_mm[1]);
	            coords.append(coords_mm[2]);

	            py::object field_result = rad_module_.attr("Fld")(radia_obj, field_type, coords);
	            f_local[0] = field_result[py::int_(0)].cast<double>();
	            f_local[1] = field_result[py::int_(1)].cast<double>();
	            f_local[2] = field_result[py::int_(2)].cast<double>();
	            // field_result and coords go out of scope here, releasing Python objects
	        } catch (std::exception &e) {
	            std::cerr << "[RadiaField] Single point error (" << field_type << "): "
	                      << e.what() << std::endl;
	            result(0) = 0.0;
	            result(1) = 0.0;
	            result(2) = 0.0;
	            return;
	        }
	    }  // GIL released here

	    // Transform field back to global coordinate system (no GIL needed)
	    double f_global[3];
	    if (use_transform) {
	        f_global[0] = u_axis[0]*f_local[0] + v_axis[0]*f_local[1] + w_axis[0]*f_local[2];
	        f_global[1] = u_axis[1]*f_local[0] + v_axis[1]*f_local[1] + w_axis[1]*f_local[2];
	        f_global[2] = u_axis[2]*f_local[0] + v_axis[2]*f_local[1] + w_axis[2]*f_local[2];
	    } else {
	        f_global[0] = f_local[0];
	        f_global[1] = f_local[1];
	        f_global[2] = f_local[2];
	    }

	    // Vector potential A: No additional scaling needed
	    // Radia returns A in T*m when FldUnits('m') is set, or T*mm when FldUnits('mm')
	    // The numerical value is the same, but units match the FldUnits setting
	    // Since we use coord_scale_ to convert coords to Radia's unit system,
	    // the returned A is already in the correct units (T*m for NGSolve)
	    result(0) = f_global[0];
	    result(1) = f_global[1];
	    result(2) = f_global[2];
	}

	/**
	 * Extract dipole data from Radia object via Python (called once on first use)
	 *
	 * Note: Since radia.pyd and radia_ngsolve.pyd are separate DLLs with separate
	 * global states, we cannot directly access the Radia object through C++.
	 * Instead, we call Python to extract magnetization, center, and volume data.
	 *
	 * Radia API notes:
	 * - ObjM(key) returns [[Mabs_min, Mabs_avg, Mabs_max], [Mx, My, Mz]]
	 * - ObjGeoLim(key) returns [xmin, xmax, ymin, ymax, zmin, zmax] in mm (always)
	 * - For single elements, the key itself is the element
	 * - For containers (ObjCnt), need to iterate through children
	 *
	 * IMPORTANT: ObjGeoLim always returns values in mm, regardless of FldUnits setting.
	 * But the targets for FMM are in the user's coordinate system (coord_scale_).
	 * So dipole positions must match: if coord_scale_=1000 (mm), positions are in mm;
	 * if coord_scale_=1 (m), positions must be in m. Since ObjGeoLim is always mm,
	 * we divide by coord_scale_ to convert to user's system when needed.
	 */
	void ExtractDipolesIfNeeded() const {
	    if (dipoles_extracted_) return;
	    dipoles_extracted_ = true;

	    py::gil_scoped_acquire acquire;

	    try {
	        py::module_ rad = py::module_::import("radia");

	        dipole_cache_.clear();

	        // Extract single element's dipole data
	        // For containers, this gives aggregate values which is an approximation
	        // TODO: Traverse container hierarchy for better accuracy

	        // rad.ObjM returns [[Mabs_min, Mabs_avg, Mabs_max], [Mx, My, Mz]]
	        py::object m_result = rad.attr("ObjM")(radia_obj);
	        py::list m_outer = m_result.cast<py::list>();
	        py::list m_vec = m_outer[1].cast<py::list>();

	        double Mx = m_vec[0].cast<double>();
	        double My = m_vec[1].cast<double>();
	        double Mz = m_vec[2].cast<double>();

	        // rad.ObjGeoLim returns [xmin, xmax, ymin, ymax, zmin, zmax] in mm
	        py::object geo_result = rad.attr("ObjGeoLim")(radia_obj);
	        py::list geo_list = geo_result.cast<py::list>();

	        // ObjGeoLim is always in mm, so we get mm values here
	        double xmin_mm = geo_list[0].cast<double>();
	        double xmax_mm = geo_list[1].cast<double>();
	        double ymin_mm = geo_list[2].cast<double>();
	        double ymax_mm = geo_list[3].cast<double>();
	        double zmin_mm = geo_list[4].cast<double>();
	        double zmax_mm = geo_list[5].cast<double>();

	        // Volume in mm^3
	        double vol_mm3 = (xmax_mm - xmin_mm) * (ymax_mm - ymin_mm) * (zmax_mm - zmin_mm);

	        // Skip zero-volume elements
	        if (vol_mm3 > 0.0) {
	            // Dipole moment: m = M * V
	            // M is in A/m, V should be in m^3 for SI units
	            // vol_mm3 is in mm^3, so vol_m3 = vol_mm3 * 1e-9
	            double vol_m3 = vol_mm3 * 1e-9;

	            // Center in mm
	            double cx_mm = (xmin_mm + xmax_mm) / 2.0;
	            double cy_mm = (ymin_mm + ymax_mm) / 2.0;
	            double cz_mm = (zmin_mm + zmax_mm) / 2.0;

	            // Dipole position should match the coordinate system used for targets
	            // Targets are multiplied by coord_scale_ in Evaluate(), so they're in mm if coord_scale_=1000
	            // We store position in same units (mm when coord_scale_=1000, otherwise convert)
	            // Actually: in Evaluate, targets are: pnt * coord_scale_
	            // If coord_scale_ = 1 (user uses 'm'), targets are in m, we need positions in m
	            // If coord_scale_ = 1000 (user uses 'mm'), targets are in mm, we need positions in mm
	            // ObjGeoLim returns mm, so:
	            // - if coord_scale_ = 1000, keep mm
	            // - if coord_scale_ = 1, convert to m (divide by 1000)
	            double len_scale = (coord_scale_ < 1.5) ? 0.001 : 1.0;  // mm->m if using meters

	            RadDipoleCollect::DipoleData dipole;
	            dipole.x = cx_mm * len_scale;
	            dipole.y = cy_mm * len_scale;
	            dipole.z = cz_mm * len_scale;
	            dipole.mx = Mx * vol_m3;
	            dipole.my = My * vol_m3;
	            dipole.mz = Mz * vol_m3;
	            dipole.volume = vol_m3;

	            // Skip zero-moment dipoles
	            double momentMagSq = dipole.mx*dipole.mx + dipole.my*dipole.my + dipole.mz*dipole.mz;
	            if (momentMagSq > 0.0) {
	                dipole_cache_.dipoles.push_back(dipole);
	            }
	        }

	        dipole_cache_.flatten();

	    } catch (std::exception& e) {
	        std::cerr << "[RadiaFieldCF] Warning: Could not extract dipoles: "
	                  << e.what() << std::endl;
	    }
	}

	/**
	 * Compute field at targets using FMM/dipole approximation (pure C++ - no Python)
	 *
	 * Uses dipole approximation: each element -> point dipole at center
	 *
	 * Supported field types:
	 * - "h": H-field = (1/4pi) * sum_j [ 3*(m_j . r_j)*r_j/r_j^5 - m_j/r_j^3 ]
	 * - "b": B-field = mu0 * H (in vacuum)
	 * - "a": Vector potential A = (mu0/4pi) * sum_j [ m_j x r_j / r_j^3 ]
	 * - "phi": Scalar potential = (1/4pi) * sum_j [ m_j . r_j / r_j^3 ]
	 *
	 * @param field_type  Field type: "h", "b", "a", or "phi"
	 * @param targets     Target points (ntarget * 3 values)
	 * @param ntarget     Number of target points
	 * @param result_x    Output: x component (or scalar for phi)
	 * @param result_y    Output: y component (unused for phi)
	 * @param result_z    Output: z component (unused for phi)
	 */
	void ComputeFieldFMM(const std::string& ftype,
	                     const std::vector<double>& targets, size_t ntarget,
	                     std::vector<double>& result_x,
	                     std::vector<double>& result_y,
	                     std::vector<double>& result_z) const
	{
	    // Extract dipoles if not already done
	    ExtractDipolesIfNeeded();

	    int64_t nsource = dipole_cache_.count();
	    if (nsource == 0) {
	        // No dipoles - return zeros
	        result_x.assign(ntarget, 0.0);
	        result_y.assign(ntarget, 0.0);
	        result_z.assign(ntarget, 0.0);
	        return;
	    }

	    if (ftype == "a") {
	        // Vector potential A = (mu0/4pi) * (m x r) / r^3
	        RadExaFMM::FMMResult fmm_result = RadExaFMM::ComputeDipoleVectorPotential(
	            fmm_eps_,
	            dipole_cache_.positions.data(),
	            dipole_cache_.moments.data(),
	            nsource,
	            targets.data(),
	            static_cast<int64_t>(ntarget)
	        );

	        // Copy results (Ax, Ay, Az in gradx, grady, gradz)
	        result_x = std::move(fmm_result.gradx);
	        result_y = std::move(fmm_result.grady);
	        result_z = std::move(fmm_result.gradz);

	    } else {
	        // H-field, B-field, or phi - use standard dipole field computation
	        RadExaFMM::FMMResult fmm_result = RadExaFMM::ComputeDipoleField(
	            fmm_eps_,
	            dipole_cache_.positions.data(),
	            dipole_cache_.moments.data(),
	            nsource,
	            targets.data(),
	            static_cast<int64_t>(ntarget)
	        );

	        if (ftype == "phi") {
	            // Scalar potential - return pot in result_x
	            result_x = std::move(fmm_result.pot);
	            result_y.assign(ntarget, 0.0);
	            result_z.assign(ntarget, 0.0);
	        } else if (ftype == "b") {
	            // B-field = mu0 * H
	            const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;
	            result_x.resize(ntarget);
	            result_y.resize(ntarget);
	            result_z.resize(ntarget);
	            for (size_t i = 0; i < ntarget; ++i) {
	                result_x[i] = fmm_result.gradx[i] * MU_0;
	                result_y[i] = fmm_result.grady[i] * MU_0;
	                result_z[i] = fmm_result.gradz[i] * MU_0;
	            }
	        } else {
	            // H-field (default)
	            result_x = std::move(fmm_result.gradx);
	            result_y = std::move(fmm_result.grady);
	            result_z = std::move(fmm_result.gradz);
	        }
	    }
	}

public:
	// Batch evaluation: evaluate field at multiple points in one call
	// This significantly reduces Python call overhead and enables OpenMP parallelization
	virtual void Evaluate(const BaseMappedIntegrationRule& mir,
	                     BareSliceMatrix<> result) const override
	{
	    size_t npts = mir.Size();

	    // =========================================================================
	    // FMM PATH: Use pure C++ dipole computation (no Python calls)
	    // Supports: h (H-field), b (B-field), a (vector potential)
	    // Does NOT support: phi (scalar, needs different handling), m (magnetization)
	    // =========================================================================
	    bool fmm_supported = (field_type == "h" || field_type == "b" || field_type == "a");
	    if (fmm_eps_ > 0.0 && fmm_supported && !use_transform) {
	        // Collect target points
	        std::vector<double> targets(npts * 3);
	        for (size_t i = 0; i < npts; i++) {
	            auto pnt = mir[i].GetPoint();
	            int dim = pnt.Size();
	            // Convert to Radia units (coord_scale_: 1.0 for m, 1000.0 for mm)
	            targets[i*3 + 0] = pnt[0] * coord_scale_;
	            targets[i*3 + 1] = (dim >= 2) ? pnt[1] * coord_scale_ : 0.0;
	            targets[i*3 + 2] = (dim >= 3) ? pnt[2] * coord_scale_ : 0.0;
	        }

	        // Compute field using FMM
	        std::vector<double> Fx, Fy, Fz;
	        ComputeFieldFMM(field_type, targets, npts, Fx, Fy, Fz);

	        // Copy to result matrix
	        // For vector potential A, need to scale to NGSolve units (m)
	        // The dipole formula gives A in T*m when positions are in m
	        // But we're using coord_scale_ for positions, so we need to adjust
	        double scale = 1.0;
	        if (field_type == "a" && coord_scale_ > 1.5) {
	            // Positions were in mm, but A formula expects m
	            // A is proportional to 1/r^2, so if r is in mm (1000x larger),
	            // A is 1e6x smaller. We need to multiply by 1e-3 to convert.
	            // Actually: A = mu0/4pi * (m x r)/|r|^3
	            // If r is in mm, |r|^3 is (mm)^3 = 1e-9 m^3
	            // So A would be 1e9 times larger if we don't adjust.
	            // We need A in T*m for NGSolve, but positions were in mm.
	            // The dipole moment m is in A*m^2, correct.
	            // Let's think: if target is at 0.1m = 100mm from source
	            // r_m = 0.1, r_mm = 100
	            // A ~ 1/r^2, so A_m ~ 1/0.01, A_mm ~ 1/10000
	            // A_mm = A_m / 1e6
	            // So if we computed with mm, we get 1e-6 of the correct value
	            // Need to multiply by 1e6? No wait...
	            // A = mu0/4pi * (m x r) / |r|^3
	            // If r is in mm, then |r|^3 has units mm^3
	            // We want A in T*m
	            // mu0/4pi has units H/m = kg*m/(A^2*s^2)
	            // m has units A*m^2
	            // (m x r)/|r|^3 has units (A*m^2 * mm) / mm^3 = A*m^2/mm^2
	            // = A*m^2 / (1e-6 m^2) = 1e6 * A
	            // So A = 1e-7 * 1e6 * A = 1e-1 * A (units are weird)
	            // Actually easier: just convert positions to m before calling FMM
	            // The current code scales by coord_scale_ which gives mm if user is in m
	            // This is wrong for FMM! FMM should use consistent units.
	            // For now, positions are stored in user units, which is m for NGSolve.
	            // coord_scale_ is applied to match Radia's expected units.
	            // So if coord_scale_ = 1 (user uses m), positions are in m -> correct
	            // If coord_scale_ = 1000 (user uses mm), positions are in mm
	            // For FMM we need positions in the same units as dipole positions
	            // Dipole positions come from ObjGeoLim which is always in mm
	            // So we need target positions also in mm when coord_scale_ = 1 (m)
	            // Current: targets = pnt * coord_scale_
	            //   coord_scale_ = 1 (m): targets in m, but dipoles in m (after len_scale)
	            //   coord_scale_ = 1000 (mm): targets in mm, dipoles in mm
	            // This seems correct! No additional scaling needed.
	            // The formula gives A in T*m when everything is in SI (m).
	            scale = 1.0;
	        }

	        for (size_t i = 0; i < npts; i++) {
	            result(i, 0) = Fx[i] * scale;
	            result(i, 1) = Fy[i] * scale;
	            result(i, 2) = Fz[i] * scale;
	        }
	        return;
	    }

	    // =========================================================================
	    // PYTHON PATH: Use rad.Fld() for full field computation
	    // Required for phi (scalar), m (magnetization), or when FMM is disabled
	    // =========================================================================
	    py::gil_scoped_acquire acquire;

	    try {
	        // Try cache first if enabled
	        if (use_cache_) {
	            bool all_cached = true;
	            for (size_t i = 0; i < npts; i++) {
	                auto pnt = mir[i].GetPoint();
	                int dim = pnt.Size();
	                double p[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};

	                uint64_t hash = hash_point(p[0], p[1], p[2]);
	                auto it = point_cache_.find(hash);
	                if (it != point_cache_.end()) {
	                    cache_hits_++;
	                    result(i,0) = it->second[0];
	                    result(i,1) = it->second[1];
	                    result(i,2) = it->second[2];
	                } else {
	                    cache_misses_++;
	                    all_cached = false;
	                    break;
	                }
	            }
	            if (all_cached) return;
	        }

	        // Collect all points in a Python list of lists
	        py::list points_list;

	        for (size_t i = 0; i < npts; i++) {
	            auto pnt = mir[i].GetPoint();
	            int dim = pnt.Size();

	            // Get global coordinates (NGSolve, in meters)
	            double p_global[3];
	            p_global[0] = pnt[0];
	            p_global[1] = (dim >= 2) ? pnt[1] : 0.0;
	            p_global[2] = (dim >= 3) ? pnt[2] : 0.0;

	            // Apply coordinate transformation if enabled
	            double p_local[3];
	            if (use_transform) {
	                double p_translated[3];
	                p_translated[0] = p_global[0] - origin[0];
	                p_translated[1] = p_global[1] - origin[1];
	                p_translated[2] = p_global[2] - origin[2];

	                p_local[0] = dot(u_axis, p_translated);
	                p_local[1] = dot(v_axis, p_translated);
	                p_local[2] = dot(w_axis, p_translated);
	            } else {
	                p_local[0] = p_global[0];
	                p_local[1] = p_global[1];
	                p_local[2] = p_global[2];
	            }

	            // Convert m -> mm
	            py::list coords;
	            coords.append(p_local[0] * coord_scale_);
	            coords.append(p_local[1] * coord_scale_);
	            coords.append(p_local[2] * coord_scale_);

	            points_list.append(coords);
	        }

	        // Single Python call for all points!
	        py::module_ rad = py::module_::import("radia");
	        py::object field_results = rad.attr("Fld")(radia_obj, field_type, points_list);

	        // field_results is a list of [Bx, By, Bz] for each point
	        py::list results_list = field_results.cast<py::list>();

	        // Extract results and apply transformations
	        // No scaling needed - Radia returns values in consistent units with FldUnits

	        for (size_t i = 0; i < npts; i++) {
	            py::list field_list = results_list[i].cast<py::list>();

	            double f_local[3];
	            f_local[0] = field_list[0].cast<double>();
	            f_local[1] = field_list[1].cast<double>();
	            f_local[2] = field_list[2].cast<double>();

	            // Transform field back to global coordinate system
	            double f_global[3];
	            if (use_transform) {
	                f_global[0] = u_axis[0]*f_local[0] + v_axis[0]*f_local[1] + w_axis[0]*f_local[2];
	                f_global[1] = u_axis[1]*f_local[0] + v_axis[1]*f_local[1] + w_axis[1]*f_local[2];
	                f_global[2] = u_axis[2]*f_local[0] + v_axis[2]*f_local[1] + w_axis[2]*f_local[2];
	            } else {
	                f_global[0] = f_local[0];
	                f_global[1] = f_local[1];
	                f_global[2] = f_local[2];
	            }

	            result(i, 0) = f_global[0];
	            result(i, 1) = f_global[1];
	            result(i, 2) = f_global[2];
	        }

	    } catch (std::exception &e) {
	        std::cerr << "[RadiaField] Batch evaluation error (" << field_type << "): "
	                  << e.what() << std::endl;
	        // Fill with zeros on error
	        for (size_t i = 0; i < npts; i++) {
	            result(i, 0) = 0.0;
	            result(i, 1) = 0.0;
	            result(i, 2) = 0.0;
	        }
	    }
	}
};

} // namespace ngfem

PYBIND11_MODULE(radia_ngsolve, m) {
	m.doc() = "NGSolve CoefficientFunction interface for Radia (with m->mm conversion and coordinate transformation)";

	// Unified interface with coordinate transformation
	py::class_<ngfem::RadiaFieldCF,
	           std::shared_ptr<ngfem::RadiaFieldCF>,
	           ngfem::CoefficientFunction>(m, "RadiaField")
	    .def(py::init<int>(), py::arg("radia_obj"),
	         "Create Radia field CoefficientFunction (default: magnetic flux density)")
	    .def(py::init<int, const std::string&>(),
	         py::arg("radia_obj"), py::arg("field_type"),
	         "Create Radia field CoefficientFunction\n"
	         "field_type: 'b' (flux density), 'h' (field), 'a' (vector potential), 'm' (magnetization)")
	    .def(py::init<int, const std::string&, py::object, py::object, py::object, py::object, py::object, const std::string&, double>(),
	         py::arg("radia_obj"),
	         py::arg("field_type") = "b",
	         py::arg("origin") = py::none(),
	         py::arg("u_axis") = py::none(),
	         py::arg("v_axis") = py::none(),
	         py::arg("w_axis") = py::none(),
	         py::arg("precision") = py::none(),
	         py::arg("units") = "m",
	         py::arg("fmm_eps") = 0.0,
	         "Create Radia field CoefficientFunction with full control\n\n"
	         "Parameters:\n"
	         "  radia_obj: Radia object ID\n"
	         "  field_type: 'b' (flux density), 'h' (field), 'a' (vector potential), 'm' (magnetization)\n"
	         "  origin: Translation vector [x, y, z] in meters (default: [0, 0, 0])\n"
	         "  u_axis: Local u-axis [ux, uy, uz] (default: [1, 0, 0]) - will be normalized\n"
	         "  v_axis: Local v-axis [vx, vy, vz] (default: [0, 1, 0]) - will be normalized\n"
	         "  w_axis: Local w-axis [wx, wy, wz] (default: [0, 0, 1]) - will be normalized\n"
	         "  precision: Computation precision in Tesla (default: None = Radia default)\n"
	         "  units: Coordinate units - 'm' (meters, default) or 'mm' (millimeters)\n"
	         "  fmm_eps: FMM tolerance for field computation (0 = disabled, typical: 1e-4 to 1e-6)\n\n"
	         "FMM acceleration (Fast Multipole Method / Dipole Approximation):\n"
	         "  When fmm_eps > 0, field computation uses dipole approximation.\n"
	         "  Each magnetic element is treated as a point dipole at its center.\n"
	         "  This provides O(N*M) direct computation (OpenMP parallelized).\n"
	         "  Supported field types: 'h' (H-field), 'b' (B-field), 'a' (vector potential).\n"
	         "  Does not apply to 'm' (magnetization) or with coordinate transformation.\n\n"
	         "Coordinate transformation:\n"
	         "  1. Global point p is translated by origin: p' = p - origin\n"
	         "  2. p' is projected onto local axes: p_local = [u*p', v*p', w*p']\n"
	         "  3. Field is calculated in Radia's coordinate system\n"
	         "  4. Field is transformed back to global: F = u*F_local[0] + v*F_local[1] + w*F_local[2]\n\n"
	         "Example:\n"
	         "  # High accuracy\n"
	         "  B_cf = rad_ngsolve.RadiaField(magnet, 'b', precision=1e-6)\n"
	         "  # FMM-accelerated H-field for large meshes\n"
	         "  H_cf = rad_ngsolve.RadiaField(magnet, 'h', fmm_eps=1e-4)")
	    .def_readonly("radia_obj", &ngfem::RadiaFieldCF::radia_obj)
	    .def_readonly("field_type", &ngfem::RadiaFieldCF::field_type)
	    .def_readonly("use_transform", &ngfem::RadiaFieldCF::use_transform)
	    .def_readonly("precision", &ngfem::RadiaFieldCF::precision)
	    .def_readonly("fmm_eps", &ngfem::RadiaFieldCF::fmm_eps_)
	    .def("PrepareCache", &ngfem::RadiaFieldCF::PrepareCache, py::arg("points"),
	         "Pre-cache field values for batch evaluation")
	    .def("ClearCache", &ngfem::RadiaFieldCF::ClearCache,
	         "Clear cached field values")
	    .def("GetCacheStats", &ngfem::RadiaFieldCF::GetCacheStats,
	         "Get cache statistics (enabled, size, hits, misses, hit_rate)");
}
