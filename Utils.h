#ifndef __MN_UTILS_H__
#define __MN_UTILS_H__

#ifdef _MSC_VER
#pragma once
#endif

#define SQ(x) ((x) * (x))

#include <vector>
#include <array>
#include <utility>
#include <limits>
#include <stdexcept>

namespace MN {
	using Real = double;	// Adjust to your taste...
	using Real2 = std::pair<Real, Real>;
	using Real3 = std::array<Real, 3>;

	// Min, Max values 
	const double maxDouble = std::numeric_limits<double>::max();
	const double minDouble = -maxDouble;
	const double maxFloat = std::numeric_limits<float>::max();
	const double minFloat = -maxFloat;

	// PI-related values
	const Real PI = (Real)3.14159265358979323846;
	const Real PI05 = PI * (Real)0.5;
	const Real PI15 = PI * (Real)1.5;
	const Real PI20 = PI * (Real)2.0;

	// Domain : Represents a closed domain of [a, b]
	//          Must be initialized with [a, b] value that [ b ] is not smaller than [ a ]
	class Domain {
	protected:
		Real2 dat;
	public:
		Domain() {
			dat.first = 0;
			dat.second = 0;
		}
		// Set this domain with given values
		// [ b ] must not be smaller than [ a ]
		inline virtual void set(Real a, Real b) {
			if (a > b)
				throw(std::runtime_error("Invalid domain : First value must be smaller than second value"));
			dat.first = a;
			dat.second = b;
		}
		// [ b ] must not be smaller than [ a ]
		inline static Domain create(Real a, Real b) {
			Domain domain;
			domain.set(a, b);
			return domain;
		}
		// Find width (end - beg) of this domain
		inline Real width() const noexcept {
			return end() - beg();
		}
		// Find whether given value [ val ] is in this domain
		inline virtual bool has(Real val) const noexcept {
			return (beg() <= val && val <= end());
		}
		// Find whether given domain [ domain ] is in this domain
		inline virtual bool has(const Domain& domain) const noexcept {
			return has(domain.beg()) && has(domain.end());
		}
		// Find middle valud of this domain
		inline Real middle() const {
			return (beg() + end()) * 0.5;
		}
		// Interpolate with given value [ t ] (must be between [0, 1], 0 for [ beg ], 1 for [ end ])
		inline Real interp(Real t) const {
			return beg() * (1.0 - t) + end() * t;
		}
		
		// Begin and End value
		inline const Real& beg() const noexcept {
			return dat.first;
		}
		inline const Real& end() const noexcept {
			return dat.second;
		}
		
		// Domain Sum
		inline Domain operator+(const Domain& domain) {
			auto b = (beg() < domain.beg()) ? beg() : domain.beg();
			auto e = (end() < domain.end()) ? domain.end() : end();
			return create(b, e);
		}
		// Domain Equality
		inline virtual bool operator==(const Domain& domain) const {
			return (beg() == domain.beg() && end() == domain.end());
		}
		inline virtual bool operator!=(const Domain& domain) const {
			return !(*this == domain);
		}
	};
	class Domain2 {
	public:
		Domain a = Domain::create(0, 0);
		Domain b = Domain::create(0, 0);
	};
	class Domain3 {
	public:
		Domain a = Domain::create(0, 0);
		Domain b = Domain::create(0, 0);
		Domain c = Domain::create(0, 0);
	};

	// piDomain : Domain whose width does not exceed 2 * PI
	class piDomain : public Domain {
	private:
		piDomain() = default;
		Real2 regular;	// Regularized domain e.g. [-30, 30] ===> [330, 30] (in angle)
	public:
		piDomain() {
			dat.first = 0;
			dat.second = 0;
			regular.first = 0;
			regular.second = 0;
		}
		// Set this pi domain with given values
		// [ b ] must not be smaller than [ a ] and must be smaller than [ a + 2 * PI ]
		inline virtual void set(Real a, Real b) {
			if (b - a > PI20)
				throw(std::runtime_error("Invalid domain : Width cannot exceed PI20"));
			Domain::set(a, b);
			regular.first = regularize(a);
			regular.second = regularize(b);
		}
		// [ b ] must not be smaller than [ a ] and must be smaller than [ a + 2 * PI ]
		inline static piDomain create(Real a, Real b) {
			piDomain domain;
			domain.set(a, b);
			return domain;
		}
		// Find whether given value [ val ] is in this pi domain, in PI-sense
		// e.g. (angle) 0 is within pi domain [330, 390]
		inline virtual bool has(Real val) const noexcept {
			Real reg = regularize(val);
			if (regular.first < regular.second)
				return (regular.first <= reg && reg <= regular.second);
			else if (regular.first > regular.second)
				return (regular.first <= reg || reg <= regular.second);
			else
				return (beg() == end() ? (reg == regular.first) : true);
		}
		
		// Turn [ x ] into equivalent value that is in [ 0, PI20 ).
		inline static Real regularize(Real x) noexcept {
			int nom = (int)(x / PI20);
			if (x < 0)
				nom -= 1;
			return x - (nom * PI20);
		}
		// Get intersection domain with [ other ] domain in this domain's viewpoint
		// At most 4 intersecting subdomains [ result ] could occur, valid bits are set in [ valid ]
		inline void intersect(const piDomain& other, piDomain result[4], bool valid[4]) const {
			valid[0] = false;
			valid[1] = false;
			valid[2] = false;
			valid[3] = false;
			if (regular.first < regular.second) {
				if (other.regular.first < other.regular.second) {
					if (other.regular.first <= regular.first) {
						if (other.regular.second >= regular.second) {
							result[0] = *this;
							valid[0] = true;
							return;
						}
						else {
							result[0] = piDomain::create(regular.first, other.regular.second);
							valid[0] = true;
							return;
						}
					}
					else if (other.regular.first >= regular.first && other.regular.first <= regular.second) {
						if (other.regular.second >= regular.second) {
							result[0] = piDomain::create(other.regular.first, regular.second);
							valid[0] = true;
							return;
						}
						else {
							result[0] = other;
							valid[0] = true;
							return;
						}
					}
					else
						return;
				}
				else {
					// [ 0, other.regular.second ]
					if (regular.first > other.regular.second)
						valid[0] = false;
					else if (regular.second > other.regular.second) {
						result[0] = piDomain::create(regular.first, other.regular.second);
						valid[0] = true;
					}
					else {
						result[0] = *this;
						valid[0] = true;
					}
					// [ other.regular.first, PI20 ]
					if (regular.second < other.regular.first)
						valid[1] = false;
					else if (regular.first < other.regular.first) {
						result[1] = piDomain::create(other.regular.first, regular.second);
						valid[1] = true;
					}
					else {
						result[1] = *this;
						valid[1] = true;
					}
				}
			}
			else {
				if (other.regular.first < other.regular.second) {
					// [ 0, regular.second ]
					if (other.regular.first > regular.second)
						valid[0] = false;
					else if (other.regular.second > regular.second) {
						result[0] = piDomain::create(other.regular.first, regular.second);
						valid[0] = true;
					}
					else {
						result[0] = other;
						valid[0] = true;
					}

					// [ regular.first, PI20 ]
					if (other.regular.second < regular.first)
						valid[1] = false;
					else if (other.regular.first < regular.first) {
						result[1] = piDomain::create(regular.first, other.regular.second);
						valid[1] = true;
					}
					else {
						result[1] = other;
						valid[1] = true;
					}
				}
				else {
					// [ 0, regular.second ] - [ 0, other.regular.second ]
					if (regular.second < other.regular.second) {
						result[0] = piDomain::create(0, regular.second);
						valid[0] = true;
					}
					else {
						result[0] = piDomain::create(0, other.regular.second);
						valid[0] = true;
					}

					// [ 0, regular.second ] - [ other.regular.first, PI20 ]
					if (regular.second < other.regular.first)
						valid[1] = false;
					else {
						result[1] = piDomain::create(other.regular.first, regular.second);
						valid[1] = true;
					}

					// [ regular.first, PI20 ] - [ 0, other.regular.second ]
					if (regular.first > other.regular.second)
						valid[2] = false;
					else {
						result[2] = piDomain::create(regular.first, other.regular.second);
						valid[2] = true;
					}

					// [ regular.first, PI20 ] - [ other.regular.first, PI20 ]
					if (regular.first < other.regular.first) {
						result[3] = piDomain::create(other.regular.first, PI20);
						valid[3] = true;
					}
					else {
						result[3] = piDomain::create(regular.first, PI20);
						valid[3] = true;
					}
				}
			}
		}
		// @TODO : Get cos, sin min max value of this domain
		inline void minmaxCos(Real& min, Real& max) const noexcept {
			min = -1;
			max = 1;
		}
		inline void minmaxSin(Real& min, Real& max) const noexcept {
			min = -1;
			max = 1;
		}
	};

	// Vec2 : Array of two real values
	class Vec2 : public std::array<Real, 2> {
	public:
		Vec2() = default;
		Vec2(Real a, Real b) {
			(*this)[0] = a;
			(*this)[1] = b;
		}

		// Create zero vector
		inline static Vec2 zero() noexcept {
			return Vec2{ 0, 0 };
		}

		// Basic arithmetic operations
		inline Vec2 operator+(const Vec2& v) const noexcept {
			return Vec2{ (*this)[0] + v[0], (*this)[1] + v[1] };
		}
		inline Vec2 operator-(const Vec2& v) const noexcept {
			return Vec2{ (*this)[0] - v[0], (*this)[1] - v[1] };
		}
		inline Vec2 operator*(const Real& c) const noexcept {
			return Vec2{ (*this)[0] * c, (*this)[1] * c };
		}
		inline Vec2 operator/(const Real& c) const noexcept {
			return Vec2{ (*this)[0] / c, (*this)[1] / c };
		}

		inline void operator+=(const Vec2& v) noexcept {
			(*this)[0] += v[0];
			(*this)[1] += v[1];
		}
		inline void operator-=(const Vec2& v) noexcept {
			(*this)[0] -= v[0];
			(*this)[1] -= v[1];
		}
		inline void operator*=(const Real& c) noexcept {
			(*this)[0] *= c;
			(*this)[1] *= c;
		}
		inline void operator/=(const Real& c) noexcept {
			(*this)[0] /= c;
			(*this)[1] /= c;
		}

		// Dot product
		inline Real dot(const Vec2& v) const noexcept {
			return (*this)[0] * v[0] + (*this)[1] * v[1];
		}
		// L2 norm of this vector
		inline Real len() const noexcept {
			return sqrt(lensq());
		}
		// Square of L2 norm of this vector
		inline Real lensq() const noexcept {
			return dot(*this);
		}
		// Distance between this point and given point [ v ]
		inline Real dist(const Vec2& v) const noexcept {
			return sqrt(distsq(v));
		}
		// Square of distance between this point and given point [ v ]
		inline Real distsq(const Vec2& v) const noexcept {
			auto diff = *this - v;
			return diff.lensq();
		}
		// Normalize this vector into unit vector
		inline void normalize() {
			auto length = len();
			if(length == 0)
				throw(std::runtime_error("Cannot normalize zero vector"));
			(*this) /= length;
		}

		// Get contribution of [ basis ] vector to [ vec ] vector
		// e.g. When vec == [ 2, 3 ], for [ basis ] E1 : 2, E2 : 3
		// @unitBasis : When [ basis ] vector is unit vector, plz turn it on for acceleration
		inline static Real factorize(const Vec2& vec, const Vec2& basis, bool unitBasis = false) noexcept {
			if (unitBasis) return vec.dot(basis);
			else return vec.dot(basis) / basis.lensq();
		}

		/* Legacy */
		/*inline static void add(const Vec2& a, const Vec2& b, Vec2& res) {
			res[0] = a[0] + b[0];
			res[1] = a[1] + b[1];
		}
		inline static void sub(const Vec2& a, const Vec2& b, Vec2& res) {
			res[0] = a[0] - b[0];
			res[1] = a[1] - b[1];
		}
		inline static void mult(const Vec2& a, const Real c, Vec2& res) {
			res[0] = a[0] * c;
			res[1] = a[1] * c;
		}
		inline static void div(const Vec2& a, const Real c, Vec2& res) {
			res[0] = a[0] / c;
			res[1] = a[1] / c;
		}*/
	};

	// Vec3
	class Vec3 : public std::array<Real, 3> {
	public:
		Vec3() = default;
		Vec3(Real a, Real b, Real c) {
			(*this)[0] = a;
			(*this)[1] = b;
			(*this)[2] = c;
		}
		// Create zero vector
		inline static Vec3 zero() noexcept {
			return Vec3{ 0, 0, 0 };
		}

		// Basic arithmetic operations
		inline Vec3 operator+(const Vec3& v) const noexcept {
			return Vec3{ (*this)[0] + v[0], (*this)[1] + v[1], (*this)[2] + v[2] };
		}
		inline Vec3 operator-(const Vec3& v) const noexcept {
			return Vec3{ (*this)[0] - v[0], (*this)[1] - v[1], (*this)[2] - v[2] };
		}
		inline Vec3 operator*(const Real& c) const noexcept {
			return Vec3{ (*this)[0] * c, (*this)[1] * c, (*this)[2] * c };
		}
		inline Vec3 operator/(const Real& c) const noexcept {
			return Vec3{ (*this)[0] / c, (*this)[1] / c, (*this)[2] / c };
		}

		inline void operator+=(const Vec3& v) noexcept {
			(*this)[0] += v[0];
			(*this)[1] += v[1];
			(*this)[2] += v[2];
		}
		inline void operator-=(const Vec3& v) noexcept {
			(*this)[0] -= v[0];
			(*this)[1] -= v[1];
			(*this)[2] -= v[2];
		}
		inline void operator*=(const Real& c) noexcept {
			(*this)[0] *= c;
			(*this)[1] *= c;
			(*this)[2] *= c;
		}
		inline void operator/=(const Real& c) noexcept {
			(*this)[0] /= c;
			(*this)[1] /= c;
			(*this)[2] /= c;
		}

		// Dot & Cross product
		inline Real dot(const Vec3& v) const noexcept {
			return (*this)[0] * v[0] + (*this)[1] * v[1] + (*this)[2] * v[2];
		}
		inline Vec3 cross(const Vec3& v) const noexcept {
			return Vec3{
				(*this)[1] * v[2] - (*this)[2] * v[1],
				(*this)[2] * v[0] - (*this)[0] * v[2],
				(*this)[0] * v[1] - (*this)[1] * v[0]
			};
		}

		// Triple cross : Get determinant of 3 x 3 matrix 
		//                [ a[0] , a[1] , a[2] ]
		//                [ b[0] , b[1] , b[2] ]
		//				  [ c[0] , c[1] , c[2] ]
		inline static Real Tcross(const Vec3& a, const Vec3& b, const Vec3& c) noexcept {
			return a[0] * (b[1] * c[2] - b[2] * c[1])
				- a[1] * (b[0] * c[2] - b[2] * c[0])
				+ a[2] * (b[0] * c[1] - b[1] * c[0]);
		}

		// L2 norm of this vector
		inline Real len() const noexcept {
			return sqrt(lensq());
		}
		// Square of L2 norm of this vector
		inline Real lensq() const noexcept {
			return dot(*this);
		}
		// Distance between this point and given point [ v ]
		inline Real dist(const Vec3& v) const noexcept {
			return sqrt(distsq(v));
		}
		// Square of distance between this point and given point [ v ]
		inline Real distsq(const Vec3& v) const noexcept {
			auto diff = *this - v;
			return diff.lensq();
		}
		// Normalize this vector into unit vector
		inline void normalize() {
			auto length = len();
			if (length == 0)
				return;
				//throw(std::runtime_error("Cannot normalize zero vector"));
			(*this) /= length;
		}

		// Get contribution of [ basis ] vector to [ vec ] vector
		// e.g. When vec == [ 2, 3, 1 ], for basis E1 : 2, E2 : 3, E3 : 1
		// @unitBasis : When [ basis ] vector is unit vector, plz turn it on for acceleration
		inline static Real factorize(const Vec3& vec, const Vec3& basis, bool unitBasis = false) noexcept {
			if (unitBasis) return vec.dot(basis);
			else return vec.dot(basis) / basis.lensq();
		}
		
		/* Legacy */
		/*inline static void add(const Vec3& a, const Vec3& b, Vec3& res) {
			res[0] = a[0] + b[0];
			res[1] = a[1] + b[1];
			res[2] = a[2] + b[2];
		}
		inline static void sub(const Vec3& a, const Vec3& b, Vec3& res) {
			res[0] = a[0] - b[0];
			res[1] = a[1] - b[1];
			res[2] = a[2] - b[2];
		}
		inline static void mult(const Vec3& a, const Real c, Vec3& res) {
			res[0] = a[0] * c;
			res[1] = a[1] * c;
			res[2] = a[2] * c;
		}
		inline static void div(const Vec3& a, const Real c, Vec3& res) {
			res[0] = a[0] / c;
			res[1] = a[1] / c;
			res[2] = a[2] / c;
		}
		inline static void cross(const Vec3& a, const Vec3& b, Vec3& res) {
			res[0] = a[1] * b[2] - a[2] * b[1];
			res[1] = a[2] * b[0] - a[0] * b[2];
			res[2] = a[0] * b[1] - a[1] * b[0];
		}*/
	};

	// Mat3 : 3 x 3 matrix in row-wise order
	class Mat3 : public std::array<std::array<Real, 3>, 3> {
	public:
		Mat3() = default;
		// Create diagonal matrix with diagonal value [ diagonal ]
		Mat3(Real diagonal) {
			(*this)[0][0] = diagonal;
			(*this)[0][1] = 0;
			(*this)[0][2] = 0;

			(*this)[1][0] = 0;
			(*this)[1][1] = diagonal;
			(*this)[1][2] = 0;

			(*this)[2][0] = 0;
			(*this)[2][1] = 0;
			(*this)[2][2] = diagonal;
		}
		// Create identity matrix
		inline static Mat3 identity() noexcept {
			return Mat3(1.0);
		}
		// Create transpose of this matrix
		inline Mat3 transpose() const noexcept {
			Mat3 m;
			m[0][0] = (*this)[0][0];
			m[0][1] = (*this)[1][0];
			m[0][2] = (*this)[2][0];

			m[1][0] = (*this)[0][1];
			m[1][1] = (*this)[1][1];
			m[1][2] = (*this)[2][1];

			m[2][0] = (*this)[0][2];
			m[2][1] = (*this)[1][2];
			m[2][2] = (*this)[2][2];

			return m;
		}

		// Matrix - Matrix multiplication
		inline Mat3 operator*(const Mat3& m) const noexcept {
			Mat3 tmp;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					tmp[i][j] = 0.0;
					for (int k = 0; k < 3; k++)
						tmp[i][j] += (*this)[i][k] * m[k][j];
				}
			return tmp;
		}
		// Matrix - Vector multiplication
		inline Vec3 operator*(const Vec3& v) const noexcept {
			return Vec3{
				(*this)[0][0] * v[0] + (*this)[0][1] * v[1] + (*this)[0][2] * v[2],
				(*this)[1][0] * v[0] + (*this)[1][1] * v[1] + (*this)[1][2] * v[2],
				(*this)[2][0] * v[0] + (*this)[2][1] * v[1] + (*this)[2][2] * v[2]
			};
		}

		// Transpose this matrix and multiply given matrix
		inline Mat3 Tmult(const Mat3& m) const noexcept {
			Mat3 tmp;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					tmp[i][j] = 0.0;
					for (int k = 0; k < 3; k++)
						tmp[i][j] += (*this)[k][i] * m[k][j];
				}
			return tmp;
		}
		// Transpose this matrix and multiply given vector
		inline Vec3 Tmult(const Vec3& v) const noexcept {
			return Vec3{
				(*this)[0][0] * v[0] + (*this)[1][0] * v[1] + (*this)[2][0] * v[2],
				(*this)[0][1] * v[0] + (*this)[1][1] * v[1] + (*this)[2][1] * v[2],
				(*this)[0][2] * v[0] + (*this)[1][2] * v[1] + (*this)[2][2] * v[2]
			};
		}

		/* Legacy */
		/*inline static void mult(const Mat3& M, const Vec3& V, Vec3& res) noexcept {
			res[0] = M[0][0] * V[0] + M[0][1] * V[1] + M[0][2] * V[2];
			res[1] = M[1][0] * V[0] + M[1][1] * V[1] + M[1][2] * V[2];
			res[2] = M[2][0] * V[0] + M[2][1] * V[1] + M[2][2] * V[2];
		}
		inline static void mult(const Mat3& M1, const Mat3& M2, Mat3& res) noexcept {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					res[i][j] = 0.0;
					for (int k = 0; k < 3; k++)
						res[i][j] += M1[i][k] * M2[k][j];
				}
		}
		inline static void Tmult(const Mat3& M, const Vec3& V, Vec3& res) noexcept {
			res[0] = M[0][0] * V[0] + M[1][0] * V[1] + M[2][0] * V[2];
			res[1] = M[0][1] * V[0] + M[1][1] * V[1] + M[2][1] * V[2];
			res[2] = M[0][2] * V[0] + M[1][2] * V[1] + M[2][2] * V[2];
		}
		inline static void Tmult(const Mat3& M1, const Mat3& M2, Mat3& res) noexcept {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					res[i][j] = 0.0;
					for (int k = 0; k < 3; k++)
						res[i][j] += M1[k][i] * M2[k][j];;
				}
		}*/
	};

	// Quaternion : Simple quaternion for representing rotations
	class Quat {
	public:
		Real x;
		Real y;
		Real z;
		Real w;

		Quat(Real x = 0.0, Real y = 0.0, Real z = 0.0, Real w = 0.0) : x(x), y(y), z(z), w(w) {

		}
		inline void set(Real x = 0.0, Real y = 0.0, Real z = 0.0, Real w = 0.0) {
			this->x = x;
			this->y = y;
			this->z = z;
			this->w = w;
		}

		// Quaternion multiplication
		inline Quat operator*(const Quat& other) const {
			return
				Quat(y * other.z - z * other.y + other.w * x + w * other.x,
					z * other.x - y * other.z + other.w * y + w * other.y,
					x * other.y - y * other.x + other.w * z + w * other.z,
					w * other.w - x * other.x - y * other.y - z * other.z);
		}
		inline void operator*=(const Quat& other) {
			set(y * other.z - z * other.y + other.w * x + w * other.x,
				z * other.x - y * other.z + other.w * y + w * other.y,
				x * other.y - y * other.x + other.w * z + w * other.z,
				w * other.w - x * other.x - y * other.y - z * other.z);
		}
		// L2 norm of this quaternion
		inline Real len() const {
			return sqrt(lensq());
		}
		// Square of L2 norm of this quaternion
		inline Real lensq() const {
			return x * x + y * y + z * z + w * w;
		}
		// Conjugate of this quaternion
		inline Quat conjugate() const {
			return Quat(-x, -y, -z, w);
		}
		// Inverse of this quaternion
		inline Quat inverse() const {
			return (conjugate() * (1.0 / lensq()));
		}

		// Find quaternion for given rotation
		// @ axis : Rotation axis goes through the origin, pointing into this direction
		// @ radian : Radian angle to rotate, following right hand rule
		inline static Quat rotation(const Vec3& axis, Real radian) {
			Real
				phi = radian * 0.5,
				sinp = sin(phi),
				cosp = cos(phi);
			Vec3 naxis = axis;
			naxis.normalize();
			naxis *= sinp;
			return Quat(naxis[0], naxis[1], naxis[2], cosp);
		}
		// Convert this quaternion to 3 x 3 matrix
		inline Mat3 mat3() const {
			Mat3 mat;
			Real
				length = len(),
				s = 2.0 / length,
				xx = x * x,
				xy = x * y,
				xz = x * z,
				xw = x * w,
				yy = y * y,
				yz = y * z,
				yw = y * w,
				zz = z * z,
				zw = z * w;
			mat[0][0] = 1.0 - s * (yy + zz);
			mat[0][1] = s * (xy - zw);
			mat[0][2] = s * (xz + yw);
			mat[1][0] = s * (xy + zw);
			mat[1][1] = 1.0 - s * (xx + zz);
			mat[1][2] = s * (yz - xw);
			mat[2][0] = s * (xz - yw);
			mat[2][1] = s * (yz + xw);
			mat[2][2] = 1.0 - s * (yy + xx);
			return mat;
		}
	};

	// Transform : Represents a configuration of an object in 3 dimensional space
	//             It is comprised of a translation vector [ T ], and a rotation matrix [ R ]
	class Transform {
	public:
		using Translation = Vec3;
		using Rotation = Mat3;
		Translation T;
		Rotation	R;

		Transform() {
			clear();
		}
		// Clear transformation information
		// [ T ] becomes zero vector
		// [ R ] becomes identity matrix
		inline void clear() {
			T = Vec3::zero();
			R = Mat3::identity();
		}

		// Apply this transformation to given point
		inline Vec3 apply(const Vec3& pt) const noexcept {
			return R * pt + T;
		}
		// Apply only translation of this transformation to given point
		inline Vec3 applyT(const Vec3& pt) const noexcept {
			return pt + T;
		}
		// Apply only rotation of this transformation to given point
		inline Vec3 applyR(const Vec3& pt) const noexcept {
			return R * pt;
		}

		// Update this transform by applying [ follow ] transform to it
		inline void update(const Transform& follow) noexcept {
			// R2(R1x + T1) + T2 
			// -> (R2R1)x + (R2T1 + T2)
			R = follow.R * R;
			T = follow.R * T + follow.T;
		}
		// Update this transform by applying [ follow ] transform and save it to [ res ]
		inline void update(const Transform& follow, Transform& res) const noexcept {
			res.R = follow.R * R;
			res.T = follow.R * T + follow.T;
		}
		// Create transform with only translation vector (Rotation matrix remains identity matrix)
		inline static Transform translation(const Vec3& vec) noexcept {
			Transform tr;
			tr.T = vec;
			return tr;
		}
		// Create transform with only rotation information (Translation vector remains zero vector)
		// @a, b : Rotation axis goes through from [ a ] to [ b ] (Starting point : [ a ], Direction vector : [ b - a ])
		// @radian : Radian angle to rotate, following right hand rule
		inline static Transform rotation(const Vec3& a, const Vec3& b, Real radian) noexcept {
			Transform tr;

			Vec3 axis = b - a, a_ = a * -1.0;
			Mat3 rot = Quat::rotation(axis, radian).mat3();

			Transform follow;
			follow.R = rot;
			follow.T = a;
			tr.translate(a_);	// Have to translate [ -a ] since we rotated through axis containing origin, not [ a ]
			tr.update(follow);

			return tr;
		}
		// Translate this transform by given vector
		inline void translate(const Vec3& vec) noexcept {
			T += vec;
		}
		// Rotate this transform with given information
		// @a, b : Rotation axis goes through from [ a ] to [ b ] (Starting point : [ a ], Direction vector : [ b - a ])
		// @radian : Radian angle to rotate, following right hand rule
		inline void rotate(const Vec3& a, const Vec3& b, Real radian) {
			Transform follow = rotation(a, b, radian);
			update(follow);
		}

		// Given transform [ A ], [ B ], get transform that moves a point in local coordinates in [ A ] to local coordinates in [ B ]
		inline static Transform connect(const Transform& A, const Transform& B) noexcept {
			Transform tr;
			tr.T = B.R.Tmult(A.T - B.T);
			tr.R = B.R.Tmult(A.R);
			return tr;
		}
		// Get inverse transform of this transform
		inline Transform inverse() const noexcept {
			Transform I;
			I.R = R.transpose();
			I.T = I.R * T * -1.0;
			return I;
		}
	};

	// Binomial : Simple structure to find binomial values
	//          B(0,0)
	//          /     \
	//       B(1,0)	  B(1,1)
	//      /     \  /     \
	//   B(2,0)  B(2,1)   B(2,2)
	class Binomial {
	private:
		Binomial() = default;
		std::vector<Real> N;
	public:
		// Max degree this binomial structure has
		// e.g. If degree == 3, B(0,0) to B(3,3) is stored in this structure
		int degree = 0;

		// Create binomial structure for given degree
		inline static Binomial create(int degree) {
			Binomial bin;
			int order = degree + 1;
			bin.degree = degree;
			bin.N.resize(order * order);

			for (int n = 0; n < order; n++) {
				for (int i = 0; i <= n; i++) {
					if (i == 0 || i == n) bin.N[bin.index(n, i)] = 1.0;
					else bin.N[bin.index(n, i)] = bin.N[bin.index(n - 1, i - 1)] + bin.N[bin.index(n - 1, i)];
				}
			}
			return bin;
		}
		// Find index into [ N ] for given (n, i) 
		inline int index(int n, int i) const {
			return n * (degree + 1) + i;
		}
		// Get binomial value B(n, i)
		inline const Real& at(int n, int i) const {
			if (n < 0 || n > degree)
				throw(std::runtime_error("[ n ] is not valid for binomial value"));
			return N[index(n, i)];
		}
	};

	/* Miscellaneous operations */

	// Check if [ val ] is between [ a ] and [ b ]
	// @include : Whether or not to include [ a ] and [ b ] for test
	inline bool isbet(Real a, Real b, Real val, bool include = true) {
		Real det = (val - a) * (val - b);
		if (det > 0) return false;
		else if (det < 0) return true;
		else return include;	// det == 0
	}
	inline Real toRadian(Real angle) {
		return (angle / 180.0) * PI;
	}
	inline Real toAngle(Real radian) {
		return (radian / PI) * 180.0;
	}
}

#endif
