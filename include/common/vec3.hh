/* 3D vector classes for integer, real and complex arguments
 *
 * Uses the standard <complex> header, all functions are straightfoward
 * Function declarations included so that they are inlined
 *
 * Created on 26/10/2017
 * Copyright (C) 2017, 2018 Attila Szabó <as2372@cam.ac.uk>
 * Based largely on the Cavlib vec3 class (credits John Richer),
 * <iostream> related functions and rotation omitted
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Creative Commons Attribution License (CC-BY),
 * version 4.0 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please find a copy of the Creative Commons CC-BY License on
 * <https://creativecommons.org/licenses/by/4.0/>.
 */
 
#ifndef vec3_hh
#define vec3_hh

#include <complex>
#include <iostream>

//--------------------- INTEGER VECTOR ----------------------------------------

class vec3_int {
    // Cartesian components
    int m_x[3];
    
    public:
    // Construction from components
    vec3_int (int x, int y, int z) {
        m_x[0] = x;
        m_x[1] = y;
        m_x[2] = z;
    }
    
    // Empty constructor
    vec3_int () {
        m_x[0] = 0;
        m_x[1] = 0;
        m_x[2] = 0;
    }
    
    // Assignment
    vec3_int& operator= (const vec3_int& t) {
        m_x[0] = t.x();
        m_x[1] = t.y();
        m_x[2] = t.z();
        return *this;
    }
    
    // Copy constructor
    vec3_int (const vec3_int& t) {*this = t;}
    
    // Access components, not modifiable
    int x() const { return m_x[0]; }
    int y() const { return m_x[1]; }
    int z() const { return m_x[2]; }
    
    // Access components, modifiable
    int& operator[]( int i ) { return m_x[i]; }
    int operator[]( int i ) const { return m_x[i]; }
    
    // Length squared
    int len2() const { 
        return m_x[0] * m_x[0] + m_x[1] * m_x[1] + m_x[2] * m_x[2];
    }
  
    // Length
    double len() const {return sqrt(len2());}

    
    // += type operators
    // /= PERFORMS INTEGER DIVISION
    vec3_int& operator/= (int f) {
        m_x[0] /= f;
        m_x[1] /= f;
        m_x[2] /= f;
        return *this;
    }
  
    vec3_int& operator*= (int f) {
        m_x[0] *= f;
        m_x[1] *= f;
        m_x[2] *= f;
        return *this;
    }
  
    vec3_int& operator+= (const vec3_int& v) {
        m_x[0] += v.x();
        m_x[1] += v.y();
        m_x[2] += v.z();
        return *this;
    }
  
    vec3_int& operator-= (const vec3_int& v) {
        m_x[0] -= v.x();
        m_x[1] -= v.y();
        m_x[2] -= v.z();
        return *this;
    }

    vec3_int& wrap(int b) {
        m_x[0] = (m_x[0] %= b) >= 0 ? m_x[0] : m_x[0] + b;
        m_x[1] = (m_x[1] %= b) >= 0 ? m_x[1] : m_x[1] + b;
        m_x[2] = (m_x[2] %= b) >= 0 ? m_x[2] : m_x[2] + b;
        return *this;
    }

    vec3_int& wrap(const vec3_int& b) {
        m_x[0] = (m_x[0] %= b[0]) >= 0 ? m_x[0] : m_x[0] + b[0];
        m_x[1] = (m_x[1] %= b[1]) >= 0 ? m_x[1] : m_x[1] + b[1];
        m_x[2] = (m_x[2] %= b[2]) >= 0 ? m_x[2] : m_x[2] + b[2];
        return *this;
    }

    
    /**
     * @brief Reads from stream, expects `[ d1 d2 d2 ]`
     * 
     * @param i 
     * @return vec3_int 
     */
    static vec3_int read_from(std::istream& i) {
        int x,y,z;
        if ((i>>std::ws).peek() != '[' || !i.ignore() || !(i >> std::ws >> x) ){
            throw std::runtime_error("unexpected tokens; expected '[ <x>'");
        }
        if ( !(i >> std::ws >> y)){
            throw std::runtime_error ("unexpected tokens; expected '<x> <y>'");
        }
        if ( !(i >> std::ws >> z)){
            throw std::runtime_error ("unexpected tokens; expected '<y> <z>'");
        }
        if ((i >> std::ws).peek() != ']' || !i.ignore() ){
            throw std::runtime_error ("unexpected tokens; expected ' <z> ]'");
        }
        return vec3_int(x,y,z);
    }
}; // End of class

// Equality
inline bool operator==(const vec3_int& a, const vec3_int& b) {
    return (a.x() == b.x()) && (a.y() == b.y()) && (a.z() == b.z());
}
inline bool operator!=(const vec3_int& a, const vec3_int& b) {
    return !(a==b);
}

// Negation
inline vec3_int operator- (const vec3_int& v) {
    return vec3_int(-v.x(), -v.y(), -v.z());
}

// +/-
inline vec3_int operator+ (const vec3_int& v1, const vec3_int& v2) {
    return vec3_int(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

inline vec3_int operator- (const vec3_int& v1, const vec3_int& v2) {
    return vec3_int(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

/* Multiplication and division with integer scalar
 * NB it's sometimes useful to be able to divide with an integer as per
 * C integer division.
 * However, it is a dangerous operation, handle with care */
inline vec3_int operator* (const vec3_int& v, int f) {
    return vec3_int(f * v.x(), f * v.y(), f * v.z());
}

inline vec3_int operator* (int f, const vec3_int& v) {
    return vec3_int(f * v.x(), f * v.y(), f * v.z());
}

inline vec3_int operator/ (const vec3_int& v, int f) {
    return vec3_int(v.x() / f, v.y() / f, v.z() / f);
}

// Inner product
inline int operator% (const vec3_int& v1, const vec3_int& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

// Cross product
inline vec3_int operator* (const vec3_int& v1, const vec3_int& v2) {
    return vec3_int(v1.y() * v2.z() - v1.z() * v2.y(),
                    v1.z() * v2.x() - v1.x() * v2.z(),
                    v1.x() * v2.y() - v1.y() * v2.x());
}

// entrywise multiplication
inline vec3_int mult(const vec3_int& v1, const vec3_int& v2){
    return vec3_int(
        v1.x() * v2.x(),
        v1.y() * v2.y(),
        v1.z() * v2.z()
    );
}

//--------------------- REAL VECTOR -------------------------------------------
class vec3 {
    // Cartesian components
    double m_x[3];
    
    public:
    // Construction from components
    vec3 (double x, double y, double z) {
        m_x[0] = x;
        m_x[1] = y;
        m_x[2] = z;
    }
    
    // Empty constructor
    vec3 () {
        m_x[0] = 0.0;
        m_x[1] = 0.0;
        m_x[2] = 0.0;
    }

    // Construction from complex number (z component is 0)
    vec3 (const std::complex<double>& z) {
        m_x[0] = z.real();
        m_x[1] = z.imag();
        m_x[2] = 0.0;
    }
    
    // Assignment
    vec3& operator= (const vec3& t) {
        m_x[0] = t.x();
        m_x[1] = t.y();
        m_x[2] = t.z();
        return *this;
    }
    
    vec3& operator= (const vec3_int& t) {
        m_x[0] = t.x();
        m_x[1] = t.y();
        m_x[2] = t.z();
        return *this;
    }
    
    // Copy constructor
    vec3 (const vec3& t) {*this = t;}

    // Constructor from integer vector
    vec3 (const vec3_int& t) {*this = t;}

    // Access components, not modifiable
    double x() const { return m_x[0]; }
    double y() const { return m_x[1]; }
    double z() const { return m_x[2]; }
    
    // Access components, modifiable
    double& operator[]( int i ) { return m_x[i]; }
    double operator[]( int i ) const { return m_x[i]; }

    /* XY components as a complex number
     * Horrible hack exploiting the memory representation of vec3 and complex */
    std::complex<double>& xy() { return *((std::complex<double>*)this); }
    const std::complex<double>& xy() const {
        return *((std::complex<double>*)this);
        // return *reinterpret_cast<const std::complex<double>* >(this);
    }
    
    // Length squared
    double len2() const { 
        return m_x[0] * m_x[0] + m_x[1] * m_x[1] + m_x[2] * m_x[2];
    }
  
    // Length
    double len() const {return sqrt(len2());}
    
    // Normalisation in place
    vec3& normalise() {return *this /= len();}
    
    // += type operators
    vec3& operator/= (double f) {
        m_x[0] /= f;
        m_x[1] /= f;
        m_x[2] /= f;
        return *this;
    }
  
    vec3& operator*= (double f) {
        m_x[0] *= f;
        m_x[1] *= f;
        m_x[2] *= f;
        return *this;
    }
  
    vec3& operator+= (const vec3& v) {
        m_x[0] += v.x();
        m_x[1] += v.y();
        m_x[2] += v.z();
        return *this;
    }
  
    vec3& operator-= (const vec3& v) {
        m_x[0] -= v.x();
        m_x[1] -= v.y();
        m_x[2] -= v.z();
        return *this;
    }

    // static vec3 read_from(std::istream& i) {
    //     double x,y,z;
    //     if ((i>>std::ws).peek() != '[' || !i.ignore() || !(i >> std::ws >> x) ){
    //         throw std::runtime_error("unexpected tokens; expected '[ <x>'");
    //     }
    //     if ((i >> std::ws).peek() != ',' || !i.ignore() || !(i >> std::ws >> y)){
    //         throw std::runtime_error ("unexpected tokens; expected '<x>, <y>'");
    //     }
    //     if ((i >> std::ws).peek() != ',' || !i.ignore() || !(i >> std::ws >> z)){
    //         throw std::runtime_error ("unexpected tokens; expected '<y>, <z>'");
    //     }
    //     if ((i >> std::ws).peek() != ']' || !i.ignore() ){
    //         throw std::runtime_error ("unexpected tokens; expected ' <z> ]'");
    //     }
    //     return vec3(x,y,z);
    // }
}; // End of class

// Equality
inline bool operator==(const vec3& a, const vec3& b) {
    return (a.x() == b.x()) && (a.y() == b.y()) && (a.z() == b.z());
}
inline bool operator!=(const vec3& a, const vec3& b) {
    return !(a==b);
}

// Negation
inline vec3 operator- (const vec3& v) {
    return vec3(-v.x(), -v.y(), -v.z());
}

// +/-
inline vec3 operator+ (const vec3& v1, const vec3& v2) {
    return vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

inline vec3 operator- (const vec3& v1, const vec3& v2) {
    return vec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

// Multiplication with scalar
inline vec3 operator* (const vec3& v, double f) {
    return vec3(f * v.x(), f * v.y(), f * v.z());
}

inline vec3 operator* (double f, const vec3& v) {
    return vec3(f * v.x(), f * v.y(), f * v.z());
}

inline vec3 operator/ (const vec3& v, double f) {
    return vec3(v.x() / f, v.y() / f, v.z() / f);
}

inline vec3 operator* (const vec3_int& v, double f) {
    return vec3(f * v.x(), f * v.y(), f * v.z());
}

inline vec3 operator* (double f, const vec3_int& v) {
    return vec3(f * v.x(), f * v.y(), f * v.z());
}

inline vec3 operator/ (const vec3_int& v, double f) {
    return vec3(v.x() / f, v.y() / f, v.z() / f);
}

// Inner product
inline double operator% (const vec3& v1, const vec3& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline double operator% (const vec3_int& v1, const vec3& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline double operator% (const vec3& v1, const vec3_int& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

// Cross product
inline vec3 operator* (const vec3& v1, const vec3& v2) {
    return vec3(v1.y() * v2.z() - v1.z() * v2.y(),
                v1.z() * v2.x() - v1.x() * v2.z(),
                v1.x() * v2.y() - v1.y() * v2.x());
}

inline vec3 operator* (const vec3_int& v1, const vec3& v2) {
    return vec3(v1.y() * v2.z() - v1.z() * v2.y(),
                v1.z() * v2.x() - v1.x() * v2.z(),
                v1.x() * v2.y() - v1.y() * v2.x());
}

inline vec3 operator* (const vec3& v1, const vec3_int& v2) {
    return vec3(v1.y() * v2.z() - v1.z() * v2.y(),
                v1.z() * v2.x() - v1.x() * v2.z(),
                v1.x() * v2.y() - v1.y() * v2.x());
}

// Return a normalised vector
inline vec3 normalise (const vec3& v) {
    return v / v.len();
}

inline vec3 normalise (const vec3_int& v) {
    return v / v.len();
}

//--------------------- COMPLEX VECTOR ----------------------------------------

class vec3_cplx {
    // Cartesian components
    std::complex<double> m_x[3];
    
    public:
    // Construction from components
    vec3_cplx (std::complex<double> x, std::complex<double> y, std::complex<double> z) {
        m_x[0] = x;
        m_x[1] = y;
        m_x[2] = z;
    }
    
    // Empty constructor
    vec3_cplx () {
        m_x[0] = 0.0;
        m_x[1] = 0.0;
        m_x[2] = 0.0;
    }
    
    // Assignment
    vec3_cplx& operator= (const vec3_cplx& t) {
        m_x[0] = t.x();
        m_x[1] = t.y();
        m_x[2] = t.z();
        return *this;
    }
    vec3_cplx& operator= (const vec3& t) {
        m_x[0] = t.x();
        m_x[1] = t.y();
        m_x[2] = t.z();
        return *this;
    }
    vec3_cplx& operator= (const vec3_int& t) {
        m_x[0] = (double)(t.x());
        m_x[1] = (double)(t.y());
        m_x[2] = (double)(t.z());
        return *this;
    }
    
    // Copy constructor
    vec3_cplx (const vec3_cplx& t) {*this = t;}
    
    // Constructor out of a real vector
    vec3_cplx (const vec3& t) {*this = t;}
    
    // Constructor out of an integer vector
    vec3_cplx (const vec3_int& t) {*this = t;}
    
    // Access components, not modifiable
    std::complex<double> x() const { return m_x[0]; }
    std::complex<double> y() const { return m_x[1]; }
    std::complex<double> z() const { return m_x[2]; }
    
    // Access components, modifiable
    std::complex<double>& operator[]( int i ) { return m_x[i]; }
    std::complex<double> operator[]( int i )  const { return m_x[i]; }

    // Real & imag. parts
    vec3 real() const {
	return vec3( x().real(), y().real(), z().real() );
    }
    vec3 imag() const {
	return vec3( x().imag(), y().imag(), z().imag() );
    }
    
    // Length squared
    double len2() const { 
        return std::norm(m_x[0]) + std::norm(m_x[1]) + std::norm(m_x[2]);
    }
  
    // Length
    double len() const {return sqrt(len2());}
    
    // Normalisation in place
    vec3_cplx& normalise() {return *this /= len();}
    
    // += type operators
    vec3_cplx& operator/= (std::complex<double> f) {
        m_x[0] /= f;
        m_x[1] /= f;
        m_x[2] /= f;
        return *this;
    }
  
    vec3_cplx& operator*= (std::complex<double> f) {
        m_x[0] *= f;
        m_x[1] *= f;
        m_x[2] *= f;
        return *this;
    }
  
    vec3_cplx& operator+= (const vec3_cplx& v) {
        m_x[0] += v.x();
        m_x[1] += v.y();
        m_x[2] += v.z();
        return *this;
    }
  
    vec3_cplx& operator-= (const vec3_cplx& v) {
        m_x[0] -= v.x();
        m_x[1] -= v.y();
        m_x[2] -= v.z();
        return *this;
    }

    // static vec3_cplx read_from(std::istream& i) {
    //     std::complex<double> x,y,z;
    //     if ((i>>std::ws).peek() != '[' || !i.ignore() || !(i >> std::ws >> x) ){
    //         throw std::runtime_error("unexpected tokens; expected '[ <x>'");
    //     }
    //     if ((i >> std::ws).peek() != ',' || !i.ignore() || !(i >> std::ws >> y)){
    //         throw std::runtime_error ("unexpected tokens; expected '<x>, <y>'");
    //     }
    //     if ((i >> std::ws).peek() != ',' || !i.ignore() || !(i >> std::ws >> z)){
    //         throw std::runtime_error ("unexpected tokens; expected '<y>, <z>'");
    //     }
    //     if ((i >> std::ws).peek() != ']' || !i.ignore() ){
    //         throw std::runtime_error ("unexpected tokens; expected ' <z> ]'");
    //     }
    //     return vec3_cplx(x,y,z);
    // }
}; // End of class

// Equality
inline bool operator==(const vec3_cplx& a, const vec3_cplx& b) {
    return (a.x() == b.x()) && (a.y() == b.y()) && (a.z() == b.z());
}
inline bool operator!=(const vec3_cplx& a, const vec3_cplx& b) {
    return !(a==b);
}

// Negation
inline vec3_cplx operator- (const vec3_cplx& v) {
    return vec3_cplx(-v.x(), -v.y(), -v.z());
}

// +/-
inline vec3_cplx operator+ (const vec3_cplx& v1, const vec3_cplx& v2) {
    return vec3_cplx(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

inline vec3_cplx operator- (const vec3_cplx& v1, const vec3_cplx& v2) {
    return vec3_cplx(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

// Multiplication with scalar
inline vec3_cplx operator* (const vec3_cplx& v, std::complex<double> f) {
    return vec3_cplx(f * v.x(), f * v.y(), f * v.z());
}

inline vec3_cplx operator* (std::complex<double> f, const vec3_cplx& v) {
    return vec3_cplx(f * v.x(), f * v.y(), f * v.z());
}

inline vec3_cplx operator/ (const vec3_cplx& v, std::complex<double> f) {
    return vec3_cplx(v.x() / f, v.y() / f, v.z() / f);
}

inline vec3_cplx operator* (const vec3& v, std::complex<double> f) {
    return vec3_cplx(f * v.x(), f * v.y(), f * v.z());
}

inline vec3_cplx operator* (std::complex<double> f, const vec3& v) {
    return vec3_cplx(f * v.x(), f * v.y(), f * v.z());
}

inline vec3_cplx operator/ (const vec3& v, std::complex<double> f) {
    return vec3_cplx(v.x() / f, v.y() / f, v.z() / f);
}

inline vec3_cplx operator* (const vec3_int& v, std::complex<double> f) {
    return vec3_cplx(f * (double)v.x(),
                     f * (double)v.y(),
                     f * (double)v.z());
}

inline vec3_cplx operator* (std::complex<double> f, const vec3_int& v) {
    return vec3_cplx(f * (double)v.x(),
                     f * (double)v.y(),
                     f * (double)v.z());
}

inline vec3_cplx operator/ (const vec3_int& v, std::complex<double> f) {
    return vec3_cplx((double)v.x() / f,
                     (double)v.y() / f,
                     (double)v.z() / f);
}

// Inner product
inline std::complex<double> operator% (const vec3_cplx& v1, const vec3_cplx& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

// Inner product of real and complex vectors
inline std::complex<double> operator% (const vec3& v1, const vec3_cplx& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline std::complex<double> operator% (const vec3_cplx& v1, const vec3& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

// Cross product of complex vectors
inline vec3_cplx operator* (const vec3_cplx& v1, const vec3_cplx& v2) {
    return vec3_cplx(v1.y() * v2.z() - v1.z() * v2.y(),
       	             v1.z() * v2.x() - v1.x() * v2.z(),
                     v1.x() * v2.y() - v1.y() * v2.x());
}

// Cross product of real and complex vectors
inline vec3_cplx operator* (const vec3& v1, const vec3_cplx& v2) {
    return vec3_cplx(v1.y() * v2.z() - v1.z() * v2.y(),
       	             v1.z() * v2.x() - v1.x() * v2.z(),
                     v1.x() * v2.y() - v1.y() * v2.x());
}

inline vec3_cplx operator* (const vec3_cplx& v1, const vec3& v2) {
    return vec3_cplx(v1.y() * v2.z() - v1.z() * v2.y(),
       	             v1.z() * v2.x() - v1.x() * v2.z(),
                     v1.x() * v2.y() - v1.y() * v2.x());
}

// Return a normalised vector
inline vec3_cplx normalise (const vec3_cplx& v) {
    return v / v.len();
}


// stream operators for vanity
inline std::ostream& operator<< (std::ostream &o, const vec3_int& v){
    o << "[ " << v[0] << " " << v[1] << " "<<v[2]<<" ]";
    return o;
}
inline std::ostream& operator<< (std::ostream &o, const vec3& v){
    o << "[ " << v[0] << " " << v[1] << " "<<v[2]<<" ]";
    return o;
}
inline std::ostream& operator<< (std::ostream &o, const vec3_cplx& v){
    o << "[ " << v[0] << " " << v[1] << " "<<v[2]<<" ]";
    return o;
}


#endif
