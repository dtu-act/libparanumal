/* The MIT License (MIT)
 *
 * Copyright (c) 2014-2018 David Medina and Tim Warburton
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 */

#ifndef OCCA_TYPES_HEADER
#define OCCA_TYPES_HEADER

#include <iostream>
#include <vector>
#include <map>
#include <stdint.h>

#include <occa/defines.hpp>
#include <occa/vector.hpp>

namespace occa {
  typedef int64_t dim_t;
  typedef uint64_t udim_t;

  typedef std::vector<int>                   intVector;
  typedef std::vector<intVector>             intVecVector;

  typedef std::vector<std::string>           strVector;
  typedef strVector::iterator                strVectorIterator;
  typedef strVector::const_iterator          cStrVectorIterator;

  typedef std::map<std::string, std::string> strToStrMap;
  typedef strToStrMap::iterator              strToStrMapIterator;
  typedef strToStrMap::const_iterator        cStrToStrMapIterator;

  typedef std::map<std::string,strVector>    strToStrsMap;
  typedef strToStrsMap::iterator             strToStrsMapIterator;
  typedef strToStrsMap::const_iterator       cStrToStrsMapIterator;

  typedef std::map<std::string, bool>        strToBoolMap;
  typedef strToBoolMap::iterator             strToBoolMapIterator;
  typedef strToBoolMap::const_iterator       cStrToBoolMapIterator;

  //---[ Dim ]--------------------------
  class dim {
  public:
    int dims;
    udim_t x, y, z;

    dim();
    dim(udim_t x_);
    dim(udim_t x_, udim_t y_);
    dim(udim_t x_, udim_t y_, udim_t z_);
    dim(int dims_, udim_t x_, udim_t y_, udim_t z_);

    dim(const dim &d);

    dim& operator = (const dim &d);

    dim operator + (const dim &d) const;
    dim operator - (const dim &d) const;
    dim operator * (const dim &d) const;
    dim operator / (const dim &d) const;

    bool hasNegativeEntries();

    udim_t& operator [] (int i);
    udim_t  operator [] (int i) const;
  };
  //====================================

  //---[ Bitfield ]---------------------
  class bitfield {
  public:
    udim_t b1, b2;

    inline bitfield() :
      b1(0),
      b2(0) {}

    inline bitfield(const udim_t b2_) :
      b1(0),
      b2(b2_) {}

    inline bitfield(const udim_t b1_,
                    const udim_t b2_) :
      b1(b1_),
      b2(b2_) {}

    inline bitfield(const bitfield &bf) :
      b1(bf.b1),
      b2(bf.b2) {}

    inline bitfield& operator = (const bitfield &bf) {
      b1 = bf.b1;
      b2 = bf.b2;
      return *this;
    }

    inline static int bits() {
      return 16 * sizeof(udim_t);
    }

    inline bool operator == (const bitfield &bf) const {
      return ((b1 == bf.b1) &&
              (b2 == bf.b2));
    }

    inline bool operator != (const bitfield &bf) const {
      return ((b1 != bf.b1) ||
              (b2 != bf.b2));
    }

    inline bool operator < (const bitfield &bf) const {
      return ((b1 < bf.b1) ||
              ((b1 == bf.b1) &&
               (b2 < bf.b2)));
    }

    inline bool operator <= (const bitfield &bf) const {
      return ((b1 < bf.b1) ||
              ((b1 == bf.b1) &&
               (b2 <= bf.b2)));
    }

    inline bool operator > (const bitfield &bf) const {
      return ((b1 > bf.b1) ||
              ((b1 == bf.b1) &&
               (b2 > bf.b2)));
    }

    inline bool operator >= (const bitfield &bf) const {
      return ((b1 > bf.b1) ||
              ((b1 == bf.b1) &&
               (b2 >= bf.b2)));
    }

    inline bitfield operator | (const bitfield &bf) const {
      return bitfield(b1 | bf.b1,
                      b2 | bf.b2);
    }

    inline bitfield operator & (const bitfield &bf) const {
      return bitfield(b1 & bf.b1,
                      b2 & bf.b2);
    };

    inline bitfield operator << (const int shift) const {
      if (shift <= 0) {
        return *this;
      }
      const int bSize = 8 * sizeof(udim_t);
      if (shift > (2 * bSize)) {
        return 0;
      }

      if (shift >= bSize) {
        return bitfield(b2 << (shift - bSize),
                        0);
      }

      const udim_t carryOver = (b2 >> (bSize - shift));
      return bitfield((b1 << shift) | carryOver,
                      b2 << shift);
    };

    inline bitfield operator >> (const int shift) const {
      if (shift <= 0) {
        return *this;
      }
      const int bSize = 8 * sizeof(udim_t);
      if (shift > (2 * bSize)) {
        return 0;
      }

      if (shift >= bSize) {
        return bitfield(0,
                        b1 >> (shift - bSize));
      }

      const udim_t carryOver = (b1 << (bSize - shift));
      return bitfield(b1 >> shift,
                      (b2 >> shift) | carryOver);
    };

    inline bitfield& operator <<= (const int shift) {
      *this = (*this << shift);
      return *this;
    }

    inline bitfield& operator >>= (const int shift) {
      *this = (*this >> shift);
      return *this;
    }

    inline operator bool () const {
      return (b1 || b2);
    }
  };
  //====================================

  //---[ Type To String ]---------------
  template <class TM>
  class primitiveinfo {
  public:
    static const std::string id;
    static const std::string name;
    static const bool isUnsigned;
  };

  template <class TM>
  class typeinfo {
  public:
    static const std::string id;
    static const std::string name;
    static const bool isUnsigned;
  };

  typedef signed char    schar_t;
  typedef unsigned char  uchar_t;
  typedef unsigned short ushort_t;
  typedef unsigned int   uint_t;
  typedef unsigned long  ulong_t;

#if !(OCCA_OS & (OCCA_LINUX_OS | OCCA_MACOS_OS))
  template <> const std::string primitiveinfo<char>::id         = "c";
  template <> const std::string primitiveinfo<char>::name       = "char";
  template <> const bool        primitiveinfo<char>::isUnsigned = false;

  template <> const std::string primitiveinfo<short>::id         = "s";
  template <> const std::string primitiveinfo<short>::name       = "short";
  template <> const bool        primitiveinfo<short>::isUnsigned = false;

  template <> const std::string primitiveinfo<int>::id         = "i";
  template <> const std::string primitiveinfo<int>::name       = "int";
  template <> const bool        primitiveinfo<int>::isUnsigned = false;

  template <> const std::string primitiveinfo<long>::id         = "l";
  template <> const std::string primitiveinfo<long>::name       = "long";
  template <> const bool        primitiveinfo<long>::isUnsigned = false;

  template <> const std::string primitiveinfo<schar_t>::id         = "sc";
  template <> const std::string primitiveinfo<schar_t>::name       = "schar_t";
  template <> const bool        primitiveinfo<schar_t>::isUnsigned = false;

  template <> const std::string primitiveinfo<uchar_t>::id         = "uc";
  template <> const std::string primitiveinfo<uchar_t>::name       = "uchar_t";
  template <> const bool        primitiveinfo<uchar_t>::isUnsigned = true;

  template <> const std::string primitiveinfo<ushort_t>::id         = "us";
  template <> const std::string primitiveinfo<ushort_t>::name       = "ushort_t";
  template <> const bool        primitiveinfo<ushort_t>::isUnsigned = true;

  template <> const std::string primitiveinfo<uint_t>::id         = "ui";
  template <> const std::string primitiveinfo<uint_t>::name       = "uint_t";
  template <> const bool        primitiveinfo<uint_t>::isUnsigned = true;

  template <> const std::string primitiveinfo<ulong_t>::id         = "ul";
  template <> const std::string primitiveinfo<ulong_t>::name       = "ulong_t";
  template <> const bool        primitiveinfo<ulong_t>::isUnsigned = true;

  template <> const std::string primitiveinfo<float>::id         = "f";
  template <> const std::string primitiveinfo<float>::name       = "float";
  template <> const bool        primitiveinfo<float>::isUnsigned = false;

  template <> const std::string primitiveinfo<double>::id         = "d";
  template <> const std::string primitiveinfo<double>::name       = "double";
  template <> const bool        primitiveinfo<double>::isUnsigned = false;

  template <> const std::string typeinfo<uint8_t>::id         = "u8";
  template <> const std::string typeinfo<uint8_t>::name       = "uint8";
  template <> const bool        typeinfo<uint8_t>::isUnsigned = true;

  template <> const std::string typeinfo<uint16_t>::id         = "u16";
  template <> const std::string typeinfo<uint16_t>::name       = "uint16";
  template <> const bool        typeinfo<uint16_t>::isUnsigned = true;

  template <> const std::string typeinfo<uint32_t>::id         = "u32";
  template <> const std::string typeinfo<uint32_t>::name       = "uint32";
  template <> const bool        typeinfo<uint32_t>::isUnsigned = true;

  template <> const std::string typeinfo<uint64_t>::id         = "u64";
  template <> const std::string typeinfo<uint64_t>::name       = "uint64";
  template <> const bool        typeinfo<uint64_t>::isUnsigned = true;

  template <> const std::string typeinfo<int8_t>::id         = "i8";
  template <> const std::string typeinfo<int8_t>::name       = "int8";
  template <> const bool        typeinfo<int8_t>::isUnsigned = false;

  template <> const std::string typeinfo<int16_t>::id         = "i16";
  template <> const std::string typeinfo<int16_t>::name       = "int16";
  template <> const bool        typeinfo<int16_t>::isUnsigned = false;

  template <> const std::string typeinfo<int32_t>::id         = "i32";
  template <> const std::string typeinfo<int32_t>::name       = "int32";
  template <> const bool        typeinfo<int32_t>::isUnsigned = false;

  template <> const std::string typeinfo<int64_t>::id         = "i64";
  template <> const std::string typeinfo<int64_t>::name       = "int64";
  template <> const bool        typeinfo<int64_t>::isUnsigned = false;

  template <> const std::string typeinfo<float>::id         = "f32";
  template <> const std::string typeinfo<float>::name       = "float";
  template <> const bool        typeinfo<float>::isUnsigned = false;

  template <> const std::string typeinfo<double>::id         = "f64";
  template <> const std::string typeinfo<double>::name       = "double";
  template <> const bool        typeinfo<double>::isUnsigned = false;

  template <> const std::string typeinfo<uchar2>::id         = "vuc2";
  template <> const std::string typeinfo<uchar2>::name       = "uchar2";
  template <> const bool        typeinfo<uchar2>::isUnsigned = true;

  template <> const std::string typeinfo<uchar4>::id         = "vuc4";
  template <> const std::string typeinfo<uchar4>::name       = "uchar4";
  template <> const bool        typeinfo<uchar4>::isUnsigned = true;

  template <> const std::string typeinfo<char2>::id         = "vc2";
  template <> const std::string typeinfo<char2>::name       = "char2";
  template <> const bool        typeinfo<char2>::isUnsigned = false;

  template <> const std::string typeinfo<char4>::id         = "vc4";
  template <> const std::string typeinfo<char4>::name       = "char4";
  template <> const bool        typeinfo<char4>::isUnsigned = false;

  template <> const std::string typeinfo<ushort2>::id         = "vus2";
  template <> const std::string typeinfo<ushort2>::name       = "ushort2";
  template <> const bool        typeinfo<ushort2>::isUnsigned = true;

  template <> const std::string typeinfo<ushort4>::id         = "vus4";
  template <> const std::string typeinfo<ushort4>::name       = "ushort4";
  template <> const bool        typeinfo<ushort4>::isUnsigned = true;

  template <> const std::string typeinfo<short2>::id         = "vs2";
  template <> const std::string typeinfo<short2>::name       = "short2";
  template <> const bool        typeinfo<short2>::isUnsigned = false;

  template <> const std::string typeinfo<short4>::id         = "vs4";
  template <> const std::string typeinfo<short4>::name       = "short4";
  template <> const bool        typeinfo<short4>::isUnsigned = false;

  template <> const std::string typeinfo<uint2>::id         = "vui2";
  template <> const std::string typeinfo<uint2>::name       = "uint2";
  template <> const bool        typeinfo<uint2>::isUnsigned = true;

  template <> const std::string typeinfo<uint4>::id         = "vui4";
  template <> const std::string typeinfo<uint4>::name       = "uint4";
  template <> const bool        typeinfo<uint4>::isUnsigned = true;

  template <> const std::string typeinfo<int2>::id         = "vi2";
  template <> const std::string typeinfo<int2>::name       = "int2";
  template <> const bool        typeinfo<int2>::isUnsigned = false;

  template <> const std::string typeinfo<int4>::id         = "vi4";
  template <> const std::string typeinfo<int4>::name       = "int4";
  template <> const bool        typeinfo<int4>::isUnsigned = false;

  template <> const std::string typeinfo<ulong2>::id         = "vul2";
  template <> const std::string typeinfo<ulong2>::name       = "ulong2";
  template <> const bool        typeinfo<ulong2>::isUnsigned = true;

  template <> const std::string typeinfo<ulong4>::id         = "vul4";
  template <> const std::string typeinfo<ulong4>::name       = "ulong4";
  template <> const bool        typeinfo<ulong4>::isUnsigned = true;

  template <> const std::string typeinfo<long2>::id         = "vl2";
  template <> const std::string typeinfo<long2>::name       = "long2";
  template <> const bool        typeinfo<long2>::isUnsigned = false;

  template <> const std::string typeinfo<long4>::id         = "vl4";
  template <> const std::string typeinfo<long4>::name       = "long4";
  template <> const bool        typeinfo<long4>::isUnsigned = false;

  template <> const std::string typeinfo<float2>::id         = "vf2";
  template <> const std::string typeinfo<float2>::name       = "float2";
  template <> const bool        typeinfo<float2>::isUnsigned = false;

  template <> const std::string typeinfo<float4>::id         = "vf4";
  template <> const std::string typeinfo<float4>::name       = "float4";
  template <> const bool        typeinfo<float4>::isUnsigned = false;

  template <> const std::string typeinfo<double2>::id         = "vd2";
  template <> const std::string typeinfo<double2>::name       = "double2";
  template <> const bool        typeinfo<double2>::isUnsigned = false;

  template <> const std::string typeinfo<double4>::id         = "vd4";
  template <> const std::string typeinfo<double4>::name       = "double4";
  template <> const bool        typeinfo<double4>::isUnsigned = false;
  //====================================
#endif
}

#endif
