/* -----------------------------------------------------------------------------

  SnOB - An FFT toolkit for the symmetric group
              development version

  Copyright (C) 2006  Imre Risi Kondor


  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the

  Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor,
  Boston, MA  02110-1301, USA.

  This software is provided for educational and research purposes.
  Commercial use is prohibited.

  See the accompanying LICENSE for details

----------------------------------------------------------------------------- */

#pragma once

#include <stdarg.h>
#include <vector>

#include "Sn.hpp"

#include "Cycles.hpp"

using namespace std;

class Sn::Element : public FiniteGroup::Element {
public:
  Element(const Sn &_group) : Element(_group.n) {}

  Element(const int _n) : n(_n) {
    p = new int[n];
    pinv = new int[n];
    for (auto i = 0; i < n; i++) {
      p[i] = i + 1;
      pinv[i] = i + 1;
    }
  }

  Element(int a1, int a2, ...) {
    va_list params;
    int arg = a2;
    vector<int> v;
    v.push_back(a1);
    va_start(params, a2);
    while (arg != 0) {
      v.push_back(arg);
      arg = va_arg(params, int);
    }
    va_end(params);
    n = v.size();
    p = new int[n];
    pinv = new int[n];
    for (auto i = 0; i < n; i++)
      p[i] = v[i];
    for (auto i = 0; i < n; i++)
      pinv[p[i] - 1] = i + 1;
  }

  Element(const int _n, int *v) : n(_n) {
    p = new int[n];
    pinv = new int[n];
    for (auto i = 0; i < n; i++)
      p[i] = v[i];
    for (auto i = 0; i < n; i++)
      pinv[p[i] - 1] = i + 1;
  }

  // Creates a coset representative
  Element(const int _n, const vector<int> fixed) : n(_n) {
    p = new int[n];
    pinv = new int[n];
    for (auto i = 0; i < n; i++)
      p[i] = i + 1;
    for (unsigned j = 0; j < fixed.size(); j++) {
      for (auto i = 0; i < n; i++)
        if (p[i] == fixed[j]) {
          p[i] = p[n - 1 - j];
          p[n - 1 - j] = fixed[j];
          break;
        }
    }
    for (auto i = 0; i < n; i++)
      pinv[p[i] - 1] = i + 1;
  }

  Element(const vector<int> &factorization, const int _n) : n(_n) {
    p = new int[n];
    pinv = new int[n];
    for (auto i = 0; i < n; i++)
      p[i] = i + 1;
    for (unsigned j = 0; j < factorization.size(); j++) {
      // cout<<":"<<factorization[j];
      int t = p[n - 1 - j];
      p[n - 1 - j] = p[factorization[j] - 1];
      p[factorization[j] - 1] = t;
    }
    cout << endl;
    for (auto i = 0; i < n; i++)
      pinv[p[i] - 1] = i + 1;
  }

  Element(const Sn::Element &o) : n(o.n) {
    p = new int[n];
    for (auto i = 0; i < n; i++)
      p[i] = o.p[i];
    pinv = new int[n];
    for (auto i = 0; i < n; i++)
      pinv[i] = o.pinv[i];
  }

  ~Element() {
    delete[] p;
    delete[] pinv;
  }

  bool operator==(const Sn::Element &o) const {
    if (n != o.n)
      return false;
    for (auto i = 0; i < n; i++)
      if (p[i] != o.p[i])
        return false;
    return true;
  }

  int action(const int i) const { return pinv[i - 1]; }
  int iaction(const int i) const { return p[i - 1]; }

  vector<int> effect() const {
    vector<int> result;
    for (auto i = 0; i < n; i++)
      result.push_back(pinv[i]);
    return result;
  }
  vector<int> ieffect() const {
    vector<int> result;
    for (auto i = 0; i < n; i++)
      result.push_back(p[i]);
    return result;
  }

  Element *operator*(const Sn::Element &o) const {
    Element *result = new Element(n);
    for (auto i = 0; i < n; i++)
      result->pinv[i] = pinv[o.pinv[i] - 1];
    for (auto i = 0; i < n; i++)
      result->p[result->pinv[i] - 1] = i + 1;
    return result;
  }

  Element *inverse() { return new Element(n, pinv); }

  // Element &CcycleL(int j, int q);

  Sn::Element &CcycleL(int j, int q) {
    int t = p[q - 1];
    for (auto i = j + 1; i <= q; i++)
      p[i - 1] = p[i - 2];
    p[j - 1] = t;
    for (auto i = 0; i < n; i++)
      pinv[p[i] - 1] = i + 1;
    return *this;
  }

  Sn::Element &CcycleR(int j, int q) {
    int t = pinv[j - 1];
    for (auto i = j; i <= q - 1; i++)
      pinv[i - 1] = pinv[i + 1 - 1];
    pinv[q - 1] = t;
    for (auto i = 0; i < n; i++)
      p[pinv[i] - 1] = i + 1;
    return *this;
  }

  string str() const {
    ostringstream result;
    result << "[ ";
    for (auto i = 0; i < n; i++)
      result << p[i] << " ";
    result << "]";
    return result.str();
  }

  int n;

private:
  int *p; // p=[\sigma^{-1}(1),....,\sigma^{-1}(n)]
  int *pinv;
};
