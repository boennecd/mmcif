#ifndef SIMPLE_MAT_H
#define SIMPLE_MAT_H

#include <memory.h>
#include <algorithm>

/// a simple matrix container in column-major order
template<class T>
class simple_mat {
  std::unique_ptr<T[]> mem{};
  size_t n_rows_v{}, n_cols_v{};
  T * external{nullptr};

public:
  /// use external memory
  simple_mat
  (T * external, size_t const n_rows, size_t const n_cols):
  n_rows_v{n_rows}, n_cols_v{n_cols}, external{external} { }

  /// allocates memory for a given size matrix
  simple_mat(size_t const n_rows, size_t const n_cols):
  mem{new T[n_cols * n_rows]},
  n_rows_v{n_rows}, n_cols_v{n_cols}, external{nullptr} { }

  /// copy constructor which always copies
  simple_mat(const simple_mat &o):
  mem{new T[o.n_rows_v * o.n_cols_v]},
  n_rows_v{o.n_rows_v},
  n_cols_v{o.n_cols_v}
  {
    std::copy(o.begin(), o.end(), mem.get());
  }

  simple_mat& operator=(const simple_mat &o) {
    mem.reset(new T[o.n_rows_v * o.n_cols_v]);
    n_rows_v = o.n_rows_v;
    n_cols_v = o.n_cols_v;
    std::copy(o.begin(), o.end(), mem.get());
    return *this;
  }

  simple_mat(simple_mat &&o):
  mem{std::move(o.mem)},
  n_rows_v{o.n_rows_v},
  n_cols_v{o.n_cols_v},
  external{o.external} { }

  simple_mat& operator=(simple_mat &&o){
    n_rows_v = o.n_rows_v;
    n_cols_v = o.n_cols_v;
    external = o.external;
    mem = std::move(o.mem);
    return *this;
  }

  size_t n_cols() const {
    return n_cols_v;
  }
  size_t n_rows() const {
    return n_rows_v;
  }

  /// returns a pointer to the memory
  T * get(){
    return external ? external : mem.get();
  }
  T const * get() const {
    return external ? external : mem.get();
  }

  /// begin and end end functions
  T * begin() {
    return get();
  }
  T const * begin() const {
    return get();
  }
  T * end() {
    return get() + n_rows() * n_cols();
  }
  T const * end() const {
    return get() + n_rows() * n_cols();
  }

  /// returns the pointer to a given column
  T * col(size_t const idx){
    return get() + idx * n_rows();
  }
  T const * col(size_t const idx) const {
    return get() + idx * n_rows();
  }
};

#endif
