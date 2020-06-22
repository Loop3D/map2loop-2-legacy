#include <pybind11/pybind11.h>

int add(int i, int j) {
  int result = i + j;
  return result;
}

PYBIND11_MODULE(example, m) {
  m.doc() = "pybind11 example plugin"; // optional module docstring

  m.def("add", &add, "A function which adds two numbers");
}
