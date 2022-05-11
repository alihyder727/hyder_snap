#ifndef UTILS_UTILS_HPP_
#define UTILS_UTILS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file utils.hpp
//  \brief prototypes of functions and class definitions for utils/*.cpp files

// C headers

// C++ headers
#include <csignal>   // sigset_t POSIX C extension
#include <cstdint>   // std::int64_t
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"

void ChangeRunDir(const char *pdir);
double ran2(std::int64_t *idum);
void ShowConfig();

//----------------------------------------------------------------------------------------
//! SignalHandler
//  \brief static data and functions that implement a simple signal handling system

namespace SignalHandler {
const int nsignal = 3;
static volatile int signalflag[nsignal];
const int ITERM = 0, IINT = 1, IALRM = 2;
static sigset_t mask;
void SignalHandlerInit();
int CheckSignalFlags();
int GetSignalFlag(int s);
void SetSignalFlag(int s);
void SetWallTimeAlarm(int t);
void CancelWallTimeAlarm();
} // namespace SignalHandler

//! test file existance
bool FileExists(std::string fname);

//! test a blank line
bool IsBlankLine(char const* line);
bool IsBlankLine(std::string const& line);

//! decomment a file
std::string DecommentFile(std::string fname);

//! get number of columns in a data table
int GetNumCols(std::string fname, char c = ' ');

//! get number of rows in a data table 
int GetNumRows(std::string fname);

//! split a string to a vector
template<typename A>
std::vector<A> Vectorize(const char* cstr, const char *delimiter = " ")
{
  std::vector<A> arr;
  char str[1028], *p;
  strcpy(str, cstr);
  p = std::strtok(str, delimiter);
  while (p != NULL) {
    arr.push_back(static_cast<A>(std::stof(p)));
    p = std::strtok(NULL, delimiter);
  }
  return arr;
}

template<>
std::vector<std::string> Vectorize(const char* cstr, const char *delimiter);

//! replace a character in a string
void replaceChar(char* buf, char c_old, char c_new);

template<typename T>
void NewCArray(T** &a, int n1, int n2)
{
  a = new T* [n1];
  a[0] = new T [n1*n2];

  for (int i = 0; i < n1; ++i)
    a[i] = a[0] + i*n2;
}

template<typename T>
void FreeCArray(T **a)
{
  delete[] a[0];
  delete[] a;
}

template<typename T>
void FreeCArray2(T **a)
{
  delete[] a[0];
  delete[] a;
}

template<typename T>
void NewCArray(T*** &a, int n1, int n2, int n3)
{
  a = new T** [n1];
  a[0] = new T* [n1*n2];
  a[0][0] = new T [n1*n2*n3];

  for (int i = 0; i < n1; ++i) {
    a[i] = a[0] + i*n2;
    for (int j = 0; j < n2; ++j)
      a[i][j] = a[0][0] + i*n2*n3 + j*n3;
  }
}

template<typename T>
void FreeCArray(T ***a)
{
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template<typename T>
void FreeCArray3(T ***a)
{
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template<typename T>
void NewCArray(T**** &a, int n1, int n2, int n3, int n4)
{
  a = new T*** [n1];
  a[0] = new T** [n1*n2];
  a[0][0] = new T* [n1*n2*n3];
  a[0][0][0] = new T [n1*n2*n3*n4];

  for (int i = 0; i < n1; ++i) {
    a[i] = a[0] + i*n2;
    for (int j = 0; j < n2; ++j) {
      a[i][j] = a[0][0] + i*n2*n3 + j*n3;
      for (int k = 0; k < n3; ++k)
        a[i][j][k] = a[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;
    }
  }
}

template<typename T>
void FreeCArray(T ****a)
{
  delete[] a[0][0][0];
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template<typename T>
void FreeCArray4(T ****a)
{
  delete[] a[0][0][0];
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template <typename T> class AthenaArray;

char* StripLine(char *line);
char* NextLine(char *line, int num, FILE* stream);
void read_data_table(char const *fname, double** data, int *rows, int *cols);
void ReadDataTable(AthenaArray<Real> &data, std::string fname, char c = ' ');

#endif // UTILS_UTILS_HPP_
