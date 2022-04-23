# Simulating Non-hydrostatic Atmospheres on Planets (SNAP)
<!-- Jenkins Status Badge in Markdown (with view), unprotected, flat style -->
<!-- In general, need to be on Princeton VPN, logged into Princeton CAS, with ViewStatus access to Jenkins instance to click on unprotected Build Status Badge, but server is configured to whitelist GitHub -->
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
<!--[![Build Status](https://travis-ci.com/luminoctum/athena19-dev.svg?token=AfxC7sH2UkyrrtpsBrob&branch=dev)](https://travis-ci.com/luminoctum/athena19-dev) -->
[![Build Status](https://github.com/luminoctum/athena19-dev/actions/workflows/autotest.yml/badge.svg)](https://github.com/luminoctum/athena19-dev/actions/workflows/autotest.yml)

<!--[![Public GitHub  issues](https://img.shields.io/github/issues/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/issues)
[![Public GitHub pull requests](https://img.shields.io/github/issues-pr/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/pulls) -->

Athena++ Atmospheric Dynamics Code

## Coding convections
1. C++ Function name is a concatenated verb phrase with the verb having all lowercase letters
   and all dependent words capitalizing the first letter.
   - Exception #1: conversion function can named with either capitalized names or
     uncaptialized names to preserve symmetry. An example is *deg2rad*. It is allowed to
     use either *2* or *To* in the name of the conversion function.
   - Exception #2: if a function name is a single verb, the first letter should be
     capitalized. An example is *Vectorize*.
1. C++ Class name is a concatenated nominal phrase capitalizaing the first letter of
   each word.
1. C++ Class name should be in singular form.
1. C++ Functions or Classes should not have *underscores* in the name.
1. C function name is a lowercase verb phrase concatenated by *underscores*.
1. The first letter of a variable name must be in lowercase.
1. Variable names can contain *underscores*.
1. Private or protected variables must contain an *underscore* at the end of the name.
1. Variable name for a pointer type must start with letter *p*.
1. One dimensional array should be a C++ vector container.
1. The units of a physical variable use SI by default. Using units different from SI should be
   specified by append *underscore units* at the end of the variable name.
   - Example #1: pres\_bar means *pressure in bars*.
1. Physical constants are declared as *const Real const* under a Class or under a
   namespace.
1. Use *nullptr* for null pointer.
1. If a function returns more than one elementary variables, pass pointers as arguments.
1. If a function returns complex quantities, pass reference.
1. The returned variables are normally placed in the first few arguments in a function
  - Exception #1: In MPI calls, the send argument is always placed in the first.
1. Use 2 spaces for indentation.
