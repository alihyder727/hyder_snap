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
   and dependents having the first letter capitalized.
   - Exception #1: conversion function is named with both capitalized names or
     uncaptialized names to preserve symmetry. An example is *deg2rad*. It is allowed to
     use either *2* or *To* in the name of the conversion function.
1. C++ Class name is a concatenated nominal phrase with the all first letters
   capitalized.
1. C++ Functions or Classes cannot have *underscores* in the name.
1. C function name is a lowercase verb phrase connected by *underscores*.
1. The first letter of a variable name must be in lowercase.
1. Variable names can contain *underscores*.
1. Private or protected variable must have an *underscore* at the end of the name.
1. Variable name for a pointer type must start with letter *p*.
1. One dimensional array must use C++ vector.
1. The units of a physical variable are SI. Using units different from SI should be
   specified by append *underscore units* at the end of the variable name.
1. Physical constants are declared as *const Real const* under a Class or under a
   namespace.
