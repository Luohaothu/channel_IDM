#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#include "Geometry.h"
#include "Field.h"

#define PI 3.1415926535898
#define INFTSM 1e-8


#define WallRscl(y, r) (y < 1 ? fmin(r*y, 1.) : fmax(2-r*(2-y), 1.))