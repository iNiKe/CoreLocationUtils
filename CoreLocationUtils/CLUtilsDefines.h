//
//  CLUtilsDefines.h
//  Indoors
//
//  Created by Rotem Rubnov on 24/10/2011.
//  Copyright (c) 2011 100 grams. All rights reserved.
//

#include <math.h> // For PI

#ifndef Indoors_CLUtilsDefines_h
#define Indoors_CLUtilsDefines_h


// ARC conditions
#if !__has_feature(objc_arc)
#define Release(obj) [obj release]
#define Retain(obj) [obj retain]
#else 
#define Release(obj)
#define Retain(obj) obj
#endif

#define kFarAway                                99999999 //99,999 km
static const double kDegreesToRadians = M_PI / 180.0;
static const double kRadiansToDegrees = 180.0 / M_PI;

#endif
