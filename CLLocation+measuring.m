//
//  CLLocation+measuring.m
//  StreetWise
//
//  Created by Rotem Rubnov on 4/5/2011.
//  Copyright 2011 100 grams. All rights reserved.
//
//	Permission is hereby granted, free of charge, to any person obtaining a copy
//	of this software and associated documentation files (the "Software"), to deal
//	in the Software without restriction, including without limitation the rights
//	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	copies of the Software, and to permit persons to whom the Software is
//	furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//	THE SOFTWARE.
//

#import "CLLocation+measuring.h"
#include <math.h> // For PI




@implementation CLLocation(measuring)


- (CLLocationDistance) distanceFromCoordinate:(CLLocationCoordinate2D) fromCoord;
{
	double earthRadius = 6371.01; // Earth's radius in Kilometers
	
	// Get the difference between our two points then convert the difference into radians
	double nDLat = (fromCoord.latitude - self.coordinate.latitude) * kDegreesToRadians;  
	double nDLon = (fromCoord.longitude - self.coordinate.longitude) * kDegreesToRadians; 
	
	double fromLat =  self.coordinate.latitude * kDegreesToRadians;
	double toLat =  fromCoord.latitude * kDegreesToRadians;
	
	double nA =	pow ( sin(nDLat/2), 2 ) + cos(fromLat) * cos(toLat) * pow ( sin(nDLon/2), 2 );
	
	double nC = 2 * atan2( sqrt(nA), sqrt( 1 - nA ));
	double nD = earthRadius * nC;
	
	return nD * 1000; // Return our calculated distance in meters
}


+ (CLCoordinateRect) boundingBoxWithCenter: (CLLocationCoordinate2D)centerCoordinate radius:(CLLocationDistance)radius;
{
	CLCoordinateRect result; 
	double earthRadius = 6371.01;
	
	// angular distance in radians on a great circle
	double radDist = radius / earthRadius;
	double radLat = centerCoordinate.latitude * kDegreesToRadians; 
	double radLon = centerCoordinate.longitude * kDegreesToRadians; 
	
	double minLat = radLat - radDist;
	double maxLat = radLat + radDist;
	
	double minLon, maxLon;
	if (minLat > -M_PI/2 && maxLat < M_PI/2) {
		double deltaLon = asin(sin(radDist) / cos(radLat));
		minLon = radLon - deltaLon;
		if (minLon < -M_PI) minLon += 2 * M_PI;
		maxLon = radLon + deltaLon;
		if (maxLon > M_PI) maxLon -= 2 * M_PI;
	} else {
		// a pole is within the distance
		minLat = fmax(minLat, -M_PI/2);
		maxLat = fmin(maxLat, M_PI/2);
		minLon = -M_PI;
		maxLon = M_PI;
	}
	
	result.topLeft.latitude = minLat * kRadiansToDegrees;
	result.topLeft.longitude = minLon  * kRadiansToDegrees;
	result.bottomRight.latitude = maxLat  * kRadiansToDegrees;
	result.bottomRight.longitude = maxLon  * kRadiansToDegrees;
    
	return result; 
}

//returns a direction (in radians) between an origin (from) coordinate and a destination (to) coordinate 
- (CLLocationDirection)directionToLocation:(CLLocation*)location;
{	
	//double tc1; //angle between current location and target location
	//double dlat = location.coordinate.latitude - currLocation.coordinate.latitude;
	//convert to radians
	double fromLong = self.coordinate.longitude * kDegreesToRadians;
	double toLong = location.coordinate.longitude * kDegreesToRadians;
	double fromLat = self.coordinate.latitude * kDegreesToRadians;
	double toLat = location.coordinate.latitude * kDegreesToRadians;
	
	double dlon = toLong - fromLong; 
	double y = sin(dlon)*cos(toLat);
	double x = cos(fromLat)*sin(toLat)-sin(fromLat)*cos(toLat)*cos(dlon);
	
	return atan2(y,x);
	
}


+ (BOOL) isCoordinateValidNonZero : (CLLocationCoordinate2D) coord
{
    return (CLLocationCoordinate2DIsValid(coord) && coord.latitude >= 0.00000001 && coord.longitude > 0.00000001);
}



//========================================================================================
//
// DistancePointLine 


+ (float) magnitudeFromCoordinate: (CLLocationCoordinate2D)coord1 toCoordinate:(CLLocationCoordinate2D)coord2
{
    CLLocationCoordinate2D vector;
    
    vector.latitude = coord2.latitude - coord1.latitude;
    vector.longitude = coord2.longitude - coord1.longitude;
    
    return (float)sqrt( pow(vector.latitude,2) + pow(vector.longitude, 2));
}


- (CLLocationDistance) distanceFromLineWithStartLocation:(CLLocation*) start 
                                             endLocation:(CLLocation*) end 
                                            intersection:(CLLocation**) intersection;
{
        
    double lineMag;
    float U;
        
    lineMag = [CLLocation magnitudeFromCoordinate:end.coordinate toCoordinate:start.coordinate];
    
    U = ( ( ( self.coordinate.latitude - start.coordinate.latitude ) * ( end.coordinate.latitude - start.coordinate.latitude ) ) +
         ( ( self.coordinate.longitude - start.coordinate.longitude ) * ( end.coordinate.longitude - start.coordinate.longitude ) ) ) / pow(lineMag, 2);
    
    if( U < 0.0f || U > 1.0f )
        return -1;   // closest point does not fall within the line segment
    
    CLLocationCoordinate2D intCoord;
    intCoord.latitude = start.coordinate.latitude + U * ( end.coordinate.latitude - start.coordinate.latitude );
    intCoord.longitude = start.coordinate.longitude + U * ( end.coordinate.longitude - start.coordinate.longitude );
    *intersection = [[[CLLocation alloc] initWithLatitude:intCoord.latitude longitude:intCoord.longitude] autorelease];
        
    double distance = [self distanceFromLocation:*intersection];
    
    return distance;    
}



@end
