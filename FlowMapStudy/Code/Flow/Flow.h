#include <ostream>
#include <iostream>
#include <cmath>

#include <Point.h>

using std::ostream;
using std::cerr;
using std::endl;
using std::fabs;

#ifndef FLOW_H
#define FLOW_H

extern const int PointsToFaces[8][3];
extern const int FacesToPoints[6][4];
extern const int EdgesToFaces[12][2];
extern const int FacesToEdges[6][4];
extern const int FaceToFace[6];
extern const int BoundsToFace[6];

/**************************************************************************************************
 *	Flow
 *		Holds the in point, out point, for a flow
 *		the cellID a flow goes through
 *		the FID a flow exits
**************************************************************************************************/
class Flow
{
  public:
	Point in;
	Point out;
	long int cellID;
	int FID;

	Flow(){ cellID = -1; FID = -1; };
	void printFlow( ostream &stream );
	void setPositions( double x, double y, double z );
	void setFID( double* bbox );
};

#endif

/*
 *
Some Comments and diagrams during coding this all up

FID:
	0 - 5:   Face
	6 - 13:  Point
	14 - 25: Edge

Point_of_flow
	x, y, z, t, FID

Flow	
	in point
	out point
	cellID //Cell that flow exists inside of

Advection

Particle Location Check 1

	If particle is on corner
		goTo flow corner function - (use the flow defined at that corner and just follow it)


	If particle is on edge
		goTo Determine cellID with Euler function
	
	If particle is on face
		if i am just seeding this particle
			goTo Determine cellID with Euler
		else 
			using face to face map set cellID
			get checkFace for flows value
			if true
				goTo interpolate flow funciton
			else
				goTo eulerfunction


Determine cellID with Euler
	do a single step of euler integration
	which cell are we in - this is the cell we will check for flows
	goTo check faces for flows function

Check faces for flows
	loop over all faces adjacent to point
		get checkFace for flow value
		if true  - face contains 4 valid flows
			goTo interpolate flow function
	if no faces true
		goTo eulerfunction

Flow corner function
	move particle and add time from flow
	goTo Particle Location Check 1

Interpolate flow function
	using the four flows defined on my face interpolate output location
	check if output location is on a face
	if so
		goTo Particle Location Check 1
	else
		goTo Euler function


Euler funciton
	euler integrate particle
	check if in cell bounds
	if so 
		goto euler
	else
		step back to boundry
		goto particle location check 1


[ Point / Edge / Face ] Lables



         AXIS:
        z+
        |
		|  y+
		| /
		|/____ x+


POINT IDS 8 [0-7]:

               6 _________________________ 7
                / _____________________  /|
               / / ___________________/ / |
              / / /| |               / /  |
             / / / | |              / / . |
            / / /| | |             / / /| |
           / / / | | |            / / / | |
          / / /  | | |           / / /| | |
         / /_/___|_|_|__________/ / / | | |
     4  /________________________/ /5 | | |
        | ______________________ | |  | | |
        | | |    | | |_________| | |__| | |
        | | |  2 | |___________| | |____| | 3
        | | |   / / ___________| | |_  / /
        | | |  / / /           | | |/ / /
        | | | / / /            | | | / /
        | | |/ / /             | | |/ /
        | | | / /              | | ' /
        | | |/_/_______________| |  /
        | |____________________| | /
        |________________________|/
       0                          1    



Edge IDS 12 [0-11]:
                               
                 _________________________ 
                / ___________6_________  /|
               / / ___________________/ / |
              / / /| |               / /  |
             / / / | |              / / . |
            / 7 /| | |             / / /| |
           / / / | | |            / 5 / | |
          / / /  |11 |           / / /| | |
         / /_/___|_|_|__________/ / / |10 |
        /______________4_________/ /  | | |
        | ______________________ | |  | | |
        | | |    | | |_________| | |__| | |
        | | |    | |________2__| | |____| | 
        | | |   / / ___________| | |_  / /
        | 8 |  / / /           | 9 |/ / /
        | | | / 3 /            | | | 1 /
        | | |/ / /             | | |/ / 
        | | | / /              | | ' /
        | | |/_/_______________| |  /
        | |__________0_________| | /
        |________________________|/
                                                         





Face IDS 6 [0-5]:

                 _________________________ 
                / _____________________  /|
               / /                    / / |
              / /                    / /  |
             / /                    / / . |
            / /         5          / / /| |
           / /                    / / / | |
          / /                    / / /  | |
         / /____________________/ / /   | |
        /________________________/ /    | |
        | ______________________ | |  1 | |
        | |                    | | |    | |
        | |                    | | |    | | 
        | |                    | | |   / /
        | |                    | | |  / /
        | |          0         | | | / /
        | |                    | | |/ /
        | |                    | | ' /
        | |                    | |  /
        | |____________________| | /
        |________________________|/
                                 
                 _________________________
                / _______________________/|
               / / ___________________  | |
              / / /| |                | | |
             / / / | |                | | |
            / / /| | |                | | |
           / / / | | |                | | |
          / / /  | | |       2        | | |
         / / /   | | |                | | |
        /_/_/    | | |                | | |
        | | | 3  | | |                | | |
        | | |    | | |________________| | |
        | | |    | |__________________| | |
        | | |   / / _________________/ / /
        | | |  / / /                / / /
        | | | / / /                / / /
        | | |/ / /       4        / / /
        | | | / /                / / /
        | | |/_/________________/ / /
        | |_____________________\/ /
        |________________________|/
        
        
*/
