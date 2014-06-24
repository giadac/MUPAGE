/*#############################################################
 *#                                                           #
 *#    Author: G.Carminati                                    #
 *#    First Release: Nov 2007                                #
 *#                                                           #
 *#############################################################
 */

#ifndef GEOMETRY_h
#define GEOMETRY_h
 
struct Vec3D 
{
  double x, y, z;

  Vec3D( double x_val, double y_val, double z_val )
    : x(x_val), y(y_val), z(z_val) {};
};

#endif /* GEOMETRY_h */


