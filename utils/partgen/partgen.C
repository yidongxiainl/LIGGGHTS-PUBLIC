#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <omp.h>

using namespace std;

struct atoms_t {
  int    *  type;
  int    ** flag;
  bool   *  wall;
  double ** x;
};

int main(int argc, char **argv)
{
  // A header message
  cout << endl;
  cout << "#################################################\n";
  cout << endl;
  cout << endl;
  cout << " A Handy Code to Generate Complex Particle Shapes\n";
  cout << " Using a Multi-Sphere Approach                   \n";
  cout << endl;
  cout << " Yidong Xia \n";
  cout << endl;
  cout << "#################################################\n";
  cout << endl;

  int itype; // type of multi-sphere body

  int nxmax = 1; // max number of particles in x-direction
  int nymax = 1; // max number of particles in y-direction
  int nzmax = 1; // max number of particles in z-direction
  unsigned int nspheres = 0; // number of actual spheres

  double xlo; // lower  bound in x-axis
  double xhi; // higher bound in x-axis
  double ylo; // lower  bound in y-axis
  double yhi; // higher bound in y-axis
  double zlo; // lower  bound in z-axis
  double zhi; // higher bound in z-axis

  double dx; // unit particle spacing in x-axis
  double dy; // unit particle spacing in y-axis
  double dz; // unit particle spacing in z-axis

  double radius; // radius of base particles
  double rx; // x-radius of an ellipsoid
  double ry; // y-radius of an ellipsoid
  double rz; // z-radius of an ellipsoid

  vector<double> xcoor; // x-coordinates of spheres
  vector<double> ycoor; // y-coordinates of spheres
  vector<double> zcoor; // z-coordinates of spheres


  ifstream finpu; // input data file
  ifstream fctrl; // input control file
  ofstream fdata; // output data file
  ofstream fdump; // output LAMMPS dump file

  // strings
  string skipl;

  fctrl.open("partgen.ctrl", ios::in);
  if (!fctrl.is_open())
  {
    cout << "Error ... the file 'partgen.ctrl' does not exist! " << endl;
    exit(0);
  }

  getline(fctrl,skipl); fctrl >> itype;
  getline(fctrl,skipl);
  getline(fctrl,skipl); fctrl >> xlo >> xhi >> dx;
  getline(fctrl,skipl);
  getline(fctrl,skipl); fctrl >> ylo >> yhi >> dy;
  getline(fctrl,skipl);
  getline(fctrl,skipl); fctrl >> zlo >> zhi >> dz;
  getline(fctrl,skipl);
  getline(fctrl,skipl); fctrl >> radius >> rx >> ry >> rz;
  getline(fctrl,skipl);


  if (xlo > xhi)
  {
    cout << "Error: xlo > xhi" << endl;
    cout << "xlo = " << xlo << endl;
    cout << "xhi = " << xhi << endl;
    exit(0);
  }
  if (ylo > yhi)
  {
    cout << "Error: ylo > yhi" << endl;
    cout << "ylo = " << ylo << endl;
    cout << "yhi = " << yhi << endl;
    exit(0);
  }
  if (zlo > zhi)
  {
    cout << "Error: zlo > zhi" << endl;
    cout << "zlo = " << zlo << endl;
    cout << "zhi = " << zhi << endl;
    exit(0);
  }
  if (rx < 0.0)
  {
    cout << "Error: rx < 0.0" << endl;
    cout << "rx = " << rx << endl;
    exit(0);
  }
  if (ry < 0.0)
  {
    cout << "Error: ry < 0.0" << endl;
    cout << "ry = " << ry << endl;
    exit(0);
  }
  if (rz < 0.0)
  {
    cout << "Error: rz < 0.0" << endl;
    cout << "rz = " << rz << endl;
    exit(0);
  }

  nxmax += (int)((xhi-xlo)/dx);
  nymax += (int)((yhi-ylo)/dy);
  nzmax += (int)((zhi-zlo)/dz);

  cout << "nxmax = " << nxmax << endl;
  cout << "nymax = " << nymax << endl;
  cout << "nzmax = " << nzmax << endl;

  fdata.open("biopart.dat", ios::out);
  fdata.precision(8);
  fdata << scientific;

  switch (itype)
  {
    // line
    case 1:

      for (unsigned int j=1; j<=nymax; j++)
      {
        double xp = 0.0;
        double yp = ylo + (float)(j-1) * dy;
        double zp = 0.0;

        if ((yp+dy) > yhi || (yp-dy) < ylo) continue;

        fdata<<right<<setw(20)<<xp<<setw(20)<<yp<<setw(20)<<zp<<setw(20)<<radius<<endl;

        xcoor.push_back(xp);
        ycoor.push_back(yp);
        zcoor.push_back(zp);

        nspheres++;
      }
      break;

    // plate
    case 2:

      for (unsigned int j = 1; j <= nymax; j++)
        for (unsigned int i = 1; i <= nxmax; i++)
        {
          double xp = xlo + (float)(i-1) * dx;
          double yp = ylo + (float)(j-1) * dy;
          double zp = 0.0;

          if ((xp+dx) > xhi || (xp-dx) < xlo) continue;
          if ((yp+dy) > yhi || (yp-dy) < ylo) continue;

          fdata<<right<<setw(20)<<xp<<setw(20)<<yp<<setw(20)<<zp<<setw(20)<<radius<< endl;

          xcoor.push_back(xp);
          ycoor.push_back(yp);
          zcoor.push_back(zp);

          nspheres++;
        }
      break;

    // cylinder
    case 3:

      for (unsigned int k = 1; k <= nzmax; k++)
        for (unsigned int j = 1; j <= nymax; j++)
          for (unsigned int i = 1; i <= nxmax; i++)
          {
            double xp = xlo + (float)(i-1) * dx;
            double yp = ylo + (float)(j-1) * dy;
            double zp = zlo + (float)(k-1) * dz;

            if ((xp*xp+yp*yp) > (radius*radius)) continue;
            if ((zp+dz) > zhi || (zp-dz) < zlo) continue;

            fdata<<right<<setw(20)<<xp<<setw(20)<<yp<<setw(20)<<zp<<setw(20)<<radius<< endl;

            xcoor.push_back(xp);
            ycoor.push_back(yp);
            zcoor.push_back(zp);

            nspheres++;
          }
      break;

    // ellipsoid
    case 4:

      for (unsigned int k = 1; k <= nzmax; k++)
        for (unsigned int j = 1; j <= nymax; j++)
          for (unsigned int i = 1; i <= nxmax; i++)
          {
            double xp = xlo + (float)(i-1) * dx;
            double yp = ylo + (float)(j-1) * dy;
            double zp = zlo + (float)(k-1) * dz;

            if ((xp*xp/rx/rx+yp*yp/ry/ry+zp*zp/rz/rz) > 1.0) continue;

            fdata<<right<<setw(20)<<xp<<setw(20)<<yp<<setw(20)<<zp<<setw(20)<<radius<< endl;

            xcoor.push_back(xp);
            ycoor.push_back(yp);
            zcoor.push_back(zp);

            nspheres++;
          }

      break;

    default:
      cout << "Error: itype = " << itype << " is not defined! \n" << endl;
  }

  fctrl.close();
  fdata.close();


  fdump.open("biopart.dump", ios::out);
  fdump.precision(8);
  fdump << scientific;

  fdump << "ITEM: TIMESTEP\n";
  fdump << "0\n";
  fdump << "ITEM: NUMBER OF ATOMS\n";
  fdump << nspheres << endl;
  fdump << "ITEM: ATOMS id type x y z radius\n";
  for (unsigned int id = 1; id <= nspheres; id++)
    fdump<<right<<setw(5)<<id<<setw(3)<<itype<<setw(20)<<xcoor[id-1]<<setw(20)<<ycoor[id-1]<<setw(20)<<zcoor[id-1]<<setw(20)<<radius<< endl;

  xcoor.clear();
  ycoor.clear();
  zcoor.clear();

  fdump.close();

  //cout << "Done processing the LAMMPS restart file...\n" << endl;
  exit(0);
}
