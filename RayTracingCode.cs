using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using GH_IO;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>



  private void RunScript(Point3d Sun, List<Point3d> Sensors, int SCount, int PCount, int BCount, double Rad, List<Mesh> RoomObjects, List<Vector3d> ObjectNormals, ref object A, ref object B, ref object C, ref object D, ref object E, ref object AllSensorsResults)
  {

    // Make all sensor points from type "Point3d" --> type "Sensor"
    // //////////////////////////////////////////////////////////////////
    List<Sensor> S = new List<Sensor>();
    S = PointToSensor(Sensors, ObjectNormals, SCount);

    // Initialize list of doubles to store % of photons with light for each sensor
    // //////////////////////////////////////////////////////////////////
    List<double> luxArray = new List<double>();

    // INITIALIZE LIST FOR VISUALIZING AND DEBUGGING
    List<Line> L = new List<Line>();

    // START LOOP HERE TO CHECK LIGHT LEVEL FOR EACH SENSOR //
    for (int i = 0; i < S.Count; i++){

      // Initialize output lists
      // //////////////////////////////////////////////////////////////////
      List<List<Photon>> photonNestedList = new List<List<Photon>>();

      // Initialize empty photon array for Sensor[i] and populate with 180 photons
      // //////////////////////////////////////////////////////////////////
      List<Photon> photonsArray = new List<Photon>();
      photonsArray = SensorToPhotons(S[i], Rad);                                      // replace this line with S[i]; currently only looking at first sensor

      // for list of photons, call function to sort photons into two lists:
      // ////////  [0] those that see light, [1] those that never see light
      // //////////////////////////////////////////////////////////////////
      photonNestedList = CheckForIntersectDebug(photonsArray, RoomObjects, ObjectNormals, Rad, BCount);

      // populate luxArray with % of photons that see light.   ;
      double lightCount = photonNestedList[0].Count;
      double darkCount = photonNestedList[1].Count;
      double lightPercentage = lightCount / (lightCount + darkCount);
      luxArray.Add(lightPercentage);

      // Count number of photons for sensor S[i] that see light [0] and that never reach light [1] - THIS CODE IS FOR DEBUGGING - CAN REMOVE
      List<Photon> lightPhotons = new List<Photon>();            // initialize empty list of light photons
      lightPhotons = photonNestedList[0];                        // light photons
      List<Photon> darkPhotons = new List<Photon>();             // initialize empty list of dark photons
      darkPhotons = photonNestedList[1];                         // dark photons
      for (int j = 0; j < lightPhotons.Count();j++){
        L.Add(lightPhotons[j].Direction);
      }

    }
    
    AllSensorsResults = luxArray;

    /*
    
    // ALL REMAINING CODE IN MAIN LOOP IS FOR DEBUGGING PURPOSES ONLY
    
    
    // output readable values - list of lists of sensor results divided by face
    //int totalSensors = luxArray.Count;
    List<List<double>> luxArrayParsed = new List<List<double>>();
    for(int f = 0; f < 5; f++){                             // cycle through each face of hanging bar
      List<double> faceSensorResults = new List<double>();
      for(int s = 0; s < SCount; s++){     // add values for all sensors on face
        faceSensorResults.Add(luxArray[f * SCount + s]);
      }
      luxArrayParsed.Add(faceSensorResults);
    }

    A = luxArrayParsed[0];
    B = luxArrayParsed[1];
    C = luxArrayParsed[2];
    D = luxArrayParsed[3];
    E = luxArrayParsed[4];
    
    */
    
  }

  // <Custom additional code> 

  // ///////////////////////////////////////////////////////////////////////////////////////
  // Function to take in a list of photons (would be one list per sensor), and run loop to
  // check intersection with all meshes per photon and sort photons into buckets of either
  // photons that reached light, or photons that never saw light within limit of possible
  // bounces. Returns two lists: [0] for photons that saw light and [1] other photons
  // ///////////////////////////////////////////////////////////////////////////////////////
  List<List<Photon>> CheckForIntersectDebug(List<Photon> currentPhotonArray, List<Mesh> roomObjects, List<Vector3d> ObjectNormals, double Rad, int bCount){

    // initialize individual lists of photon "buckets" to sort. Also initialize master list of lists to hold both lists
    List<Photon> photonReachLightArray = new List<Photon>();
    List<Photon> photonUpdateBounceArray = new List<Photon>();
    List<Photon> photonTemporaryBounceArray = new List<Photon>();
    List<Photon> photonTemporary2BounceArray = new List<Photon>();
    List<List<Photon>> masterPhotonArray = new List<List<Photon>>();

    // initialize temporary array (2) to store photons for nested loop
    for(int k = 0; k < currentPhotonArray.Count; k++){
      photonTemporary2BounceArray.Add(currentPhotonArray[k]);
    }

    for(int b = 0; b < bCount; b++){                                          // check all photons that can still bounce for intersections until max bounces are reached
      for (int i = 0; i < photonTemporary2BounceArray.Count(); i++){          // for each photon, check if it intersects with any object in room

        // initialize point and point array for intersection
        Point3d[] intersectionPointArray = new Point3d[0];
        Point3d intersection = new Point3d();

        // initialize check if photon[i] intersects any mesh
        bool intersectTrue = false;

        for (int m = 0; m < roomObjects.Count(); m++){                        // check all meshes [m]
          for(int r = 0; r < roomObjects[m].Faces.Count; r++){                // check all faces [r] for mesh [m]

            // initialize mesh face index and check mesh faces for intersection
            int[] faceIndex = new int[r];
            intersectionPointArray = Rhino.Geometry.Intersect.Intersection.MeshLine(roomObjects[m], photonTemporary2BounceArray[i].Direction, out faceIndex);

            // if yes intersection, reflect photon and send to "bounce" bucket
            if (intersectionPointArray.Count() > 0){
              intersectTrue = true;
              intersection = intersectionPointArray[0];                                   // get intersection point
              Vector3d v = new Vector3d(photonTemporary2BounceArray[i].Direction.Direction);       // create vector to represent current direction vector of photon
              Vector3d reflectedVector = new Vector3d();                                  // initialize direction vector of reflected photon
              reflectedVector = Reflect(ObjectNormals[m], v);                             // call function to get direction vector of reflected photon

              double x = intersection.X + (0.0001 * reflectedVector.X);
              double y = intersection.Y + (0.0001 * reflectedVector.Y);
              double z = intersection.Z + (0.0001 * reflectedVector.Z);
              Point3d offsetIntersection = new Point3d(x, y, z);


              Line reflectedLine = new Line(offsetIntersection, reflectedVector, Rad);          // initialize and populate direction line of reflected photon
              Photon reflectedPhoton = new Photon(offsetIntersection, reflectedLine, photonTemporary2BounceArray[i].Bounce + 1);   // initialize and populate new photon; update bounce count
              photonTemporaryBounceArray.Add(reflectedPhoton);                               // add reflected photon to "bounced" bucket

              goto nextPhoton;
            }
          } // next face
        } // next mesh

        // if no intersection, send to "light" bucket
        if (intersectTrue == false){
          photonReachLightArray.Add(photonTemporary2BounceArray[i]);
        }

        nextPhoton:;

      } // next photon

      // reset temporary lists for next bounce itteration
      photonTemporary2BounceArray.Clear();
      for(int j = 0; j < photonTemporaryBounceArray.Count; j++){
        photonTemporary2BounceArray.Add(photonTemporaryBounceArray[j]);
      }
      photonTemporaryBounceArray.Clear();

    } // next bounce

    // pass final list of "bounce" photons from photonTemporary2BounceArray to photonUpdateBounceArray
    for(int j = 0; j < photonTemporary2BounceArray.Count; j++){
      photonUpdateBounceArray.Add(photonTemporary2BounceArray[j]);
    }

    // populate master list with each "bucket" array
    masterPhotonArray.Add(photonReachLightArray);
    masterPhotonArray.Add(photonUpdateBounceArray);
    return masterPhotonArray;
  }

  // ///////////////////////////////////////////////////////////////////////////////////////
  // Function to take in points and normals and make types sensors
  // ///////////////////////////////////////////////////////////////////////////////////////
  List<Sensor> PointToSensor(List < Point3d > points, List<Vector3d> normals, int sCount) {
    List<Sensor> templist = new List<Sensor>();

    for(int c = 0; c < (points.Count / sCount); c++){                    // uses sensors.Count/sCount to do one loop per face
      for(int i = 0; i < sCount; i++){                           // cycle through all sensors on that face (all haves have sCount number of sensors)
        double x = points[(c * sCount) + i].X + (0.05 * normals[c].X);
        double y = points[(c * sCount) + i].Y + (0.05 * normals[c].Y);
        double z = points[(c * sCount) + i].Z + (0.05 * normals[c].Z);
        Point3d offsetStart = new Point3d(x, y, z);
        Sensor s = new Sensor(offsetStart, normals[c]);
        templist.Add(s);
      }
    }

    return templist;

  }

  // ///////////////////////////////////////////////////////////////////////////////////////
  // Function to take in type sensor and create 180 photons in a hemisphere as list
  // ///////////////////////////////////////////////////////////////////////////////////////
  List<Photon> SensorToPhotons(Sensor spt, double L){      // c = number photons per count or pCount
    //List<Photon> SensorToPhotons(Sensor spt, double L){
    // p = sensor point
    // n = normal vector
    // c is the number of photons (hard setting 180)
    // L is the length hardset to large value so it will always intersect

    // Make an array of objects Photons called photonsArray (will return this at the end of the function
    List<Photon> photonsArray = new List<Photon>();

    // Initalize the values to generate photons and sensors
    Point3d initPoint = new Point3d();
    initPoint = spt.Position;
    Vector3d initNorm = new Vector3d();
    initNorm = spt.Normal;
    Matrix rotateMat = new Matrix(3, 3);

    // Write a check for the normal vector to define the rotation matrix used
    if (initNorm.X == 1.0) {
      rotateMat[0, 0] = 0;
      rotateMat[0, 1] = 0;
      rotateMat[0, 2] = 1;
      rotateMat[1, 0] = 0;
      rotateMat[1, 1] = 1;
      rotateMat[1, 2] = 0;
      rotateMat[2, 0] = -1;
      rotateMat[2, 1] = 0;
      rotateMat[2, 2] = 0;
    }
    else if (initNorm.X == -1.0) {
      rotateMat[0, 0] = 0;
      rotateMat[0, 1] = 0;
      rotateMat[0, 2] = -1;
      rotateMat[1, 0] = 0;
      rotateMat[1, 1] = 1;
      rotateMat[1, 2] = 0;
      rotateMat[2, 0] = 1;
      rotateMat[2, 1] = 0;
      rotateMat[2, 2] = 0;
    }

    else if (initNorm.Y == 1.0) {
      rotateMat[0, 0] = 1;
      rotateMat[0, 1] = 0;
      rotateMat[0, 2] = 0;
      rotateMat[1, 0] = 0;
      rotateMat[1, 1] = 0;
      rotateMat[1, 2] = 1;
      rotateMat[2, 0] = 0;
      rotateMat[2, 1] = -1;
      rotateMat[2, 2] = 0;
    }

    else if (initNorm.Y == -1.0) {
      rotateMat[0, 0] = 1;
      rotateMat[0, 1] = 0;
      rotateMat[0, 2] = 0;
      rotateMat[1, 0] = 0;
      rotateMat[1, 1] = 0;
      rotateMat[1, 2] = -1;
      rotateMat[2, 0] = 0;
      rotateMat[2, 1] = 1;
      rotateMat[2, 2] = 0;
    }

    else if (initNorm.Z == 1.0) {
      rotateMat[0, 0] = 1;
      rotateMat[0, 1] = 0;
      rotateMat[0, 2] = 0;
      rotateMat[1, 0] = 0;
      rotateMat[1, 1] = 1;
      rotateMat[1, 2] = 0;
      rotateMat[2, 0] = 0;
      rotateMat[2, 1] = 0;
      rotateMat[2, 2] = 1;
    }
    else if (initNorm.Z == -1.0) {
      rotateMat[0, 0] = 1;
      rotateMat[0, 1] = 0;
      rotateMat[0, 2] = 0;
      rotateMat[1, 0] = 0;
      rotateMat[1, 1] = -1;
      rotateMat[1, 2] = 0;
      rotateMat[2, 0] = 0;
      rotateMat[2, 1] = 0;
      rotateMat[2, 2] = -1;
    }


    for (int thetaIter = 0; thetaIter < 10; thetaIter++){
      double theta = thetaIter * 0.1571;
      for (int phiIter = 0; phiIter < 18; phiIter++ ){
        double phi = phiIter * 0.3491;
        double xPhotonVec = Math.Cos(phi) * Math.Cos(theta);
        double yPhotonVec = Math.Sin(phi) * Math.Cos(theta);
        double zPhotonVec = Math.Sin(theta);
        Matrix matVec = new Matrix(3, 1);
        Matrix matVecRotate = new Matrix(3, 1);

        matVec[0, 0] = xPhotonVec;
        matVec[1, 0] = yPhotonVec;
        matVec[2, 0] = zPhotonVec;
        matVecRotate = rotateMat * matVec;
        // For each photon set the original bounce count to zero
        int bounceCountInit = 0;

        Vector3d photonVec = new Vector3d(matVecRotate[0, 0], matVecRotate[1, 0], matVecRotate[2, 0]);
        //Vector3d photonVec = new Vector3d(xPhotonVec, yPhotonVec, zPhotonVec);

        Line photonLine = new Line(initPoint, photonVec, L);
        // Need to make a line not a vector

        Photon populatePhoton = new Photon(initPoint, photonLine, bounceCountInit);
        photonsArray.Add(populatePhoton);

      }
    }

    return photonsArray;

  }

  // //////////////////////////////////////////////////////////////////////////////
  // Reflection Function - n = surface normal, v = photon's current direction vector
  // returns: new photon direction vector /////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////////
  Vector3d Reflect(Vector3d n, Vector3d v) {
    Vector3d v2 = new Vector3d();
    v2 = v - (2 * n * (v.X * n.X + v.Y * n.Y + v.Z * n.Z));
    return v2;
  }

  // //////////////////////////////////////////////////////////////////////////////
  // Class Definition /////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////////

  public struct Photon
  {
    public Point3d Position;
    public Line Direction;
    public int Bounce;
    public Photon(Point3d pos, Line l, int b)
    {
      Position = pos;
      Direction = l;
      Bounce = b;
    }
  }

  public struct Sensor
  {
    public Point3d Position;
    public Vector3d Normal;
    public Sensor(Point3d pos, Vector3d n)
    {
      Position = pos;
      Normal = n;
    }
  }



  // </Custom additional code> 
}
