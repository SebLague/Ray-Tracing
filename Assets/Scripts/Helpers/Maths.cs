using UnityEngine;
using System;
using System.Collections.Generic;
using static System.Math;
using static System.MathF;

namespace Seb.Helpers
{
	// ---- Version 0.6 [31/Aug/2025] ----
	public static class Maths
	{
		/*
		 * OVERVIEW:
		 * #Constants
		 * #Basics
		 * #Easing
		 * #Ray3D
		 * #Ray2D
		 * #ShapesAndLines2D
		 * #ClosestPoint
		 * #RNGUtils
		 * #StatsAndCounting
		 * #SphericalGeometry
		 * #Sampling
		 * #Layout
		 * #Miscellaneous
		 * #Structures
		 */

		#region #Constants

		public const float PI = 3.1415926f;
		public const float TAU = 2 * PI;
		public const float Epsilon = 1.175494351E-38f;

		#endregion

		#region #Basics

		public static float Min(float a, float b) => a < b ? a : b;
		public static float Max(float a, float b) => a > b ? a : b;
		public static float Mod(float value, float max) => (value % max + max) % max;
		public static int Mod(int value, int max) => (value % max + max) % max;

		#endregion

		#region #Easing

		public enum EaseType
		{
			Linear,
			QuadIn,
			QuadOut,
			QuadInOut,
			CubeIn,
			CubeOut,
			CubeInOut,
			ElasticOut
		}

		public static float EaseLinear(float t) => Clamp01(t);
		public static float EaseQuadIn(float t) => Square(Clamp01(t));
		public static float EaseQuadOut(float t) => 1 - Square(1 - Clamp01(t));
		public static float EaseQuadInOut(float t) => 3 * Square(Clamp01(t)) - 2 * Cube(Clamp01(t));

		public static float EaseCubeIn(float t) => Cube(Clamp01(t));
		public static float EaseCubeOut(float t) => 1 - Cube(1 - Clamp01(t));

		public static float EaseCubeInOut(float t)
		{
			t = Clamp01(t);
			int r = (int)System.Math.Round(t);
			return 4 * Cube(t) * (1 - r) + (1 - 4 * Cube(1 - t)) * r;
		}

		public static float EaseElasticOut(float t)
		{
			t = Clamp01(t);
			if (t is 0 or 1)
			{
				return t;
			}

			float a = System.MathF.Pow(2, -10 * t);
			float b = System.MathF.Sin((t - .075f) * (2 * System.MathF.PI) / 0.3f);
			return a * b + 1;
		}

		public static float Ease(EaseType type, float t)
		{
			return type switch
			{
				EaseType.Linear => EaseLinear(t),
				EaseType.QuadIn => EaseQuadIn(t),
				EaseType.QuadOut => EaseQuadOut(t),
				EaseType.QuadInOut => EaseQuadInOut(t),
				EaseType.CubeIn => EaseCubeIn(t),
				EaseType.CubeOut => EaseCubeOut(t),
				EaseType.CubeInOut => EaseCubeInOut(t),
				EaseType.ElasticOut => EaseElasticOut(t),
				_ => t
			};
		}

		public static float Clamp01(float t) => Mathf.Clamp(t, 0, 1);
		public static float Square(float x) => x * x;
		public static float Cube(float x) => x * x * x;
		public static float Quart(float x) => x * x * x * x;
		public static float Abs(float x) => System.Math.Abs(x);

		#endregion

		#region #Ray3D

		public static RaySphereResult RaySphere(Vector3 rayOrigin, Vector3 rayDir, Vector3 centre, float radius)
		{
			Vector3 offset = rayOrigin - centre;
			const float a = 1;
			float b = 2 * Vector3.Dot(offset, rayDir);
			float c = Vector3.Dot(offset, offset) - radius * radius;
			float d = b * b - 4 * c; // Discriminant from quadratic formula

			// Number of intersections: 0 when d < 0; 1 when d = 0; 2 when d > 0
			if (d > 0)
			{
				float s = Mathf.Sqrt(d);
				float dstToSphereNear = Mathf.Max(0, -b - s) / (2 * a);
				float dstToSphereFar = (-b + s) / (2 * a);

				// Ignore intersections that occur behind the ray
				if (dstToSphereFar >= 0)
				{
					return new RaySphereResult()
					{
						intersects = true,
						dstToSphere = dstToSphereNear,
						dstThroughSphere = dstToSphereFar - dstToSphereNear
					};
				}
			}

			// Ray did not intersect sphere
			return new RaySphereResult()
			{
				intersects = false,
				dstToSphere = Mathf.Infinity,
				dstThroughSphere = 0
			};
		}

		// Thanks to https://tavianator.com/2011/ray_box.html
		public static (bool hit, bool isInside, float dst) RayBoundingBox(Vector3 rayOrigin, Vector3 rayDir, Vector3 boxMin, Vector3 boxMax)
		{
			float invDirX = rayDir.x == 0 ? float.PositiveInfinity : 1 / rayDir.x;
			float invDirY = rayDir.y == 0 ? float.PositiveInfinity : 1 / rayDir.y;
			float invDirZ = rayDir.z == 0 ? float.PositiveInfinity : 1 / rayDir.z;

			float tx1 = (boxMin.x - rayOrigin.x) * invDirX;
			float tx2 = (boxMax.x - rayOrigin.x) * invDirX;
			float tmin = Mathf.Min(tx1, tx2);
			float tmax = Mathf.Max(tx1, tx2);

			float ty1 = (boxMin.y - rayOrigin.y) * invDirY;
			float ty2 = (boxMax.y - rayOrigin.y) * invDirY;
			tmin = Mathf.Max(tmin, Mathf.Min(ty1, ty2));
			tmax = Mathf.Min(tmax, Mathf.Max(ty1, ty2));

			float tz1 = (boxMin.z - rayOrigin.z) * invDirZ;
			float tz2 = (boxMax.z - rayOrigin.z) * invDirZ;
			tmin = Mathf.Max(tmin, Mathf.Min(tz1, tz2));
			tmax = Mathf.Min(tmax, Mathf.Max(tz1, tz2));

			bool hit = tmax >= tmin && tmax > 0;
			bool isInside = tmin <= 0;
			float dst = isInside ? tmax : tmin;
			return (hit, isInside, dst);
		}

		// Calculate the intersection of a ray with a triangle using Moller-Trumbore algorithm
		// Thanks to https://stackoverflow.com/a/42752998
		public static (bool hit, float dst) RayTriangle(Vector3 rayOrigin, Vector3 rayDir, Vector3 triA, Vector3 triB, Vector3 triC)
		{
			Vector3 edgeAB = triB - triA;
			Vector3 edgeAC = triC - triA;
			Vector3 ao = rayOrigin - triA;
			Vector3 dao = Vector3.Cross(ao, rayDir);
			Vector3 normalVector = Vector3.Cross(edgeAB, edgeAC);

			float determinant = -Vector3.Dot(rayDir, normalVector);
			float invDet = 1 / determinant;

			// Calculate dst to triangle & barycentric coordinates of intersection point
			float dst = Vector3.Dot(ao, normalVector) * invDet;
			float u = Vector3.Dot(edgeAC, dao) * invDet;
			float v = -Vector3.Dot(edgeAB, dao) * invDet;
			float w = 1 - u - v;

			// Initialize hit info
			bool hit = determinant >= 1E-6 && dst >= 0 && u >= 0 && v >= 0 && w >= 0;
			return (hit, dst);
		}

		#endregion

		#region #Ray2D

		public static (bool hit, Vector2 near, Vector2 far) RayCircle(Vector2 origin, Vector2 dir, Vector2 centre, float radius)
		{
			// Thanks to https://www.bluebill.net/circle_ray_intersection.html
			Vector2 u = centre - origin;
			float offsetDotRay = Vector2.Dot(u, dir);
			Vector2 u1 = offsetDotRay * dir;
			Vector2 u2 = u - u1;
			float desc = radius * radius - u2.sqrMagnitude;

			bool inside = u.sqrMagnitude <= radius * radius;
			bool behind = offsetDotRay < 0 && !inside;
			if (desc < 0 || behind) return (false, Vector2.zero, Vector2.zero);

			float m = Mathf.Sqrt(desc);
			Vector2 near = origin + u1 + m * dir * (inside ? 1 : -1);
			Vector2 far = origin + u1 + m * dir * (inside ? -1 : 1);

			return (true, near, far);
		}


		/// <summary> Test if ray intersects infinite line, and get point of intersection </summary>
		public static (bool intersects, Vector2 point) RayIntersectsLine(Vector2 rayOrigin, Vector2 rayDir, Vector2 lineA, Vector2 lineB)
		{
			Vector2 lineOffset = lineA - lineB;
			float d = Determinant(-rayDir, lineOffset);
			// Check if parallel
			if (ApproximatelyEqual(d, 0))
			{
				return (false, Vector2.zero);
			}

			float n = Determinant(rayOrigin - lineA, lineOffset);
			float t = n / d;
			Vector2 intersectionPoint = rayOrigin + rayDir * t;
			bool intersectsInFrontOfRay = t >= 0;
			return (intersectsInFrontOfRay, intersectionPoint);
		}

		/// <summary> Test if ray intersects line segment, and get point of intersection </summary>
		public static (bool intersects, Vector2 point) RayIntersectsLineSegment(Vector2 rayOrigin, Vector2 rayDir, Vector2 lineA, Vector2 lineB)
		{
			Vector2 ab = lineA - lineB;
			Vector2 abPerp = new Vector2(-ab.y, ab.x);
			float rayDotABPerp = Vector2.Dot(rayDir, abPerp);

			if (rayDotABPerp == 0) return (false, rayOrigin);

			float dst = Vector2.Dot(lineA - rayOrigin, abPerp) / rayDotABPerp;
			Vector2 intersectionPoint = rayOrigin + rayDir * dst;
			bool isOnSegment = Vector2.Dot(lineA - intersectionPoint, lineB - intersectionPoint) <= 0;

			return (dst >= 0 && isOnSegment, intersectionPoint);
		}

		#endregion

		#region #ShapesAndLines2D

		/// <summary> Returns the area of a triangle in 3D space. </summary>
		public static float TriangleArea(Vector3 a, Vector3 b, Vector3 c)
		{
			// Thanks to https://math.stackexchange.com/a/1951650
			Vector3 ortho = Vector3.Cross(c - a, b - a);
			float parallogramArea = ortho.magnitude;
			return parallogramArea * 0.5f;
		}

		public static bool TriangleIsClockwiseTest(Vector2 a, Vector2 b, Vector2 c)
		{
			return Vector2.Dot(b - a, Perpendicular(c - b)) > 0;
		}

		public static bool PolygonIsClockwiseTest(Vector2[] points)
		{
			float extremeX = float.MinValue;
			float extremeY = float.MinValue;
			int extremeIndex = -1;

			for (int i = 0; i < points.Length; i++)
			{
				Vector2 p = points[i];
				if (p.x > extremeX || (p.x == extremeX && p.y > extremeY))
				{
					extremeX = p.x;
					extremeY = p.y;
					extremeIndex = i;
				}
			}

			Vector2 a = points[extremeIndex];
			Vector2 b = points[(extremeIndex + 1) % points.Length];
			Vector2 c = points[(extremeIndex - 1 + points.Length) % points.Length];

			return TriangleIsClockwiseTest(a, b, c);
		}

		public static Vector2 Perpendicular(Vector2 v) => new Vector2(-v.y, v.x);


		/// <summary>
		/// Returns the signed area of a polygon (negative if clockwise)
		/// Note: does not matter whether endpoints are duplicate
		/// </summary>
		public static float PolygonAreaSigned(Vector2[] points)
		{
			// Ignore last point if it is a duplicate of the first
			int numPoints = (points[^1] - points[0]).sqrMagnitude < 0.00001f ? points.Length - 1 : points.Length;

			float area = 0;
			for (int i = 0; i < numPoints; i++)
			{
				Vector2 a = points[i];
				Vector2 b = points[(i + 1) % points.Length];
				// We can calculate the polygon area by summing signed areas of triangles between each edge and a fixed point (here
				// chosen as the origin to simplify calculations). This is the same idea as rendering glyphs with even-odd rule, since
				// the sign of the triangle area cancels out overlapping parts when coming around the other side of the polygon.
				// Note also that the factor of 0.5 in triangle area is deferred to the end. 
				area += (a.x + b.x) * (b.y - a.y);
			}

			return area * 0.5f;
		}


		/// <summary>
		/// Returns the centroid (centre of mass) of a polygon
		/// Note: order of points does not matter, nor whether endpoints are duplicate
		/// </summary>
		public static Vector2 PolygonCentreOfMass(Vector2[] points)
		{
			// Ignore last point if it is a duplicate of the first
			int numPoints = (points[^1] - points[0]).sqrMagnitude < 0.00001f ? points.Length - 1 : points.Length;

			float xSum = 0;
			float ySum = 0;
			float area = 0;

			for (int i = 0; i < numPoints; i++)
			{
				Vector2 a = points[i];
				Vector2 b = points[(i + 1) % points.Length];
				area += (a.x + b.x) * (b.y - a.y) / 2;

				float x = (a.x + b.x) * (a.x * b.y - b.x * a.y);
				float y = (a.y + b.y) * (a.x * b.y - b.x * a.y);
				xSum += x;
				ySum += y;
			}

			return new Vector2(xSum, ySum) / (6 * area);
		}

		/// <summary>
		/// Returns the signed area of a triangle in 2D space.
		/// The sign depends on whether the points are given in clockwise (negative) or counter-clockwise (positive) order.
		/// </summary>
		public static float TriangleAreaSigned2D(Vector2 a, Vector2 b, Vector2 c)
		{
			return 0.5f * ((a.x - b.x) * (a.y - c.y) + (a.y - b.y) * (c.x - a.x));
			/*  // ---- Unoptimized version ----
			 *  // Consider AC as the base, and calculate line perpendicular to it (of same length).
			 *	// Then, take adjacent edge such as AB and project it onto that to get height * base
			 *	Vector2 baseEdge = c - a;
			 *	Vector2 basePerp = new Vector2(-baseEdge.y, baseEdge.x);
			 *	return 0.5f * Vector2.Dot(basePerp, a - b);
			 */
		}

		public static bool PointInCircle2D(Vector2 point, Vector2 circleCentre, float circleRadius)
		{
			return (point - circleCentre).sqrMagnitude <= circleRadius * circleRadius;
		}

		/// <summary> Test if point p is inside the triangle (a, b, c) </summary>
		public static bool TriangleContainsPoint(Vector2 p, Vector2 a, Vector2 b, Vector2 c)
		{
			// Thanks to https://stackoverflow.com/a/14382692
			float area = TriangleAreaSigned2D(a, b, c);
			float s = (a.y * c.x - a.x * c.y + (c.y - a.y) * p.x + (a.x - c.x) * p.y) * Mathf.Sign(area);
			float t = (a.x * b.y - a.y * b.x + (a.y - b.y) * p.x + (b.x - a.x) * p.y) * Mathf.Sign(area);
			return s >= 0 && t >= 0 && s + t < 2 * Mathf.Abs(area);
		}

		/// <summary> Determines whether the given 2D triangle is wound in a clockwise order</summary>
		public static bool TriangleIsClockwise(Vector2 a, Vector2 b, Vector2 c)
		{
			return TriangleAreaSigned2D(a, b, c) < 0;
		}

		/// <summary>
		/// Test if a 2D polygon contains the given point.
		/// Points can be ordered clockwise or counterclockwise.
		/// Note: function doesn't care if last point in polygon is duplicate of first point or not
		/// </summary>
		public static bool PolygonContainsPoint(Vector2 p, Vector2[] points)
		{
			// Thanks to Dan Sunday
			int windingNumber = 0;

			// Ignore last point if it is a duplicate of the first
			int numPoints = (points[^1] - points[0]).sqrMagnitude < 0.00001f ? points.Length - 1 : points.Length;

			for (int i = 0; i < numPoints; i++)
			{
				Vector2 a = points[i];
				Vector2 b = points[(i + 1) % points.Length];

				if (a.y <= p.y)
				{
					if (b.y > p.y && PointIsOnLeftSideOfLine(p, a, b))
					{
						windingNumber++;
					}
				}
				else if (b.y <= p.y && !PointIsOnLeftSideOfLine(p, a, b))
				{
					windingNumber--;
				}
			}

			return windingNumber != 0;

			// Calculate which side of line AB point P is on
			bool PointIsOnLeftSideOfLine(Vector2 p, Vector2 a, Vector2 b)
			{
				return (b.x - a.x) * (p.y - a.y) - (p.x - a.x) * (b.y - a.y) > 0;
			}
		}

		public static bool QuadContainsPoint(Vector2 point, Vector2 boxCentre, Vector2 boxSize)
		{
			float ox = Mathf.Abs(point.x - boxCentre.x);
			float oy = Mathf.Abs(point.y - boxCentre.y);
			return ox < boxSize.x / 2 && oy < boxSize.y / 2;
		}

		public static bool QuadsOverlap(Vector2 centreA, Vector2 sizeA, Vector2 centreB, Vector2 sizeB)
		{
			float leftA = centreA.x - sizeA.x / 2;
			float rightA = centreA.x + sizeA.x / 2;
			float topA = centreA.y + sizeA.y / 2;
			float bottomA = centreA.y - sizeA.y / 2;

			float leftB = centreB.x - sizeB.x / 2;
			float rightB = centreB.x + sizeB.x / 2;
			float topB = centreB.y + sizeB.y / 2;
			float bottomB = centreB.y - sizeB.y / 2;

			return leftA <= rightB && rightA >= leftB && topA >= bottomB && bottomA <= topB;
		}


		/// <summary> Get point on the line segment (a1, a2) that's closest to the given point (p) </summary>
		public static Vector3 ClosestPointOnLineSegment(Vector3 p, Vector3 a1, Vector3 a2)
		{
			Vector3 lineDelta = a2 - a1;
			Vector3 pointDelta = p - a1;
			float sqrLineLength = lineDelta.sqrMagnitude;

			if (sqrLineLength == 0)
				return a1;

			float t = Mathf.Clamp01(Vector3.Dot(pointDelta, lineDelta) / sqrLineLength);
			return a1 + lineDelta * t;
		}

		/// <summary> Calculates smallest distance from given point to the line segment (a1, a2)</summary>
		public static float DistanceToLineSegment(Vector3 p, Vector3 a1, Vector3 a2)
		{
			Vector3 closestPoint = ClosestPointOnLineSegment(p, a1, a2);
			return (p - closestPoint).magnitude;
		}

		/// <summary> Calculates smallest distance from given point to the line segment (a1, a2)</summary>
		public static float SqrDistanceToLineSegment(Vector3 p, Vector3 a1, Vector3 a2)
		{
			return (p - ClosestPointOnLineSegment(p, a1, a2)).sqrMagnitude;
		}

		/// <summary> Test if two infinite 2D lines intersect (true unless parallel), and get point of intersection </summary>
		public static (bool intersects, Vector2 point) LineIntersectsLine(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2)
		{
			float d = (a1.x - a2.x) * (b1.y - b2.y) - (a1.y - a2.y) * (b1.x - b2.x);
			// Check if parallel
			if (ApproximatelyEqual(d, 0))
			{
				return (false, Vector2.zero);
			}

			float n = (a1.x - b1.x) * (b1.y - b2.y) - (a1.y - b1.y) * (b1.x - b2.x);
			float t = n / d;
			Vector2 intersectionPoint = a1 + (a2 - a1) * t;
			return (true, intersectionPoint);
		}

		static float Determinant(Vector2 a, Vector2 b)
		{
			return a.x * b.y - a.y * b.x;
		}

		public static (bool intersects, Vector2 point) LineSegIntersectsLineSegTest(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2)
		{
			float d = (a1.x - a2.x) * (b1.y - b2.y) - (a1.y - a2.y) * (b1.x - b2.x);
			// Check if parallel
			if (d == 0) return (false, Vector2.zero);

			float n = (a1.x - b1.x) * (b1.y - b2.y) - (a1.y - b1.y) * (b1.x - b2.x);
			float t = n / d;
			Vector2 intersectionPoint = a1 + (a2 - a1) * t;

			bool onSegA = Vector2.Dot(a1 - intersectionPoint, a2 - intersectionPoint) <= 0;
			bool onSegB = Vector2.Dot(b1 - intersectionPoint, b2 - intersectionPoint) <= 0;
			return (onSegA && onSegB, intersectionPoint);
		}

		public static bool LineSegmentsIntersect(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2)
		{
			float d = (b2.x - b1.x) * (a1.y - a2.y) - (a1.x - a2.x) * (b2.y - b1.y);
			if (d == 0)
				return false;
			float t = ((b1.y - b2.y) * (a1.x - b1.x) + (b2.x - b1.x) * (a1.y - b1.y)) / d;
			float u = ((a1.y - a2.y) * (a1.x - b1.x) + (a2.x - a1.x) * (a1.y - b1.y)) / d;

			return t >= 0 && t <= 1 && u >= 0 && u <= 1;
		}


		/// <summary> Returns -1 or +1 depending on which side point p is of the line (a1, a2). Returns 0 if on line. </summary>
		public static int SideOfLine(Vector2 p, Vector2 a, Vector2 b)
		{
			float det = Determinant(b - a, p - a);
			return System.Math.Sign(det);
		}

		/// <summary> Test if points p1 and p2 are on the same side of the line (a1, a2) </summary>
		public static bool PointOnSameSideOfLine(Vector2 p1, Vector2 p2, Vector2 a1, Vector2 a2)
		{
			return SideOfLine(p1, a1, a2) == SideOfLine(p2, a1, a2);
		}

		#endregion

		#region #ClosestPoint

		/// <summary>
		/// Given an infinite plane defined by some point that the plane passes through, as well as the normal vector of the plane,
		/// this function returns the nearest point on that plane to the given point p.
		/// </summary>
		public static Vector3 ClosestPointOnPlane(Vector3 anyPointOnPlane, Vector3 planeNormal, Vector3 p)
		{
			float signedDstToPlane = Vector3.Dot(anyPointOnPlane - p, planeNormal);
			Vector3 closestPointOnPlane = p + planeNormal * signedDstToPlane;
			return closestPointOnPlane;
		}

		public static Vector3 ClosestPointOnBox(Vector3 boxSize, Vector3 p)
		{
			Vector3 result = p;
			Vector3 halfSize = boxSize / 2;

			// Handle outside box
			if (Mathf.Abs(p.x) > halfSize.x) result.x = halfSize.x * Mathf.Sign(p.x);
			if (Mathf.Abs(p.y) > halfSize.y) result.y = halfSize.y * Mathf.Sign(p.y);
			if (Mathf.Abs(p.z) > halfSize.z) result.z = halfSize.z * Mathf.Sign(p.z);

			// Handle inside box
			Vector3 o = halfSize - Maths.Abs(p);
			bool nearestX = o.x < o.y && o.x < o.z;
			bool nearestY = !nearestX && o.y < o.z;
			if (nearestX) result.x = halfSize.x * Mathf.Sign(p.x);
			else if (nearestY) result.y = halfSize.y * Mathf.Sign(p.y);
			else result.z = halfSize.z * Mathf.Sign(p.z);

			return result;
		}

		// Calculate normal vector of nearest face to p on box at origin
		public static Vector3 ClosestBoxFaceNormal(Vector3 boxSize, Vector3 p)
		{
			Vector3 o = (boxSize * 0.5f - Maths.Abs(p));
			bool nearestX = o.x < o.y && o.x < o.z;
			bool nearestY = !nearestX && o.y < o.z;

			if (nearestX) return Vector3.right * Mathf.Sign(p.x);
			if (nearestY) return Vector3.up * Mathf.Sign(p.y);
			return Vector3.forward * Mathf.Sign(p.z);
		}

		#endregion

		#region #RNGUtils

		/// <summary> Random point inside of circle (uniform distribution) </summary>
		public static Vector2 RandomPointInCircle(System.Random rng)
		{
			Vector2 pointOnCircle = RandomPointOnCircle(rng);
			float r = Mathf.Sqrt((float)rng.NextDouble());
			return pointOnCircle * r;
		}

		/// <summary> Random point on circumference of circle </summary>
		public static Vector2 RandomPointOnCircle(System.Random rng)
		{
			float angle = (float)rng.NextDouble() * 2 * PI;
			float x = Mathf.Cos(angle);
			float y = Mathf.Sin(angle);
			return new Vector2(x, y);
		}

		/// <summary> Random point on surface of sphere (i.e. random direction) </summary>
		public static Vector3 RandomPointOnSphere(System.Random rng)
		{
			float x = RandomNormal(rng, 0, 1);
			float y = RandomNormal(rng, 0, 1);
			float z = RandomNormal(rng, 0, 1);
			return new Vector3(x, y, z).normalized;
		}

		/// <summary> Random point inside a triangle (with uniform distribution). </summary>
		public static Vector3 RandomPointInTriangle(Vector3 a, Vector3 b, Vector3 c, System.Random rng)
		{
			double randA = rng.NextDouble();
			double randB = rng.NextDouble();
			if (randA + randB > 1)
			{
				randA = 1 - randA;
				randB = 1 - randB;
			}

			return a + (b - a) * (float)randA + (c - a) * (float)randB;
		}

		/// <summary>
		/// Pick random index, weighted by the weights array.
		/// For example, if the array contains {1, 6, 3}...
		/// The possible indices would be (0, 1, 2)
		/// and the probabilities for these would be (1/10, 6/10, 3/10)
		/// </summary>
		public static int WeightedRandomIndex(System.Random rng, float[] weights)
		{
			float weightSum = 0;
			for (int i = 0; i < weights.Length; i++)
			{
				weightSum += weights[i];
			}

			float randomValue = (float)rng.NextDouble() * weightSum;
			float cumul = 0;

			for (int i = 0; i < weights.Length; i++)
			{
				cumul += weights[i];
				if (randomValue < cumul)
				{
					return i;
				}
			}

			return weights.Length - 1;
		}

		/// <summary> Randomly shuffles the elements of the given array </summary>
		public static void ShuffleArray<T>(T[] array, System.Random rng)
		{
			// wikipedia.org/wiki/Fisher–Yates_shuffle#The_modern_algorithm
			for (int i = 0; i < array.Length - 1; i++)
			{
				int randomIndex = rng.Next(i, array.Length);
				(array[randomIndex], array[i]) = (array[i], array[randomIndex]); // Swap
			}
		}

		/// <summary> Randomly shuffles the elements of the given list </summary>
		public static void ShuffleList<T>(IList<T> list, System.Random rng)
		{
			// wikipedia.org/wiki/Fisher–Yates_shuffle#The_modern_algorithm
			for (int i = 0; i < list.Count - 1; i++)
			{
				int randomIndex = rng.Next(i, list.Count);
				(list[randomIndex], list[i]) = (list[i], list[randomIndex]); // Swap
			}
		}

		/// <summary> Populates array with indices from 0 up to array.Length-1, but with the order shuffled randomly </summary>
		public static void PopulateWithUniqueRandomIndices(int[] array, System.Random rng)
		{
			for (int i = 0; i < array.Length; i++)
			{
				array[i] = i;
			}

			ShuffleArray(array, rng);
		}

		/// <summary> Create an array with indices from 0 up to Count-1, but with the order shuffled randomly </summary>
		public static int[] CreateUniqueRandomIndices(int count, System.Random rng)
		{
			int[] array = new int[count];
			PopulateWithUniqueRandomIndices(array, rng);
			return array;
		}

		#endregion

		#region #StatsAndCounting

		/// <summary>
		/// Returns a random value with normal distribution.
		/// The mean determines the 'centre' of the distribution. 
		/// The standardDeviation controls the spread of the distribution (i.e. how likely it is to get values that are far from the mean).
		/// See https://www.desmos.com/calculator/0dnzmd0x0h for example.
		/// </summary>
		public static float RandomNormal(System.Random rng, float mean = 0, float standardDeviation = 1)
		{
			// Thanks to https://stackoverflow.com/a/6178290
			float theta = 2 * Mathf.PI * (float)rng.NextDouble();
			float rho = Mathf.Sqrt(-2 * Mathf.Log((float)rng.NextDouble()));
			float scale = standardDeviation * rho;
			return mean + scale * Mathf.Cos(theta);
		}

		/// <summary>
		/// Get min max float value from function.
		/// Example usage: Maths.GetMinMax((i) => particles[i].velocity.magnitude, particles.Count);
		/// </summary>
		public static (float min, float max, float mean) GetMinMaxMean(Func<int, float> GetValue, int count)
		{
			float min = float.MaxValue;
			float max = float.MinValue;
			float sum = 0;

			for (int i = 0; i < count; i++)
			{
				float val = GetValue(i);
				min = Math.Min(min, val);
				max = Math.Max(max, val);
				sum += val;
			}

			float mean = sum / count;
			return (min, max, mean);
		}

		public static (float min, float max, float mean) GetMinMaxMean(ReadOnlySpan<float> values)
		{
			float min = float.MaxValue;
			float max = float.MinValue;
			float sum = 0;

			foreach (var val in values)
			{
				min = Math.Min(min, val);
				max = Math.Max(max, val);
				sum += val;
			}

			float mean = sum / values.Length;
			return (min, max, mean);
		}

		#endregion

		#region #VectorUtils

		// Thanks to https://math.stackexchange.com/a/4112622
		// Calculates arbitrary normalized vector that is perpendicular to the given direction
		public static Vector3 CalculateOrthonormal(Vector3 dir)
		{
			float a = Mathf.Sign((Mathf.Sign(dir.x) + 0.5f) * (Mathf.Sign(dir.z) + 0.5f));
			float b = Mathf.Sign((Mathf.Sign(dir.y) + 0.5f) * (Mathf.Sign(dir.z) + 0.5f));
			return new Vector3(a * dir.z, b * dir.z, -a * dir.x - b * dir.y).normalized;
		}

		public static Vector3 Refract(Vector3 inDir, Vector3 normal, float iorA, float iorB)
		{
			// Thanks to https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
			float refractRatio = iorA / iorB;
			float cosAngleIn = -Vector3.Dot(inDir, normal);
			float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
			if (sinSqrAngleOfRefraction > 1) return Vector3.zero; // Ray is fully reflected, no refraction occurs

			Vector3 refractDir = refractRatio * inDir + (refractRatio * cosAngleIn - Mathf.Sqrt(1 - sinSqrAngleOfRefraction)) * normal;
			return refractDir;
		}

		public static Vector3 Reflect(Vector3 inDir, Vector3 normal)
		{
			return inDir - 2 * Vector3.Dot(inDir, normal) * normal;
		}

		public static float Fresnel(Vector3 inDir, Vector3 normal, float iorA, float iorB)
		{
			// Thanks to https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
			float refractRatio = iorA / iorB;
			float cosAngleIn = -Vector3.Dot(inDir, normal);
			float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
			if (sinSqrAngleOfRefraction > 1) return 1; // Ray is fully reflected, no refraction occurs

			float cosAngleOfRefraction = Mathf.Sqrt(1 - sinSqrAngleOfRefraction);
			float rPerp = (iorA * cosAngleIn - iorB * cosAngleOfRefraction) / (iorA * cosAngleIn + iorB * cosAngleOfRefraction);
			float rPar = (iorB * cosAngleIn - iorA * cosAngleOfRefraction) / (iorB * cosAngleIn + iorA * cosAngleOfRefraction);

			return (rPerp * rPerp + rPar * rPar) / 2;
		}

		public static Vector2 Rotate2D(Vector2 p, float angle)
		{
			float sinAngle = MathF.Sin(angle);
			float cosAngle = MathF.Cos(angle);
			Vector2 iHat = new Vector2(cosAngle, sinAngle);
			Vector2 jHat = new Vector2(-sinAngle, cosAngle);
			return iHat * p.x + jHat * p.y;
		}

		public static Vector2 RotateAroundPoint2D(Vector2 p, Vector2 anchor, float angle)
		{
			return Rotate2D(p - anchor, angle) + anchor;
		}

		public static Vector2 Abs(Vector2 vec) => new(MathF.Abs(vec.x), MathF.Abs(vec.y));
		public static Vector3 Abs(Vector3 vec) => new(MathF.Abs(vec.x), MathF.Abs(vec.y), MathF.Abs(vec.z));
		public static Vector4 Abs(Vector4 vec) => new(MathF.Abs(vec.x), MathF.Abs(vec.y), MathF.Abs(vec.z), MathF.Abs(vec.w));

		#endregion

		#region #SphericalGeometry

		/// <summary>
		/// Returns the length of the shortest arc between two points on the surface of a unit sphere.
		/// </summary>
		public static float ArcLengthBetweenPointsOnUnitSphere(Vector3 a, Vector3 b)
		{
			// Thanks to https://www.movable-type.co.uk/scripts/latlong-vectors.html
			return Mathf.Atan2(Vector3.Cross(a, b).magnitude, Vector3.Dot(a, b));
			// Note: The following (simpler) approach works too, but is less precise for small angles:
			// return Mathf.Acos(Vector3.Dot(a, b));
		}

		/// <summary>
		/// Returns the length of the shortest arc between two points on the surface of a sphere with the specified radius.
		/// </summary>
		public static float ArcLengthBetweenPointsOnSphere(Vector3 a, Vector3 b, float sphereRadius)
		{
			return ArcLengthBetweenPointsOnUnitSphere(a.normalized, b.normalized) * sphereRadius;
		}

		#endregion

		#region #Sampling

		public static float SampleBilinear(float[,] map, float u, float v, bool wrap = false)
		{
			int width = map.GetLength(0);
			int height = map.GetLength(1);
			u -= 0.5f / width;
			v -= 0.5f / height;

			if (wrap)
			{
				u = Mod(u, 1);
				v = Mod(v, 1);
			}

			float px = u * (width);
			float py = v * (height);

			int x = (int)px;
			int y = (int)py;

			float xFrac = Mathf.Clamp01(px - x);
			float yFrac = Mathf.Clamp01(py - y);

			float bottomLeft = Sample(x, y);
			float bottomRight = Sample(x + 1, y);
			float topLeft = Sample(x, y + 1);
			float topRight = Sample(x + 1, y + 1);

			float interpolatedTop = Mathf.Lerp(topLeft, topRight, xFrac);
			float interpolatedBottom = Mathf.Lerp(bottomLeft, bottomRight, xFrac);
			return Mathf.Lerp(interpolatedBottom, interpolatedTop, yFrac);

			float Sample(int x, int y) => map[HandleBoundary(x, width), HandleBoundary(y, height)];

			int HandleBoundary(int value, int dim)
			{
				return wrap ? Mod(value, dim) : Clamp(value, 0, dim - 1);
			}
		}

		#endregion

		#region #Layout

		public static (Vector2 size, Vector2 centre) HorizontalLayout(int numElements, int elementIndex, Vector2 centre, Vector2 size, float spacing)
		{
			float spaceTotal = (numElements - 1) * spacing;
			float elementWidth = (size.x - spaceTotal) / numElements;
			float posX = centre.x - size.x / 2 + elementWidth / 2 + (spacing + elementWidth) * elementIndex;
			return (new Vector2(elementWidth, size.y), new Vector2(posX, centre.y));
		}

		#endregion

		#region #Miscellaneous

		public static float AbsoluteMax(float a, float b)
		{
			return Mathf.Abs(a) > Mathf.Abs(b) ? a : b;
		}

		public static Vector2Int CalculateSpawnCountPerAxisBox2D(Vector2 size, float spawnDensity)
		{
			float area = size.x * size.y;
			int targetTotal = Mathf.CeilToInt(area * spawnDensity);

			float lenSum = size.x + size.y;
			Vector2 t = size / lenSum;
			float m = Mathf.Sqrt(targetTotal / (t.x * t.y));
			int nx = Mathf.CeilToInt(t.x * m);
			int ny = Mathf.CeilToInt(t.y * m);

			return new Vector2Int(nx, ny);
		}

		public static Vector3Int CalculateSpawnCountPerAxisBox3D(Vector3 size, int spawnDensity)
		{
			float volume = size.x * size.y * size.z;
			int targetParticleCount = Mathf.CeilToInt(volume * spawnDensity);

			float lenSum = size.x + size.y + size.z;
			Vector3 t = size / lenSum;

			float m = MathF.Cbrt(targetParticleCount / (t.x * t.y * t.z));

			int nx = Mathf.CeilToInt(t.x * m);
			int ny = Mathf.CeilToInt(t.y * m);
			int nz = Mathf.CeilToInt(t.z * m);

			return new Vector3Int(nx, ny, nz);
		}

		public static (Vector2 centre, Vector2 size) BoundingBox(Vector2[] points)
		{
			if (points.Length == 0)
			{
				return (Vector2.zero, Vector2.zero);
			}

			Vector2 min = points[0];
			Vector2 max = points[0];
			for (int i = 1; i < points.Length; i++)
			{
				Vector2 p = points[i];
				min = new Vector2(Min(min.x, p.x), Min(min.y, p.y));
				max = new Vector2(Max(max.x, p.x), Max(max.y, p.y));
			}

			Vector2 centre = (min + max) / 2;
			Vector2 size = max - min;
			return (centre, size);
		}

		// Calculate point on sphere given longitude and latitude (in radians), and the radius of the sphere
		public static Vector3 CoordinateToSpherePoint(float latitude, float longitude, float radius = 1)
		{
			float y = Mathf.Sin(latitude);
			float r = Mathf.Cos(latitude); // radius of 2d circle cut through sphere at 'y'
			float x = Mathf.Sin(longitude) * r;
			float z = -Mathf.Cos(longitude) * r;

			return new Vector3(x, y, z) * radius;
		}

		public static (float longitude, float latitude) PointToCoordinate(Vector3 pointOnUnitSphere)
		{
			float latitude = Mathf.Asin(pointOnUnitSphere.y);
			float a = pointOnUnitSphere.x;
			float b = -pointOnUnitSphere.z;

			float longitude = Mathf.Atan2(a, b);
			return (longitude, latitude);
		}


		/// <summary>
		/// Rotates the point around the axis by the given angle (in radians)
		/// </summary>
		public static Vector3 RotateAroundAxis(Vector3 point, Vector3 axis, float angle)
		{
			return RotateAroundAxis(point, axis, Mathf.Sin(angle), Mathf.Cos(angle));
		}

		/// <summary>
		/// Rotates given vector by the rotation that aligns startDir with endDir
		/// </summary>
		public static Vector3 RotateBetweenDirections(Vector3 vector, Vector3 startDir, Vector3 endDir)
		{
			Vector3 rotationAxis = Vector3.Cross(startDir, endDir);
			float sinAngle = rotationAxis.magnitude;
			float cosAngle = Vector3.Dot(startDir, endDir);

			return RotateAroundAxis(vector, rotationAxis.normalized, sinAngle, cosAngle);
			// Note: this achieves the same as doing: 
			// return Quaternion.FromToRotation(originalDir, newDir) * point;
		}

		static Vector3 RotateAroundAxis(Vector3 point, Vector3 axis, float sinAngle, float cosAngle)
		{
			// Rodrigues' rotation formula: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
			return point * cosAngle + Vector3.Cross(axis, point) * sinAngle + axis * Vector3.Dot(axis, point) * (1 - cosAngle);
		}

		public static float SampleData(float[] data, double t)
		{
			if (data.Length == 0) return float.NaN;
			if (data.Length == 1) return data[0];

			// Sample
			int indexA = (int)(t * (data.Length - 1));
			int indexB = Math.Min(data.Length - 1, indexA + 1);
			float valA = data[indexA];
			float valB = data[indexB];

			// Interpolate
			double interval = 1.0 / (data.Length - 1);
			double intervalA = indexA * interval;
			double intervalB = indexB * interval;
			if (intervalA >= intervalB) return valB;
			double intervalT = (t - intervalA) / (intervalB - intervalA);
			return (float)(valA + (valB - valA) * intervalT);
		}

		// Calculates the real roots of a cubic equation (i.e. values of x satisfying: ax^3 + bx^2 + cx + d = 0)
		public static float[] RealCubicRoots(float a, float b, float c, float d)
		{
			float p = (3 * a * c - b * b) / (3 * a * a);
			float q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
			float tToX = -b / (3 * a);
			float rootTest = 4 * p * p * p + 27 * q * q;
			bool onlyOneRealRoot = (rootTest > 0 && p < 0) || p > 0;

			if (onlyOneRealRoot && p > 0)
			{
				float asinh = Asinh(3 * q / (2 * p) * Sqrt(3 / p));
				float root0 = -2 * Sqrt(p / 3) * Sinh(1 / 3f * asinh) + tToX;
				return new float[] { root0, float.NaN, float.NaN };
			}
			else if (onlyOneRealRoot)
			{
				float acosh = Acosh(-3 * Math.Abs(q) / (2 * p) * Sqrt(-3 / p));
				float root0 = -2 * Math.Sign(q) * Sqrt(-p / 3) * Cosh(1 / 3f * acosh) + tToX;
				return new float[] { root0, float.NaN, float.NaN };
			}
			else return new float[] { CalcRoot(0), CalcRoot(1), CalcRoot(2) };

			float CalcRoot(int k) // where k is 0, 1, or 2
			{
				float acos = Acos(3 * q / (2 * p) * Sqrt(-3 / p));
				float cos = Cos(1 / 3f * acos - 2 * PI * k / 3);
				return 2 * Sqrt(-p / 3) * cos + tToX;
			}
		}

		public static bool ApproximatelyEqual(float a, float b) => System.Math.Abs(a - b) < Epsilon;


		public static int IntLog2(int value)
		{
			int result = 0;
			while (value > 1)
			{
				value >>= 1;
				result++;
			}

			return result;
		}

		public static int TwosComplement(uint unsignedValue, int numBits)
		{
			if (numBits < 32)
			{
				uint unsignedRange = 1u << numBits;
				uint firstNegativeValue = unsignedRange >> 1;

				if (unsignedValue >= firstNegativeValue)
				{
					return (int)(unsignedValue - unsignedRange);
				}
			}

			return (int)unsignedValue;
		}

		/// <summary>
		/// Returns n points distributed reasonably evenly on a sphere.
		/// Uses fibonacci spiral technique.
		/// </summary>
		public static Vector3[] GetPointsOnSphereSurface(int numPoints, float radius = 1)
		{
			// Thanks to https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/44164075#44164075
			Vector3[] points = new Vector3[numPoints];
			const double goldenRatio = 1.618033988749894; // (1 + sqrt(5)) / 2
			const double angleIncrement = System.Math.PI * 2 * goldenRatio;

			System.Threading.Tasks.Parallel.For(0, numPoints, i =>
			{
				double t = (double)i / numPoints;
				double inclination = System.Math.Acos(1 - 2 * t);
				double azimuth = angleIncrement * i;

				double x = Sin(inclination) * Cos(azimuth);
				double y = Sin(inclination) * Sin(azimuth);
				double z = Cos(inclination);
				points[i] = new Vector3((float)x, (float)y, (float)z) * radius;
			});
			return points;
		}
		
		// Moment of inertia (as well as mass, and centre of mass) of convex polygon
		public static (float moi, float mass, Vector2 com) PolygonMomentOfInertia(Vector2[] polygonPoints, float density)
		{
			// Thanks to https://physics.stackexchange.com/a/709017
			float polygonArea = 0;
			Vector2 triCentreSum = Vector2.zero;
			float moiAboutOrigin = 0;

			for (int i = 2; i < polygonPoints.Length; i++)
			{
				Vector2 a = polygonPoints[i];
				Vector2 b = polygonPoints[0];
				Vector2 c = polygonPoints[i - 1];

				float triArea = Mathf.Abs(Maths.TriangleArea(a, b, c));
				polygonArea += triArea;
				triCentreSum += (a + b + c) / 3 * triArea;

				moiAboutOrigin += triArea * (Dot(a, a) + Dot(b, b) + Dot(c, c) + Dot(a, b) + Dot(b, c) + Dot(c, a)) / 6;
			}

			Vector2 com = triCentreSum / polygonArea;
			float mass = density * Mathf.Abs(polygonArea);

			moiAboutOrigin *= density;
			float moiAboutCOM = moiAboutOrigin - mass * Dot(com, com);
			return (moiAboutCOM, mass, com);

			static float Dot(Vector2 a, Vector2 b) => Vector2.Dot(a, b);
		}

		#endregion

		#region #Structures

		public struct RaySphereResult
		{
			public bool intersects;
			public float dstToSphere;
			public float dstThroughSphere;
		}

		#endregion
	}
}