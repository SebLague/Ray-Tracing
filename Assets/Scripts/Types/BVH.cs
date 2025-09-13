using UnityEngine;
using System;
using System.Collections.Generic;
using System.Text;

namespace Seb.AccelerationStructures
{
	// ---- Version 0.1 [13/Sep/2025] ---- 
	public class BVH
	{
		public enum Quality
		{
			Low, // Split along longest axis
			High, // Surface area heuristic
			Disabled // No splitting (all triangles inside single bounding box)
		}

		public readonly Triangle[] Triangles;
		public readonly Node[] Nodes;
		public BuildStats stats;

		readonly NodeList nodeList;
		readonly BVHTriangle[] buildTriangles;
		readonly Quality quality;

		public BVH(Vector3[] verts, int[] indices, Vector3[] normals, Quality quality = Quality.High)
		{
			// Start recording stats
			var sw = System.Diagnostics.Stopwatch.StartNew();
			stats = new BuildStats(quality);

			// Construct BVH
			this.quality = quality;
			nodeList = new();
			buildTriangles = new BVHTriangle[indices.Length / 3];

			float xMin = float.MaxValue;
			float xMax = float.MinValue;
			float yMin = float.MaxValue;
			float yMax = float.MinValue;
			float zMin = float.MaxValue;
			float zMax = float.MinValue;

			for (int i = 0; i < indices.Length; i += 3)
			{
				Vector3 a = verts[indices[i + 0]];
				Vector3 b = verts[indices[i + 1]];
				Vector3 c = verts[indices[i + 2]];

				BVHTriangle tri = new(a, b, c, i);
				buildTriangles[i / 3] = tri;

				if (tri.MinX < xMin) xMin = tri.MinX;
				if (tri.MinY < yMin) yMin = tri.MinY;
				if (tri.MinZ < zMin) zMin = tri.MinZ;
				if (tri.MaxX > xMax) xMax = tri.MaxX;
				if (tri.MaxY > yMax) yMax = tri.MaxY;
				if (tri.MaxZ > zMax) zMax = tri.MaxZ;
			}

			nodeList.Add(new Node(xMin, yMin, zMin, xMax, yMax, zMax, -1, -1));
			if (quality == Quality.Disabled)
			{
				nodeList.Nodes[0].StartIndex = 0;
				nodeList.Nodes[0].TriangleCount = buildTriangles.Length;
			}
			else Split(0, 0, buildTriangles.Length);

			Triangles = new Triangle[buildTriangles.Length];
			for (int i = 0; i < buildTriangles.Length; i++)
			{
				BVHTriangle buildTri = buildTriangles[i];
				Vector3 a = verts[indices[buildTri.Index + 0]];
				Vector3 b = verts[indices[buildTri.Index + 1]];
				Vector3 c = verts[indices[buildTri.Index + 2]];
				Vector3 norm_a = normals[indices[buildTri.Index + 0]];
				Vector3 norm_b = normals[indices[buildTri.Index + 1]];
				Vector3 norm_c = normals[indices[buildTri.Index + 2]];
				Triangles[i] = new Triangle(a, b, c, norm_a, norm_b, norm_c);
			}

			Nodes = nodeList.Nodes.AsSpan(0, nodeList.NodeCount).ToArray();

			// Finish recording stats
			sw.Stop();
			stats.TimeMs = (int)sw.ElapsedMilliseconds;
		}

		void Split(int parentIndex, int triGlobalStart, int triNum, int depth = 0)
		{
			const int MaxDepth = 32;
			ref Node parent = ref nodeList.Nodes[parentIndex];

			float sizeX = parent.MaxX - parent.MinX;
			float sizeY = parent.MaxY - parent.MinY;
			float sizeZ = parent.MaxZ - parent.MinZ;
			float parentCost = NodeCost(sizeX, sizeY, sizeZ, triNum);

			(int splitAxis, float splitPos, float cost) = ChooseSplit(parent, triGlobalStart, triNum);

			if (cost < parentCost && depth < MaxDepth)
			{
				float xMinLeft = float.MaxValue;
				float xMaxLeft = float.MinValue;
				float yMinLeft = float.MaxValue;
				float yMaxLeft = float.MinValue;
				float zMinLeft = float.MaxValue;
				float zMaxLeft = float.MinValue;

				float xMinRight = float.MaxValue;
				float xMaxRight = float.MinValue;
				float yMinRight = float.MaxValue;
				float yMaxRight = float.MinValue;
				float zMinRight = float.MaxValue;
				float zMaxRight = float.MinValue;
				int numOnLeft = 0;

				for (int i = triGlobalStart; i < triGlobalStart + triNum; i++)
				{
					BVHTriangle tri = buildTriangles[i];

					float c = splitAxis switch
					{
						0 => tri.CentreX,
						1 => tri.CentreY,
						_ => tri.CentreZ
					};

					if (c < splitPos)
					{
						if (tri.MinX < xMinLeft) xMinLeft = tri.MinX;
						if (tri.MinY < yMinLeft) yMinLeft = tri.MinY;
						if (tri.MinZ < zMinLeft) zMinLeft = tri.MinZ;
						if (tri.MaxX > xMaxLeft) xMaxLeft = tri.MaxX;
						if (tri.MaxY > yMaxLeft) yMaxLeft = tri.MaxY;
						if (tri.MaxZ > zMaxLeft) zMaxLeft = tri.MaxZ;

						BVHTriangle swap = buildTriangles[triGlobalStart + numOnLeft];
						buildTriangles[triGlobalStart + numOnLeft] = tri;
						buildTriangles[i] = swap;
						numOnLeft++;
					}
					else
					{
						if (tri.MinX < xMinRight) xMinRight = tri.MinX;
						if (tri.MinY < yMinRight) yMinRight = tri.MinY;
						if (tri.MinZ < zMinRight) zMinRight = tri.MinZ;
						if (tri.MaxX > xMaxRight) xMaxRight = tri.MaxX;
						if (tri.MaxY > yMaxRight) yMaxRight = tri.MaxY;
						if (tri.MaxZ > zMaxRight) zMaxRight = tri.MaxZ;
					}
				}

				int numOnRight = triNum - numOnLeft;
				int triStartLeft = triGlobalStart + 0;
				int triStartRight = triGlobalStart + numOnLeft;

				// Split parent into two children
				Node childLeft = new(xMinLeft, yMinLeft, zMinLeft, xMaxLeft, yMaxLeft, zMaxLeft, triStartLeft, 0);
				Node childRight = new(xMinRight, yMinRight, zMinRight, xMaxRight, yMaxRight, zMaxRight, triStartRight, 0);
				int childIndexLeft = nodeList.Add(childLeft);
				int childIndexRight = nodeList.Add(childRight);

				// Update parent
				parent.StartIndex = childIndexLeft;
				nodeList.Nodes[parentIndex] = parent;
				stats.RecordNode(depth, false);

				// Recursively split children
				Split(childIndexLeft, triGlobalStart, numOnLeft, depth + 1);
				Split(childIndexRight, triGlobalStart + numOnLeft, numOnRight, depth + 1);
			}
			else
			{
				// Parent is actually leaf, assign all triangles to it
				parent.StartIndex = triGlobalStart;
				parent.TriangleCount = triNum;
				nodeList.Nodes[parentIndex] = parent;
				stats.RecordNode(depth, true, triNum);
			}
		}

		(int axis, float pos, float cost) ChooseSplit(Node node, int start, int count)
		{
			if (count <= 1) return (0, 0, float.PositiveInfinity);

			float sizeX = node.MaxX - node.MinX;
			float sizeY = node.MaxY - node.MinY;
			float sizeZ = node.MaxZ - node.MinZ;

			if (quality == Quality.Low)
			{
				int largestAxisIndex = sizeX > sizeY && sizeX > sizeZ ? 0 : sizeY > sizeZ ? 1 : 2;
				float pos = largestAxisIndex switch
				{
					0 => node.MinX + sizeX * 0.5f,
					1 => node.MinY + sizeY * 0.5f,
					_ => node.MinZ + sizeZ * 0.5f
				};

				return (largestAxisIndex, pos, EvaluateSplit(largestAxisIndex, pos, start, count));
			}

			float bestSplitPos = 0;
			int bestSplitAxis = 0;
			int maxSplitTests = count < 10 ? 3 : 5;
			float maxAxis = Mathf.Max(sizeX, sizeY, sizeZ);
			float bestCost = float.MaxValue;

			// Estimate best split pos
			for (int axis = 0; axis < 3; axis++)
			{
				float axisSize;
				float axisMin;

				switch (axis)
				{
					case 0:
						axisSize = sizeX;
						axisMin = node.MinX;
						break;
					case 1:
						axisSize = sizeY;
						axisMin = node.MinY;
						break;
					default:
						axisSize = sizeZ;
						axisMin = node.MinZ;
						break;
				}

				int numSplitTests = Mathf.CeilToInt(axisSize / maxAxis * maxSplitTests);
				numSplitTests = Math.Clamp(numSplitTests, 1, maxSplitTests);

				for (int i = 0; i < numSplitTests; i++)
				{
					float splitT = (i + 1) / (numSplitTests + 1f);
					float splitPos = axisMin + axisSize * splitT;
					float cost = EvaluateSplit(axis, splitPos, start, count);
					if (cost < bestCost)
					{
						bestCost = cost;
						bestSplitPos = splitPos;
						bestSplitAxis = axis;
					}
				}
			}

			return (bestSplitAxis, bestSplitPos, bestCost);
		}


		float EvaluateSplit(int splitAxis, float splitPos, int start, int count)
		{
			int numOnLeft = 0;
			int numOnRight = 0;

			float xMinLeft = float.MaxValue;
			float xMaxLeft = float.MinValue;
			float yMinLeft = float.MaxValue;
			float yMaxLeft = float.MinValue;
			float zMinLeft = float.MaxValue;
			float zMaxLeft = float.MinValue;

			float xMinRight = float.MaxValue;
			float xMaxRight = float.MinValue;
			float yMinRight = float.MaxValue;
			float yMaxRight = float.MinValue;
			float zMinRight = float.MaxValue;
			float zMaxRight = float.MinValue;

			int end = start + count;

			for (int i = start; i < end; i++)
			{
				ref BVHTriangle tri = ref buildTriangles[i];
				float c = splitAxis switch
				{
					0 => tri.CentreX,
					1 => tri.CentreY,
					_ => tri.CentreZ
				};

				if (c < splitPos)
				{
					if (tri.MinX < xMinLeft) xMinLeft = tri.MinX;
					if (tri.MinY < yMinLeft) yMinLeft = tri.MinY;
					if (tri.MinZ < zMinLeft) zMinLeft = tri.MinZ;
					if (tri.MaxX > xMaxLeft) xMaxLeft = tri.MaxX;
					if (tri.MaxY > yMaxLeft) yMaxLeft = tri.MaxY;
					if (tri.MaxZ > zMaxLeft) zMaxLeft = tri.MaxZ;

					numOnLeft++;
				}
				else
				{
					if (tri.MinX < xMinRight) xMinRight = tri.MinX;
					if (tri.MinY < yMinRight) yMinRight = tri.MinY;
					if (tri.MinZ < zMinRight) zMinRight = tri.MinZ;
					if (tri.MaxX > xMaxRight) xMaxRight = tri.MaxX;
					if (tri.MaxY > yMaxRight) yMaxRight = tri.MaxY;
					if (tri.MaxZ > zMaxRight) zMaxRight = tri.MaxZ;
					numOnRight++;
				}
			}


			float costA = NodeCost(xMaxLeft - xMinLeft, yMaxLeft - yMinLeft, zMaxLeft - zMinLeft, numOnLeft);
			float costB = NodeCost(xMaxRight - xMinRight, yMaxRight - yMinRight, zMaxRight - zMinRight, numOnRight);
			return costA + costB;
		}

		static float NodeCost(float x, float y, float z, int numTriangles)
		{
			if (numTriangles == 0) return 0;
			float area = x * y + x * z + y * z;
			return area * numTriangles;
		}
		
		// ---- Traversal ---
		public (bool hit, float dst, Vector3 pos, bool backface, Vector3 normal) Search(Vector3 rayOrigin, Vector3 rayDir)
		{
			Stack<Node> stack = new();
			stack.Push(nodeList.Nodes[0]);

			float minDst = float.MaxValue;
			int hitTriangleIndex = -1;
			Vector3 hitPoint = Vector3.zero;
			bool backface = false;
			Vector3 normal = Vector3.zero;

			while (stack.Count > 0)
			{
				Node node = stack.Pop();

				if (node.TriangleCount > 0)
				{
					for (int i = 0; i < node.TriangleCount; i++)
					{
						int triIndex = node.StartIndex + i;
						Triangle tri = Triangles[triIndex];
						var hitInfo = RayTriangle(rayOrigin, rayDir, tri);
						if (hitInfo.hit)
						{
							if (hitInfo.dst < minDst)
							{
								minDst = hitInfo.dst;
								hitTriangleIndex = triIndex;
								hitPoint = rayOrigin + rayDir * minDst;
								backface = hitInfo.backface;
								normal = hitInfo.normal;
							}
						}
					}
				}
				else
				{
					Node childA = nodeList.Nodes[node.StartIndex];
					Node childB = nodeList.Nodes[node.StartIndex + 1];

					float dstA = RayBoundingBox(rayOrigin, rayDir, childA).dst;
					float dstB = RayBoundingBox(rayOrigin, rayDir, childB).dst;

					if (dstA > dstB)
					{
						if (dstA < minDst) stack.Push(childA);
						if (dstB < minDst) stack.Push(childB);
					}
					else
					{
						if (dstB < minDst) stack.Push(childB);
						if (dstA < minDst) stack.Push(childA);
					}
				}
			}

			return (hitTriangleIndex != -1, minDst, hitPoint, backface, normal.normalized);
		}

		public static (bool hit, float dst) RayBoundingBox(Vector3 rayOrigin, Vector3 rayDir, Node node)
		{
			float invDirX = rayDir.x == 0 ? float.PositiveInfinity : 1 / rayDir.x;
			float invDirY = rayDir.y == 0 ? float.PositiveInfinity : 1 / rayDir.y;
			float invDirZ = rayDir.z == 0 ? float.PositiveInfinity : 1 / rayDir.z;

			float tx1 = (node.MinX - rayOrigin.x) * invDirX;
			float tx2 = (node.MaxX - rayOrigin.x) * invDirX;
			float tmin = Mathf.Min(tx1, tx2);
			float tmax = Mathf.Max(tx1, tx2);

			float ty1 = (node.MinY - rayOrigin.y) * invDirY;
			float ty2 = (node.MaxY - rayOrigin.y) * invDirY;
			tmin = Mathf.Max(tmin, Mathf.Min(ty1, ty2));
			tmax = Mathf.Min(tmax, Mathf.Max(ty1, ty2));

			float tz1 = (node.MinZ - rayOrigin.z) * invDirZ;
			float tz2 = (node.MaxZ - rayOrigin.z) * invDirZ;
			tmin = Mathf.Max(tmin, Mathf.Min(tz1, tz2));
			tmax = Mathf.Min(tmax, Mathf.Max(tz1, tz2));

			bool hit = tmax >= tmin && tmax > 0;
			float dst = tmin > 0 ? tmin : tmax;
			if (!hit) dst = Mathf.Infinity;
			return (hit, dst);
		}

		static (bool hit, float dst, bool backface, Vector3 normal) RayTriangle(Vector3 rayOrigin, Vector3 rayDir, Triangle tri)
		{
			Vector3 edgeAB = tri.B - tri.A;
			Vector3 edgeAC = tri.C - tri.A;
			Vector3 ao = rayOrigin - tri.A;
			Vector3 dao = Vector3.Cross(ao, rayDir);
			Vector3 faceVector = Vector3.Cross(edgeAB, edgeAC);

			float determinant = -Vector3.Dot(rayDir, faceVector);
			float invDet = 1 / determinant;

			// Calculate dst to triangle & barycentric coordinates of intersection point
			float dst = Vector3.Dot(ao, faceVector) * invDet;
			float u = Vector3.Dot(edgeAC, dao) * invDet;
			float v = -Vector3.Dot(edgeAB, dao) * invDet;
			float w = 1 - u - v;

			// Initialize hit info
			bool backface = determinant < 0;
			bool hit = Mathf.Abs(determinant) >= 1E-8 && dst >= 0 && u >= 0 && v >= 0 && w >= 0;
			return (hit, dst, backface, faceVector);
		}
		
		// ---- Structures and classes ----

		public struct Node
		{
			public float MinX;
			public float MinY;
			public float MinZ;
			public float MaxX;
			public float MaxY;
			public float MaxZ;

			// Index of first child (if triangle count is negative) otherwise index of first triangle
			public int StartIndex;
			public int TriangleCount;


			public Node(float minX, float minY, float minZ, float maxX, float maxY, float maxZ, int startIndex, int triCount)
			{
				MinX = minX;
				MinY = minY;
				MinZ = minZ;
				MaxX = maxX;
				MaxY = maxY;
				MaxZ = maxZ;
				StartIndex = startIndex;
				TriangleCount = triCount;
			}
		}

		public readonly struct BVHTriangle
		{
			public readonly float CentreX;
			public readonly float CentreY;
			public readonly float CentreZ;
			public readonly float MinX;
			public readonly float MinY;
			public readonly float MinZ;
			public readonly float MaxX;
			public readonly float MaxY;
			public readonly float MaxZ;

			public readonly int Index;

			public BVHTriangle(Vector3 a, Vector3 b, Vector3 c, int index)
			{
				float ax = a.x;
				float ay = a.y;
				float az = a.z;
				float bx = b.x;
				float by = b.y;
				float bz = b.z;
				float cx = c.x;
				float cy = c.y;
				float cz = c.z;

				CentreX = (ax + bx + cx) / 3;
				CentreY = (ay + by + cy) / 3;
				CentreZ = (az + bz + cz) / 3;
				MinX = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx);
				MinY = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy);
				MinZ = az < bz ? (az < cz ? az : cz) : (bz < cz ? bz : cz);
				MaxX = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx);
				MaxY = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy);
				MaxZ = az > bz ? (az > cz ? az : cz) : (bz > cz ? bz : cz);
				Index = index;
			}
		}

		public class NodeList
		{
			public Node[] Nodes = new Node[256];
			int Index;

			public int Add(Node node)
			{
				if (Index >= Nodes.Length)
				{
					Array.Resize(ref Nodes, Nodes.Length * 2);
				}

				int nodeIndex = Index;
				Nodes[Index++] = node;
				return nodeIndex;
			}

			public int NodeCount => Index;
		}

		[System.Serializable]
		public class BuildStats
		{
			public int TimeMs;
			public int TriangleCount;
			public int TotalNodeCount;
			public int LeafNodeCount;

			public int LeafDepthMax;
			public int LeafDepthMin = int.MaxValue;
			public int LeafDepthSum;

			public int LeafMaxTriCount;
			public int LeafMinTriCount = int.MaxValue;
			public Quality Quality;

			public BuildStats(Quality quality)
			{
				this.Quality = quality;
			}

			public void RecordNode(int depth, bool isLeaf, int triCount = 0)
			{
				TotalNodeCount++;

				if (isLeaf)
				{
					LeafNodeCount++;
					LeafDepthSum += depth;
					LeafDepthMin = Mathf.Min(LeafDepthMin, depth);
					LeafDepthMax = Mathf.Max(LeafDepthMax, depth);
					TriangleCount += triCount;

					LeafMaxTriCount = Mathf.Max(LeafMaxTriCount, triCount);
					LeafMinTriCount = Mathf.Min(LeafMinTriCount, triCount);
				}
			}


			public override string ToString()
			{
				var sb = new StringBuilder();
				sb.AppendLine($"Time (BVH): {TimeMs} ms (quality = {Quality})");
				sb.AppendLine($"Triangles: {TriangleCount}");
				sb.AppendLine($"Node Count: {TotalNodeCount}");
				sb.AppendLine($"Leaf Count: {LeafNodeCount}");
				sb.AppendLine($"Leaf Depth:");
				sb.AppendLine($" - Min: {LeafDepthMin}");
				sb.AppendLine($" - Max: {LeafDepthMax}");
				sb.AppendLine($" - Mean: {LeafDepthSum / (float)LeafNodeCount:0.####}");
				sb.AppendLine($"Leaf Tris:");
				sb.AppendLine($" - Min: {LeafMinTriCount}");
				sb.AppendLine($" - Max: {LeafMaxTriCount}");
				sb.AppendLine($" - Mean: {TriangleCount / (float)LeafNodeCount:0.####}");
				// sb.AppendLine($"Leaf Nodes:")

				return sb.ToString();
			}
		}


		public readonly struct Triangle
		{
			public readonly Vector3 A;
			public readonly Vector3 B;
			public readonly Vector3 C;

			public readonly Vector3 normalA;
			public readonly Vector3 normalB;
			public readonly Vector3 normalC;

			public Triangle(Vector3 posA, Vector3 posB, Vector3 posC, Vector3 normalA, Vector3 normalB, Vector3 normalC)
			{
				this.A = posA;
				this.B = posB;
				this.C = posC;
				this.normalA = normalA;
				this.normalB = normalB;
				this.normalC = normalC;
			}
		}
	}
}